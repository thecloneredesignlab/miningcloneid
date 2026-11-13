#!/usr/bin/env python3
"""Owned-job submission and coherent integer accounting without aggregate caps."""
import argparse
import contextlib
import csv
import fcntl
import hashlib
import io
import os
import re
import subprocess
import tempfile
import uuid
from pathlib import Path

SLURM=Path('/cm/shared/apps/slurm/current/bin')
FIELDS=['task','role','mode','release','job_id','comment','state','exit_code','cores','wall_seconds','reserved_core_seconds','actual_core_seconds','intent_utc','receipt','spec_sha256','dispatch_sha256','policy','root_id','root_family','source_task','selector','density_index','kind']
TERMINAL={'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'}

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''): h.update(chunk)
    return h.hexdigest()
def read_table(p):
    with Path(p).open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def atomic_text(p,s):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.pending_',dir=p.parent)
    try:
        with os.fdopen(fd,'w') as f: f.write(s); f.flush(); os.fsync(f.fileno())
        os.replace(tmp,p)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)
def table(p,rows,fields):
    f=io.StringIO(); w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows); atomic_text(p,f.getvalue())
def utc():
    from datetime import datetime,timezone
    return datetime.now(timezone.utc).isoformat()
def verify_manifest(root,path):
    for row in read_table(path):
        if sha(root/row['path'])!=row['sha256']: raise ValueError('Source/control fingerprint mismatch: '+row['path'])
def verify_payload(p):
    p=Path(p)
    if not (p/'COMPLETE').exists() or (p/'COMPLETE').read_text().strip()!=sha(p/'payload_manifest.tsv'): raise ValueError('Incomplete task payload')
    verify_manifest(p,p/'payload_manifest.tsv')
def audit_formats(root):
    errors=[]
    for p in Path(root).rglob('*'):
        if not p.is_file(): continue
        if p.suffix.lower() in {'.json','.jsonl','.ndjson','.js','.mjs','.cjs','.java','.jar','.class'}: errors.append(str(p))
        if p.suffix.lower() in {'.csv','.tsv'}:
            with p.open(newline='') as f:
                for row in csv.reader(f,delimiter='\t' if p.suffix=='.tsv' else ','):
                    if any(re.match(r'^\s*(?:\{\s*"|\[\s*[\{"0-9-])',v) for v in row): errors.append(str(p))
    if errors: raise ValueError('Forbidden file/transport: '+repr(errors[:8]))
class State:
    def __init__(self,study): self.study=Path(study)
    @contextlib.contextmanager
    def locked(self):
        (self.study/'state').mkdir(exist_ok=True)
        with (self.study/'state/LOCK').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            yield
            fcntl.flock(f,fcntl.LOCK_UN)
    def load(self):
        p=self.study/'state/CURRENT'
        if not p.exists(): return []
        g=p.read_text().strip()
        if not re.fullmatch('g_[0-9a-f]+',g): raise ValueError('Invalid state generation')
        d=self.study/'state'/g; verify_manifest(d,d/'manifest.tsv'); return read_table(d/'jobs.tsv')
    def save(self,rows):
        g='g_'+uuid.uuid4().hex; d=self.study/'state'/g; d.mkdir()
        table(d/'jobs.tsv',rows,FIELDS); table(d/'manifest.tsv',[dict(path='jobs.tsv',sha256=sha(d/'jobs.tsv'))],['path','sha256'])
        atomic_text(self.study/'state/CURRENT',g+'\n')
def closed(j): return j['state'].split('+')[0].split()[0] in TERMINAL
def charge(j):
    return int(j['actual_core_seconds']) if closed(j) and j['actual_core_seconds'] else int(j['reserved_core_seconds'])

def kv(path): return {r['key']:r['value'] for r in read_table(path)}
def typed(path,values):
    rows=[]
    for k,v in values.items():
        typ='boolean' if isinstance(v,bool) else 'integer' if isinstance(v,int) else 'double' if isinstance(v,float) else 'string'
        value=str(v).lower() if isinstance(v,bool) else format(v,'.17g') if isinstance(v,float) else str(v)
        rows.append(dict(key=k,type=typ,value=value,missing='FALSE'))
    table(path,rows,['key','type','value','missing'])
def slurm_time(seconds): return '{:02d}:{:02d}:{:02d}'.format(seconds//3600,(seconds%3600)//60,seconds%60)
def reconcile(study,rows):
    ids=[j['job_id'] for j in rows if j['job_id']]
    if not ids:return
    result=subprocess.run([str(SLURM/'sacct'),'-n','-P','-j',','.join(ids),'--format=JobIDRaw,State,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS'],check=True,universal_newlines=True,stdout=subprocess.PIPE)
    atomic_text(study/'tables/raw_accounting.txt',result.stdout)
    raw={r.split('|')[0]:r.split('|') for r in result.stdout.splitlines() if r}
    for j in rows:
        r=raw.get(j['job_id'])
        if not r:continue
        j.update(state=r[1],exit_code=r[2])
        if closed(j):
            if int(r[3])!=int(j['cores']):raise ValueError('Allocated CPU identity changed')
            j['actual_core_seconds']=str(int(r[3])*int(r[4]))
            if int(j['actual_core_seconds'])>int(j['reserved_core_seconds']):raise ValueError('Job wall reservation exceeded')
def snapshot(study):
    s=State(study)
    with s.locked(): rows=s.load();reconcile(study,rows);s.save(rows)
    table(study/'tables/jobs.tsv',rows,FIELDS)
    used=sum(map(charge,rows));actual=sum(int(j['actual_core_seconds']) for j in rows if closed(j) and j['actual_core_seconds'])
    active=[j for j in rows if not closed(j)]
    atomic_text(study/'STATUS.md','# Improvement-driven in vivo continuation\n\nUpdated '+utc()+'.\n\n'+str(len(rows))+' owned parents; '+str(len(active))+' active/uncertain.\n\nNew actual '+str(actual)+' core-seconds; actual plus live '+str(used)+'. Historical lifetime301846/formal271355remain charged. No aggregate core-hour or attempt ceiling. Culture protected; pinned SIF. Profiles require both baselines.\n\n'+'\n'.join(j['task']+' '+j['job_id']+' '+j['state'] for j in active)+'\n')
    print('new_actual_ch',actual/3600,'lifetime_used_plus_live_ch',(301846+used)/3600,'formal_used_plus_live_ch',(271355+used)/3600)
    for j in active:print('ACTIVE',j['task'],j['job_id'],j['state'])
    return rows
def gate(study,name):
    p=study/'tables'/name
    if not p.exists():raise ValueError('Required gate absent: '+name)
    z=kv(p)
    if z.get('pass')!='true' or z['input_manifest_sha256']!=sha(study/'manifests/input_manifest.tsv'):raise ValueError('Current source gate does not pass')
    verify_manifest(study,study/z['evidence_manifest'])
    return z
def submit(study,task,**options):
    if not re.fullmatch('[A-Za-z0-9_]+',task):raise ValueError('Invalid owned task')
    defaults=dict(mode='local',role='baseline_local',policy='production',root_id='',root_family='archive',selector='import',import_id='',source_task='',density_index=0,seed=53001,init_mode='hybrid',sigma_frac=.15,cores=1,wall_seconds=7200,margin_seconds=600,max_calls=20000,kind='baseline',release='cont01')
    if (study/'tables/execution_release.tsv').exists():defaults['release']=kv(study/'tables/execution_release.tsv')['release']
    defaults.update(options);spec=defaults;release=spec['release']
    state=State(study)
    with state.locked():
        rows=state.load();reconcile(study,rows)
        if any(j['task']==task for j in rows):raise ValueError('Existing task; reconcile without duplicate dispatch')
        if any(not j['job_id'] for j in rows):raise ValueError('Uncertain receipt must be reconciled')
        reviewed=set()
        if (study/'tables/reviewed_failures.tsv').exists():reviewed={r['task'] for r in read_table(study/'tables/reviewed_failures.tsv')}
        if any(closed(j) and (j['state']!='COMPLETED' or j['exit_code']!='0:0') and j['task'] not in reviewed for j in rows):raise ValueError('Unreviewed failure blocks new work')
        verify_manifest(study,study/'manifests/input_manifest.tsv');verify_manifest(study/'releases'/release,study/'releases'/release/'manifest.tsv')
        if spec['kind']!='engineering':
            g=gate(study,'engineering_gate.tsv')
            if g['release']!=release:raise ValueError('Wrong engineering release')
        if spec['kind']=='profile':
            gate(study,'baseline_gate.tsv')
            if any(j['kind']=='baseline' and not closed(j) for j in rows):raise ValueError('Baseline accounting still active')
        if spec['source_task']:
            parent=next((j for j in rows if j['task']==spec['source_task']),None)
            if parent is None or parent['state']!='COMPLETED' or parent['exit_code']!='0:0':raise ValueError('Parent not closed')
            root=study/('fits' if parent['mode'] in ('local','global') else 'tests')/parent['task'];verify_payload(root)
            if spec['root_family']!='profile' and (parent['root_id']!=spec['root_id'] or parent['policy']!=spec['policy']):raise ValueError('Cross-root or cross-policy start rejected')
        if spec['max_calls']<1 or spec['wall_seconds']<1 or spec['cores'] not in (1,8):raise ValueError('Invalid per-window allowance')
        # There is intentionally no total/phase/role/attempt cap here.
        spec['input_manifest_sha256']=sha(study/'manifests/input_manifest.tsv');spec['release_manifest_sha256']=sha(study/'releases'/release/'manifest.tsv')
        dispatch=study/'dispatch'/task;dispatch.mkdir()
        if spec['selector']=='selection':
            name=spec.pop('selection_name');p=study/'tables'/(name+'.rds');(dispatch/'selection.rds').write_bytes(p.read_bytes())
        typed(dispatch/'spec.tsv',spec);table(dispatch/'manifest.tsv',[dict(path=p.name,sha256=sha(p)) for p in sorted(dispatch.iterdir())],['path','sha256'])
        comment='ivc:'+study.name+':'+task
        j=dict(task=task,role=spec['role'],mode=spec['mode'],release=release,job_id='',comment=comment,state='INTENT',exit_code='',cores=str(spec['cores']),wall_seconds=str(spec['wall_seconds']),reserved_core_seconds=str(spec['cores']*spec['wall_seconds']),actual_core_seconds='',intent_utc=utc(),receipt='',spec_sha256=sha(dispatch/'spec.tsv'),dispatch_sha256=sha(dispatch/'manifest.tsv'))
        for k in ('policy','root_id','root_family','source_task','selector','density_index','kind'):j[k]=str(spec[k])
        rows.append(j);state.save(rows)
        cmd=[str(SLURM/'sbatch'),'--parsable','--partition=red','--qos=xxlarge','--no-requeue','--nodes=1','--ntasks=1','--cpus-per-task='+j['cores'],'--mem='+('16G' if spec['cores']==8 else '8G'),'--time='+slurm_time(spec['wall_seconds']),'--job-name=ivc_'+task,'--comment='+comment,'--chdir='+str(study),'--output='+str(study/'logs'/(task+'_%j.log')),str(study/'releases'/release/'batch.sh'),str(study),spec['mode'],task,release]
        atomic_text(study/'logs'/(task+'_submission.txt'),'\n'.join(cmd)+'\n')
        result=subprocess.run(cmd,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
        atomic_text(study/'logs'/(task+'_receipt.txt'),result.stdout+'\n'+result.stderr)
        receipt=result.stdout.strip().split(';')[0]
        if result.returncode or not receipt.isdigit():raise ValueError('Uncertain submission retained; no retry')
        j.update(job_id=receipt,state='SUBMITTED',receipt=str(study/'logs'/(task+'_receipt.txt')));state.save(rows)
        print(task,receipt,flush=True)

def advance(study,wave):
    rows=snapshot(study);by={j['task']:j for j in rows}
    roots=read_table(study/'manifests/roots.tsv')
    if wave=='engineering':
        tasks=[('e_parallel',dict(mode='engineering_parallel',cores=8,root_id='engineering_reference',root_family='reference',selector='reference'))]
        for policy in ('production','data_only'):
            root=next(r for r in roots if r['policy']==policy and r['kind']=='anchor')
            tasks.append(('e_replay_'+policy,dict(mode='replay',policy=policy,root_id=root['root_id'],root_family='archive',selector='import',import_id=root['import_id'])))
            for fixed in (False,True):tasks.append(('e_'+policy+('_fixed' if fixed else '_free'),dict(mode='engineering_serial',policy=policy,root_id=root['root_id'],root_family='archive',selector='import',import_id=root['import_id'],density_index=4 if fixed else 0)))
        version=kv(study/'tables/execution_release.tsv')['release'] if (study/'tables/execution_release.tsv').exists() else 'cont01'
        for task,opts in tasks:
            task=task+'_'+version
            if task not in by:submit(study,task,kind='engineering',role='engineering',wall_seconds=1800,margin_seconds=300,max_calls=500,**opts)
    elif wave=='initial':
        for r in roots:
            task='b_'+r['import_id']+'_w0001'
            if task not in by:submit(study,task,policy=r['policy'],root_id=r['root_id'],root_family=r['root_family'],import_id=r['import_id'])
    elif wave=='checks':
        for j in rows:
            if j['kind']!='baseline' or j['mode']!='local' or j['state']!='COMPLETED' or j['exit_code']!='0:0':continue
            prefix=j['task'];fresh=prefix+'_fresh';direction=prefix+'_directions'
            if fresh not in by:submit(study,fresh,mode='fresh',role='fresh',policy=j['policy'],root_id=j['root_id'],root_family=j['root_family'],selector='task',source_task=prefix,wall_seconds=900,margin_seconds=120,max_calls=500,kind=j['kind'],density_index=int(j['density_index']))
            elif by[fresh]['state']=='COMPLETED' and direction not in by:submit(study,direction,mode='directions',role='directions',policy=j['policy'],root_id=j['root_id'],root_family=j['root_family'],selector='task',source_task=prefix,wall_seconds=900,margin_seconds=120,max_calls=500,kind=j['kind'],density_index=int(j['density_index']))
    snapshot(study)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('wave',choices=('status','engineering','initial','checks'));a=p.parse_args();advance(a.study.resolve(),a.wave)
