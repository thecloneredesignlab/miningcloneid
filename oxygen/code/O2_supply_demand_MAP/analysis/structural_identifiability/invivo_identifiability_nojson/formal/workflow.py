#!/usr/bin/env python3
"""Owned-job submission, coherent TSV state and integer formal-phase allocation limits."""
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
FIELDS=['task','role','mode','release','job_id','comment','state','exit_code','cores','wall_seconds','reserved_core_seconds','actual_core_seconds','intent_utc','receipt','spec_sha256','dispatch_sha256']
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
def prior_debit(study):
    rows=read_table(study/'prior/jobs.tsv')
    records={r.split('|')[0]:r.split('|') for r in (study/'prior/raw_accounting.txt').read_text().splitlines() if r}
    if len(rows)!=29 or len({r['job_id'] for r in rows})!=29 or not all(closed(r) for r in rows): raise ValueError('Prior ledger is not closed and unique')
    debit=sum(int(records[r['job_id']][3])*int(records[r['job_id']][4]) for r in rows)
    if debit!=30491: raise ValueError('Prior allocated debit mismatch')
    return debit
def slurm_time(seconds):
    if seconds<=0 or seconds%60: raise ValueError('Wall time must be positive whole minutes')
    return '{:02d}:{:02d}:{:02d}'.format(seconds//3600,(seconds%3600)//60,seconds%60)
def resolve_matches(rows,candidates):
    for j in rows:
        if j['job_id']: continue
        matches=sorted(set(job_id for job_id,comment in candidates if comment==j['comment'] and job_id.isdigit()))
        if len(matches)!=1: raise ValueError('Uncertain submission has '+str(len(matches))+' matches; no retry is allowed')
        j.update(job_id=matches[0],state='RECONCILED')
    return rows
def check_budget(rows,budget,role,prior=30491):
    b=next(x for x in budget if x['role']==role);request=int(b['cores'])*int(b['wall_seconds'])
    same=[j for j in rows if j['role']==role]
    if len(same)>=int(b['attempts']): raise ValueError('Role attempt cap reached')
    if sum(map(charge,same))+request>int(b['cap_core_seconds']): raise ValueError('Role allocation cap reached')
    total=sum(map(charge,rows))+request
    if total>648000: raise ValueError('Formal 180 core-hour cap reached')
    if prior+total>691200: raise ValueError('Cumulative 192 core-hour cap reached')
    return b,request
def reconcile(study,rows):
    owned=[x['job_id'] for x in rows if x['job_id']]
    if not owned: return
    result=subprocess.run([str(SLURM/'sacct'),'-n','-P','-j',','.join(owned),'--format=JobIDRaw,State,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS'],check=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    atomic_text(study/'tables/raw_accounting.txt',result.stdout)
    records={r.split('|')[0]:r.split('|') for r in result.stdout.splitlines() if r}
    for j in rows:
        if not j['job_id']: continue
        r=records.get(j['job_id'])
        if not r: continue
        j['state']=r[1]; j['exit_code']=r[2]
        if closed(j):
            if int(r[3])!=int(j['cores']): raise ValueError('Unexpected allocated core count')
            actual=int(r[3])*int(r[4])
            if actual>int(j['reserved_core_seconds']): raise ValueError('Actual allocation exceeded reservation')
            j['actual_core_seconds']=str(actual)
def gate(study,row,release,rows):
    mode=row['mode']
    if mode in {'engineering_parallel','engineering_serial','reference'}: return
    required=['engineering_gate.tsv']
    if row['role'].startswith('density_') or row['role'].startswith('reserve_density_'): required.append('baseline_gate.tsv')
    for name in required:
        file=study/'tables'/name
        specific=study/'tables'/('engineering_gate_'+release+'.tsv')
        if name=='engineering_gate.tsv' and specific.exists(): file=specific
        values={r['key']:r['value'] for r in read_table(file)}
        if values.get('pass')!='true' or values.get('release')!=release or values.get('input_manifest_sha256')!=sha(study/'manifests/input_manifest.tsv') or values.get('release_manifest_sha256')!=sha(study/'releases'/release/'manifest.tsv'): raise ValueError('Required current-release gate failed: '+name)
        if release=='formal02' and name=='engineering_gate.tsv':
            if values.get('parent_engineering_gate_sha256')!=sha(study/'tables/engineering_gate.tsv'): raise ValueError('Parent engineering proof changed')
            verify_manifest(study,study/'tables/release_bridge_formal02.tsv')
        if name=='baseline_gate.tsv': verify_manifest(study,study/'tables/baseline_gate_payloads.tsv')
    bytask={j['task']:j for j in rows}
    for field in ('source_task','required_check'):
        task=row[field]
        if not task: continue
        j=bytask.get(task)
        if not j or j['state']!='COMPLETED' or j['exit_code']!='0:0': raise ValueError('Dependency has not closed successfully: '+task)
        verify_payload(study/('tests' if j['mode'] in {'fresh','directions','profile_fresh','reference','engineering_parallel','engineering_serial'} else 'fits')/task)

def validate_override(plan,overrides):
    if not overrides: return {}
    if plan['role']!='reserve_baseline_global' or plan['mode']!='global': raise ValueError('Override is outside the reserve-global role')
    if set(overrides)!={'init_mode','start_seed'}: raise ValueError('Only family and archived center may be selected')
    mode=overrides['init_mode'];seed=str(overrides['start_seed'])
    if mode not in ('hybrid','uniform') or not seed.isdigit(): raise ValueError('Invalid initialization selector')
    if mode=='uniform' and int(seed)!=0: raise ValueError('Uniform initialization cannot insert an endpoint')
    if mode=='hybrid' and int(seed) not in (25,74,366,292,392): raise ValueError('Unverified archived center')
    return dict(override_init_mode=mode,override_start_seed=seed)

def submit(study,task,release='formal01',overrides=None):
    if not re.fullmatch('[A-Za-z0-9_]+',task) or not re.fullmatch('[A-Za-z0-9_]+',release): raise ValueError('Invalid task/release identity')
    s=State(study)
    with s.locked():
        rows=s.load();reconcile(study,rows)
        plan=next((r for r in read_table(study/'manifests/task_plan.tsv') if r['task']==task),None)
        if plan is None: raise ValueError('Task is outside the frozen registry')
        override_fields=validate_override(plan,overrides)
        if override_fields and release!='formal02': raise ValueError('Explicit initialization selection requires the versioned dispatch release')
        existing=next((j for j in rows if j['task']==task),None)
        if existing:
            if not existing['job_id']: raise ValueError('Task submission uncertain; reconcile rather than resubmit')
            if existing['mode']!=plan['mode'] or existing['role']!=plan['role'] or existing['release']!=release: raise ValueError('Existing task identity differs')
            if override_fields:
                dispatched={r['key']:r['value'] for r in read_table(study/'dispatch'/task/'context.tsv')}
                if any(dispatched.get(k)!=v for k,v in override_fields.items()): raise ValueError('Existing initialization selection differs')
            s.save(rows);print('Existing task; no duplicate:',task,existing['job_id'],flush=True);return
        if any(not j['job_id'] for j in rows): raise ValueError('Uncertain receipt blocks new submissions')
        if any(closed(j) and (j['state']!='COMPLETED' or j['exit_code']!='0:0') for j in rows): raise ValueError('Unexpected failed allocation blocks dependent submissions')
        verify_manifest(study,study/'manifests/input_manifest.tsv')
        verify_manifest(study/'releases'/release,study/'releases'/release/'manifest.tsv')
        for directory in ('manifests','source','reference','prior'): audit_formats(study/directory)
        audit_formats(study/'releases'/release);gate(study,plan,release,rows)
        budget,request=check_budget(rows,read_table(study/'manifests/budget.tsv'),plan['role'],prior_debit(study))
        dispatch=study/'dispatch'/task;dispatch.mkdir(parents=True,exist_ok=False)
        context=dict(task=task,role=plan['role'],mode=plan['mode'],policy=plan['policy'],release=release,
            input_manifest_sha256=sha(study/'manifests/input_manifest.tsv'),release_manifest_sha256=sha(study/'releases'/release/'manifest.tsv'),
            task_specs_sha256=sha(study/'manifests/task_specs.rds'),selection_origin='')
        context.update(override_fields)
        if plan['selector'] in {'baseline_best','baseline_probe','density_selected'}:
            name=('baseline_'+plan['policy'] if plan['selector']=='baseline_best' else 'baseline_'+plan['policy']+'_probe_best' if plan['selector']=='baseline_probe' else 'density_'+plan['policy']+'_'+plan['density_index']+'_selected')+'.rds'
            selection=study/'tables'/name
            (dispatch/'selection.rds').write_bytes(selection.read_bytes());context['selection_origin']='tables/'+name
        table(dispatch/'context.tsv',[dict(key=k,type='string',value=v,missing='FALSE') for k,v in context.items()],['key','type','value','missing'])
        table(dispatch/'manifest.tsv',[dict(path=p.name,sha256=sha(p)) for p in sorted(dispatch.iterdir())],['path','sha256'])
        comment='ivi:'+study.name+':'+task
        j=dict(task=task,role=plan['role'],mode=plan['mode'],release=release,job_id='',comment=comment,state='INTENT',exit_code='',cores=budget['cores'],wall_seconds=budget['wall_seconds'],reserved_core_seconds=str(request),actual_core_seconds='',intent_utc=utc(),receipt='',spec_sha256=sha(study/'manifests/task_specs.rds'),dispatch_sha256=sha(dispatch/'manifest.tsv'))
        rows.append(j);s.save(rows)
        cmd=[str(SLURM/'sbatch'),'--parsable','--partition=red','--qos=xxlarge','--no-requeue','--nodes=1','--ntasks=1','--cpus-per-task='+budget['cores'],'--mem='+('16G' if int(budget['cores'])>1 else '8G'),'--time='+slurm_time(int(budget['wall_seconds'])),'--job-name=ivi_f010554Z_'+task,'--comment='+comment,'--chdir='+str(study),'--output='+str(study/'logs'/(task+'_%j.log')),str(study/'releases'/release/'batch.sh'),str(study),plan['mode'],task,release]
        atomic_text(study/'logs'/(task+'_submission.txt'),'\n'.join(cmd)+'\n')
        result=subprocess.run(cmd,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
        atomic_text(study/'logs'/(task+'_receipt.txt'),result.stdout+'\n'+result.stderr)
        receipt=result.stdout.strip().split(';')[0]
        if result.returncode or not receipt.isdigit(): raise RuntimeError('Submission uncertain; reservation retained, no automatic retry')
        j.update(job_id=receipt,state='SUBMITTED',receipt=str(study/'logs'/(task+'_receipt.txt')));s.save(rows)
        print(task,receipt,str(request)+' core-seconds',flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['submit','status','audit','resolve']);p.add_argument('study',type=Path);p.add_argument('--task');p.add_argument('--release',default='formal01');a=p.parse_args();study=a.study.resolve()
    if a.command=='submit': submit(study,a.task,a.release)
    elif a.command=='audit': audit_formats(study);print('Operational format audit passed')
    elif a.command=='resolve':
        s=State(study)
        with s.locked():
            rows=s.load();candidates=[]
            for j in rows:
                if j['job_id']: continue
                name='ivi_f010554Z_'+j['task']
                for cmd in ([str(SLURM/'squeue'),'-u','4482173','--name='+name,'-h','-o','%i|%k'],[str(SLURM/'sacct'),'-u','4482173','-S',j['intent_utc'][:10],'-n','-P','--name='+name,'--format=JobIDRaw,Comment%200']):
                    r=subprocess.run(cmd,check=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    candidates.extend(tuple(line.rstrip('|').split('|',1)) for line in r.stdout.splitlines() if '|' in line.rstrip('|'))
            resolve_matches(rows,candidates);s.save(rows)
    else:
        s=State(study)
        with s.locked(): rows=s.load();reconcile(study,rows);s.save(rows)
        table(study/'tables/jobs.tsv',rows,FIELDS)
        actual=sum(int(j['actual_core_seconds']) for j in rows if closed(j) and j['actual_core_seconds'])
        print('Formal actual core-hours:',actual/3600,'formal used-plus-committed:',sum(map(charge,rows))/3600,'cumulative used-plus-committed:',(prior_debit(study)+sum(map(charge,rows)))/3600)
        for j in rows: print(j['task'],j['job_id'],j['state'],j['exit_code'],j['actual_core_seconds'])
if __name__=='__main__': main()
