#!/usr/bin/env python3
"""Transactional TSV-only study control. All analysis is native R/C++.
No application concurrency throttle. Scheduler closure and complete payloads are
separate gates. Exact native DE continuation is deliberately not implemented.
"""
import argparse
import contextlib
import csv
import datetime
import fcntl
import hashlib
import io
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import uuid

SLURM = Path('/cm/shared/apps/slurm/current/bin')
JOB_FIELDS = ['slot','task','mode','role','phase','job_id','token','state','exit_code','cpus','wall_seconds','reserved_ch','actual_ch','submitted_utc']
TASK_FIELDS = ['task_id','mode','phase','origin','family','case','source_fit','seed','alternative','min_generations','max_generations']
TERMINAL = {'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'}

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()
def atomic_text(path,text):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.pending_',dir=path.parent)
    try:
        with os.fdopen(fd,'w',encoding='utf-8',newline='') as f: f.write(text); f.flush(); os.fsync(f.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)
def write_table(path,rows,fields):
    f=io.StringIO(newline=''); w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n',extrasaction='raise'); w.writeheader(); w.writerows(rows); atomic_text(path,f.getvalue())
def read_table(path):
    with open(path,newline='',encoding='utf-8') as f: return list(csv.DictReader(f,delimiter='\t'))
def write_kv(path,values):
    rows=[]
    for k,v in values.items():
        typ='boolean' if isinstance(v,bool) else 'integer' if isinstance(v,int) else 'double' if isinstance(v,float) else 'string'
        s='NA' if v is None else ('true' if v else 'false') if isinstance(v,bool) else format(v,'.17g') if isinstance(v,float) else str(v)
        rows.append(dict(key=k,type=typ,value=s))
    write_table(path,rows,['key','type','value'])
def read_kv(path):
    ans={}
    for r in read_table(path):
        v=r['value']; t=r['type']
        if r['key'] in ans: raise ValueError('Duplicate typed key')
        ans[r['key']]=None if v=='NA' else {'boolean':lambda: {'true':True,'false':False}[v],'integer':lambda:int(v),'double':lambda:float(v),'string':lambda:v}[t]()
    return ans

def audit(study):
    errors=[]
    for p in Path(study).rglob('*'):
        if not p.is_file(): continue
        if p.suffix.lower() in {'.json','.jsonl','.ndjson','.js','.mjs','.cjs','.java','.jar','.class'}: errors.append((str(p),'forbidden_file'))
        if p.suffix.lower() in {'.tsv','.csv'}:
            try:
                with p.open(newline='',encoding='utf-8') as f:
                    for row in csv.reader(f,delimiter='\t' if p.suffix=='.tsv' else ','):
                        for v in row:
                            if re.match(r'^\s*(?:\{\s*"|\[\s*[\{"0-9-])',v): errors.append((str(p),'embedded_serialization'))
            except (UnicodeError,csv.Error): errors.append((str(p),'invalid_text'))
    if errors: raise ValueError('Format audit failed: '+repr(errors[:8]))

def verify_code(study):
    for r in read_table(study/'manifests/control_manifest.tsv'):
        if sha(study/r['path'])!=r['sha256']: raise ValueError('Stale study control')
    for r in read_table(study/'manifests/code_manifest.tsv'):
        if sha(study/'implementation'/r['path'])!=r['sha256']: raise ValueError('Stale implementation')

class State:
    def __init__(self,study): self.study=Path(study)
    @contextlib.contextmanager
    def locked(self):
        self.study.joinpath('state').mkdir(exist_ok=True)
        with open(self.study/'state/LOCK','a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            yield
            fcntl.flock(f,fcntl.LOCK_UN)
    def load(self):
        ptr=self.study/'state/CURRENT'
        if not ptr.exists(): return {'jobs':[],'tasks':[],'events':[]}
        generation=ptr.read_text().strip()
        if not re.fullmatch(r'g_[0-9a-f]+',generation): raise ValueError('Invalid state pointer')
        d=self.study/'state'/generation
        for r in read_table(d/'manifest.tsv'):
            if sha(d/r['path'])!=r['sha256']: raise ValueError('State generation corrupted')
        return {n:read_table(d/(n+'.tsv')) for n in ['jobs','tasks','events']}
    def save(self,data):
        name='g_'+uuid.uuid4().hex; tmp=self.study/'state'/('.'+name); tmp.mkdir()
        for n,fields in [('jobs',JOB_FIELDS),('tasks',TASK_FIELDS),('events',['utc','event','detail'])]: write_table(tmp/(n+'.tsv'),data[n],fields)
        write_table(tmp/'manifest.tsv',[{'path':n+'.tsv','sha256':sha(tmp/(n+'.tsv'))} for n in ['jobs','tasks','events']],['path','sha256'])
        os.rename(tmp,self.study/'state'/name); atomic_text(self.study/'state/CURRENT',name+'\n')
    def event(self,data,event,detail): data['events'].append(dict(utc=now(),event=event,detail=detail))

def task_record(tid,mode,phase,origin,family='',source=None,seed=0,alternative=False,case='unrestricted',minimum=300,maximum=1000):
    return dict(task_id=tid,mode=mode,phase=phase,origin=origin,family=family,case=case,source_fit=source,
                seed=seed,alternative=alternative,min_generations=minimum,max_generations=maximum)
def add_task(state,data,t):
    if any(r['task_id']==t['task_id'] for r in data['tasks']): return
    path=state.study/'tasks'/(t['task_id']+'.tsv')
    if path.exists():
        if read_kv(path)!=t: raise ValueError('Task identity collision')
    else: write_kv(path,t)
    data['tasks'].append({k:('NA' if v is None else str(v)) for k,v in t.items()})
    state.save(data)

def job_ok(j): return j is not None and j['state']=='COMPLETED' and j['exit_code']=='0:0'
def is_closed(j): return j['state'].split('+')[0] in TERMINAL
def get_job(data,tid,mode): return next((j for j in data['jobs'] if j['task']==tid and j['mode']==mode),None)
def resources(role):
    return {'global':(8,43200,16),'local':(1,21600,8),'screen':(1,3600,8),'check':(1,1800,8),'validate':(8,10800,16),'controller':(1,43200,4)}[role]
def committed(data,role=None):
    return sum(float(j['actual_ch']) if is_closed(j) and j['actual_ch'] else float(j['reserved_ch']) for j in data['jobs'] if role is None or (j['role'] not in {'global','local'} if role=='nonfit' else j['role']==role))
def budget_check(data,role,cpus,seconds):
    request=cpus*seconds/3600
    if committed(data)+request>1440+1e-9: raise ValueError('Total allocation ceiling reached')
    if role not in {'global','local'} and committed(data,'nonfit')+request>120+1e-9: raise ValueError('Non-fit allocation ceiling reached')
    if role=='local' and sum(j['role']=='local' for j in data['jobs'])>=28: raise ValueError('Local attempt cap reached')
    if role=='global' and sum(j['role']=='global' for j in data['jobs'])>=12: raise ValueError('Global attempt cap reached')
    return request

def submit(state,data,project,tid,mode,role,phase):
    old=get_job(data,tid,mode)
    if old: return old
    cpus,seconds,mem=resources(role); request=budget_check(data,role,cpus,seconds)
    token='bnj_'+uuid.uuid4().hex[:18]
    row=dict(slot=uuid.uuid4().hex,task=tid,mode=mode,role=role,phase=phase,job_id='',token=token,state='INTENT',exit_code='',cpus=str(cpus),wall_seconds=str(seconds),reserved_ch=format(request,'.17g'),actual_ch='',submitted_utc=now())
    data['jobs'].append(row); state.event(data,'submit_intent',token); state.save(data)
    hours=seconds//3600; minutes=(seconds%3600)//60
    command=[str(SLURM/'sbatch'),'--parsable','--no-requeue','--partition=red','--qos=xxlarge',f'--cpus-per-task={cpus}',f'--mem={mem}G',f'--time={hours:02d}:{minutes:02d}:00',f'--job-name={token}',f'--output={state.study}/logs/{token}_%j.out',f'--error={state.study}/logs/{token}_%j.err',str(state.study/'implementation/batch.sh'),str(project),str(state.study),mode,tid]
    write_table(state.study/'submissions'/(token+'_arguments.tsv'),[{'order':i,'argument':a} for i,a in enumerate(command)],['order','argument'])
    try:
        p=subprocess.run(command,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
        atomic_text(state.study/'submissions'/(token+'_stdout.txt'),p.stdout)
        atomic_text(state.study/'submissions'/(token+'_stderr.txt'),p.stderr)
        match=re.fullmatch(r'\s*(\d+)(?:;[^\s]+)?\s*',p.stdout)
        if p.returncode or not match: raise RuntimeError('Ambiguous submission response; reconcile token before retry')
        row['job_id']=match.group(1); row['state']='SUBMITTED'; state.event(data,'submission_receipt',token+' '+row['job_id']); state.save(data)
    except BaseException:
        row['state']='AMBIGUOUS'; state.save(data); raise
    return row

def refresh(state,data):
    uncertain=[j for j in data['jobs'] if not j['job_id'] and not is_closed(j)]
    for j in uncertain:
        p=subprocess.run([str(SLURM/'sacct'),'-X','-nP','--starttime=2026-11-11',f'--name={j["token"]}','--format=JobIDRaw,JobName%80'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,check=True)
        found=[r.split('|')[0] for r in p.stdout.splitlines() if r.split('|')[1]==j['token']]
        if len(found)!=1: raise RuntimeError('Unresolved ambiguous submission: '+j['token'])
        j['job_id']=found[0]; j['state']='SUBMITTED'
    ids=[j['job_id'] for j in data['jobs'] if not is_closed(j)]
    if ids:
        p=subprocess.run([str(SLURM/'sacct'),'-X','-nP','--starttime=2026-11-11','-j',','.join(ids),'--format=JobIDRaw,State,ExitCode,AllocCPUS,ElapsedRaw,CPUTimeRAW'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,check=True)
        records={r.split('|')[0]:r.split('|') for r in p.stdout.splitlines() if '|' in r}
        for j in data['jobs']:
            if j['job_id'] not in records: continue
            r=records[j['job_id']]; status=r[1].split()[0].rstrip('+')
            if is_closed(j): continue
            j['state']=status; j['exit_code']=r[2]
            if status in TERMINAL:
                actual=float(r[5])/3600 if r[5] else int(r[3])*int(r[4])/3600
                j['actual_ch']=format(actual,'.17g'); state.event(data,'scheduler_closed',j['job_id']+' '+status+' '+r[2])
                if status != 'COMPLETED' or r[2] != '0:0':
                    if j['role'] in {'global','local'} and status == 'FAILED':
                        state.save(data); raise ValueError('Unexpected fitting execution failure; dependent work held: '+j['job_id'])
                    error_logs=[]
                    for log in (state.study/'logs').glob(j['token']+'_*'):
                        text=log.read_text(errors='replace'); error_logs.append(text)
                        if re.search(r'checkpoint_write_failed|Missing immutable improvement|Stale or incomplete improvement|Wrong actual CPU allocation|Global workers remained active|Stale (implementation|study|checkpoint)|fingerprint.*(mismatch|changed)|Immutable culture|Wrong provided population|Retained sample contract|score replay mismatch',text,re.I):
                            state.save(data); raise ValueError('Integrity failure in closed job '+j['job_id'])
                    if j['role']=='screen' and status=='FAILED' and 'Initialization unresolved; prospective origin attempt consumed' not in ''.join(error_logs):
                        state.save(data); raise ValueError('Unexpected screening execution failure; dependent work held: '+j['job_id'])
        atomic_text(state.study/'tables/latest_parent_accounting.txt',p.stdout)
    for j in data['jobs']:
        step=state.study/'tables/accounting_steps'/(j['job_id']+'.txt')
        if is_closed(j) and j['job_id'].isdigit() and not step.exists():
            raw=subprocess.run([str(SLURM/'sacct'),'-nP','--starttime=2026-11-11','-j',j['job_id'],'--format=JobIDRaw,State,ExitCode,AllocCPUS,ElapsedRaw,CPUTimeRAW,MaxRSS'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,check=True)
            atomic_text(step,raw.stdout)
    state.save(data)

def accepted(data,study,tid):
    j=get_job(data,tid,'check'); d=study/'fits'/tid
    if not job_ok(j) or not (d/'COMPLETE').exists(): return None
    check=read_kv(d/'artifact_validation.tsv'); result=read_kv(d/'result.tsv')
    if not check['pass'] or check['result_sha256']!=sha(d/'result.tsv') or check['accepted_fit_sha256']!=sha(d/'accepted_fit.rds') or check['payload_manifest_sha256']!=sha(d/'payload_manifest.tsv'): raise ValueError('Stale accepted artifact')
    return result

def launch_ready(state,data,project,t):
    tid=t['task_id']; mode=t['mode']; j=get_job(data,tid,mode)
    if j is None: submit(state,data,project,tid,mode,mode,t['phase']); return False
    if not is_closed(j): return False
    if not job_ok(j): return True
    c=get_job(data,tid,'check')
    if c is None: submit(state,data,project,tid,'check','check',t['phase']); return False
    if not is_closed(c): return False
    if not job_ok(c): raise ValueError('Independent artifact validation failed: '+tid)
    if accepted(data,state.study,tid) is None: raise ValueError('Missing complete validated payload')
    return True

def get_task(data,tid): return read_kv(Path(data['_study'])/'tasks'/(tid+'.tsv'))
def fit_source(study,project,tid): return str((study/'fits'/tid/'accepted_fit.rds').relative_to(project))
def primary_roots(data): return [read_kv(Path(data['_study'])/'tasks'/(r['task_id']+'.tsv')) for r in data['tasks'] if r['mode']=='global' and r['phase']=='primary']

def candidates(data,study,project):
    rows=[]
    for t in data['tasks']:
        if t['phase']=='fixture' or t['mode'] not in {'local','global'}: continue
        r=accepted(data,study,t['task_id'])
        if r:
            r.update(source_fit=fit_source(study,project,t['task_id']),phase=t['phase']); rows.append(r)
    # Immutable old minimum remains eligible for Q_star, independently replayed in Phase A.
    settings=read_kv(study/'settings.tsv')
    rows.append(dict(task_id='historical_reference',origin='review_035',family='',mode='reference',objective=read_kv(study/'tests/replays/review_035_rep1/replay.tsv')['objective'] if (study/'tests/replays/review_035_rep1/replay.tsv').exists() else 3.85522952512099,local_stable=False,
        source_fit=settings['review_relative_path']+'/fits/review_035/accepted_fit.rds',phase='historical'))
    return sorted(rows,key=lambda r:(r['objective'],r['task_id']))

def gate(data,study,project):
    cs=candidates(data,study,project); q=cs[0]['objective']; roots={}
    for r in cs:
        if r['family'] in {'A','B'} and r['mode']=='local' and r['origin'] not in roots: roots[r['origin']]=r
    good=[r for r in roots.values() if r['objective']-q<=.001 and r['local_stable']]
    repeated=[]
    for r in cs:
        if r['phase'] in {'consistency','reserve_local'} and r['mode']=='local' and r['local_stable']:
            t=read_kv(study/'tasks'/(r['task_id']+'.tsv'))
            source=next((c for c in cs if c['source_fit']==t['source_fit']),None)
            if source and source['objective']-r['objective']<1e-4 and r['objective']-q<1e-8: repeated.append(r)
    passed={r['family'] for r in good}=={'A','B'} and len(good)>=2 and bool(repeated)
    return dict(pass_gate=passed,q_star=q,eligible_origins=len(good),families=';'.join(sorted({r['family'] for r in good})),strong_repeat=bool(repeated)),cs,roots

def export_progress(state,data,project,status):
    g,cs,roots=gate(data,state.study,project)
    write_kv(state.study/'tables/recovery_gate.tsv',dict(g,status=status,updated_utc=now(),allocated_used_or_committed=committed(data),nonfit_used_or_committed=committed(data,'nonfit')))
    rows=[]
    for r in cs:
        rows.append(dict({k:r.get(k,'') for k in ['task_id','origin','family','mode','phase','objective','local_stable','source_fit']}, gap=r['objective']-g['q_star']))
    write_table(state.study/'tables/recovery_candidates.tsv',rows,['task_id','origin','family','mode','phase','objective','local_stable','source_fit','gap'])
    write_table(state.study/'tables/jobs.tsv',data['jobs'],JOB_FIELDS)
    atomic_text(state.study/'STATUS.md',f'# Baseline recovery status\n\nStatus: {status}. Updated: {now()}.\n\nLowest verified Q: {g["q_star"]:.17g}. Prospective complete-chain gate: {g["pass_gate"]}. Eligible new origins: {g["eligible_origins"]}. Families: {g["families"] or "none"}.\n\nUsed or committed allocations: {committed(data):.6f} core-hours; hard ceiling: 1440. Non-fit used or committed: {committed(data,"nonfit"):.6f}; cap: 120. Historical global-only recovery remains failed. No profiles or tumor fits are authorized.\n')
    return g,cs,roots

def drive(state,data,project):
    data['_study']=str(state.study)
    refresh(state,data); verify_code(state.study); audit(state.study)
    v=get_job(data,'validation','validate')
    if v is None: submit(state,data,project,'validation','validate','validate','validation'); return 'validation'
    if not is_closed(v): return 'validation'
    if not job_ok(v) or not (state.study/'tests/global_validation.tsv').exists(): raise ValueError('Phase-A global validation failed; recovery submissions blocked')
    gv=read_kv(state.study/'tests/global_validation.tsv')
    if not gv['pass'] or gv['code_sha256']!=sha(state.study/'manifests/code_manifest.tsv'): raise ValueError('Stale global validation')
    lvjob=get_job(data,'validation_local','validate_local')
    if lvjob is None:
        submit(state,data,project,'validation_local','validate_local','check','validation'); return 'separate_one_core_local_validation'
    if not is_closed(lvjob): return 'separate_one_core_local_validation'
    if not job_ok(lvjob) or not (state.study/'tests/local_validation.tsv').exists(): raise ValueError('Phase-A one-core local validation failed')
    lv=read_kv(state.study/'tests/local_validation.tsv')
    if not lv['pass'] or lv['code_sha256']!=gv['code_sha256']: raise ValueError('Stale local validation')
    if not (state.study/'tests/validation.tsv').exists():
        audit(state.study)
        write_kv(state.study/'tests/validation.tsv',{'pass':True,'code_sha256':gv['code_sha256'],'global_job_id':v['job_id'],'local_job_id':lvjob['job_id'],'fresh_process_replays':12,'legacy_endpoint_replays':gv.get('legacy_endpoint_replays',0),'append_only_improvement_archives':True,'actual_eight_worker_fixtures':3,'actual_one_core_local_fixtures':3,'completed_utc':now()})
    vg=read_kv(state.study/'tests/validation.tsv')
    if not vg['pass'] or vg['code_sha256']!=gv['code_sha256']: raise ValueError('Invalid Phase-A gate')
    settings=read_kv(state.study/'settings.tsv')
    hist=['review_001','review_002','review_038','review_039','review_035','review_006']
    finished=[]
    for i,h in enumerate(hist):
        t=task_record('history_'+h,'local','diagnostic',h,source=settings['review_relative_path']+'/fits/'+h+'/accepted_fit.rds',seed=9001+i)
        add_task(state,data,t); finished.append(launch_ready(state,data,project,t))
    if not all(finished): return 'historical_endpoint_diagnostics'
    if not (state.study/'manifests/method_freeze.tsv').exists():
        # Prospective method is frozen before any prospective screening/population evaluation.
        write_kv(state.study/'manifests/method_freeze.tsv',dict(code_sha256=sha(state.study/'manifests/code_manifest.tsv'),settings_sha256=sha(state.study/'settings.tsv'),frozen_utc=now(),diagnostics_closed=True))
    seeds=read_table(state.study/'origins.tsv')
    roots=[r for r in seeds if r['stage']=='primary']
    all_done=True
    for origin in roots:
        t=task_record(origin['origin']+'_global','global','primary',origin['origin'],origin['family'],seed=int(origin['seed']))
        add_task(state,data,t)
        s=get_job(data,t['task_id'],'screen')
        if s is None: submit(state,data,project,t['task_id'],'screen','screen','primary'); all_done=False; continue
        if not is_closed(s): all_done=False; continue
        if not job_ok(s): continue
        if not launch_ready(state,data,project,t): all_done=False; continue
        gr=accepted(data,state.study,t['task_id'])
        if gr:
            l=task_record(origin['origin']+'_local','local','primary',origin['origin'],origin['family'],source=fit_source(state.study,project,t['task_id']),seed=int(origin['seed']))
            add_task(state,data,l)
            if not launch_ready(state,data,project,l): all_done=False
    if not all_done: return 'prospective_primary_chains'
    # Deterministic distinct source selection, frozen once primary chains close.
    if not (state.study/'manifests/consistency_sources.tsv').exists():
        cs=candidates(data,state.study,project); selected=[]
        for pool in [cs,[r for r in cs if r['family']=='A' and r['mode']=='local'],[r for r in cs if r['family']=='B' and r['mode']=='local'],[r for r in cs if r['family'] in {'A','B'} and r['mode']=='local']]:
            choose=next((r for r in pool if r['origin'] not in {x['origin'] for x in selected}),None)
            if choose: selected.append(choose)
        for r in cs:
            if len(selected)>=4: break
            if r['origin'] not in {x['origin'] for x in selected}: selected.append(r)
        write_table(state.study/'manifests/consistency_sources.tsv',[{k:r.get(k,'') for k in ['source_fit','origin','family']} for r in selected],['source_fit','origin','family'])
    done=[]
    for i,s in enumerate(read_table(state.study/'manifests/consistency_sources.tsv')):
        t=task_record('consistency_'+str(i+1),'local','consistency',s['origin'],s['family'],source=s['source_fit'],alternative=True)
        add_task(state,data,t); done.append(launch_ready(state,data,project,t))
    if not all(done): return 'primary_consistency_checks'
    g,cs,rootresults=gate(data,state.study,project)
    if g['pass_gate']: return 'complete_pass'
    # One reserve decision only; no method tuning or threshold change.
    decision=state.study/'manifests/reserve_origins.tsv'
    if not decision.exists():
        near={r['family'] for r in cs if r['family'] in {'A','B'} and r['objective']-g['q_star']<=.001}
        missing={'A','B'}-near
        allocation={'A':2,'B':2} if len(missing)==2 else {next(iter(missing)):4} if missing else {}
        selected=[]
        for fam,count in sorted(allocation.items()): selected.extend([r for r in seeds if r['stage']=='reserve' and r['family']==fam][:count])
        write_table(decision,selected,['origin','family','stage','seed'])
    reserved=read_table(decision); closed=[]
    for origin in reserved:
        t=task_record(origin['origin']+'_global','global','reserve',origin['origin'],origin['family'],seed=int(origin['seed']))
        add_task(state,data,t); s=get_job(data,t['task_id'],'screen')
        if s is None: submit(state,data,project,t['task_id'],'screen','screen','reserve'); closed.append(False); continue
        if not is_closed(s): closed.append(False); continue
        if not job_ok(s): closed.append(True); continue
        if not launch_ready(state,data,project,t): closed.append(False); continue
        r=accepted(data,state.study,t['task_id'])
        if r:
            l=task_record(origin['origin']+'_local','local','reserve',origin['origin'],origin['family'],source=fit_source(state.study,project,t['task_id']),seed=int(origin['seed']))
            add_task(state,data,l); closed.append(launch_ready(state,data,project,l))
        else: closed.append(True)
    if not all(closed): return 'reserve_chains'
    # At most six targeted repairs/repeats, selected serially because each depends
    # on its own updated endpoint. Other ready chains were submitted together.
    for i in range(6):
        g,cs,rootresults=gate(data,state.study,project)
        if g['pass_gate']: return 'complete_pass'
        tid='reserve_check_'+str(i+1); existing=next((r for r in data['tasks'] if r['task_id']==tid),None)
        if existing: t=read_kv(state.study/'tasks'/(tid+'.tsv'))
        else:
            previous=[r for r in data['tasks'] if r['phase']=='reserve_local']
            prospective={}
            for r in cs:
                if r['family'] in {'A','B'} and r['origin'] not in prospective: prospective[r['origin']]=r
            unstable=[r for r in prospective.values() if r['objective']-g['q_star']<=.001 and (r['mode']!='local' or not r['local_stable'])]
            target=sorted(unstable,key=lambda r:(sum(p['origin']==r['origin'] for p in previous),r['objective']))[0] if unstable else cs[0]
            t=task_record(tid,'local','reserve_local',target['origin'],target['family'],source=target['source_fit'],alternative=True)
            add_task(state,data,t)
        if not launch_ready(state,data,project,t): return 'reserve_local_checks'
    g,cs,roots=gate(data,state.study,project)
    return 'complete_pass' if g['pass_gate'] else 'complete_unresolved'

def validate(project,study):
    verify_code(study); audit(study)
    runner=study/'implementation/run.sh'
    def run(mode,tid): subprocess.run(['bash',str(runner),str(project),str(study),mode,tid],check=True)
    subprocess.run([sys.executable,str(study/'implementation/test_workflow.py')],check=True)
    run('unit','unit')
    for ident in ['review_035','review_009','review_014','review_001']:
        for rep in range(1,4): run('replay',ident+'_rep'+str(rep))
    legacy_replays=0
    for ident in ['review_001','review_002','review_038','review_039','review_035','review_006']:
        tid='history_'+ident
        if (study/'fits'/tid/'COMPLETE').exists():
            run('replay',tid+'_rep_migration'); legacy_replays+=1
    for case,source in [('unrestricted','review_035'),('induction_zero','review_009'),('survival_zero','review_014')]:
        tid='fixture_'+case+'_global'; run('global',tid); run('check',tid)
    audit(study)
    write_kv(study/'tests/global_validation.tsv',{'pass':True,'code_sha256':sha(study/'manifests/code_manifest.tsv'),'fresh_process_replays':12,'legacy_endpoint_replays':legacy_replays,'actual_eight_worker_fixtures':3,'completed_utc':now()})

def validate_local(project,study):
    verify_code(study); audit(study)
    for case in ['unrestricted','induction_zero','survival_zero']:
        tid='fixture_'+case+'_local'
        for mode in ['local','check']: subprocess.run(['bash',str(study/'implementation/run.sh'),str(project),str(study),mode,tid],check=True)
    audit(study)
    write_kv(study/'tests/local_validation.tsv',{'pass':True,'code_sha256':sha(study/'manifests/code_manifest.tsv'),'actual_one_core_local_fixtures':3,'completed_utc':now()})

def controller(project,study,once=False):
    state=State(study); start=time.monotonic()
    while True:
        with state.locked():
            data=state.load()
            try:
                status=drive(state,data,project); export_progress(state,data,project,status)
                if status.startswith('complete_'):
                    report=get_job(data,'final_report','report')
                    if report is None:
                        submit(state,data,project,'final_report','report','check','report')
                        status='final_static_report'
                    elif not is_closed(report): status='final_static_report'
                    elif not job_ok(report) or not (study/'reports/COMPLETE').exists(): raise ValueError('Final report failed')
                if status.startswith('complete_'):
                    write_kv(study/'tables/final_status.tsv',dict(status=status,completed_utc=now()))
                    atomic_text(study/'CLOSURE.md',(study/'STATUS.md').read_text()+'\nThe numerical evidence is in tables/recovery_candidates.tsv, tables/recovery_gate.tsv, tables/jobs.tsv and each independently validated fit directory. Native global iteration limits, local stability and independent-origin recovery are separate quantities. These diagnostics do not provide likelihood confidence intervals or structural-identifiability proofs.\n')
                    return
                if once: return
                if time.monotonic()-start>42300:
                    submit(state,data,project,'controller_'+uuid.uuid4().hex[:8],'controller','controller','control'); return
            except BaseException as exc:
                state.event(data,'dependent_submissions_stopped',str(exc)); state.save(data)
                atomic_text(study/'BLOCKED.txt',now()+' '+str(exc)+'\n'); export_progress(state,data,project,'blocked_integrity_or_budget')
                raise
        time.sleep(30)

def main():
    p=argparse.ArgumentParser(); p.add_argument('mode',choices=['controller','validate','validate_local','start','audit']); p.add_argument('project',type=Path); p.add_argument('study',type=Path); p.add_argument('--once',action='store_true'); a=p.parse_args()
    if a.mode=='audit': audit(a.study); verify_code(a.study)
    elif a.mode=='validate': validate(a.project,a.study)
    elif a.mode=='validate_local': validate_local(a.project,a.study)
    elif a.mode=='controller': controller(a.project,a.study,a.once)
    elif a.mode=='start':
        state=State(a.study)
        with state.locked():
            data=state.load(); submit(state,data,a.project,'validation','validate','validate','validation')
            if not any(j['role']=='controller' and not is_closed(j) for j in data['jobs']):
                number=sum(j['role']=='controller' for j in data['jobs'])+1
                submit(state,data,a.project,'controller_start_'+str(number),'controller','controller','control')
if __name__=='__main__': main()
