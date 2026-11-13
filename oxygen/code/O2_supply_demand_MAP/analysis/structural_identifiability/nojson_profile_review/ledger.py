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
TASK_FIELDS = ['task_id','mode','phase','origin','family','case','point','target','value','source_fit','fallback_source','seed','alternative','min_generations','max_generations']
RESERVATION_FIELDS = ['task','role','phase','reserved_ch']
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
        if not ptr.exists(): return {'jobs':[],'tasks':[],'events':[],'reservations':[]}
        generation=ptr.read_text().strip()
        if not re.fullmatch(r'g_[0-9a-f]+',generation): raise ValueError('Invalid state pointer')
        d=self.study/'state'/generation
        for r in read_table(d/'manifest.tsv'):
            if sha(d/r['path'])!=r['sha256']: raise ValueError('State generation corrupted')
        return {n:read_table(d/(n+'.tsv')) for n in ['jobs','tasks','events','reservations']}
    def save(self,data):
        name='g_'+uuid.uuid4().hex; tmp=self.study/'state'/('.'+name); tmp.mkdir()
        for n,fields in [('jobs',JOB_FIELDS),('tasks',TASK_FIELDS),('events',['utc','event','detail']),('reservations',RESERVATION_FIELDS)]: write_table(tmp/(n+'.tsv'),data[n],fields)
        write_table(tmp/'manifest.tsv',[{'path':n+'.tsv','sha256':sha(tmp/(n+'.tsv'))} for n in ['jobs','tasks','events','reservations']],['path','sha256'])
        os.rename(tmp,self.study/'state'/name); atomic_text(self.study/'state/CURRENT',name+'\n')
    def event(self,data,event,detail): data['events'].append(dict(utc=now(),event=event,detail=detail))

def add_task(state,data,t):
    if any(r['task_id']==t['task_id'] for r in data['tasks']):
        if read_kv(state.study/'tasks'/(t['task_id']+'.tsv'))!=t: raise ValueError('Task identity collision')
        return
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
    return {'global':(8,43200,16),'local':(1,21600,8),'screen':(1,3600,8),'check':(1,1800,8),'report':(1,1800,8),'validate':(8,10800,16),'controller':(1,43200,4)}[role]
def committed(data,role=None):
    actual=sum(float(j['actual_ch']) if is_closed(j) and j['actual_ch'] else float(j['reserved_ch']) for j in data['jobs'] if role is None or (j['role'] not in {'global','local'} if role=='nonfit' else j['role']==role))
    return actual+sum(float(r['reserved_ch']) for r in data.get('reservations',[]) if role is None or role==r['role'])
def budget_check(data,role,cpus,seconds):
    request=cpus*seconds/3600
    if committed(data)+request>4000+1e-9: raise ValueError('Total allocation ceiling reached')
    # Keep one core-hour inside the approved non-fit cap for final closure.
    nonfit_limit=120 if role=='report' else 119
    if role not in {'global','local'} and committed(data,'nonfit')+request>nonfit_limit+1e-9: raise ValueError('Non-fit allocation ceiling reached')
    if role=='local' and sum(j['role']=='local' for j in data['jobs'])+len(data.get('reservations',[]))>=70: raise ValueError('Local attempt cap reached')
    return request

def submit(state,data,project,tid,mode,role,phase):
    old=get_job(data,tid,mode)
    if old: return old
    cpus,seconds,mem=resources(role); request=budget_check(data,role,cpus,seconds)
    token='cpr_'+uuid.uuid4().hex[:18]
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
    j=get_job(data,tid,'check'); directory=study/'fits'/tid
    if not job_ok(j) or not (directory/'COMPLETE').exists(): return None
    receipt=read_kv(directory/'artifact_validation.tsv');result=read_kv(directory/'result.tsv')
    if not receipt['pass'] or receipt['result_sha256']!=sha(directory/'result.tsv') or receipt['payload_manifest_sha256']!=sha(directory/'payload_manifest.tsv'): raise ValueError('Stale accepted artifact')
    if result.get('skipped',False):
        if receipt['accepted_fit_sha256'] is not None or not (directory/'TRANSFER_SKIPPED').exists(): raise ValueError('Invalid skipped transfer')
    elif receipt['accepted_fit_sha256']!=sha(directory/'accepted_fit.rds'): raise ValueError('Stale accepted fit')
    return result
