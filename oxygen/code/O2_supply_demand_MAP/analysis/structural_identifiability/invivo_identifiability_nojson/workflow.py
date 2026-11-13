#!/usr/bin/env python3
"""Owned-job submission, coherent TSV state and strict pilot allocation limits."""
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
FIELDS=['task','role','mode','release','job_id','comment','state','exit_code','cores','wall_seconds','reserved_ch','actual_ch','intent_utc','receipt']
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
def charge(j): return float(j['actual_ch']) if closed(j) and j['actual_ch'] else float(j['reserved_ch'])
def resolve_matches(rows,candidates):
    for j in rows:
        if j['job_id']: continue
        matches=sorted(set(job_id for job_id,comment in candidates if comment==j['comment'] and job_id.isdigit()))
        if len(matches)!=1: raise ValueError('Uncertain submission has '+str(len(matches))+' matches; no retry is allowed')
        j.update(job_id=matches[0],state='RECONCILED')
    return rows
def check_budget(rows,budget,role):
    b=next(x for x in budget if x['role']==role); request=int(b['cores'])*int(b['wall_seconds'])/3600
    same=[x for x in rows if x['role']==role]
    if len(same)>=int(b['attempts']): raise ValueError('Role attempt cap reached')
    if sum(map(charge,same))+request>float(b['cap_core_hours'])+1e-9: raise ValueError('Role allocation cap reached')
    if sum(map(charge,rows))+request>192+1e-9: raise ValueError('Total allocation cap reached')
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
            actual=int(r[3])*int(r[4])/3600
            if actual>float(j['reserved_ch'])+1e-6: raise ValueError('Actual allocation exceeded reservation')
            j['actual_ch']=format(actual,'.17g')
def gate(study,mode,task):
    if mode in {'replay','replay_suite'}: return
    required=['replay_gate.tsv']
    if mode in {'global','local'}: required.append('engineering_gate.tsv')
    for name in required:
        rows=read_table(study/'tables'/name)
        if next(r['value'] for r in rows if r['key']=='pass')!='true': raise ValueError('Required gate has not passed: '+name)
    if mode=='local': verify_payload(study/'fits'/task.replace('local_','global_'))
def submit(study,role,mode,task,release):
    if not re.fullmatch('[A-Za-z0-9_]+',task) or not re.fullmatch('[A-Za-z0-9_]+',release): raise ValueError('Invalid task/release identity')
    s=State(study)
    with s.locked():
        rows=s.load(); reconcile(study,rows)
        existing=next((j for j in rows if j['task']==task),None)
        if existing:
            if not existing['job_id']: raise ValueError('Task submission is uncertain; reconcile rather than resubmit')
            s.save(rows); print('Existing task; no duplicate submission:',task,existing['job_id']); return
        if any(not j['job_id'] for j in rows): raise ValueError('Uncertain submission must be reconciled before any further submission')
        verify_manifest(study,study/'manifests/input_manifest.tsv')
        verify_manifest(study/'releases'/release,study/'releases'/release/'manifest.tsv')
        audit_formats(study); gate(study,mode,task)
        budget,request=check_budget(rows,read_table(study/'manifests/budget.tsv'),role)
        comment='ivi:'+study.name+':'+task
        j=dict(task=task,role=role,mode=mode,release=release,job_id='',comment=comment,state='INTENT',exit_code='',cores=budget['cores'],wall_seconds=budget['wall_seconds'],reserved_ch=format(request,'.17g'),actual_ch='',intent_utc=utc(),receipt='')
        rows.append(j); s.save(rows)
        hours=int(budget['wall_seconds'])//3600
        cmd=[str(SLURM/'sbatch'),'--parsable','--partition=red','--qos=xxlarge','--no-requeue','--nodes=1','--ntasks=1','--cpus-per-task='+budget['cores'],'--mem='+('16G' if int(budget['cores'])>1 else '8G'),'--time='+str(hours)+':00:00','--job-name=ivi_202133Z_'+task,'--comment='+comment,'--chdir='+str(study),'--output='+str(study/'logs'/(task+'_%j.log')),str(study/'releases'/release/'batch.sh'),str(study),mode,task,release]
        atomic_text(study/'logs'/(task+'_submission.txt'),'\n'.join(cmd)+'\n')
        result=subprocess.run(cmd,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
        atomic_text(study/'logs'/(task+'_receipt.txt'),result.stdout+'\n'+result.stderr)
        receipt=result.stdout.strip().split(';')[0]
        if result.returncode or not receipt.isdigit():
            s.save(rows); raise RuntimeError('Submission uncertain; reservation retained and no automatic retry')
        j.update(job_id=receipt,state='SUBMITTED',receipt=str(study/'logs'/(task+'_receipt.txt'))); s.save(rows)
        print(task,receipt,request,flush=True)
def main():
    p=argparse.ArgumentParser(); p.add_argument('command',choices=['submit','status','audit','resolve']); p.add_argument('study',type=Path); p.add_argument('--role'); p.add_argument('--mode'); p.add_argument('--task'); p.add_argument('--release',default='replay01'); a=p.parse_args()
    study=a.study.resolve()
    if a.command=='submit': submit(study,a.role,a.mode,a.task,a.release)
    elif a.command=='audit': audit_formats(study); print('Format audit passed')
    elif a.command=='resolve':
        s=State(study)
        with s.locked():
            rows=s.load(); candidates=[]
            result=subprocess.run([str(SLURM/'squeue'),'-u','4482173','-h','-o','%i|%k'],check=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            candidates.extend(tuple(line.split('|',1)) for line in result.stdout.splitlines() if '|' in line)
            for j in rows:
                if j['job_id']: continue
                result=subprocess.run([str(SLURM/'sacct'),'-u','4482173','-S',j['intent_utc'][:10],'-n','-P','--name=ivi_202133Z_'+j['task'],'--format=JobIDRaw,Comment%200'],check=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                candidates.extend(tuple(line.rstrip('|').split('|',1)) for line in result.stdout.splitlines() if '|' in line.rstrip('|'))
            resolve_matches(rows,candidates); s.save(rows)
    else:
        s=State(study)
        with s.locked():
            rows=s.load(); reconcile(study,rows); s.save(rows)
        table(study/'tables/jobs.tsv',rows,FIELDS)
        print('Used plus committed core-hours:',sum(map(charge,rows)))
        for j in rows: print(j['task'],j['job_id'],j['state'],j['exit_code'],j['actual_ch'])
if __name__=='__main__': main()
