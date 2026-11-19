"""Explicit engineering or A400 dispatch with cumulative N20 accounting."""
import fcntl
import math
from pathlib import Path
import statistics
import subprocess
import sys

from common import read, table, verify, sha

CAP = 184*3600
HISTORY = 119
OLD_CALLS = 2
TERMINAL = {'COMPLETED','FAILED','TIMEOUT','CANCELLED','NODE_FAIL','OUT_OF_MEMORY',
            'PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED','SPECIAL_EXIT'}
FIELDS = ('action','job_id','state','exit_code','settled_core_seconds','reserved_core_seconds')
PHASES = {'engineering':('E',),'benchmark':('N_A_import','N_A_bench')}

def site_time():
    raw=subprocess.check_output(['sacctmgr','-nP','show','qos','xxlarge','format=MaxWall'],text=True)
    text=raw.strip().rstrip('|')
    days,sep,clock=text.partition('-')
    if not sep:clock,days=days,'0'
    h,m,s=map(int,clock.split(':'))
    return int(days)*86400+3600*h+60*m+s

def count(path):
    return sum(r['event']=='START' for r in read(path)) if path.is_file() else 0

class Review:
    def __init__(self,root):
        self.root=root
        self.jobs_path=root/'accounting/jobs.tsv'
        self.jobs=read(self.jobs_path) if self.jobs_path.exists() else []
    def save(self):
        table(self.jobs_path,self.jobs,FIELDS)
    def refresh(self):
        if not self.jobs:return
        ids=[r['job_id'] for r in self.jobs]
        if any(not s.isdigit() for s in ids):
            raise RuntimeError('Unreconciled native submission receipt')
        x=subprocess.check_output(['sacct','-nP','-j',','.join(ids),
            '--format=JobIDRaw,State%40,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,QOS'],text=True)
        (self.root/'accounting/slurm_snapshot.txt').write_text(x)
        raw={r[0]:r for line in x.splitlines() if len(r:=line.split('|'))>=7}
        for job in self.jobs:
            jid=job['job_id']
            if jid not in raw:continue
            r=raw[jid]
            state=r[1].split()[0]
            if r[6]!='xxlarge':
                raise RuntimeError('Unexpected job QoS')
            job.update(state=state,exit_code=r[2])
            if state in TERMINAL:
                cost=max(int(y[3] or 0)*int(y[4] or 0) for k,y in raw.items() if k==jid or k.startswith(jid+'.'))
                job.update(settled_core_seconds=str(cost),reserved_core_seconds='0')
        self.save()
    def budget(self):
        actual=sum(int(j['settled_core_seconds']) for j in self.jobs)
        pending=sum(int(j['reserved_core_seconds']) for j in self.jobs)
        e=count(self.root/'verification/E/forward_attempts.tsv')
        a=count(self.root/'verification/N_A_bench/evaluations.tsv')
        if actual+pending+HISTORY>CAP or e>3 or a>400:
            raise RuntimeError('Scope or cumulative cap exceeded')
        result=dict(historical_core_seconds=HISTORY,supplement_core_seconds=actual,
          cumulative_core_seconds=HISTORY+actual,live_reserved_core_seconds=pending,cumulative_cap_core_seconds=CAP,
          historical_forwards=OLD_CALLS,engineering_forwards=e,A_benchmark_forwards=a)
        table(self.root/'accounting/cumulative_budget.tsv',[result])
        return result

    def done(self,action,gate=True):
        j=next((r for r in self.jobs if r['action']==action),None)
        if not j or (j['state'],j['exit_code'])!=('COMPLETED','0:0'):return False
        if not gate:return True
        for phase in PHASES[action]:
            p=self.root/'verification'/phase
            manifest=p/'phase_manifest.tsv'
            if not manifest.exists():return False
            for m in read(manifest):
                path=p/m['path']
                if path.stat().st_size!=int(m['bytes']) or sha(path)!=m['sha256']:
                    raise RuntimeError('Phase receipt changed: '+m['path'])
            app,=read(p/'application_result.tsv')
            policy,=read(p/'policy_gate.tsv')
            codes,=read(p/'phase_exit_codes.tsv')
            if app['status']!='PASS' or app['job']!=j['job_id'] or policy['status']!='PASS' or any(int(x) for x in codes.values()):
                raise RuntimeError('Application/source/policy admission failure: '+phase)
            if (p/'read_denials.tsv').exists() and (p/'read_denials.tsv').read_text().strip():
                raise RuntimeError('Unexpected read denial')
        return True

    def submit(self,action):
        self.refresh()
        verify(self.root)
        if action not in PHASES or any(j['action']==action for j in self.jobs):
            raise RuntimeError('One allocation per approved action')
        if action=='benchmark' and not self.done('engineering'):
            raise RuntimeError('Engineering verification has not passed')
        b=self.budget()
        left=CAP-b['cumulative_core_seconds']-b['live_reserved_core_seconds']
        envelope=(4*3600-HISTORY) if action=='engineering' else 28*3600
        cap=min(site_time(),left,envelope)//60*60
        if cap<60:raise RuntimeError('No remaining allocation budget')
        j=dict(action=action,job_id='submission_pending',state='RESERVED',exit_code='',
               settled_core_seconds='0',reserved_core_seconds=str(cap))
        self.jobs.append(j)
        self.save()
        args=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1',
             '--cpus-per-task=1','--hint=nomultithread','--mem=16G','--no-requeue','--time='+str(cap//60),
             '--job-name=n20score_'+action,'--output='+str(self.root/'logs'/(action+'_%j.txt')),
             str(self.root/'release/launch_scope.sh'),str(self.root),action]
        jid=subprocess.check_output(args,text=True).strip().split(';')[0]
        if not jid.isdigit():raise RuntimeError('Uncertain sbatch ID; preserve reservation')
        j.update(job_id=jid,state='SUBMITTED')
        self.save()
        self.budget()
        print('SUBMITTED',action,jid)

    def cost(self):
        if not self.done('benchmark'):raise RuntimeError('Benchmark has no passing terminal receipt')
        p=self.root/'verification/N_A_bench'
        t=read(p/'callback_timings.tsv')
        f=sorted(float(r['seconds']) for r in t if r['feasible']=='TRUE')
        values=[float(r['seconds']) for r in t]
        starts=count(p/'evaluations.tsv')
        if len(t)!=400 or starts!=400 or len(f)<200 or any(not math.isfinite(x) or x<0 for x in values):
            raise RuntimeError('Not exactly 400 complete benchmark callbacks with verified initial feasibility')
        j=next(r for r in self.jobs if r['action']=='benchmark')
        overhead=max(0,int(j['settled_core_seconds'])-sum(values))
        q95=f[math.ceil(.95*len(f))-1]
        c=max(1.0,1.25*statistics.mean(values),q95)
        forecast=math.ceil(84000*c+3600+18*max(0,overhead-200))
        status='PASS' if forecast<=28*3600 and 12200*c+overhead<site_time() else 'HOLD'
        records=[s.split('|') for s in (self.root/'accounting/slurm_snapshot.txt').read_text().splitlines()]
        rss=[r[5] for r in records if len(r)>=7 and r[0].split('.')[0]==j['job_id'] and r[5]]
        if not rss:raise RuntimeError('Benchmark memory accounting unavailable')
        table(self.root/'accounting/benchmark_N.tsv',[dict(status=status,job_id=j['job_id'],
            callbacks=400,feasible_callbacks=len(f),observed_mean_seconds=statistics.mean(values),feasible_p95_seconds=q95,
            admission_seconds_per_callback=c,allocated_core_seconds=j['settled_core_seconds'],overhead_seconds=overhead,
            forecast_core_seconds=forecast,forecast_core_hours=forecast/3600,peak_memory=';'.join(rss),
            continuation_released=False,other_blind_chains_released=False)])
        print('BENCHMARK',status,forecast)

def main(root,verb,*args):
    root=Path(root).resolve()
    with (root/'accounting/owner.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        r=Review(root)
        if verb=='submit':r.submit(args[0])
        elif verb in ('status','cost'):
            r.refresh()
            print(r.budget())
            for j in r.jobs:
                print(*(j[k] for k in FIELDS),r.done(j['action']),sep='\t')
            if verb=='cost':r.cost()
        else:raise RuntimeError('Unregistered action')

if __name__=='__main__':
    main(*sys.argv[1:])
