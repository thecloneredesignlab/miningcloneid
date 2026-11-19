"""Single finite owner. Metadata only on login host; all scientific work uses xxlarge."""
from pathlib import Path
import fcntl
import math
import statistics
import subprocess
import sys
import time
from common import *
TERMINAL={'COMPLETED','FAILED','TIMEOUT','CANCELLED','NODE_FAIL','OUT_OF_MEMORY','PREEMPTED','BOOT_FAIL','DEADLINE'}
CHAIN=['screen','bench','de1b','de2','de3','de4','de5','de6','local','pattern','verify1','verify2']
FIELDS=['job_id','task','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds']
def true(v):return str(v).lower()=='true'
def cost_admission(rows):
    values=[float(r['seconds']) for r in rows if true(r['feasible'])]
    if len(values)<8:raise ValueError('Fewer than eight representative feasible callbacks')
    # The mean includes every active callback, not only feasible outcomes.
    active=[float(r['seconds']) for r in rows]
    p95=sorted(values)[math.ceil(.95*len(values))-1]
    cost=max(1.,1.25*statistics.mean(active),p95)
    return cost,84000*cost+3600

def result_quality(rows,kind):
    if kind=='D':return sum(true(r['prediction_pass']) for r in rows)>=2
    if not rows:return False
    best=min(rows,key=lambda r:float(r['score']));q=float(best['score'])
    return (sum(float(r['score'])<=q+1e-4 for r in rows)>=2 and true(best['local_pass']) and true(best['complete_poll'])
            and all(float(r['best_found'])>=q-1e-12 for r in rows)
            and best['key']==best['best_found_key'])

class Controller:
 def __init__(self,study):
    self.s=Path(study).resolve();self.jobs=read(self.s/'accounting/jobs.tsv');self.site=int(read(self.s/'accounting/site_limits.tsv')[0]['max_seconds'])
    self.used=int(kv(self.s/'accounting/retry_ledger.tsv')['used'])
 def save(self):table(self.s/'accounting/jobs.tsv',self.jobs,FIELDS)
 def stop(self,status,detail):
    table(self.s/'accounting/stop_record.tsv',[dict(status=status,detail=detail)])
    (self.s/'controller/STOP').write_text(status+'\n'+detail+'\n')
    print(status,detail,flush=True)
 def refresh(self):
    if not self.jobs:return
    ids=','.join(j['job_id'] for j in self.jobs);assert all(j['job_id'].isdigit() for j in self.jobs),'Uncertain submission requires reconciliation'
    raw=subprocess.check_output(['sacct','-n','-P','-j',ids,'--format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,QOS'],text=True)
    (self.s/'accounting/slurm_snapshot.txt').write_text(raw)
    records={v[0]:v for line in raw.splitlines() if len(v:=line.split('|'))>=7}
    for j in self.jobs:
      if j['job_id'] not in records:continue
      r=records[j['job_id']];j.update(state=r[1],exit_code=r[2],allocated_cpus=r[3],elapsed_seconds=r[4])
      if r[1].split()[0] in TERMINAL:
        parts=[v for k,v in records.items() if k==j['job_id'] or k.startswith(j['job_id']+'.')]
        j.update(actual_core_seconds=max(int(v[3] or 0)*int(v[4] or 0) for v in parts),reserved_core_seconds=0)
    self.save()
    actual=sum(int(j['actual_core_seconds']) for j in self.jobs);live=sum(int(j['reserved_core_seconds']) for j in self.jobs)
    table(self.s/'accounting/budget_ledger.tsv',[dict(cap_core_seconds=CAP,settled_core_seconds=actual,live_reserved_core_seconds=live,retry_hold=0 if self.used else RETRY,collection_hold=SETTLEMENT,remaining_unreserved=CAP-actual-live)])
    if actual+live>CAP:raise RuntimeError('Accounting ceiling violation')
 def job(self,t):return next((j for j in self.jobs if j['task']==t),None)
 def ready(self,t):
    j=self.job(t)
    if not j or j['state']!='COMPLETED':return False
    if j['exit_code']!='0:0':raise RuntimeError('Nonzero terminal exit: '+t)
    p=self.s/'verification'/t
    if not (p/'phase_manifest.tsv').exists():return False
    for r in read(p/'phase_manifest.tsv'):
      q=p/r['path']
      if q.stat().st_size!=int(r['bytes']) or sha(q)!=r['sha256']:raise RuntimeError('Phase receipt mismatch '+str(q))
    for name in ('gate.tsv','policy_gate.tsv','application_result.tsv'):
      if read(p/name)[0]['status']!='PASS':raise RuntimeError('Independent gate failed '+t+'/'+name)
    exits=read(p/'phase_exit_codes.tsv')[0]
    if any(int(v)!=0 for v in exits.values()):raise RuntimeError('Phase exits failed '+t)
    if t.startswith(('D_','N_')) and (p/'read_denials.tsv').exists() and (p/'read_denials.tsv').stat().st_size:raise RuntimeError('Unexpected blinded-read denial '+t)
    return True
 def submit(self,t):
    if self.job(t):return
    verify(self.s)
    actual=sum(int(j['actual_core_seconds']) for j in self.jobs);live=sum(int(j['reserved_core_seconds']) for j in self.jobs)
    if t.startswith(('D_','N_')):
      origin=t.rsplit('_',1)[0];spent=sum(int(j['actual_core_seconds']) for j in self.jobs if j['task'].startswith(origin+'_'))
      purpose=100800-spent
    elif t=='I02_R1':
      purpose=RETRY
    else:
      purpose=28800-sum(int(j['actual_core_seconds']) for j in self.jobs if j['task'] in ('I01','I02','I03'))
    seconds=min(self.site,purpose,CAP-actual-live-SETTLEMENT-(0 if self.used else RETRY))//60*60
    if seconds<=0:raise RuntimeError('Budget admission stopped '+t)
    row=dict(job_id='submission_pending',task=t,state='RESERVED',exit_code='',allocated_cpus='',elapsed_seconds='',actual_core_seconds=0,reserved_core_seconds=seconds)
    self.jobs.append(row);self.save()
    cmd=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1','--cpus-per-task=1','--hint=nomultithread','--mem=16G','--no-requeue','--time='+str(seconds//60),'--job-name=cr19_'+t,'--output='+str(self.s/'logs'/('%s_%%j.txt'%t)),str(self.s/'release/runner.sh'),str(self.s),t]
    jid=subprocess.check_output(cmd,text=True).strip().split(';')[0]
    if not jid.isdigit():raise RuntimeError('Ambiguous submission '+repr(jid))
    row.update(job_id=jid,state='SUBMITTED');self.save();print('SUBMITTED',t,jid,'reserved_seconds',seconds,flush=True)
 def benchmark(self,kind):
    path=self.s/'accounting'/('cost_admission_'+kind+'.tsv')
    if path.exists():return read(path)[0]['status']=='PASS'
    rows=[]
    for stage in ('screen','bench'):rows+=read(self.s/'verification'/(kind+'_A_'+stage)/'callback_timings.tsv')
    cost,forecast=cost_admission(rows);status='PASS' if forecast<=100800 else 'RESOURCE_ADMISSION_HOLD'
    actual=sum(int(j['actual_core_seconds']) for j in self.jobs)
    # Admit primaries only if all remaining registered primary allocations fit;
    # later reserve admission is separate and cannot spend reserved retry/closeout.
    future_primary=(4 if kind=='N' else 8)*forecast
    if actual+future_primary+SETTLEMENT+(0 if self.used else RETRY)>CAP:status='RESOURCE_ADMISSION_HOLD'
    table(path,[dict(status=status,c_admission_seconds=cost,forecast_origin_seconds=forecast,origin_planning_seconds=100800,active_callbacks=len(rows),reason='measured full callbacks including journal and improvements')])
    if status!='PASS':self.stop(status,kind+' benchmark does not admit the registered complete search; controls were not weakened')
    return status=='PASS'
 def quality(self,origin):
    for stage in ('verify2','verify1'):
      t=origin+'_'+stage
      if self.ready(t):return read(self.s/'verification'/t/'quality.tsv')[0]
    return None
 def advance_origin(self,origin,benchmark=False):
    for stage in CHAIN:
      t=origin+'_'+stage
      if stage=='verify2' and self.ready(origin+'_verify1'):
        q=read(self.s/'verification'/(origin+'_verify1')/'quality.tsv')[0]
        if true(q['local_pass']) and q['key']==q['best_found_key']:return True
      if not self.job(t):self.submit(t);return False
      if not self.ready(t):return False
      if stage=='screen':
        # Scientific screen failure is finite nonattainment, not a technical retry.
        r=kv(self.s/'verification'/t/'result.tsv')
        if r['status']=='SCIENTIFIC_SCREEN_INSUFFICIENT':return True
      if stage=='bench' and benchmark and not self.benchmark(origin[0]):return False
    return True
 def case(self,kind):
    a_done=self.advance_origin(kind+'_A',benchmark=True)
    if (self.s/'controller/STOP').exists():return False
    if not (self.s/'accounting'/('cost_admission_'+kind+'.tsv')).exists():
      if a_done:self.stop('SCIENTIFIC_SCREEN_INSUFFICIENT',kind+' first origin did not permit the benchmark')
      return False
    done=[a_done]+[self.advance_origin(kind+'_'+letter) for letter in 'BCD']
    if not all(done):return False
    prim=[self.quality(kind+'_'+letter) for letter in 'ABCD'];rows=[r for r in prim if r]
    if result_quality(rows,kind):self.lock_case(kind,rows);return True
    ledger=read(self.s/'accounting/reserve_ledger.tsv');reg=next(r for r in ledger if r['case']==kind)
    if reg['rounds_used']=='0':
      actual=sum(int(j['actual_core_seconds']) for j in self.jobs);live=sum(int(j['reserved_core_seconds']) for j in self.jobs)
      if actual+live+201600+SETTLEMENT+(0 if self.used else RETRY)>CAP:self.stop('RESERVE_RESOURCE_HOLD',kind+' two complete reserve origins cannot be admitted');return False
      reg.update(rounds_used=1,origin_E='released',origin_F='released');table(self.s/'accounting/reserve_ledger.tsv',ledger)
    ends=[self.advance_origin(kind+'_'+letter) for letter in 'EF']
    if not all(ends):return False
    rows=[self.quality(kind+'_'+letter) for letter in 'ABCDEF'];rows=[r for r in rows if r]
    if result_quality(rows,kind):self.lock_case(kind,rows);return True
    self.stop('DET_OPTIMIZER_UNRESOLVED' if kind=='D' else 'NOISY_ESTIMATE_PROVISIONAL_OPTIMIZER_UNRESOLVED','Finite origins and one reserve round exhausted without the registered gate')
    return False
 def lock_case(self,kind,rows):
    p=self.s/'reports'/('selection_lock_'+kind+'.tsv')
    if not p.exists():table(p,sorted(rows,key=lambda r:float(r['score'])))
 def tick(self):
    self.refresh()
    for j in self.jobs:
      if j['state'].split()[0] in TERMINAL and (j['state']!='COMPLETED' or j['exit_code']!='0:0'):
        if j['task']=='I02' and self.used==1 and (self.s/'accounting/technical_retry.tsv').exists():continue
        self.stop('TECHNICAL_FAILURE_REVIEW_REQUIRED',j['task']+' job '+j['job_id']+' '+j['state']+' '+j['exit_code']+'; no implicit retry');return
    for t in ('I01','I02_R1','I03'):
      if not self.job(t):self.submit(t);return
      if not self.ready(t):return
    if not self.case('D'):return
    if not self.case('N'):return
    (self.s/'controller/FITTING_FINISHED').write_text('D19 and N20 passed their registered finite gates. Selection locked before truth disclosure.\n')
    print('FITTING_FINISHED',flush=True)

def main():
 s=Path(sys.argv[1]).resolve();lock=(s/'controller/owner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 (s/'controller/ownership.tsv').write_text('pid\tstudy\n'+str(__import__('os').getpid())+'\t'+str(s)+'\n')
 c=Controller(s)
 while not any((s/'controller'/p).exists() for p in ('STOP','FITTING_FINISHED')):
  try:c.tick()
  except Exception as e:c.stop('CONTROLLER_RECONCILIATION_HOLD',str(e));break
  time.sleep(20) # metadata polling interval, never a research timeout
if __name__=='__main__':main()
