"""Single metadata owner for the new32-hour study. No model/response operations."""
from pathlib import Path
import fcntl,math,os,statistics,subprocess,sys,time
from common import *
TERMINAL={'COMPLETED','FAILED','TIMEOUT','CANCELLED','NODE_FAIL','OUT_OF_MEMORY','PREEMPTED','BOOT_FAIL','DEADLINE'}
FIELDS=['job_id','task','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds']
CHAIN=['bench','de1b','de2','de3','de4','de5','de6','fill']
def truth(v):return str(v).lower()=='true'
def vtask(t):return t.endswith(('_check1','_check2','_robust'))
def estimate(rows,feasible):
 if len(feasible)<8:raise ValueError('Eight complete feasible engineering cost samples required')
 actual=[float(r['seconds']) for r in rows];samples=[float(r['seconds']) for r in feasible]
 c=max(1.,1.25*statistics.mean(actual or samples),sorted(samples)[math.ceil(.95*len(samples))-1]);return c,12000*c+1800
class Owner:
 def __init__(self,s):
  self.s=Path(s).resolve();self.jobs=read(self.s/'accounting/jobs.tsv');self.site=int(read(self.s/'accounting/site_limits.tsv')[0]['max_seconds']);self.cache={}
 def save(self):table(self.s/'accounting/jobs.tsv',self.jobs,FIELDS)
 def retry_map(self):
  p=self.s/'accounting/retry_authorization.tsv';return read(p) if p.exists() else []
 def mapped(self,t):return next((r['retry_task'] for r in self.retry_map() if r['original_task']==t),t)
 def job(self,t):return next((j for j in self.jobs if j['task']==self.mapped(t)),None)
 def halt(self,status,detail):
  table(self.s/'accounting/stop_record.tsv',[dict(status=status,detail=detail)]);(self.s/'controller/STOP').write_text(status+'\n'+detail+'\n');print(status,detail,flush=True)
 def refresh(self):
  if not self.jobs:return
  assert all(j['job_id'].isdigit() for j in self.jobs),'Uncertain dispatch must be reconciled'
  raw=subprocess.check_output(['sacct','-n','-P','-j',','.join(j['job_id'] for j in self.jobs),'--format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,QOS'],text=True)
  (self.s/'accounting/slurm_snapshot.txt').write_text(raw);records={v[0]:v for line in raw.splitlines() if len(v:=line.split('|'))>=7}
  for j in self.jobs:
   if j['job_id'] not in records:continue
   r=records[j['job_id']];j.update(state=r[1],exit_code=r[2],allocated_cpus=r[3],elapsed_seconds=r[4])
   if r[1].split()[0] in TERMINAL:
    parts=[v for k,v in records.items() if k==j['job_id'] or k.startswith(j['job_id']+'.')];j.update(actual_core_seconds=max(int(v[3] or 0)*int(v[4] or 0) for v in parts),reserved_core_seconds=0)
  self.save();a=sum(int(j['actual_core_seconds']) for j in self.jobs);live=sum(int(j['reserved_core_seconds']) for j in self.jobs)
  table(self.s/'accounting/budget_ledger.tsv',[dict(cap_core_seconds=CAP,settled_core_seconds=a,live_reserved_core_seconds=live,technical_retry_hold=self.retry_hold(),collection_hold=SETTLEMENT,validation_hold=self.validation_hold(),uncommitted=CAP-a-live)])
  assert a+live<=CAP,'Allocated budget ceiling exceeded'
 def retry_hold(self):return RETRY if kv(self.s/'accounting/retry_ledger.tsv')['used']=='0'else 0
 def validation_hold(self):
  charged=sum(int(j['actual_core_seconds'])+int(j['reserved_core_seconds']) for j in self.jobs if vtask(j['task']));return max(0,14400-charged)
 def ready(self,t):
  j=self.job(t)
  if not j or j['state']!='COMPLETED':return False
  assert j['exit_code']=='0:0','Terminal exit not0:0 '+t
  task=j['task'];p=self.s/'verification'/task;manifest=p/'phase_manifest.tsv'
  if not manifest.exists():return False
  key=sha(manifest)
  if self.cache.get(task)==key:return True
  for r in read(manifest):
   q=p/r['path'];assert q.stat().st_size==int(r['bytes']) and sha(q)==r['sha256'],'Bad phase receipt '+str(q)
  for name in ('application_result.tsv','gate.tsv','policy_gate.tsv'):assert read(p/name)[0]['status']=='PASS','Gate failed '+t+'/'+name
  assert all(int(v)==0 for v in read(p/'phase_exit_codes.tsv')[0].values()),'Independent phase exits failed '+t
  if task!='P1C' and (p/'read_denials.tsv').exists():assert not (p/'read_denials.tsv').read_text().strip(),'Unexpected read denial '+task
  self.cache[task]=key;return True
 def submit(self,t):
  t=self.mapped(t)
  if self.job(t):return
  verify(self.s);actual=sum(int(j['actual_core_seconds']) for j in self.jobs);live=sum(int(j['reserved_core_seconds']) for j in self.jobs)
  if t.endswith('_R1'):purpose=RETRY-sum(int(j['actual_core_seconds']) for j in self.jobs if j['task'].endswith('_R1'))
  elif t in ('P0','P1','P1C'):purpose=14400-sum(int(j['actual_core_seconds']) for j in self.jobs if j['task'] in ('P0','P1','P1C'))
  elif t.startswith('P2_'):purpose=1800
  elif vtask(t):
   origin=t[:3];purpose=3600-sum(int(j['actual_core_seconds']) for j in self.jobs if j['task'].startswith(origin+'_') and vtask(j['task']))
  else:
   origin=t[:3];p=self.s/'accounting/cost_admission.tsv';share=float(read(p)[0]['forecast_origin_seconds']) if p.exists() else 14400
   purpose=share-sum(int(j['actual_core_seconds']) for j in self.jobs if j['task'].startswith(origin+'_') and not vtask(j['task']))
  # Validation allocations transfer, rather than double-charge, their protected pool.
  protect=SETTLEMENT+self.retry_hold()+(0 if vtask(t) else self.validation_hold())
  available=CAP-actual-live-protect;seconds=int(min(self.site,purpose,available)//60*60)
  if seconds<=0:raise RuntimeError('Budget admission exhausted '+t)
  row=dict(job_id='submission_pending',task=t,state='RESERVED',exit_code='',allocated_cpus='',elapsed_seconds='',actual_core_seconds=0,reserved_core_seconds=seconds);self.jobs.append(row);self.save()
  args=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1','--cpus-per-task=1','--hint=nomultithread','--mem=16G','--no-requeue','--time='+str(seconds//60),'--job-name=fdiag_'+t,'--output='+str(self.s/'logs'/(t+'_%j.txt')),str(self.s/'release/runner.sh'),str(self.s),t]
  jid=subprocess.check_output(args,text=True).strip().split(';')[0];assert jid.isdigit(),'Ambiguous submission';row.update(job_id=jid,state='SUBMITTED');self.save();print('SUBMITTED',t,jid,'reserved_core_seconds',seconds,flush=True)
 def admission(self):
  feasible=[r for r in read(self.s/'verification'/self.mapped('P1')/'complete_cost_samples.tsv') if truth(r['feasible'])];rows=[]
  for j in self.jobs:
   if j['task'].startswith('F_') and not vtask(j['task']) and self.ready(j['task']):
    data=read(self.s/'verification'/j['task']/'callback_timings.tsv');rows+=data;feasible += [r for r in data if truth(r['feasible'])]
  c,forecast=estimate(rows,feasible);actual=sum(int(j['actual_core_seconds']) for j in self.jobs)
  # Full-chain conservative forecast, not just the next block; completed families
  # need no further construction reservation. Live shares already bound allocations.
  remaining=sum(not self.bank_path('F_'+letter) for letter in 'ABCD')
  fits=actual+remaining*forecast+14400+SETTLEMENT+self.retry_hold()<=CAP
  table(self.s/'accounting/cost_admission.tsv',[dict(status='PASS' if fits else 'RESOURCE_ADMISSION_HOLD',c_admission_seconds=c,forecast_origin_seconds=forecast,remaining_families=remaining,benchmark_callbacks=len(rows),complete_feasible_cost_samples=len(feasible),full_origin_calls=12000)])
  if not fits:self.halt('RESOURCE_ADMISSION_HOLD','Registered complete constructor forecast does not fit independent32-hour package')
  return fits
 def bank_path(self,origin):
  for stage in CHAIN:
   t=origin+'_'+stage
   if self.ready(t):
    p=self.s/'verification'/self.mapped(t)/'result.tsv'
    if kv(p)['status']=='BANK_READY':return self.mapped(t)
  return None
 def construction(self,origin):
  if self.bank_path(origin):return True
  for stage in CHAIN:
   t=origin+'_'+stage
   if not self.job(t):self.submit(t);return False
   if not self.ready(t):return False
   status=kv(self.s/'verification'/self.mapped(t)/'result.tsv')['status']
   if status=='BANK_READY':return True
   if status=='FEASIBILITY_CONSTRUCTION_UNRESOLVED':return True
  return True
 def validate(self,origin):
  path=self.bank_path(origin)
  if not path:return True
  record=self.s/'controller'/(origin+'_bank_path.txt')
  if not record.exists():record.write_text(path+'\n')
  assert record.read_text().strip()==path,'Bank pointer changed after selection'
  for stage in ('check1','check2','robust'):
   t=origin+'_'+stage
   if not self.job(t):self.submit(t);return False
   if not self.ready(t):return False
  return True
 def diagnosis_summary(self):
  p=self.s/'reports/diagnosis_gate.tsv'
  if p.exists():return
  failures=[];ledgers=[];refs=[]
  for i in range(1,5):
   base=self.s/'verification'/self.mapped('P2_'+str(i));failures+=read(base/'candidate_failure_ledger.tsv');ledgers+=read(base/'passage_supply_ledger.tsv');refs+=read(base/'reference_replay_checks.tsv')
  assert len(failures)==2000 and len(refs)==64 and len({r['candidate'] for r in failures})==1000
  assert all(r['status']=='PASS' for r in refs)
  table(self.s/'reports/candidate_failure_ledger.tsv',failures);table(self.s/'reports/passage_supply_ledger.tsv',ledgers);table(self.s/'reports/reference_replay_checks.tsv',refs)
  from collections import Counter
  counts=Counter((r['cohort'],r['parent'],r['child'],r['feasible']) for r in failures)
  table(self.s/'reports/cohort_failure_summary.tsv',[dict(cohort=k[0],parent=k[1],child=k[2],feasible=k[3],candidates=n) for k,n in sorted(counts.items())])
  table(p,[dict(status='PASS',historical_candidates=1000,cohort_records=2000,reference_checks=64,interpretation='First-failure records; later edges remain censored NOT_REACHED')])
 def tick(self):
  self.refresh()
  retries=self.retry_map()
  for j in self.jobs:
   if j['state'].split()[0] in TERMINAL and (j['state']!='COMPLETED' or j['exit_code']!='0:0'):
    if any(r['original_task']==j['task'] and r['original_job']==j['job_id'] for r in retries):continue
    self.halt('TECHNICAL_REVIEW_REQUIRED',j['task']+' '+j['job_id']+' '+j['state']+' '+j['exit_code']+'; no implicit retry');return
  for t in ('P0','P1','P1C'):
   if not self.job(t):self.submit(t);return
   if not self.ready(t):return
  for i in range(1,5):
   if not self.job('P2_'+str(i)):self.submit('P2_'+str(i))
  if not all(self.ready('P2_'+str(i)) for i in range(1,5)):return
  self.diagnosis_summary()
  if not self.job('F_A_bench'):self.submit('F_A_bench');return
  if not self.ready('F_A_bench'):return
  if not self.admission():return
  finished=[self.construction('F_'+letter) for letter in 'ABCD']
  checked=[self.validate('F_'+letter) for letter in 'ABCD']
  if all(finished) and all(checked):
   paths=[self.bank_path('F_'+letter) for letter in 'ABCD'];status='BANKS_READY_FOR_D19_REVIEW' if all(paths)else'FEASIBILITY_CONSTRUCTION_UNRESOLVED'
   table(self.s/'reports/bank_release_gate.tsv',[dict(status=status,validated_banks=sum(bool(p) for p in paths),D19_fit=0,N20=0,profiles=0,tumor=0)])
   (self.s/'controller/FINISHED').write_text(status+'\n');print(status,flush=True)
def main():
 s=Path(sys.argv[1]).resolve();lock=(s/'controller/owner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 table(s/'controller/ownership.tsv',[dict(pid=os.getpid(),study=str(s),scope='metadata-only owner; no old tasks/controllers')]);owner=Owner(s)
 while not any((s/'controller'/f).exists() for f in ('STOP','FINISHED')):
  try:owner.tick()
  except Exception as e:owner.halt('CONTROLLER_RECONCILIATION_HOLD',str(e));break
  time.sleep(20)
if __name__=='__main__':main()
