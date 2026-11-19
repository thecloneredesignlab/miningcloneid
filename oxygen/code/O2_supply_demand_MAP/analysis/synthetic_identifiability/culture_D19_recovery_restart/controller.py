"""Single resident metadata controller for the authorized200-core-hour D19 study."""
from pathlib import Path
import fcntl,math,os,statistics,subprocess,sys,time
from common import *
TERMINAL={'COMPLETED','FAILED','TIMEOUT','CANCELLED','NODE_FAIL','OUT_OF_MEMORY','PREEMPTED','BOOT_FAIL','DEADLINE'}
FIELDS=['job_id','task','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds']
DCHAIN=['import','bench','de1b','de2','de3','de4','de5','de6','local','pattern','verify1','replay1a','replay1b','verify2','replay2a','replay2b']
FCHAIN=['bench','de1b','de2','de3','de4','de5','de6','fill']
def truth(v):return str(v).lower()=='true'
def cost(rows):
 feasible=[float(r['seconds']) for r in rows if truth(r['feasible'])]
 assert len(feasible)>=8,'Eight complete feasible callback costs required'
 return max(1.,1.25*statistics.mean(float(r['seconds']) for r in rows),sorted(feasible)[math.ceil(.95*len(feasible))-1])
class Owner:
 def __init__(self,s):
  self.s=Path(s).resolve();self.jobs=read(self.s/'accounting/jobs.tsv');self.site=int(read(self.s/'accounting/site_limits.tsv')[0]['max_seconds']);self.cache={}
 def save(self):table(self.s/'accounting/jobs.tsv',self.jobs,FIELDS)
 def halt(self,status,detail):
  table(self.s/'accounting/stop_record.tsv',[dict(status=status,detail=detail)]);(self.s/'controller/STOP').write_text(status+'\n'+detail+'\n');print(status,detail,flush=True)
 def job(self,t):return next((j for j in self.jobs if j['task']==t),None)
 def used(self,predicate=lambda j:True):return sum(int(j['actual_core_seconds'])+int(j['reserved_core_seconds']) for j in self.jobs if predicate(j))
 def refresh(self):
  if not self.jobs:return
  assert all(j['job_id'].isdigit() for j in self.jobs),'Ambiguous dispatch requires reconciliation; never submit a duplicate'
  raw=subprocess.check_output(['sacct','-n','-P','-j',','.join(j['job_id'] for j in self.jobs),'--format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,QOS'],text=True)
  (self.s/'accounting/slurm_snapshot.txt').write_text(raw);records={v[0]:v for line in raw.splitlines() if len(v:=line.split('|'))>=7}
  for j in self.jobs:
   if j['job_id'] not in records:continue
   r=records[j['job_id']];j.update(state=r[1],exit_code=r[2],allocated_cpus=r[3],elapsed_seconds=r[4])
   if r[1].split()[0] in TERMINAL:
    parts=[v for k,v in records.items() if k==j['job_id'] or k.startswith(j['job_id']+'.')]
    j.update(actual_core_seconds=max(int(v[3] or 0)*int(v[4] or 0) for v in parts),reserved_core_seconds=0)
  self.save();self.budget()
 def future(self):
  # Complete-work forecasts protect uncommitted work without double-counting live reservations.
  costpath=self.s/'accounting/cost_admission.tsv';forecast=float(read(costpath)[0]['forecast_primary_seconds']) if costpath.exists() else 100800
  future=0.;items=[]
  locked=(self.s/'controller/SELECTION_LOCKED').exists()
  for letter in 'ABCD':
   o='D_'+letter;done=self.ready(o+'_replay2b');remaining=0 if done or locked else max(0,forecast-self.used(lambda j:j['task'].startswith(o+'_')))
   items.append((o,remaining));future+=remaining
  triggered=(self.s/'accounting/reserve_release.tsv').exists()
  if not triggered and not locked:items.append(('protected_E_F',216000));future+=216000
  elif not locked:
   for letter in 'EF':
    done=(self.s/'controller'/('reserve_'+letter+'_terminal.txt')).exists()
    charge=self.used(lambda j:j['task'].startswith(('D_'+letter+'_','F_'+letter+'_')))
    remaining=0 if done else max(0,108000-charge);future+=remaining;items.append(('reserve_'+letter,remaining))
  engineering_left=0 if self.ready('T03F') and self.ready('T02D') and self.ready('T01') else max(0,28800-self.used(lambda j:j['task'].startswith('T')))
  future+=engineering_left;items.append(('engineering',engineering_left))
  retry=28800 if kv(self.s/'accounting/retry_ledger.tsv')['used']=='0'else 0
  assessment=max(0,28800-self.used(lambda j:j['task']=='C0'));collection=max(0,14400-self.used(lambda j:j['task']=='C1'))
  future+=retry+assessment+collection;items += [('technical_retry',retry),('assessment',assessment),('collection',collection)]
  return future,items
 def budget(self):
  settled=sum(int(j['actual_core_seconds']) for j in self.jobs);live=sum(int(j['reserved_core_seconds']) for j in self.jobs)
  future,items=self.future();table(self.s/'accounting/future_obligations.tsv',[dict(purpose=k,uncommitted_core_seconds=v) for k,v in items])
  table(self.s/'accounting/budget_ledger.tsv',[dict(cap_core_seconds=CAP,settled_core_seconds=settled,live_reserved_core_seconds=live,uncommitted_protected_future_core_seconds=future,unencumbered_core_seconds=CAP-settled-live-future)])
  assert settled+live<=CAP,'Allocated budget cap exceeded'
  # Complete work, not just the next submission, must remain affordable.
  assert settled+live+future<=CAP+1e-6,'Complete-work obligations exceed independent budget'
 def ready(self,t):
  j=self.job(t)
  if not j or j['state']!='COMPLETED':return False
  assert j['exit_code']=='0:0','Nonzero terminal exit '+t
  p=self.s/'verification'/t;manifest=p/'phase_manifest.tsv'
  if not manifest.exists():return False
  h=sha(manifest)
  if self.cache.get(t)==h:return True
  for r in read(manifest):
   q=p/r['path'];assert q.stat().st_size==int(r['bytes']) and sha(q)==r['sha256'],'Invalid receipt '+str(q)
  for name in ('application_result.tsv','gate.tsv','policy_gate.tsv'):assert read(p/name)[0]['status']=='PASS','Failed '+t+'/'+name
  assert all(int(v)==0 for v in read(p/'phase_exit_codes.tsv')[0].values()),'Independent exits failed '+t
  if t not in ('T02D','T03F') and (p/'read_denials.tsv').exists():assert not (p/'read_denials.tsv').read_text().strip(),'Unexpected read denial '+t
  self.cache[t]=h;return True
 def submit(self,t):
  if self.job(t):return
  verify(self.s);self.budget();actual=self.used()
  if t=='C0':envelope=28800-self.used(lambda j:j['task']=='C0')
  elif t=='C1':envelope=14400-self.used(lambda j:j['task']=='C1')
  elif t.startswith('T'):envelope=28800-self.used(lambda j:j['task'].startswith('T'))
  elif t[:3] in ('D_A','D_B','D_C','D_D'):
   path=self.s/'accounting/cost_admission.tsv';share=float(read(path)[0]['forecast_primary_seconds']) if path.exists() else 100800
   envelope=share-self.used(lambda j:j['task'].startswith(t[:3]+'_'))
  else:
   letter=t[2];envelope=108000-self.used(lambda j:j['task'].startswith(('D_'+letter+'_','F_'+letter+'_')))
  # These reservations are site/budget-derived; no application stopwatch exists.
  seconds=int(min(self.site,envelope,CAP-actual)//60*60);assert seconds>0,'Resource envelope exhausted '+t
  row=dict(job_id='submission_pending',task=t,state='RESERVED',exit_code='',allocated_cpus='',elapsed_seconds='',actual_core_seconds=0,reserved_core_seconds=seconds);self.jobs.append(row);self.save()
  try:self.budget()
  except Exception:raise RuntimeError('Full remaining-work resource admission failed before '+t)
  runner=self.s/('assessment_release/assessment_runner.sh' if t=='C0' else 'assessment_release/audit_runner.sh' if t=='C1' else 'release/runner.sh')
  if t in ('C0','C1'):
   for r in read(self.s/'assessment_release/manifest.tsv'):
    p=self.s/'assessment_release'/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],'Post-selection source mismatch'
  args=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1','--cpus-per-task=1','--hint=nomultithread','--mem=16G','--no-requeue','--time='+str(seconds//60),'--job-name=d19r_'+t,'--output='+str(self.s/'logs'/(t+'_%j.txt')),str(runner),str(self.s),t]
  jid=subprocess.check_output(args,text=True).strip().split(';')[0];assert jid.isdigit(),'Uncertain sbatch result';row.update(job_id=jid,state='SUBMITTED');self.save();print('SUBMITTED',t,jid,'reserved_core_seconds',seconds,flush=True)
 def admission(self):
  rows=[]
  for j in self.jobs:
   if j['task'].startswith('D_') and self.ready(j['task']):
    p=self.s/'verification'/j['task']/'callback_timings.tsv'
    if p.exists():rows+=read(p)
  c=cost(rows);forecast=84000*c+3600
  # Planned28h per primary; unused completed work can support a documented forecast.
  # Reforecast through the same nonoverlapping ledger used for submissions.
  table(self.s/'accounting/cost_admission.tsv',[dict(status='REFORECAST',c_D_core_seconds=c,forecast_primary_seconds=forecast)])
  fits=self.used()+self.future()[0]<=CAP
  table(self.s/'accounting/cost_admission.tsv',[dict(status='PASS' if fits else'RESOURCE_ADMISSION_HOLD',c_D_core_seconds=c,forecast_primary_seconds=forecast,callbacks=len(rows),feasible_callbacks=sum(truth(r['feasible']) for r in rows),full_primary_calls=84000,complete_work_admitted=fits)])
  if not fits:self.halt('RESOURCE_ADMISSION_HOLD','Complete registered remaining work plus reserve/quality/collection/retry protection is unaffordable')
  return fits
 def chain(self,o):
  for stage in DCHAIN:
   t=o+'_'+stage
   if not self.job(t):self.submit(t);return False
   if not self.ready(t):return False
  return True
 def bank_path(self,o):
  for stage in FCHAIN:
   t=o+'_'+stage
   if self.ready(t) and kv(self.s/'verification'/t/'result.tsv')['status']=='BANK_READY':return t
  return None
 def reserve(self,l):
  marker=self.s/'controller'/('reserve_'+l+'_terminal.txt')
  if marker.exists():return True
  o='F_'+l;path=self.bank_path(o)
  if not path:
   for stage in FCHAIN:
    t=o+'_'+stage
    if not self.job(t):self.submit(t);return False
    if not self.ready(t):return False
    status=kv(self.s/'verification'/t/'result.tsv')['status']
    if status=='BANK_READY':path=t;break
    if status=='FEASIBILITY_CONSTRUCTION_UNRESOLVED':marker.write_text('BANK_NOT_ADMITTED\n');return True
  pointer=self.s/'controller'/(o+'_bank_path.txt')
  if not pointer.exists():pointer.write_text(path+'\n')
  assert pointer.read_text().strip()==path,'Reserve bank pointer changed'
  for stage in ['check1','check2','robust']:
   t=o+'_'+stage
   if not self.job(t):self.submit(t);return False
   if not self.ready(t):return False
  # Full validation and construction costs feed complete family admission.
  rows=[]
  for stage in FCHAIN:
   phase=self.s/'verification'/(o+'_'+stage)/'callback_timings.tsv'
   if phase.exists() and self.ready(o+'_'+stage):rows+=read(phase)
  cF=cost(rows)
  v=[]
  for stage in ['check1','check2','robust']:v+=read(self.s/'verification'/(o+'_'+stage)/'callback_timings.tsv')
  cV=cost(v);cD=float(read(self.s/'accounting/cost_admission.tsv')[0]['c_D_core_seconds']);forecast=84000*cD+12000*cF+688*cV+3600
  table(self.s/'accounting'/('reserve_'+l+'_admission.tsv'),[dict(c_F=cF,c_V=cV,c_D=cD,forecast_family_core_seconds=forecast,status='PASS' if forecast<=108000 else'RESOURCE_ADMISSION_HOLD')])
  if forecast>108000:self.halt('RESOURCE_ADMISSION_HOLD','Complete reserve family forecast exceeds protected share; reforecast required');return False
  if self.chain('D_'+l):marker.write_text('D19_CHAIN_COMPLETE\n');return True
  return False
 def classifications(self,letters):
  rows=[]
  for l in letters:
   t='D_'+l+'_replay2b'
   if self.ready(t):rows+=read(self.s/'verification'/t/'quality.tsv')
  table(self.s/'reports/verified_candidates.tsv',rows or [dict(status='NO_VERIFIED_CANDIDATE')])
  return {r['origin'] for r in rows if truth(r.get('prediction_pass',False))},rows
 def close_selection(self,letters):
  passes,rows=self.classifications(letters);best=[]
  for l in letters:
   p=self.s/'verification'/('D_'+l+'_replay2b')/'best_identity.tsv'
   if p.exists():
    # RDS is not loaded by metadata owner. Final best identity was exported by R.
    q=read(self.s/'verification'/('D_'+l+'_replay2b')/'best_identity.tsv');best+=q
  assert best,'No completed bounded origin'
  lowest=min(best,key=lambda r:float(r['D']));q=next((r for r in rows if r['key']==lowest['key'] and r['origin']==lowest['origin']),None)
  qualified=q is not None and truth(q['verified']) and truth(q['local_pass'])
  status='DET_OPTIMIZER_UNRESOLVED' if len(passes)<2 else'D19_PREDICTION_RECOVERY_PASS' if qualified else'D19_PREDICTION_RECOVERED_LOCAL_UNRESOLVED'
  table(self.s/'reports/selection_lock.tsv',[dict(status=status,independent_prediction_origins=len(passes),best_found_origin=lowest['origin'],best_found_key=lowest['key'],best_found_D=lowest['D'],best_found_replay_verified=q is not None,best_found_local_pass=qualified,N20=0,profiles=0,tumor=0)])
  (self.s/'controller/SELECTION_LOCKED').write_text(status+'\n');print('SELECTION_LOCKED',status,flush=True)
 def finish_report(self):
  lock=read(self.s/'reports/selection_lock.tsv')[0];used=sum(int(j['actual_core_seconds']) for j in self.jobs)
  assert used<=CAP and all(j['state']=='COMPLETED' and j['exit_code']=='0:0' for j in self.jobs),'Final allocated/numerical review required'
  assessments=read(self.s/'verification/C0/assessment_gate.tsv');audits=read(self.s/'verification/C1/audit_gate.tsv')
  assert assessments[0]['status']=='PASS' and audits[0]['status']=='PASS'
  rows=read(self.s/'reports/verified_candidates.tsv');counts=[]
  for l in 'ABCDEF':
   p=self.s/'verification'/('D_'+l+'_replay2b')/'quality.tsv'
   if p.exists():
    q=read(p);b=min(q,key=lambda r:float(r['D']));g=[r for r in q if truth(r['prediction_pass'])]
    counts.append(dict(origin='D_'+l,best_verified_D=b['D'],prediction_pass=bool(g),best_verified_local_pass=b['local_pass'],gate_key=min(g,key=lambda r:float(r['D']))['key'] if g else 'NONE',admission='BOUNDED_D19_CHAIN_COMPLETE'))
   elif l in 'EF' and (self.s/'accounting/reserve_release.tsv').exists():
    counts.append(dict(origin='D_'+l,best_verified_D='NA',prediction_pass=False,best_verified_local_pass=False,gate_key='NONE',admission=(self.s/'controller'/('reserve_'+l+'_terminal.txt')).read_text().strip()))
  table(self.s/'reports/origin_verdicts.tsv',counts)
  external=sum(int(r['full_forward_attempts']) for r in read(self.s/'verification/T01/attempt_receipt.tsv'))
  assert external<=512
  attempts=[]
  for l in 'ABCDEF':
   o='D_'+l;n=0
   for j in self.jobs:
    if j['task'].startswith(o+'_'):
     p=self.s/'verification'/j['task']/'evaluations.tsv'
     if p.exists():n+=sum(r['event']=='START' for r in read(p))
   if n:assert n<=84000;attempts.append(dict(origin=o,kind='D19_fit',full_bundle_attempts=n,cap=84000))
  for l in 'EF':
   for kind in ['constructor','validation']:
    n=0
    for j in self.jobs:
     if not j['task'].startswith('F_'+l+'_'):continue
     validation=j['task'].endswith(('_check1','_check2','_robust'))
     if validation!=(kind=='validation'):continue
     p=self.s/'verification'/j['task']/'evaluations.tsv'
     if p.exists():n+=sum(r['event']=='START' for r in read(p))
    if n:assert n<=(688 if kind=='validation'else 12000);attempts.append(dict(origin='F_'+l,kind=kind,full_bundle_attempts=n,cap=688 if kind=='validation'else 12000))
  attempts.append(dict(origin='T01',kind='external_engineering',full_bundle_attempts=external,cap=512));assert sum(r['full_bundle_attempts'] for r in attempts)<=529888
  table(self.s/'accounting/candidate_attempts.tsv',attempts)
  table(self.s/'accounting/final_budget.tsv',[dict(cap_core_seconds=CAP,settled_core_seconds=used,settled_core_hours=used/3600,live_reserved_core_seconds=0,unspent_core_seconds=CAP-used,technical_retries_used=kv(self.s/'accounting/retry_ledger.tsv')['used'],scientific_reserve_rounds_used=int((self.s/'accounting/reserve_release.tsv').exists()),external_full_forwards=external)])
  lines=['# D19 deterministic recovery restart report','',lock['status'],'',
   'This is a conditional culture prediction-recovery engineering study. It is not a structural-identifiability proof, a global-optimum certificate, a confidence interval or a noisy-data recovery result.','',
   'All14 mechanisms and6 nuisance roles were retained. The19 unknown coordinates include the four root-distribution parameters and sigma_kary; sigma_growth was the known generating control and was injected directly from the original D19 RDS.','',
   '| Origin | Best verified D | Prediction gate | Best verified local gate |','|---|---:|---|---|']
  lines += ['| {origin} | {best_verified_D} | {prediction_pass} | {best_verified_local_pass} |'.format(**r) for r in counts]
  lines += ['', 'Objective-selected lowest found point: '+lock['best_found_origin']+' / '+lock['best_found_key']+'; D='+lock['best_found_D']+'. Its own replay/local flags, rather than a worse checked point, determine the optimizer-quality label.','',
   'Sigma_kary remained free over0.01-10. Its raw coordinate, observation-kernel difference, latent/observed distributions and root-width trade-offs are reported separately. A local binary64-flat kernel does not establish a whole-domain structural symmetry.','',
   'The post-selection assessor opened full truth only after selection/quality locking. See verification/C0/all20_parameter_recovery.tsv, sigma_kary_kernel_diagnostics.tsv, function_error_summary.tsv, the complete growth/karyotype tables and D19_recovery_diagnostics.pdf. No raw1% parameter-error gate was imposed.','',
   'Actual allocated use: '+format(used/3600,'.6f')+' of200 core-hours. Scientific reserve rounds: '+str(int((self.s/'accounting/reserve_release.tsv').exists()))+' of1. Technical retries: '+kv(self.s/'accounting/retry_ledger.tsv')['used']+' of1. External engineering full forwards: '+str(external)+' of512. All original and duplicate/interrupted fitting calls remain charged in origin states/journals.','',
   'All recorded jobs used xxlarge. No application concurrency cap or discretionary elapsed-time limit was introduced. Scheduler reservations came from site/budget envelopes. Full raw traces/native checkpoints and the authoritative inventory remain on HPC; local collection is a separately verified thin subset.','',
   'The12 unresolved structural cases and outstanding experimental SOP/flow/live-dead questions retain their prior status. No N20, profiles, tumor fitting, new observations or additional generating regime was launched.']
  (self.s/'REPORT.md').write_text('\n'.join(lines)+'\n')
  table(self.s/'reports/final_gate.tsv',[dict(status=lock['status'],source_numerical_read_audit='PASS',postselection_assessment='PASS',jobs=len(self.jobs),active_owned_jobs=0,N20=0,profiles=0,tumor=0)])
  table(self.s/'accounting/budget_ledger.tsv',[dict(cap_core_seconds=CAP,settled_core_seconds=used,live_reserved_core_seconds=0,uncommitted_protected_future_core_seconds=0,unencumbered_core_seconds=CAP-used,status='CLOSED')])
  table(self.s/'accounting/future_obligations.tsv',[dict(purpose='CLOSED; unused authority is not released for new work',uncommitted_core_seconds=0)])
  (self.s/'controller/FINISHED').write_text(lock['status']+'\n');print('FINISHED',lock['status'],flush=True)
 def tick(self):
  self.refresh()
  for j in self.jobs:
   if j['state'].split()[0] in TERMINAL and (j['state']!='COMPLETED' or j['exit_code']!='0:0'):
    self.halt('TECHNICAL_REVIEW_REQUIRED',j['task']+' '+j['job_id']+' '+j['state']+' '+j['exit_code']+'; no implicit restart');return
  if (self.s/'controller/SELECTION_LOCKED').exists():
   for t in ['C0','C1']:
    if not self.job(t):self.submit(t);return
    if not self.ready(t):return
   self.finish_report();return
  for t in ['T00','T01','T02D','T03F']:
   if not self.job(t):self.submit(t);return
   if not self.ready(t):return
  for t in ['D_A_import','D_A_bench']:
   if not self.job(t):self.submit(t);return
   if not self.ready(t):return
  if not self.admission():return
  complete=[self.chain('D_'+l) for l in 'ABCD']
  if not all(complete):return
  passes,rows=self.classifications('ABCD')
  if len(passes)>=2:self.close_selection('ABCD');return
  release=self.s/'accounting/reserve_release.tsv'
  if not release.exists():table(release,[dict(round=1,families='E,F',reason='fewer than two primary verified prediction passes',extra_rounds=0)])
  complete=[self.reserve(l) for l in 'EF']
  if all(complete):self.close_selection('ABCDEF')
def main():
 s=Path(sys.argv[1]).resolve();lock=(s/'controller/owner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 table(s/'controller/ownership.tsv',[dict(pid=os.getpid(),scope='single new D19 metadata-only owner; old work untouched',study=str(s))]);owner=Owner(s)
 while not any((s/'controller'/f).exists() for f in ['STOP','FINISHED']):
  try:owner.tick()
  except Exception as e:owner.halt('CONTROLLER_RECONCILIATION_HOLD',str(e));break
  time.sleep(20)
if __name__=='__main__':main()
