"""Independent finite 96-core-hour reduction audit. TSV only; old studies read-only."""
from pathlib import Path
import argparse,csv,datetime,fcntl,hashlib,os,subprocess,tempfile
SLURM=Path('/cm/shared/apps/slurm/current/bin')
CAPS={k:int(v*3600) for k,v in dict(engineering=8,time=12,aggregation=12,within=16,tail=24,history=12,report=4,reserve=8).items()}
CALL_CAPS={'native':640,'reference':2000,'history':3}
FIELDS=['task','phase','case_id','release','job_id','state','exit_code','cpus','wall_seconds','reserved_core_seconds','actual_core_seconds','retry_of','created_utc']
TERMINAL={'COMPLETED','FAILED','TIMEOUT','CANCELLED','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE'}
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(path):
 with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for data in iter(lambda:f.read(1048576),b''):h.update(data)
 return h.hexdigest()
def table(path,rows,fields):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
 with tempfile.NamedTemporaryFile('w',dir=p.parent,delete=False,newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows);f.flush();os.fsync(f.fileno());temp=f.name
 os.replace(temp,p)
def kv(path):return {r['key']:r['value'] for r in read(path)}
def typed(path,data):table(path,[dict(key=k,value=str(v).lower() if isinstance(v,bool) else str(v)) for k,v in data.items()],['key','value'])
def position(rows):
 out=[];total=0
 for phase,cap in CAPS.items():
  actual=sum(int(r['actual_core_seconds']) for r in rows if r['phase']==phase and r['state'] in TERMINAL)
  live=sum(int(r['reserved_core_seconds']) for r in rows if r['phase']==phase and r['state'] not in TERMINAL)
  if actual+live>cap:raise ValueError('Phase budget exceeded: '+phase)
  total+=actual+live;out.append(dict(phase=phase,actual=actual,committed=live,cap=cap,remaining=cap-actual-live))
 if total>96*3600:raise ValueError('Independent 96-core-hour ceiling exceeded')
 return out
def verify(study,release):
 for r in read(release/'manifest.tsv'):
  if sha(release/r['path'])!=r['sha256']:raise ValueError('Frozen release changed: '+r['path'])
 for r in read(study/'inputs/manifest.tsv'):
  if sha(study/'inputs'/r['path'])!=r['sha256']:raise ValueError('Frozen input changed: '+r['path'])
 project=Path(kv(study/'study.tsv')['project'])
 for r in read(study/'inputs/external_sources.tsv'):
  if sha(project/r['path'])!=r['sha256']:raise ValueError('Read-only source changed: '+r['path'])
def reserve_call(study,kind,label,job_id):
 if kind not in CALL_CAPS:raise ValueError('Unregistered call type')
 if (study/'COMPLETE').exists():raise ValueError('Study closed')
 with (study/'state/CALL_LOCK').open('a') as f:
  fcntl.flock(f,fcntl.LOCK_EX)
  p=study/'accounting/call_ledger.tsv';rows=read(p) if p.exists() else []
  count=sum(r['kind']==kind for r in rows)
  if count>=CALL_CAPS[kind]:raise ValueError('Cumulative call ceiling reached: '+kind)
  rows.append(dict(slot=len(rows)+1,kind=kind,kind_slot=count+1,label=label,job_id=job_id,utc=utc(),status='consumed_before_call'))
  table(p,rows,['slot','kind','kind_slot','label','job_id','utc','status']);return count+1
def reconcile(study):
 rows=read(study/'accounting/jobs.tsv')
 if (study/'COMPLETE').exists():return rows
 for r in rows:
  if r['state'] in TERMINAL:continue
  if not r['job_id']:raise ValueError('Unresolved submission intent')
  output=subprocess.run([str(SLURM/'sacct'),'-n','-P','-X','-j',r['job_id'],'--format=JobIDRaw,State%30,ExitCode,AllocCPUS,ElapsedRaw'],check=True,stdout=subprocess.PIPE,universal_newlines=True).stdout
  record=[line.split('|') for line in output.splitlines() if line.split('|')[0]==r['job_id']]
  if record:
   v=record[0];r['state']=v[1].split()[0];r['exit_code']=v[2]
   if r['state'] in TERMINAL:
    r['actual_core_seconds']=str(int(v[3])*int(v[4]))
    if int(r['actual_core_seconds'])>int(r['reserved_core_seconds']):raise ValueError('Allocation reservation overrun')
 table(study/'accounting/jobs.tsv',rows,FIELDS);table(study/'accounting/budget_ledger.tsv',position(rows),['phase','actual','committed','cap','remaining']);return rows
def payload_ok(path):
 if not (path/'COMPLETE').exists():return False
 for r in read(path/'payload_manifest.tsv'):
  if sha(path/r['path'])!=r['sha256']:raise ValueError('Result payload mismatch')
 return (path/'COMPLETE').read_text().strip()==sha(path/'payload_manifest.tsv')
def admission(study,phase):
 if phase=='engineering':
  if kv(study/'inputs/local_preflight.tsv').get('pass')!='true':raise ValueError('Local tests not passed')
 elif phase=='time':
  for name in ('pilot_c01','pilot_c23'):
   p=study/'results'/name
   if not payload_ok(p) or kv(p/'result.tsv').get('reference_gate')!='true':raise ValueError('Pilot/reference gate not passed')
 elif phase in ('aggregation','tail'):
  for c in read(study/'inputs/cases.tsv'):
   p=study/'results'/('time_'+c['case_id'])
   if not payload_ok(p) or kv(p/'result.tsv').get('time_gate')!='true':raise ValueError('Time/reference gate incomplete')
 elif phase=='within':
  for c in read(study/'inputs/cases.tsv'):
   if not payload_ok(study/'results'/('aggregation_'+c['case_id'])):raise ValueError('Aggregation evidence incomplete')
 elif phase=='history':
  if kv(study/'reports/candidate_selection.tsv').get('release_history')!='true':raise ValueError('Candidate/history gate not passed')
def submit(study,task,phase,case_id,release,hours,memory,retry_of=''):
 if (study/'COMPLETE').exists():raise ValueError('Study closed')
 rows=reconcile(study);verify(study,study/'releases'/release)
 if any(r['task']==task for r in rows):raise ValueError('Task already registered')
 failed=[r for r in rows if r['state'] in TERMINAL and r['state']!='COMPLETED']
 if phase=='reserve':
  if any(r['phase']=='reserve' for r in rows):raise ValueError('Single review/retry allowance exhausted')
  if retry_of not in {r['task'] for r in failed}:raise ValueError('Retry must cite a failed owned task')
  if not (study/'reports'/('review_'+retry_of+'.md')).exists():raise ValueError('Failure review missing')
 elif failed:
  reviewed={r['task'] for r in read(study/'accounting/reviewed_failures.tsv')} if (study/'accounting/reviewed_failures.tsv').exists() else set()
  if any(r['task'] not in reviewed for r in failed):raise ValueError('Unreviewed failure; stop dispatch')
 run_phase=phase
 if phase=='reserve':run_phase=next(r['phase'] for r in rows if r['task']==retry_of)
 admission(study,run_phase)
 seconds=int(hours*3600);r=dict(task=task,phase=phase,case_id=case_id,release=release,job_id='',state='INTENT',exit_code='',cpus='1',wall_seconds=seconds,reserved_core_seconds=seconds,actual_core_seconds='',retry_of=retry_of,created_utc=utc())
 position(rows+[r]);table(study/'accounting/jobs.tsv',rows+[r],FIELDS)
 out=subprocess.run([str(SLURM/'sbatch'),'--parsable','--partition=red','--qos=xxlarge','--no-requeue','--cpus-per-task=1','--mem='+memory,'--time='+str(int(hours*60)),'--job-name=rr96_'+task,'--comment=rr96:'+study.name+':'+task,'--chdir='+str(study),'--output='+str(study/'logs'/('%x_%j.out')),'--error='+str(study/'logs'/('%x_%j.err')),str(study/'releases'/release/'batch.sh'),str(study),release,task,run_phase,case_id],check=True,stdout=subprocess.PIPE,universal_newlines=True)
 r['job_id']=out.stdout.strip().split(';')[0];r['state']='SUBMITTED';table(study/'accounting/jobs.tsv',rows+[r],FIELDS);print(task,r['job_id'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('action',choices=['verify','submit','reconcile','call']);p.add_argument('--release',default='r01');p.add_argument('--task');p.add_argument('--phase',choices=CAPS);p.add_argument('--case-id',default='none');p.add_argument('--hours',type=float,default=.5);p.add_argument('--memory',default='16G');p.add_argument('--retry-of',default='');p.add_argument('--kind',choices=CALL_CAPS);a=p.parse_args();s=a.study.resolve()
 if a.action=='verify':verify(s,s/'releases'/a.release)
 elif a.action=='call':print(reserve_call(s,a.kind,a.task,os.environ.get('SLURM_JOB_ID','local_fixture')))
 else:
  with (s/'state/JOB_LOCK').open('a') as f:
   fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
   if a.action=='reconcile':print(position(reconcile(s)))
   else:submit(s,a.task,a.phase,a.case_id,a.release,a.hours,a.memory,a.retry_of)
