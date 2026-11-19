"""Finite, isolated Step 04 allocations; TSV state, no automatic resubmission."""
from pathlib import Path
import argparse,csv,datetime,fcntl,hashlib,os,subprocess,tempfile
SLURM=Path('/cm/shared/apps/slurm/current/bin')
CAPS={'tool':8*3600,'reduction':24*3600,'structural':120*3600,'report':8*3600,'infrastructure':8*3600}
FIELDS=['task','phase','mode','release','job_id','state','exit_code','cpus','wall_seconds','reservation','actual','created_utc']
TERMINAL={'COMPLETED','FAILED','TIMEOUT','CANCELLED','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE'}
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,rows,fields):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with tempfile.NamedTemporaryFile('w',dir=p.parent,delete=False,newline='') as f:
  writer=csv.DictWriter(f,fields,delimiter='\t');writer.writeheader();writer.writerows(rows);f.flush();os.fsync(f.fileno());tmp=f.name
 os.replace(tmp,p)
def verify(release):
 for r in read(release/'manifest.tsv'):
  if sha(release/r['path'])!=r['sha256']:raise ValueError('Frozen release mismatch: '+r['path'])
def position(rows):
 total=0;out=[]
 for phase,cap in CAPS.items():
  actual=sum(int(r['actual']) for r in rows if r['phase']==phase and r['state'] in TERMINAL)
  live=sum(int(r['reservation']) for r in rows if r['phase']==phase and r['state'] not in TERMINAL)
  if actual+live>cap:raise ValueError('Phase ceiling exceeded: '+phase)
  total+=actual+live;out.append(dict(phase=phase,actual_core_seconds=actual,committed_core_seconds=live,ceiling_core_seconds=cap,remaining_core_seconds=cap-actual-live))
 if total>200*3600:raise ValueError('Independent 200-core-hour ceiling exceeded')
 return out
def reconcile(study):
 p=study/'accounting/jobs.tsv';rows=read(p)
 for r in rows:
  if r['state'] in TERMINAL:continue
  if not r['job_id']:raise ValueError('Ambiguous submission intent; manual review required')
  result=subprocess.run([str(SLURM/'sacct'),'-n','-P','-X','-j',r['job_id'],'--format=JobIDRaw,State%30,ExitCode,AllocCPUS,ElapsedRaw'],check=True,stdout=subprocess.PIPE,universal_newlines=True)
  found=[line.split('|') for line in result.stdout.splitlines() if line.split('|')[0]==r['job_id']]
  if found:
   x=found[0];r['state']=x[1].split()[0];r['exit_code']=x[2]
   if r['state'] in TERMINAL:
    r['actual']=str(int(x[3])*int(x[4]))
    if int(r['actual'])>int(r['reservation']):raise ValueError('Allocation overrun')
 table(p,rows,FIELDS);table(study/'accounting/budget_ledger.tsv',position(rows),['phase','actual_core_seconds','committed_core_seconds','ceiling_core_seconds','remaining_core_seconds'])
 return rows
def submit(study,task,phase,mode,release,hours,memory):
 rows=reconcile(study);verify(study/'releases'/release)
 expected={'tool_setup':'tool','tool_postfreeze':'tool','reduction':'reduction','exact_analysis':'structural','julia_global':'structural','report':'report'}
 if mode not in expected or (phase!=expected[mode] and phase!='infrastructure'):raise ValueError('Unregistered mode/phase')
 if mode in ('reduction','exact_analysis','julia_global'):
  gate={r['key']:r['value'] for r in read(study/'inputs/tool_gate.tsv')}
  if gate.get('pass')!='true':raise ValueError('Tool gate not passed')
 if mode in ('exact_analysis','julia_global'):
  gate={r['key']:r['value'] for r in read(study/'results/reduction_native_rate_exact/result.tsv')}
  if gate.get('source_and_kernel_integrity')!='true':raise ValueError('Source/kernel validation not passed')
  if not (study/'results/tool_postfreeze/CONSTRUCTOR_PASS').exists():raise ValueError('Fresh frozen-tool constructor check not passed')
 if any(r['task']==task for r in rows):raise ValueError('Task already registered')
 if any(r['state'] in TERMINAL and r['state']!='COMPLETED' for r in rows) and phase!='infrastructure':
  reviewed=read(study/'accounting/reviewed_failures.tsv') if (study/'accounting/reviewed_failures.tsv').exists() else []
  if any(r['state'] in TERMINAL and r['state']!='COMPLETED' and r['task'] not in {x['task'] for x in reviewed} for r in rows):raise ValueError('Unreviewed allocation failure')
 if phase=='structural' and sum(r['phase']=='structural' for r in rows)>=24:raise ValueError('Finite symbolic-task count reached')
 if phase=='infrastructure' and any(r['phase']=='infrastructure' for r in rows):raise ValueError('Infrastructure retry already consumed')
 seconds=int(hours*3600);receipt=dict(task=task,phase=phase,mode=mode,release=release,job_id='',state='INTENT',exit_code='',cpus='1',wall_seconds=seconds,reservation=seconds,actual='',created_utc=utc())
 position(rows+[receipt]);table(study/'accounting/jobs.tsv',rows+[receipt],FIELDS)
 result=subprocess.run([str(SLURM/'sbatch'),'--parsable','--partition=red','--qos=xxlarge','--no-requeue','--cpus-per-task=1','--mem='+memory,'--time='+str(int(hours*60)),'--job-name=si4_'+task,'--chdir='+str(study),'--output='+str(study/'logs'/('%x_%j.out')),'--error='+str(study/'logs'/('%x_%j.err')),str(study/'releases'/release/'batch.sh'),str(study),release,task,mode],check=True,stdout=subprocess.PIPE,universal_newlines=True)
 receipt['job_id']=result.stdout.strip().split(';')[0];receipt['state']='SUBMITTED';table(study/'accounting/jobs.tsv',rows+[receipt],FIELDS);print(task,receipt['job_id'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('action',choices=['reconcile','submit','verify','native']);p.add_argument('--task');p.add_argument('--phase',choices=CAPS);p.add_argument('--mode');p.add_argument('--release',default='r01');p.add_argument('--hours',type=float,default=2);p.add_argument('--memory',default='64G');a=p.parse_args();s=a.study.resolve()
 if a.action=='verify':verify(s/'releases'/a.release)
 else:
  with (s/'state/LOCK').open('a') as f:
   fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
   if a.action=='native':
    p=s/'state/native_calls.tsv';rows=read(p) if p.exists() else []
    if len(rows)>=48:raise ValueError('Cumulative native-call allowance exhausted')
    rows.append(dict(slot=len(rows)+1,task=a.task,status='consumed_before_call',source='',sha256=''));table(p,rows,['slot','task','status','source','sha256']);print(len(rows))
   elif a.action=='submit':submit(s,a.task,a.phase,a.mode,a.release,a.hours,a.memory)
   else:print(position(reconcile(s)))
