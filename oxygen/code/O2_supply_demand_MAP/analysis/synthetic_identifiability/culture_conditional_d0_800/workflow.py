"""Finite TSV-only dispatch and allocated-parent accounting for conditional D0."""
import contextlib,csv,datetime,fcntl,hashlib,io,os,subprocess,sys,tempfile
from pathlib import Path
sys.dont_write_bytecode=True
IMAGE_SHA='dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372'
SLURM=Path('/cm/shared/apps/slurm/current/bin')
CAPS={k:v*3600 for k,v in dict(engineering=24,positive=6,primary=420,reserve=210,report=12,infrastructure=128).items()}
PLAN={'prepare':('engineering',1,3600),'replay':('engineering',1,3600),'replay_final':('engineering',1,3600),'parallel':('engineering',8,1800),'engineering_check':('engineering',1,7200),'positive_local':('positive',1,18000),'positive_check':('positive',1,3600),'report':('report',1,10800)}
PLAN['local_api_check']=('engineering',1,1800)
PLAN['report_resume01']=('report',1,10800)
for origin in ('A','B','C','D','E','F'):
    phase='primary' if origin in 'ABCD' else 'reserve'
    for step,cores,wall in (('screen',1,3600),('global',8,43200),('globalcheck',1,3600),('local',1,21600),('check',1,3600)):
        PLAN[step+'_'+origin]=(phase,cores,wall)
FIELDS=['task','phase','mode','release','job_id','state','exit_code','cores','wall_seconds','reserved_core_seconds','actual_core_seconds','intent_utc','dispatch_sha256','comment']
TERMINAL={'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'}
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def text(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);fd,tmp=tempfile.mkstemp(prefix='.pending_',dir=str(p.parent))
    try:
        with os.fdopen(fd,'w') as f:f.write(value);f.flush();os.fsync(f.fileno())
        os.replace(tmp,str(p))
    finally:
        if os.path.exists(tmp):os.unlink(tmp)
def table(p,rows,fields):
    s=io.StringIO();w=csv.DictWriter(s,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows);text(p,s.getvalue())
def kv(p):return {r['key']:r['value'] for r in read(p)}
def typed(p,x):table(p,[dict(key=k,value=str(v).lower() if isinstance(v,bool) else str(v)) for k,v in x.items()],['key','value'])
@contextlib.contextmanager
def lock(s,name='LOCK'):
    p=Path(s)/'state'/name;p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('a') as f:fcntl.flock(f,fcntl.LOCK_EX);yield
def jobs(s):
    p=Path(s)/'accounting/jobs.tsv';return read(p) if p.exists() else []
def closed(j):return j['state'].split('+')[0].split()[0] in TERMINAL
def success(j):return j['state']=='COMPLETED' and j['exit_code']=='0:0'
def verify(root,manifest):
    for r in read(manifest):
        if sha(Path(root)/r['path'])!=r['sha256']:raise ValueError('Fingerprint mismatch: '+r['path'])
def verify_inputs(s):
    s=Path(s);verify(Path(kv(s/'study.tsv')['project']),s/'manifests/external_manifest.tsv')
    if (s/'truth/manifest.tsv').exists():verify(s,s/'truth/manifest.tsv')
def payload(root):
    root=Path(root)
    if (root/'COMPLETE').read_text().strip()!=sha(root/'payload_manifest.tsv'):raise ValueError('Unsealed payload')
    verify(root,root/'payload_manifest.tsv')
def position(rows):
    if any(j['phase'] not in CAPS for j in rows):raise ValueError('Unknown accounting phase')
    out={}
    for phase,cap in CAPS.items():
        actual=sum(int(j['actual_core_seconds']) for j in rows if j['phase']==phase and closed(j))
        live=sum(int(j['reserved_core_seconds']) for j in rows if j['phase']==phase and not closed(j))
        out[phase]=dict(actual=actual,live=live,cap=cap,remaining=cap-actual-live)
        if actual+live>cap:raise ValueError('Phase cap exceeded: '+phase)
    if sum(x['actual']+x['live'] for x in out.values())>800*3600:raise ValueError('Whole package cap exceeded')
    for origin in 'ABCDEF':
        v=[j for j in rows if j['task'].endswith('_'+origin) and j['phase']!='infrastructure']
        cost=sum(int(j['actual_core_seconds']) if closed(j) else int(j['reserved_core_seconds']) for j in v)
        if cost>105*3600:raise ValueError('Origin cap exceeded: '+origin)
    return out
def reconcile(s):
    s=Path(s)
    with lock(s):
        rows=jobs(s)
        if rows:
            if any(not j['job_id'] for j in rows):raise ValueError('Ambiguous submission; reconcile comment before any resubmission')
            ids=','.join(j['job_id'] for j in rows)
            r=subprocess.run([str(SLURM/'sacct'),'-n','-X','-P','-j',ids,'--format=JobIDRaw,State%30,ExitCode,AllocCPUS,ElapsedRaw'],check=True,universal_newlines=True,stdout=subprocess.PIPE)
            text(s/'accounting/raw_parent_accounting.txt',r.stdout)
            raw={l.split('|')[0]:l.split('|') for l in r.stdout.splitlines() if l.strip()}
            for j in rows:
                if j['job_id'] not in raw:continue
                z=raw[j['job_id']];j.update(state=z[1],exit_code=z[2])
                if closed(j):
                    if int(z[3])!=int(j['cores']):raise ValueError('Allocated CPU mismatch')
                    j['actual_core_seconds']=str(int(z[3])*int(z[4]))
                    if int(j['actual_core_seconds'])>int(j['reserved_core_seconds']):raise ValueError('Reservation overrun')
            table(s/'accounting/jobs.tsv',rows,FIELDS)
        pos=position(rows);table(s/'accounting/budget_ledger.tsv',[dict(phase=k,**v) for k,v in pos.items()],['phase','actual','live','cap','remaining'])
        return rows
def submit(s,task,mode=None,release=None,retry=False):
    s=Path(s);mode=mode or task
    if release is None:release=kv(s/'tables/execution_release.tsv')['release'] if (s/'tables/execution_release.tsv').exists() else 'r01'
    with lock(s):
        rows=jobs(s)
        if any(j['task']==task for j in rows):return
        if (s/'COMPLETE').exists():raise ValueError('Study closed')
        if mode not in PLAN:raise ValueError('Undeclared task '+mode)
        phase,cores,wall=PLAN[mode]
        if retry:
            if len([j for j in rows if j['phase']=='infrastructure'])>=2:raise ValueError('Infrastructure retry allowance exhausted')
            phase='infrastructure'
        elif task!=mode:raise ValueError('Only explicit bounded infrastructure retries may rename tasks')
        if mode.startswith(('screen_','global_','globalcheck_','local_','check_','positive_')):
            gate=kv(s/'tables/engineering_gate.tsv')
            if gate['pass']!='true' or gate['cost_fits']!='true':raise ValueError('Engineering/cost release held')
            evidence=s/'tables/engineering_evidence_r02.tsv'
            verify(s,evidence if evidence.exists() else s/'tables/engineering_evidence.tsv')
        if mode.startswith('local_') and mode!='local_api_check':
            api=s/'engineering/local_api_check';payload(api)
            if kv(api/'result.tsv')['pass']!='true':raise ValueError('Local stopping API gate held')
        if mode.startswith(('screen_','global_','local_')) and mode[-1] in 'BCD':
            if kv(s/'tables/benchmark_gate.tsv')['pass']!='true':raise ValueError('First complete blind chain gate held')
        if position(rows)[phase]['remaining']<cores*wall:raise ValueError('Whole allocation does not fit remaining phase budget')
        rel=s/'releases'/release;verify(rel,rel/'manifest.tsv');verify_inputs(s)
        d=s/'dispatch'/task;d.mkdir(parents=True,exist_ok=False)
        typed(d/'spec.tsv',dict(task=task,mode=mode,phase=phase,cores=cores,wall_seconds=wall,release=release,release_manifest_sha256=sha(rel/'manifest.tsv')))
        comment='d0800:'+s.name+':'+task
        j=dict(task=task,phase=phase,mode=mode,release=release,job_id='',state='INTENT',exit_code='',cores=str(cores),wall_seconds=str(wall),reserved_core_seconds=str(cores*wall),actual_core_seconds='',intent_utc=utc(),dispatch_sha256=sha(d/'spec.tsv'),comment=comment)
        rows.append(j);position(rows);table(s/'accounting/jobs.tsv',rows,FIELDS)
        cmd=[str(SLURM/'sbatch'),'--parsable','--partition=red','--qos=xxlarge','--no-requeue','--nodes=1','--ntasks=1','--cpus-per-task='+str(cores),'--mem='+('64G' if cores==8 else '16G'),'--time={:02d}:{:02d}:{:02d}'.format(wall//3600,(wall%3600)//60,wall%60),'--job-name=d0800_'+task,'--comment='+comment,'--chdir='+str(s),'--output='+str(s/'logs'/(task+'_%j.log')),str(rel/'batch.sh'),str(s),task]
        text(d/'submission.txt','\n'.join(cmd)+'\n')
        r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        text(d/'receipt.txt',r.stdout);text(d/'submission_stderr.txt',r.stderr)
        if r.returncode or not r.stdout.split(';')[0].strip().isdigit():raise ValueError('Ambiguous sbatch receipt; dispatch held')
        j.update(job_id=r.stdout.split(';')[0].strip(),state='SUBMITTED');table(s/'accounting/jobs.tsv',rows,FIELDS)
        print(task,j['job_id'],flush=True)
def reserve_forward(s,task,job,pid):
    s=Path(s)
    with lock(s,'forward_LOCK'):
        p=s/'accounting/engineering_forward_slots.tsv';rows=read(p) if p.exists() else []
        if len(rows)>=128:raise ValueError('Engineering forward-call ceiling exhausted')
        row=dict(ordinal=len(rows)+1,task=task,job_id=job,pid=pid,utc=utc(),status='consumed_before_forward')
        rows.append(row);table(p,rows,list(row));print(row['ordinal'])
def initial_ordinal(s,task):
    s=Path(s);p=s/'fits'/task/'initial_counter.tsv'
    with lock(s,task+'_initial_LOCK'):
        n=int(kv(p)['count']) if p.exists() else 0
        if n>=256:print(0);return
        n+=1;typed(p,dict(count=n))
        if n==256:text(s/'fits'/task/'initial_counter_COMPLETE','256 reserved inputs\n')
        print(n)
if __name__=='__main__':
    mode=sys.argv[1];s=Path(sys.argv[2]).resolve()
    if mode=='counter':reserve_forward(s,*sys.argv[3:6])
    elif mode=='initial':initial_ordinal(s,sys.argv[3])
    elif mode=='verify':verify_inputs(s);rel=Path(sys.argv[3]);verify(rel,rel/'manifest.tsv');print('Frozen input and source fingerprints verified')
    elif mode=='snapshot':print(position(reconcile(s)))
    elif mode=='submit':submit(s,sys.argv[3])
    else:raise ValueError(mode)
