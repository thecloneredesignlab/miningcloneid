"""Single-owner, integer-accounted dispatch for the approved culture D0 pilot."""
import argparse,contextlib,csv,datetime,fcntl,hashlib,io,os,subprocess,sys,tempfile,uuid
from pathlib import Path
sys.dont_write_bytecode=True
IMAGE_SHA='dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372'
SLURM=Path('/cm/shared/apps/slurm/current/bin')
CAPS={'engineering':28800,'fit_A':208800,'fit_B':208800,'verification':7200,'buffer':7200}
PLAN={'prepare':('engineering',1,3600),'parallel':('engineering',8,1800),
      'engineering_check':('engineering',1,10800),
      'global_A':('fit_A',8,25200),'local_A':('fit_A',1,7200),
      'global_B':('fit_B',8,25200),'local_B':('fit_B',1,7200),
      'check_A':('verification',1,3600),'check_B':('verification',1,3600)}
FIELDS=['task','phase','mode','release','job_id','state','exit_code','cores','wall_seconds',
        'reserved_core_seconds','actual_core_seconds','intent_utc','dispatch_sha256','comment']
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
    s=Path(s);verify(s,s/'manifests/input_manifest.tsv')
    verify(Path(kv(s/'study.tsv')['project']),s/'manifests/external_manifest.tsv')
    if (s/'truth/manifest.tsv').exists():verify(s,s/'truth/manifest.tsv')
    extension=s/'authorizations/attempt_extension_01'
    if extension.exists():verify(extension,extension/'manifest.tsv')
    resume=s/'authorizations/local_A_resume_01'
    if resume.exists():verify(resume,resume/'manifest.tsv')
def payload(root):
    root=Path(root)
    if (root/'COMPLETE').read_text().strip()!=sha(root/'payload_manifest.tsv'):raise ValueError('Unsealed payload')
    verify(root,root/'payload_manifest.tsv')
def position(rows):
    out={}
    for phase,cap in CAPS.items():
        actual=sum(int(j['actual_core_seconds']) for j in rows if j['phase']==phase and closed(j))
        live=sum(int(j['reserved_core_seconds']) for j in rows if j['phase']==phase and not closed(j))
        out[phase]=dict(actual=actual,live=live,cap=cap,remaining=cap-actual-live)
        if actual+live>cap:raise ValueError('Phase cap exceeded: '+phase)
    if sum(x['actual']+x['live'] for x in out.values())>128*3600:raise ValueError('Whole package cap exceeded')
    return out
def extra_attempt_config(s,rows,mode):
    s=Path(s);extension=s/'authorizations/attempt_extension_01'
    if not extension.exists():raise ValueError('No explicit extra-attempt authorization')
    verify(extension,extension/'manifest.tsv');a=kv(extension/'approval.tsv')
    if a.get('authorized')!='true' or a.get('task')!='engineering_extra_01' or mode!='prepare':
        raise ValueError('Wrong extra-attempt scope')
    if (a.get('cores'),a.get('wall_seconds'),a.get('total_core_hour_cap'),a.get('engineering_core_hour_cap'),a.get('forward_cap'))!=('1','3600','128','8','64'):
        raise ValueError('Extra attempt cannot alter approved caps')
    original=[j for j in rows if j['mode']=='prepare']
    if len(original)!=2 or any(not closed(j) for j in original) or sum(int(j['actual_core_seconds']) for j in original)!=23:
        raise ValueError('Prior attempts/costs must remain attached')
    return 'buffer',1,3600
def local_resume_config(s,rows,mode):
    s=Path(s);d=s/'authorizations/local_A_resume_01'
    if not d.exists():raise ValueError('No explicit local continuation authorization')
    verify(d,d/'manifest.tsv');a=kv(d/'approval.tsv')
    if a.get('authorized')!='true' or a.get('task')!='local_A_resume_01' or mode!='local_A':
        raise ValueError('Wrong local continuation scope')
    prior=[j for j in rows if j['task']=='local_A']
    if len(prior)!=1 or prior[0].get('job_id')!='21902503' or prior[0]['state']!='FAILED' or prior[0]['actual_core_seconds']!='22':
        raise ValueError('Original local failure/cost must remain attached')
    calls=[r for p in (s/'fits/local_A').glob('evaluations_*.tsv') for r in read(p)]
    if len(calls)!=2:raise ValueError('Original local callbacks must remain attached')
    if kv(s/'fits/local_A/start_identity.tsv').get('difference')!='0':raise ValueError('Original local start not restored')
    wanted=dict(cores='1',wall_seconds='7178',remaining_callbacks='11181',remaining_soft_seconds='6878',
        original_stage_core_seconds='7200',original_stage_callbacks='11183',total_core_hour_cap='128',
        engineering_core_hour_cap='8',forward_cap='64')
    if any(a.get(k)!=v for k,v in wanted.items()):raise ValueError('Continuation cannot expand stage/call/package caps')
    return 'fit_A',1,7178
def reconcile(s):
    s=Path(s)
    with lock(s):
        rows=jobs(s)
        if rows:
            ids=','.join(j['job_id'] for j in rows if j['job_id'])
            if any(not j['job_id'] for j in rows):raise ValueError('Ambiguous submission: reconcile receipt, never resubmit')
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
        pos=position(rows);table(s/'accounting/phase_accounting.tsv',[dict(phase=k,**v) for k,v in pos.items()],['phase','actual','live','cap','remaining'])
        return rows
def submit(s,task,mode=None,release=None):
    s=Path(s);mode=mode or task
    if release is None:
        p=s/'tables/execution_release.tsv';release=kv(p)['release'] if p.exists() else 'r01'
    with lock(s):
        rows=jobs(s)
        if any(j['task']==task for j in rows):return
        if (s/'tables/final_accounting.tsv').exists():raise ValueError('Study closed')
        if task in PLAN:
            phase,cores,wall=PLAN[task]
        elif task=='engineering_retry':
            phase,cores,wall='buffer',1,7200
            if mode not in ('prepare','engineering_check'):raise ValueError('Only one bounded engineering recheck')
        elif task=='engineering_extra_01':phase,cores,wall=extra_attempt_config(s,rows,mode)
        elif task=='local_A_resume_01':phase,cores,wall=local_resume_config(s,rows,mode)
        else:raise ValueError('No declared allowance for '+task)
        if task.startswith(('global_','local_')):
            gate=kv(s/'tables/engineering_gate.tsv')
            if gate['pass']!='true' or gate['cost_fits']!='true':raise ValueError('Engineering/cost release held')
            verify(s,s/'tables/engineering_evidence.tsv')
        if task.startswith('local_'):
            parent='global_'+mode[-1];j=next(j for j in rows if j['task']==parent)
            if not success(j):raise ValueError('Local requires a successful own global parent')
            payload(s/'fits'/parent)
        p=position(rows)[phase]
        if p['remaining']<cores*wall:raise ValueError('No complete allocation fits the phase')
        if phase in ('engineering','buffer'):
            first=sum(v['actual']+v['live'] for k,v in position(rows).items() if k in ('engineering','buffer'))
            if first+cores*wall>8*3600:raise ValueError('First engineering window including fault costs exceeds8core-hours')
        rel=s/'releases'/release;verify(rel,rel/'manifest.tsv');verify_inputs(s)
        d=s/'dispatch'/task;d.mkdir(parents=True,exist_ok=False)
        spec=dict(task=task,mode=mode,phase=phase,cores=cores,wall_seconds=wall,
            release=release,release_manifest_sha256=sha(rel/'manifest.tsv'),input_manifest_sha256=sha(s/'manifests/input_manifest.tsv'))
        if task=='engineering_extra_01':spec['extra_authorization_manifest_sha256']=sha(s/'authorizations/attempt_extension_01/manifest.tsv')
        if task=='local_A_resume_01':spec['local_continuation_authorization_manifest_sha256']=sha(s/'authorizations/local_A_resume_01/manifest.tsv')
        typed(d/'spec.tsv',spec)
        comment='cd0:'+s.name+':'+task
        j=dict(task=task,phase=phase,mode=mode,release=release,job_id='',state='INTENT',exit_code='',cores=str(cores),wall_seconds=str(wall),reserved_core_seconds=str(cores*wall),actual_core_seconds='',intent_utc=utc(),dispatch_sha256=sha(d/'spec.tsv'),comment=comment)
        rows.append(j);position(rows);table(s/'accounting/jobs.tsv',rows,FIELDS)
        cmd=[str(SLURM/'sbatch'),'--parsable','--partition=red','--qos=xxlarge','--no-requeue','--nodes=1','--ntasks=1','--cpus-per-task='+str(cores),'--mem='+('64G' if cores==8 else '16G'),'--time={:02d}:{:02d}:{:02d}'.format(wall//3600,(wall%3600)//60,wall%60),'--job-name=cd0_'+task,'--comment='+comment,'--chdir='+str(s),'--output='+str(s/'logs'/(task+'_%j.log')),str(rel/'batch.sh'),str(s),task]
        text(d/'submission.txt','\n'.join(cmd)+'\n')
        r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        text(d/'receipt.txt',r.stdout);text(d/'submission_stderr.txt',r.stderr)
        if r.returncode or not r.stdout.split(';')[0].strip().isdigit():raise ValueError('Ambiguous sbatch receipt; fitting held')
        j.update(job_id=r.stdout.split(';')[0].strip(),state='SUBMITTED');table(s/'accounting/jobs.tsv',rows,FIELDS)
        print(task,j['job_id'],flush=True)

def reserve_forward(s,task,job,pid):
    s=Path(s)
    with lock(s,'forward_LOCK'):
        p=s/'accounting/engineering_forward_slots.tsv';rows=read(p) if p.exists() else []
        if len(rows)>=64:raise ValueError('64 full-forward slots consumed; no reset/retry')
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
    elif mode=='verify':
        verify_inputs(s);rel=Path(sys.argv[3]);verify(rel,rel/'manifest.tsv');print('Source, input and generated-object fingerprints verified')
    elif mode=='snapshot':print(position(reconcile(s)))
    elif mode=='submit':submit(s,sys.argv[3],sys.argv[4] if len(sys.argv)>4 else None)
    else:raise ValueError(mode)
