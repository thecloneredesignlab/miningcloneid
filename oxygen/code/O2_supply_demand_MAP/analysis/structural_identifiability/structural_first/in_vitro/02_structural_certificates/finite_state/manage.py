"""Single foreground submission owner; no background or legacy controller."""
from pathlib import Path
import shlex
import subprocess
import sys
from protocol import read,table,sha,kv,reservation,CAP

study=Path(sys.argv[1]).resolve();action=sys.argv[2]
remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+study.name
ssh=['ssh','-o','LogLevel=ERROR','-o','BatchMode=yes','4482173@red.moffitt.org']
def command(body):
    r=subprocess.run(ssh+['bash -lc '+shlex.quote(body)],capture_output=True,text=True)
    if r.returncode:raise RuntimeError(r.stderr or r.stdout)
    return r.stdout
def duration(value):
    if value.upper() in ('UNLIMITED','INFINITE'):return CAP
    days=0
    if '-' in value:d,value=value.split('-',1);days=int(d)
    parts=list(map(int,value.split(':')))
    if len(parts)==3:h,m,s=parts
    elif len(parts)==2:h,m,s=0,parts[0],parts[1]
    else:raise ValueError('Unrecognized site duration')
    return days*86400+h*3600+m*60+s
path=study/'accounting/jobs.tsv';jobs=read(path)
if action=='constraints':
    qos=command('sacctmgr -n -P show qos xxlarge format=Name,MaxWall')
    (study/'accounting/qos_constraints.txt').write_text(qos)
    rows=[x.split('|') for x in qos.splitlines() if x.startswith('xxlarge|')];assert len(rows)==1
    table(study/'accounting/site_limits.tsv',[dict(qos='xxlarge',max_seconds=duration(rows[0][1]))])
    (study/'accounting/partition_constraints.txt').write_text(command('scontrol show partition'))
    print('SITE_LIMIT_RECORDED',rows[0][1])
elif action=='submit':
    assert not (study/'COMPLETE').exists() and not (study/'accounting/stop_record.tsv').exists(),'Closed or stopped package cannot dispatch additional work'
    task,release,runner,program=sys.argv[3:7]
    assert task not in {j['task'] for j in jobs},'No implicit task retry or duplicate dispatch'
    assert all(j['job_id'].isdigit() for j in jobs),'Uncertain prior submission must be reconciled'
    if task=='B0':assert all(r['passed']=='True' for r in read(study/'reports/static_B0_tests.tsv'))
    else:
        assert all(r['passed']=='True' for r in read(study/'reports'/('static_'+release+'_tests.tsv')))
        assert read(study/'verification/B0/semantic_gate.tsv')[0]['status']=='PASS'
        assert read(study/'verification/B0/policy_gate.tsv')[0]['status']=='PASS'
    if task in ('B2','B2_review','B6') or task.startswith(('B01_','B02_','B03_','B04_','B05_','B06_')):
        assert read(study/'verification/B1/validation_gate.tsv')[0]['status']=='PASS'
    if task=='B6':
        accepted=kv(study/'accounting/operator_acceptance.tsv')['phase']
        assert read(study/'verification'/accepted/'operator_gate.tsv')[0]['status']=='PASS'
        for j in jobs:
            if j['task']=='B2':continue  # Failed original and saved-result review remain separate.
            assert j['state']=='COMPLETED' and j['exit_code']=='0:0'
            assert read(study/'verification'/j['task']/'policy_gate.tsv')[0]['status']=='PASS'
    if task=='B2_review':
        assert (study/'accounting/technical_retry.tsv').exists()
        assert len(read(study/'verification/B2_native/kernel_calls.tsv'))==9
        assert runner=='job_python.sh' and program=='compare_b2.py','Native kernels cannot be repeated'
    if task.startswith(('B01_','B02_','B03_','B04_','B05_','B06_')):
        case,witness=task.split('_');assert witness in ('W1','W2','W3')
        accepted=kv(study/'accounting/operator_acceptance.tsv')['phase']
        assert read(study/'verification'/accepted/'operator_gate.tsv')[0]['status']=='PASS'
        assert read(study/'verification'/accepted/'policy_gate.tsv')[0]['status']=='PASS'
        if case in ('B05','B06'):
            prior=('B03_' if case=='B05' else 'B02_')+witness
            j=next(j for j in jobs if j['task']==prior)
            assert j['state']=='COMPLETED' and j['exit_code']=='0:0'
            assert read(study/'verification'/prior/'policy_gate.tsv')[0]['status']=='PASS'
            assert read(study/'verification'/prior/'application_result.tsv')[0]['status']=='PASS'
            assert read(study/'verification'/prior/'case_result.tsv')[0]['status']!='CERTIFIED_FULL_LOCAL'
    manifest=study/('source_manifest_'+release+'.tsv')
    assert manifest.exists()
    for r in read(manifest):assert sha(study/r['path'])==r['sha256'],r['path']
    settled=sum(int(j['actual_core_seconds']) for j in jobs)
    live=sum(int(j['reserved_core_seconds']) for j in jobs)
    used_retry=int(kv(study/'accounting/technical_retry.tsv').get('used',0)) if (study/'accounting/technical_retry.tsv').exists() else 0
    seconds=reservation(settled,live,int(read(study/'accounting/site_limits.tsv')[0]['max_seconds']),used_retry)
    row=dict(job_id='submission_pending',task=task,release=release,state='RESERVED',exit_code='',allocated_cpus='',elapsed_seconds='',actual_core_seconds=0,reserved_core_seconds=seconds)
    jobs.append(row);table(path,jobs)
    args=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1','--cpus-per-task=1',
          '--hint=nomultithread','--mem=16G','--time='+str(seconds//60),'--no-requeue','--job-name=stageB_'+task,
          '--output='+remote+'/logs/'+task+'_%j.txt',remote+'/releases/'+release+'/'+runner,remote,release,task,program]
    result=command(' '.join(shlex.quote(a) for a in args)).strip().split(';')[0];assert result.isdigit(),result
    row.update(job_id=result,state='SUBMITTED');table(path,jobs)
    print('SUBMITTED',task,result,'xxlarge; budget/site-derived core-second reservation',seconds)
elif action=='refresh':
    assert jobs and all(j['job_id'].isdigit() for j in jobs)
    ids=','.join(j['job_id'] for j in jobs)
    raw=command('sacct -n -P -j '+ids+' --format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,QOS,AllocTRES')
    (study/'accounting/slurm_snapshot.txt').write_text(raw)
    rec={l.split('|')[0]:l.split('|') for l in raw.splitlines() if '|' in l}
    terminal=('COMPLETED','FAILED','TIMEOUT','CANCELLED','NODE_FAIL','OUT_OF_MEMORY','PREEMPTED','BOOT_FAIL','DEADLINE')
    for j in jobs:
        if j['job_id'] not in rec:continue
        v=rec[j['job_id']];j.update(state=v[1],exit_code=v[2],allocated_cpus=v[3],elapsed_seconds=v[4])
        if v[1].split()[0] in terminal:
            matches=[r for k,r in rec.items() if k==j['job_id'] or k.startswith(j['job_id']+'.')]
            j['actual_core_seconds']=max(int(r[3] or 0)*int(r[4] or 0) for r in matches);j['reserved_core_seconds']=0
    table(path,jobs)
    (study/'accounting/owned_queue.txt').write_text(command('squeue -h -j '+ids+' -o '+shlex.quote('%i|%T')))
    for folder in ('verification','logs'):
        subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','-e','ssh -o LogLevel=ERROR -o BatchMode=yes',
                        '4482173@red.moffitt.org:'+remote+'/'+folder+'/',str(study/folder)+'/'],check=True)
    print(jobs)
else:raise ValueError('Unknown foreground action')
