"""Single foreground owner of only the new Step5 allocations and resource ledger."""
from pathlib import Path
import shlex
import subprocess
import sys
from common import *

study=Path(sys.argv[1]).resolve();action=sys.argv[2]
remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+study.name
ssh=['ssh','-o','LogLevel=ERROR','-o','BatchMode=yes','4482173@red.moffitt.org']
def command(body):
    r=subprocess.run(ssh+['bash -lc '+shlex.quote(body)],capture_output=True,text=True)
    if r.returncode:raise RuntimeError(r.stderr or r.stdout)
    return r.stdout
path=study/'accounting/jobs.tsv';jobs=read(path)
if action=='constraints':
    qos=command('sacctmgr -n -P show qos xxlarge format=Name,MaxWall')
    (study/'accounting/qos_constraints.txt').write_text(qos)
    value=next(line.split('|')[1] for line in qos.splitlines() if line.startswith('xxlarge|'))
    h,m,s=map(int,value.split(':'));seconds=3600*h+60*m+s
    table(study/'accounting/site_limits.tsv',[dict(qos='xxlarge',max_seconds=seconds)])
    (study/'accounting/original_terminal_snapshot.txt').write_text(command('sacct -n -P -j 21971421,21971431,21971432,21971436 --format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw'))
    (study/'accounting/partition_snapshot.txt').write_text(command('scontrol show partition'))
    print('CURRENT_XXLARGE_SITE_LIMIT',value)
elif action=='submit':
    task=sys.argv[3];assert task in ('M01','O01','R01','R02','O02','R03','G01')
    assert not (study/'CLOSED').exists() and not (study/'accounting/stop_record.tsv').exists()
    assert all(j['job_id'].isdigit() for j in jobs),'Reconcile uncertain submission'
    assert task not in {j['task'] for j in jobs},'No duplicate dispatch or implicit retry'
    assert all(r['passed']=='True' for r in read(study/'reports/static_tests.tsv'))
    prereq={'M01':[], 'O01':['M01'], 'R01':['M01'], 'R02':['R01'], 'O02':['R01'], 'R03':['R02','O02'], 'G01':['R03','O02','O01']}
    for prior in prereq[task]:
        j=next(j for j in jobs if j['task']==prior)
        assert j['state']=='COMPLETED' and j['exit_code']=='0:0'
        assert read(study/'verification'/prior/'application_result.tsv')[0]['status']=='PASS'
        assert read(study/'verification'/prior/'policy_gate.tsv')[0]['status']=='PASS'
        assert read(study/'verification'/prior/'gate.tsv')[0]['status'].startswith('PASS')
        for r in read(study/'verification'/prior/'phase_manifest.tsv'):
            pp=study/'verification'/prior/r['path'];assert sha(pp)==r['sha256']
    verify(study)
    if task in ('O02','G01'):
        assert all(r['passed']=='True' for r in read(study/'reports/additional_validation_tests.tsv'))
        for r in read(study/'validation_source_manifest.tsv'):
            pp=study/r['path'];assert sha(pp)==r['sha256']
    settled=sum(int(j['actual_core_seconds']) for j in jobs);live=sum(int(j['reserved_core_seconds']) for j in jobs)
    seconds=reservation(settled,live,int(read(study/'accounting/site_limits.tsv')[0]['max_seconds']),int(kv(study/'accounting/retry_ledger.tsv')['used']))
    row=dict(job_id='submission_pending',task=task,state='RESERVED',exit_code='',allocated_cpus='',elapsed_seconds='',actual_core_seconds=0,reserved_core_seconds=seconds)
    jobs.append(row);table(path,jobs)
    args=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1','--cpus-per-task=1','--hint=nomultithread','--mem=16G',
        '--time='+str(seconds//60),'--no-requeue','--job-name=step5_' +task,'--output='+remote+'/logs/'+task+'_%j.txt',remote+('/validation_release/runner.sh' if task in ('O02','G01') else '/release/runner.sh'),remote,task]
    result=command(' '.join(shlex.quote(a) for a in args)).strip().split(';')[0];assert result.isdigit(),result
    row.update(job_id=result,state='SUBMITTED');table(path,jobs)
    print('SUBMITTED',task,result,'xxlarge; site/budget-derived reservation',seconds)
elif action=='refresh':
    assert jobs and all(j['job_id'].isdigit() for j in jobs)
    ids=','.join(j['job_id'] for j in jobs)
    raw=command('sacct -n -P -j '+ids+' --format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,QOS,AllocTRES')
    (study/'accounting/slurm_snapshot.txt').write_text(raw)
    rec={line.split('|')[0]:line.split('|') for line in raw.splitlines() if '|' in line}
    terminal=('COMPLETED','FAILED','TIMEOUT','CANCELLED','NODE_FAIL','OUT_OF_MEMORY','PREEMPTED','BOOT_FAIL','DEADLINE')
    for j in jobs:
        if j['job_id'] not in rec:continue
        r=rec[j['job_id']];j.update(state=r[1],exit_code=r[2],allocated_cpus=r[3],elapsed_seconds=r[4])
        if r[1].split()[0] in terminal:
            matches=[v for k,v in rec.items() if k==j['job_id'] or k.startswith(j['job_id']+'.')]
            j['actual_core_seconds']=max(int(v[3] or 0)*int(v[4] or 0) for v in matches);j['reserved_core_seconds']=0
    table(path,jobs)
    queue=command('squeue -h -u 4482173 -o '+shlex.quote('%i|%T'))
    owned={j['job_id'] for j in jobs}
    (study/'accounting/owned_queue.txt').write_text(''.join(line+'\n' for line in queue.splitlines() if line.split('|')[0] in owned))
    for folder in ('verification','logs'):
        subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','-e','ssh -o LogLevel=ERROR -o BatchMode=yes',
            '4482173@red.moffitt.org:'+remote+'/'+folder+'/',str(study/folder)+'/'],check=True)
    print(jobs)
else:raise ValueError('Unknown action')
