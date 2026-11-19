"""Foreground metadata admission/collection with aggregate resource accounting."""
from pathlib import Path
import shlex
import subprocess
import sys
from protocol import read,table,allocation_seconds,CAP

study=Path(sys.argv[1]).resolve();action=sys.argv[2]
remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+study.name
ssh=['ssh','-o','LogLevel=ERROR','-o','BatchMode=yes','4482173@red.moffitt.org']
def command(body):
    result=subprocess.run(ssh+['bash -lc '+shlex.quote(body)],capture_output=True,text=True)
    if result.returncode:raise RuntimeError(result.stderr or result.stdout)
    return result.stdout

def duration(value):
    if value.upper() in ('UNLIMITED','INFINITE'):return CAP
    days=0
    if '-' in value:
        d,value=value.split('-',1);days=int(d)
    parts=list(map(int,value.split(':')))
    if len(parts)==3:h,m,s=parts
    elif len(parts)==2:h,m,s=0,*parts
    else:raise ValueError('Unrecognized scheduler duration')
    return days*86400+h*3600+m*60+s

path=study/'accounting/jobs.tsv';jobs=read(path)
if action=='constraints':
    qos=command('sacctmgr -n -P show qos xxlarge format=Name,MaxWall')
    (study/'accounting/qos_constraints.txt').write_text(qos)
    selected=[line.split('|') for line in qos.splitlines() if line.startswith('xxlarge|')]
    assert len(selected)==1 and selected[0][1]
    partition=command('scontrol show partition')
    (study/'accounting/partition_constraints.txt').write_text(partition)
    table(study/'accounting/site_limits.tsv',[dict(qos='xxlarge',max_seconds=duration(selected[0][1]),scope='QoS maximum; partition/account limits may be stricter at scheduler admission')])
    print('QOS_CONSTRAINTS_RECORDED',selected[0][1])
elif action=='submit':
    assert not jobs,'Submission uncertainty or prior attempt must be reconciled before any retry'
    assert all(r['passed']=='True' for r in read(study/'reports/static_tests.tsv'))
    site=int(read(study/'accounting/site_limits.tsv')[0]['max_seconds'])
    auth={r['key']:r['value'] for r in read(study/'inputs/authorization.tsv')}
    already_used_retry=int(auth['technical_retries_used'])
    seconds=allocation_seconds(0,retry=bool(already_used_retry),site_seconds=site)
    row=dict(job_id='submission_pending',attempt=0,state='RESERVED',exit_code='',allocated_cpus='',elapsed_seconds='',actual_core_seconds=0,reserved_core_seconds=seconds)
    table(path,[row])
    args=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1','--cpus-per-task=1',
          '--hint=nomultithread','--mem=8G','--time='+str(seconds//60),'--no-requeue',
          '--job-name=initial_observer_A','--output='+remote+'/logs/slurm_%j.txt',remote+'/releases/r01/job.sh',remote]
    output=command(' '.join(shlex.quote(a) for a in args)).strip()
    result=output.split(';')[0];assert result.isdigit(),output
    row.update(job_id=result,state='SUBMITTED');table(path,[row])
    table(study/'accounting/admission.tsv',[dict(cap_core_seconds=CAP,settled=0,live_reservation=seconds,retry_reserve=0 if already_used_retry else 3600,settlement_reserve=3600,technical_retries_used=already_used_retry,discretionary_timeout='none')])
    print('SUBMITTED',result,'xxlarge; budget-derived reservation',seconds,'core-seconds; no per-case timer')
elif action=='refresh':
    assert jobs and all(r['job_id'].isdigit() for r in jobs)
    ids=','.join(r['job_id'] for r in jobs)
    snapshot=command('sacct -n -P -j '+ids+' --format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,QOS,AllocTRES')
    (study/'accounting/slurm_snapshot.txt').write_text(snapshot)
    records={line.split('|')[0]:line.split('|') for line in snapshot.splitlines() if '|' in line}
    terminal=('COMPLETED','FAILED','TIMEOUT','NODE_FAIL','OUT_OF_MEMORY','CANCELLED','PREEMPTED','BOOT_FAIL','DEADLINE')
    for row in jobs:
        if row['job_id'] not in records:continue
        v=records[row['job_id']];row.update(state=v[1],exit_code=v[2],allocated_cpus=v[3],elapsed_seconds=v[4])
        if v[1].split()[0] in terminal:
            matches=[b for a,b in records.items() if a==row['job_id'] or a.startswith(row['job_id']+'.')]
            row['actual_core_seconds']=max(int(b[3] or 0)*int(b[4] or 0) for b in matches)
            row['reserved_core_seconds']=0
    table(path,jobs)
    queue=command('squeue -h -j '+ids+' -o '+shlex.quote('%i|%T'))
    (study/'accounting/owned_queue.txt').write_text(queue)
    for folder in ('verification','logs'):
        subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','-e','ssh -o LogLevel=ERROR -o BatchMode=yes',
                        '4482173@red.moffitt.org:'+remote+'/'+folder+'/',str(study/folder)+'/'],check=True)
    print(jobs)
else:raise ValueError('Unknown foreground management action')
