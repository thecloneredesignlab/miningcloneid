"""Budget-derived resource admission, not a discretionary case time cutoff."""
from pathlib import Path
import shlex
import subprocess
import sys
from protocol import read,table

study=Path(sys.argv[1]).resolve();action=sys.argv[2]
remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+study.name
ssh=['ssh','-o','LogLevel=ERROR','-o','BatchMode=yes','4482173@red.moffitt.org']
def command(text):
    r=subprocess.run(ssh+['bash -lc '+shlex.quote(text)],capture_output=True,text=True);assert r.returncode==0,r.stderr;return r.stdout
path=study/'accounting/jobs.tsv';jobs=read(path)
if action=='submit':
    assert not jobs and all(r['passed']=='True' for r in read(study/'reports/static_tests.tsv'))
    auth={r['key']:r['value'] for r in read(study/'inputs/authorization.tsv')}
    minutes=int(auth['budget_derived_allocation_minutes']);reserve=minutes*60
    assert reserve+int(auth['previous_charged_core_seconds'])<=int(auth['cap_total_core_seconds'])
    row=dict(job_id='submission_pending',state='RESERVED',exit_code='',allocated_cpus='',elapsed_seconds='',actual_core_seconds=0,reserved_core_seconds=reserve)
    table(path,[row])
    args=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1','--cpus-per-task=1','--hint=nomultithread','--mem=8G','--time='+str(minutes),'--no-requeue','--job-name=C04_exact_certificate','--output='+remote+'/logs/slurm_%j.txt',remote+'/releases/r01/job.sh',remote]
    result=command(' '.join(shlex.quote(a) for a in args)).strip().split(';')[0];assert result.isdigit()
    row.update(job_id=result,state='SUBMITTED');table(path,[row])
    print('SUBMITTED',result,'xxlarge; no case timeout; whole remaining core-hour budget reserved, rounded down to scheduler minutes')
elif action=='refresh':
    assert len(jobs)==1;row=jobs[0]
    text=command('sacct -n -P -j '+row['job_id']+' --format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,QOS,AllocTRES')
    (study/'accounting/slurm_snapshot.txt').write_text(text)
    records={line.split('|')[0]:line.split('|') for line in text.splitlines() if '|' in line}
    if row['job_id'] in records:
        v=records[row['job_id']];row.update(state=v[1],exit_code=v[2],allocated_cpus=v[3],elapsed_seconds=v[4])
        if v[1] in ('COMPLETED','FAILED','TIMEOUT','NODE_FAIL','OUT_OF_MEMORY','CANCELLED'):
            row['actual_core_seconds']=max(int(b[3] or 0)*int(b[4] or 0) for a,b in records.items() if a==row['job_id'] or a.startswith(row['job_id']+'.'));row['reserved_core_seconds']=0
    table(path,[row])
    for folder in ('verification','logs'):
        r=subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','-e','ssh -o LogLevel=ERROR -o BatchMode=yes','4482173@red.moffitt.org:'+remote+'/'+folder+'/',str(study/folder)+'/'],capture_output=True,text=True);assert r.returncode==0,r.stderr
    print(row)
else:raise ValueError('Unknown metadata action')
