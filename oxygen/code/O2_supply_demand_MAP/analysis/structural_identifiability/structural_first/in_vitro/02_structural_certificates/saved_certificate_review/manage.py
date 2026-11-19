"""Foreground metadata admission/accounting; never performs scientific verification."""
from pathlib import Path
import shlex
import subprocess
import sys
from base_contract import read,table,sha,policy_findings

study=Path(sys.argv[1]).resolve();action=sys.argv[2]
remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+study.name
ssh=['ssh','-o','LogLevel=ERROR','-o','BatchMode=yes','-o','ConnectTimeout=15','4482173@red.moffitt.org']
path=study/'accounting/jobs.tsv';jobs=read(path)
fields=['task','job_id','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds','retry']
def cmd(command):
    r=subprocess.run(ssh+['bash -lc '+shlex.quote(command)],capture_output=True,text=True)
    assert r.returncode==0,r.stderr
    return r.stdout
def save():table(path,jobs,fields)
def gate():
    r=next((r for r in reversed(jobs) if r['task']=='review_C03'),None)
    assert r and r['state']=='COMPLETED' and r['exit_code']=='0:0','Saved-review scheduler gate'
    folder=study/'review/C03';assert not policy_findings(folder)
    receipt=read(folder/'result.tsv')[0]
    assert receipt['status']=='PASS' and receipt['run_id']==study.name and receipt['new_native_C03_calls']=='0'
    assert receipt['input_sha256']==sha(study/'saved/C03/native_app/models/case_model.py')
    assert all(r['solver_calls']=='0' for r in read(folder/'future_interface_gate.tsv'))
    assert read(folder/'policy_gate.tsv')[0]['status']=='PASS'
    return True
if action=='submit':
    task=sys.argv[3];assert task in ('review_C03','native_C04','native_C05')
    assert len(jobs)<4 and not any(r['task']==task for r in jobs)
    if task!='review_C03':gate()
    else:assert not jobs and all(r['passed']=='True' for r in read(study/'reports/static_tests.tsv'))
    assert sum(int(r['actual_core_seconds'])+int(r['reserved_core_seconds']) for r in jobs)+1800<=7200
    row=dict(task=task,job_id='submission_pending',state='RESERVED',exit_code='',allocated_cpus='',elapsed_seconds='',actual_core_seconds=0,reserved_core_seconds=1800,retry=False)
    jobs.append(row);save()
    args=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1','--cpus-per-task=1','--hint=nomultithread','--mem=8G','--time=00:30:00','--no-requeue','--job-name=savedctrl_'+task,'--output='+remote+'/logs/slurm_%j.txt',remote+'/releases/r01/job.sh',remote,task]
    result=cmd(' '.join(shlex.quote(s) for s in args)).strip().split(';')[0];assert result.isdigit()
    row.update(job_id=result,state='SUBMITTED');save();print('SUBMITTED',task,result,'xxlarge; reserve 1800 core-seconds')
elif action=='refresh':
    ids=','.join(r['job_id'] for r in jobs);assert ids
    snapshot=cmd('sacct -n -P -j '+ids+' --format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,QOS')
    (study/'accounting/slurm_snapshot.txt').write_text(snapshot)
    records={line.split('|')[0]:line.split('|') for line in snapshot.splitlines() if '|' in line}
    for r in jobs:
        if r['job_id'] not in records:continue
        v=records[r['job_id']];r.update(state=v[1],exit_code=v[2],allocated_cpus=v[3],elapsed_seconds=v[4])
        if v[1] in ('COMPLETED','FAILED','TIMEOUT','NODE_FAIL','OUT_OF_MEMORY','CANCELLED'):
            r['actual_core_seconds']=max(int(b[3] or 0)*int(b[4] or 0) for a,b in records.items() if a==r['job_id'] or a.startswith(r['job_id']+'.'))
            r['reserved_core_seconds']=0
    save()
    for folder in ('review','controls','validation_retry','logs'):
        r=subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','-e','ssh -o LogLevel=ERROR -o BatchMode=yes','4482173@red.moffitt.org:'+remote+'/'+folder+'/',str(study/folder)+'/'],text=True,capture_output=True)
        assert r.returncode==0,r.stderr
    print('\n'.join(r['task']+' '+r['job_id']+' '+r['state']+' '+r['exit_code'] for r in jobs))
    print('Conservative actual core-seconds:',sum(int(r['actual_core_seconds']) for r in jobs))
elif action=='gate':assert gate();print('SAVED_C03_GATE_PASSED_WITH_ZERO_NATIVE_CALLS')
else:raise ValueError('Unknown metadata action')
