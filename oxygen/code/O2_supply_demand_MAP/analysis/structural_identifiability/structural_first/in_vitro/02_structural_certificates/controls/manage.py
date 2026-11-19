"""Foreground admission/accounting helper for the finite rerun. No scientific calls."""
from pathlib import Path
import shlex
import subprocess
import sys
from contract import CASES,read,table,sha,receipt_valid,policy_findings

study=Path(sys.argv[1]).resolve();action=sys.argv[2]
remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+study.name
ssh=['ssh','-o','LogLevel=ERROR','-o','BatchMode=yes','-o','ConnectTimeout=15','4482173@red.moffitt.org']
jobs_path=study/'accounting/jobs.tsv';jobs=read(jobs_path)
fields=['task','job_id','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds','retry']

def command(text):
    result=subprocess.run(ssh+['bash -lc '+shlex.quote(text)],text=True,capture_output=True)
    if result.returncode:raise RuntimeError(result.stderr)
    return result.stdout

def save():table(jobs_path,jobs,fields)

def gate(task):
    row=next((r for r in reversed(jobs) if r['task']==task),None)
    assert row and row['state']=='COMPLETED' and row['exit_code']=='0:0','Scheduler prerequisite not passed: '+task
    cases=CASES if task=='preflight' else (task,)
    for case in cases:
        folder=study/(('preflight_retry' if str(row['retry'])=='True' else 'preflight') if task=='preflight' else 'controls')/case
        assert not policy_findings(folder),'Policy gate failed'
        result=read(folder/'result.tsv');assert len(result)==1
        assert receipt_valid(result[0],case,sha(study/'fixtures'/case/'model.py'),study.name),'Invalid/stale result receipt'
        assert read(folder/'policy_gate.tsv')[0]['status']=='PASS'
        assert (folder/'application_status.tsv').read_text().strip()=='application_exit_code\t0'
    return True

if action in ('submit','retry_preflight'):
    task=sys.argv[3];assert task in CASES+('preflight',)
    retry=action=='retry_preflight'
    if retry:
        assert task=='preflight' and len(jobs)==1 and jobs[0]['state']=='FAILED'
        assert not any(str(r['retry'])=='True' for r in jobs)
        assert (study/'reports/technical_retry_registration.md').exists()
    else:assert not any(r['task']==task for r in jobs),'Duplicate task requires separately documented retry handling'
    if task=='C02':gate('preflight')
    elif task=='C03':gate('C02')
    elif task in ('C04','C05'):gate('C03')
    else:assert (retry or not jobs) and all(r['passed']=='True' for r in read(study/'reports/static_input_checks.tsv'))
    assert sum(int(r['actual_core_seconds'])+int(r['reserved_core_seconds']) for r in jobs)+1800<=10800
    row=dict(task=task,job_id='submission_pending',state='RESERVED',exit_code='',allocated_cpus='',elapsed_seconds='',actual_core_seconds=0,reserved_core_seconds=1800,retry=retry)
    jobs.append(row);save()
    release='r02' if (study/'releases/r02/job.sh').exists() else 'r01'
    args=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1','--cpus-per-task=1','--hint=nomultithread','--mem=8G','--time=00:30:00','--no-requeue','--job-name=fourctrl_'+task,'--output='+remote+'/logs/slurm_%j.txt',remote+'/releases/'+release+'/job.sh',remote,task]
    reply=command(' '.join(shlex.quote(a) for a in args)).strip()
    jid=reply.split(';')[0];assert jid.isdigit(),reply
    row.update(job_id=jid,state='SUBMITTED');save()
    print('SUBMITTED',task,jid,'xxlarge; reserved 1800 core-seconds; no concurrency throttle')
elif action=='refresh':
    ids=[r['job_id'] for r in jobs];assert ids
    snapshot=command('sacct -n -P -j '+','.join(ids)+' --format=JobID,State,ExitCode,AllocCPUS,ElapsedRaw,ReqTRES,AllocTRES')
    (study/'accounting/slurm_snapshot.txt').write_text(snapshot)
    records={line.split('|')[0]:line.split('|') for line in snapshot.splitlines() if '|' in line}
    for row in jobs:
        data=records.get(row['job_id'])
        if data:
            row.update(state=data[1],exit_code=data[2],allocated_cpus=data[3],elapsed_seconds=data[4])
            if data[1] in ('COMPLETED','FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED','NODE_FAIL'):
                row['actual_core_seconds']=max(int(v[3] or 0)*int(v[4] or 0) for k,v in records.items() if k==row['job_id'] or k.startswith(row['job_id']+'.'))
                row['reserved_core_seconds']=0
    save()
    for folder in ('preflight','preflight_retry','controls','logs'):
        if folder=='preflight_retry' and not (study/'releases/r02/job.sh').exists():continue
        result=subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','-e','ssh -o LogLevel=ERROR -o BatchMode=yes','4482173@red.moffitt.org:'+remote+'/'+folder+'/',str(study/folder)+'/'],capture_output=True,text=True)
        assert result.returncode==0,result.stderr
    print('\n'.join(r['task']+' '+r['job_id']+' '+r['state']+' '+r['exit_code'] for r in jobs))
    print('Conservative actual core-seconds:',sum(int(r['actual_core_seconds']) for r in jobs))
elif action=='gate':
    assert gate(sys.argv[3]);print('ACCEPTANCE_GATE_PASSED',sys.argv[3])
else:raise ValueError('Unknown foreground action')
