#!/usr/bin/env python3
"""Recover reporting after the canonical/global report module-name collision.

Preserves all original model/analysis jobs, cached blocks, and frozen plan.
"""
import argparse
import csv
from pathlib import Path
import shutil
import subprocess

from common import load,save,sha,now,locked
from submit import FILES,RELATIVE,output


def recover(repo,root,revision):
    if (root/'reporting_recovery_receipt.json').exists():
        print((root/'reporting_recovery_receipt.json').read_text());return
    p=load(root/'plan.json');original=load(root/'launch_receipt.json')
    jobs={r['stage']:r['job_id'] for r in original['jobs']}
    if p['scope']!='global' or p['revision']!=original['revision']:
        raise ValueError('Original batch identity differs')
    overlay=root/'slurm'/('reporting_overlay_'+revision[:9]);overlay.mkdir(exist_ok=True)
    files=list(FILES)+['reporting_entry.py','run_reporting_overlay.sbatch','recover_reporting.py']
    for name in files:
        data=subprocess.check_output(['git','-C',str(repo),'show',revision+':'+RELATIVE+name])
        target=overlay/name
        if target.exists() and target.read_bytes()!=data:raise ValueError('Overlay differs')
        if not target.exists():target.write_bytes(data);target.chmod(0o444)
    # The scientific estimator/model/report sources are the original certified files.
    for name in ['report.py','plot_helpers.py','common.py','evaluate_efast.R']:
        if sha(overlay/name)!=p['controller_sha256'][name]:
            raise ValueError('Reporting repair unexpectedly changes mathematics: '+name)
    save(overlay/'manifest.json',dict(revision=revision,base_plan_sha256=sha(root/'plan.json'),scope='reporting_only',
        reason='qualified frozen-report import avoids canonical report.py name collision',
        files={name:sha(overlay/name) for name in files}))
    ledger=root/'slurm/reporting_recovery_jobs.json'
    records=load(ledger) if ledger.exists() else []
    new={r['stage']:r['job_id'] for r in records}
    stages=[('sobol_recovered','sobol_report',0,[jobs['prepare']])]
    previous=None
    for n in p['resolutions']:
        name='report_'+str(n)+'_recovered'
        stages.append((name,'report',n,[jobs['analyze_'+str(n)]]+([previous] if previous else [])))
        previous=name
    stages.append(('final_recovered','final_report',0,['sobol_recovered']+['report_'+str(n)+'_recovered' for n in p['resolutions']]))
    for stage,role,n,deps in stages:
        if stage in new:continue
        dependencies=[new.get(x,x) for x in deps]
        name='gef_fix_'+stage
        args=['sbatch','--parsable','--hold','--qos=xxlarge','--time=12:00:00','--ntasks=1','--cpus-per-task=1',
              '--mem=16G','--job-name='+name,'--chdir='+str(root),'--export=ALL',
              '--dependency=afterok:'+':'.join(dependencies),
              '--output='+str(root/'slurm/logs'/(stage+'_%A_%a.out')),
              '--error='+str(root/'slurm/logs'/(stage+'_%A_%a.err')),
              str(overlay/'run_reporting_overlay.sbatch'),str(root),role,str(n),str(overlay)]
        found=set(output(['squeue','-h','-u','4482173','-n',name,'-o','%A']).splitlines())
        if len(found)>1:raise ValueError('Duplicate reporting recovery job')
        job=found.pop() if found else output(args).split(';')[0]
        if not job.isdigit():raise ValueError('Invalid reporting receipt')
        if 'WorkDir='+str(root) not in output(['scontrol','show','job','-o',job]):raise ValueError('Unrelated recovery job')
        records.append(dict(stage=stage,role=role,N=n,job_id=job,dependencies=dependencies,sbatch_args=args,submitted_utc=now()))
        save(ledger,records);new[stage]=job
        print('REPORTING_RECOVERY',stage,job,flush=True)
    monitor=root/'time_extension_reporting'
    for d in ['logs','bin','source_scripts']:(monitor/d).mkdir(parents=True,exist_ok=True)
    for name,target in [('time_extension.py','monitor.py'),('run_time_extension.sbatch','run_monitor.sbatch'),('s_add_job_hours_proxy.sh','bin/s-add-job-hours')]:
        shutil.copyfile(str(overlay/name),str(monitor/target))
    (monitor/'bin/s-add-job-hours').chmod(0o755)
    for name in ['submit_hourly_add_hours.sh','add_hours_to_running_array_jobs.sh']:
        shutil.copyfile(str(root/'time_extension/source_scripts'/name),str(monitor/'source_scripts'/name))
    save(monitor/'targets.json',dict(target_jobs=[r['job_id'] for r in records],base_plan_sha256=sha(root/'plan.json'),overlay_sha256=sha(overlay/'manifest.json')))
    subprocess.run(['python3',str(monitor/'monitor.py'),'submit'],check=True)
    with (monitor/'submission_jobs.tsv').open() as handle:monitors=list(csv.DictReader(handle,delimiter='\t'))
    if {r['target_job_id'] for r in monitors}!={r['job_id'] for r in records}:raise ValueError('Recovery monitoring incomplete')
    cancelled=[]
    for stage in ['report_'+str(n) for n in p['resolutions']]+['final']:
        job=jobs[stage]
        state=output(['squeue','-h','-j',job,'-o','%T'])
        if state:
            if set(state.splitlines())!={'PENDING'}:raise ValueError('Refusing to cancel an active original reporting job')
            subprocess.run(['scancel',job],check=True);cancelled.append(job)
    for row in records:subprocess.run(['scontrol','release',row['job_id']],check=True)
    save(root/'reporting_recovery_receipt.json',dict(status='submitted',revision=revision,base_revision=p['revision'],
        overlay=str(overlay),overlay_manifest_sha256=sha(overlay/'manifest.json'),reason='report module namespace collision',
        original_failed_sobol_report=jobs['sobol'],cancelled_pending_reporting_jobs=cancelled,
        preserved_model_jobs=[jobs['model_'+str(n)] for n in p['resolutions']],
        preserved_analysis_jobs=[jobs['analyze_'+str(n)] for n in p['resolutions']],
        model_jobs_resubmitted=0,jobs=records,monitors=monitors,submitted_utc=now()))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=Path('/share/lab_crd/taoli/Project/soft_couping_eFAST'))
    parser.add_argument('--revision',required=True)
    args=parser.parse_args();root=args.repo/'oxygen/results/eFAST/global_finalize_20261112'
    with locked(root/'reporting_recovery.lock'):recover(args.repo,root,args.revision)
