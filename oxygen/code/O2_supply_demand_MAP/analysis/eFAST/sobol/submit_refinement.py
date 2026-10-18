#!/usr/bin/env python3
"""Idempotent, monitored submission of the approved 2026-10-18 Sobol tails.

Runs with RED's host Python 3.6. Numerical preparation and all computation are
Slurm jobs in the original SIF; the login node only freezes/submits manifests.
"""
import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from refinement_run import load, save, sha, now, task_layout

FILES=('refinement_run.py','submit_refinement.py','run_refinement.sbatch',
       'test_refinement.py','retry_container_startup.sh','time_extension.py',
       'run_time_extension.sbatch','s_add_job_hours_proxy.sh')


def output(args):
    return subprocess.check_output(args,universal_newlines=True).strip()


def stage_specs(selection):
    models, analyses=task_layout(selection)
    stages=[('prepare','prepare','',0,0,'4G',[]),
            ('audit_global','audit','global',0,1,'4G',['prepare']),
            ('audit_neighborhood','audit','neighborhood10pct',0,500,'4G',['prepare']),
            ('smoke','smoke','global',8192,0,'4G',['audit_global'])]
    for n in (8192,16384,32768):
        key='g%d'%n; model='model_'+key
        stages.append((model,'model','global',n,models[key],'4G',['smoke']))
        deps=['model_g%d'%v for v in (8192,16384,32768) if v<=n]
        if n>8192: deps.append('analyze_g%d'%(n//2))
        stages.append(('analyze_'+key,'analyze','global',n,analyses[key],'16G',deps))
    stages += [('report_global','report','global',0,0,'16G',['analyze_g8192','analyze_g16384','analyze_g32768']),
               ('model_n8192','model','neighborhood10pct',8192,models['n8192'],'4G',['audit_neighborhood','smoke']),
               ('analyze_n8192','analyze','neighborhood10pct',8192,analyses['n8192'],'8G',['model_n8192']),
               ('report_neighborhood','report','neighborhood10pct',0,0,'16G',['analyze_n8192'])]
    return stages


def freeze(repo,root,revision):
    selection=load(root/'selection.json')
    if sha(root/'targets.tsv.gz')!=selection['targets_sha256']: raise ValueError('Frozen targets changed')
    source=repo/'oxygen/results/eFAST/sobol/full_N4096'
    if sha(source/'slurm/submission_plan.json')!=selection['source_plan_sha256']:
        raise ValueError('Source design does not match selected local results')
    if sha(source/'status.json')!=selection['source_status_sha256']:
        raise ValueError('Source completion changed')
    plan_path=root/'plan.json'
    if plan_path.exists():
        p=load(plan_path)
        if p['revision']!=revision or p['selection_sha256']!=sha(root/'selection.json'):
            raise ValueError('Existing refinement belongs to another code/selection')
        for name,expected in p['controller_sha256'].items():
            if sha(root/'slurm/code_snapshot'/name)!=expected: raise ValueError('Frozen controller changed')
        return p
    relative='oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/'
    snapshot=root/'slurm/code_snapshot'; snapshot.mkdir(parents=True,exist_ok=True)
    for name in FILES:
        data=subprocess.check_output(['git','-C',str(repo),'show',revision+':'+relative+name])
        dest=snapshot/name
        if dest.exists() and dest.read_bytes()!=data: raise ValueError('Refusing to overwrite frozen code')
        dest.write_bytes(data); dest.chmod(0o444)
    for d in ('slurm/logs','audits','failures'): (root/d).mkdir(parents=True,exist_ok=True)
    # Copy the exact existing automatic-extension implementation and source helpers.
    monitor=root/'time_extension'
    for d in ('bin','logs','source_scripts'): (monitor/d).mkdir(parents=True,exist_ok=True)
    for name,target in [('time_extension.py','monitor.py'),('run_time_extension.sbatch','run_monitor.sbatch'),
                        ('s_add_job_hours_proxy.sh','bin/s-add-job-hours')]:
        shutil.copyfile(str(snapshot/name),str(monitor/target))
    (monitor/'bin/s-add-job-hours').chmod(0o755)
    for name in ('submit_hourly_add_hours.sh','add_hours_to_running_array_jobs.sh'):
        src=Path('/home/4482173')/name
        if not src.is_file(): raise ValueError('Missing authorized home extension helper: '+str(src))
        shutil.copyfile(str(src),str(monitor/'source_scripts'/name))
    parent=root.parent
    p=dict(status='frozen',revision=revision,created_utc=now(),source_root=str(source),
       source_plan_sha256=selection['source_plan_sha256'],selection_sha256=sha(root/'selection.json'),
       global_root=str(parent/'global_refinement'),neighborhood_root=str(parent/'neighborhood_refinement'),
       cases=selection['cases'],global_N=32768,neighborhood_N=8192,n_repeats=5,
       n_target_seeds=497,n_target_seed_oxygen=55089,
       controller_sha256={name:sha(snapshot/name) for name in FILES},
       qos='xxlarge',initial_time='12:00:00',specified_node=None,array_concurrency_limit=None,
       convergence_policy='diagnostic_only; preserve all original endpoints',
       original_model_evaluations_to_repeat=0)
    save(plan_path,p)
    plan_path.chmod(0o444)
    return p


def sbatch_args(root,stage,dependencies):
    name,role,scope,n,count,memory,_=stage
    exports='ALL,SOBOL_REFINEMENT_ROOT={},SOBOL_REFINEMENT_ROLE={}'.format(root,role)
    if scope: exports+=',SOBOL_REFINEMENT_SCOPE='+scope
    if n: exports+=',SOBOL_REFINEMENT_N='+str(n)
    args=['sbatch','--parsable','--hold','--qos=xxlarge','--time=12:00:00','--ntasks=1',
          '--cpus-per-task=1','--mem='+memory,'--job-name=sref_'+name,
          '--comment=sref:'+hashlib.sha256(str(root).encode()).hexdigest()[:16],
          '--chdir='+str(root),'--export='+exports,
          '--output='+str(root/'slurm/logs'/(name+'_%A_%a.out')),
          '--error='+str(root/'slurm/logs'/(name+'_%A_%a.err'))]
    if count: args.append('--array=1-'+str(count))
    if dependencies: args.append('--dependency=afterok:'+':'.join(dependencies))
    return args+[str(root/'slurm/code_snapshot/run_refinement.sbatch')]


def submit(repo,root,revision):
    p=freeze(repo,root,revision)
    if (root/'launch_receipt.json').exists():
        print(json.dumps(load(root/'launch_receipt.json'),indent=2),flush=True)
        return
    ledger=root/'slurm/jobs.json'; rows=load(ledger) if ledger.exists() else []
    jobs={row['stage']:row['job_id'] for row in rows}
    selection=load(root/'selection.json')
    for stage in stage_specs(selection):
        name=stage[0]
        if name in jobs: continue
        # Recover a submit/ledger interruption without duplicating a held job.
        found=output(['squeue','-h','-u','4482173','-n','sref_'+name,'-o','%A']).splitlines()
        found=sorted(set(x.strip() for x in found if x.strip()))
        if len(found)>1: raise ValueError('Multiple jobs for stage '+name)
        if found:
            details=output(['scontrol','show','job','-o',found[0]])
            if ('WorkDir='+str(root)) not in details or ('Command='+str(root/'slurm/code_snapshot/run_refinement.sbatch')) not in details:
                raise ValueError('Same-name job belongs to a different task; refusing adoption')
        args=sbatch_args(root,stage,[jobs[d] for d in stage[-1]])
        job=found[0] if found else output(args).split(';')[0]
        if not job.isdigit(): raise ValueError('Invalid submitted job id')
        row=dict(stage=name,job_id=job,role=stage[1],scope=stage[2],N=stage[3],
                 array_count=stage[4],memory=stage[5],dependencies=[jobs[d] for d in stage[-1]],
                 submitted_utc=now(),sbatch_args=args)
        rows.append(row); save(ledger,rows); jobs[name]=job
        print('SUBMITTED',name,job,flush=True)
    monitor=root/'time_extension'; registry=monitor/'targets.json'
    targets=dict(target_jobs=[r['job_id'] for r in rows],refinement_plan_sha256=sha(root/'plan.json'))
    if registry.exists() and load(registry)!=targets: raise ValueError('Monitor registry changed')
    if not registry.exists(): save(registry,targets)
    subprocess.run(['python3',str(monitor/'monitor.py'),'submit'],check=True)
    monitors=read_monitors(monitor/'submission_jobs.tsv')
    if {r['target_job_id'] for r in monitors}!=set(targets['target_jobs']): raise ValueError('Missing extension monitor')
    # All jobs were held until every stage and its monitor had a durable receipt.
    for row in rows:
        state=output(['squeue','-h','-j',row['job_id'],'-o','%T'])
        if state and set(state.splitlines())=={'PENDING'}:
            subprocess.run(['scontrol','release',row['job_id']],check=True)
    save(root/'launch_receipt.json',dict(status='submitted',revision=revision,jobs=rows,monitors=monitors,
         n_model_tasks=75600,n_analysis_tasks=2875,n_repeats=5,n_target_seeds=497,n_target_seed_oxygen=55089,
         initial_qos='xxlarge',initial_time='12:00:00',specified_node=None,array_concurrency_limit=None,
         extension_threshold_hours=6,extension_increment_hours=6,submitted_utc=now()))
    print(json.dumps(dict(status='submitted',jobs=jobs,monitors=monitors),indent=2),flush=True)


def read_monitors(path):
    with path.open() as h: return list(csv.DictReader(h,delimiter='\t'))


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--repo',type=Path,default=Path('/share/lab_crd/taoli/Project/soft_couping_eFAST'))
    a.add_argument('--revision',required=True)
    args=a.parse_args(); root=args.repo/'oxygen/results/eFAST/sobol/refinement_20261018'
    root.mkdir(parents=True,exist_ok=True)
    with (root/'submission.lock').open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        submit(args.repo,root,args.revision)


if __name__=='__main__': main()
