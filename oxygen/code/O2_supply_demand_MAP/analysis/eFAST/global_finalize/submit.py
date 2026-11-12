#!/usr/bin/env python3
"""Freeze and submit one global eFAST batch, without an upper N restriction."""
import argparse
import csv
import fcntl
import hashlib
from pathlib import Path
import shutil
import subprocess

from common import (INITIAL_NS, M, ACTIVE, OUTPUTS, now, load, save, sha,
                    valid_n, model_count, resolutions, efast_directory, next_n)

FILES = ('common.py', 'workflow.py', 'report.py', 'submit.py', 'test_workflow.py', 'collect.py', 'finish_report.py',
         'evaluate_efast.R', 'run.sbatch', 'plot_helpers.py', 'time_extension.py',
         'run_time_extension.sbatch', 's_add_job_hours_proxy.sh', 'retry_global_startup.sh')
RELATIVE = 'oxygen/code/O2_supply_demand_MAP/analysis/eFAST/global_finalize/'


def output(args):
    return subprocess.check_output(args, universal_newlines=True).strip()


def stages(ns, refinement=False):
    rows = [('prepare', 'prepare', 0, 0, '8G', []),
            ('smoke', 'smoke', ns[0], 0, '4G', ['prepare'])]
    if not refinement:
        rows.append(('sobol', 'sobol_report', 0, 0, '8G', ['prepare']))
    previous = None
    for n in ns:
        key = str(n)
        model_deps = [previous] if refinement and previous else ['smoke']
        rows += [('model_'+key, 'model', n, model_count(n), '4G', model_deps),
                 ('analyze_'+key, 'analyze', n, 70, '8G', ['model_'+key]),
                 ('report_'+key, 'report', n, 0, '8G', ['analyze_'+key] + ([previous] if previous else []))]
        previous = 'report_'+key
    rows.append(('final', 'final_report', 0, 0, '16G',
                 ([] if refinement else ['sobol'])+['report_'+str(n) for n in ns]))
    return rows


def freeze(repo, root, revision, ns, previous_root=None):
    ns = sorted(set(map(valid_n, ns)))
    if not ns:
        raise ValueError('A batch needs at least one N')
    if (root/'plan.json').exists():
        p = load(root/'plan.json')
        if (p['revision'] != revision or p['resolutions'] != ns or
                p.get('previous_root') != (str(previous_root) if previous_root else None)):
            raise ValueError('Resume the immutable batch; use a separate root for a new batch')
        return p
    initial = repo/'oxygen/results/eFAST/sobol/full_N4096'
    sobol = repo/'oxygen/results/eFAST/sobol/global_convergence_20261021'
    old = load(initial/'slurm/submission_plan.json')
    done = load(sobol/'N131072/completion.json')
    if done['status'] != 'complete' or done['N'] != 131072 or done['n_repeats'] != 5 or done['n_oxygen'] != 201:
        raise ValueError('Require certified completed global Sobol N131072')
    if sha(sobol/'N131072/indices_summary.tsv.gz') != done['summary_sha256']:
        raise ValueError('Sobol summary differs from its completion receipt')
    snapshot = root/'slurm/code_snapshot'
    snapshot.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        data = subprocess.check_output(['git', '-C', str(repo), 'show', revision+':'+RELATIVE+name])
        target = snapshot/name
        if target.exists() and target.read_bytes() != data:
            raise ValueError('Refusing to replace a frozen controller')
        if not target.exists():
            target.write_bytes(data); target.chmod(0o444)
    inputs = root/'inputs'
    inputs.mkdir(exist_ok=True)
    source = {}

    def copy(src, rel):
        target = root/rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and sha(target) != sha(src):
            raise ValueError('Existing copied input differs: '+rel)
        if not target.exists():
            shutil.copyfile(str(src), str(target))
        source[rel] = sha(target)

    for name in ['indices_summary.tsv.gz', 'indices_repeats.npz', 'completion.json', 'classification.tsv', 'row_order.tsv']:
        copy(sobol/'N131072'/name, 'inputs/sobol_N131072/'+name)
    copy(initial/'global/full/global/parameter_bounds.tsv', 'inputs/parameter_bounds.tsv')
    copy(repo/'oxygen/results/eFAST/parameter_ranges.tsv', 'inputs/legacy_parameter_ranges.tsv')
    copy(repo/'oxygen/results/eFAST/indices.tsv', 'inputs/legacy_efast_indices.tsv')
    for name in ['figure4_parameter_groups_source.tsv', 'figure4_parameter_group_palette_source.tsv', 'figure4b_o2_classification_source.tsv']:
        copy(sobol/'inputs'/name, 'inputs/'+name)
    for name in ['model_equations.md', 'parameter_equation_map.tsv', 'parameter_pair_mechanism_map.tsv']:
        copy(initial/'model_structure'/name, 'inputs/model_structure/'+name)
    copy(initial/'response_surface_inventory.tsv', 'inputs/response_surface_inventory.tsv')
    for suffix in ['pdf', 'png']:
        copy(repo/'oxygen/results/eFAST/figures'/('figure4b_spearman_direction.'+suffix),
             'inputs/figure4b_spearman_direction.'+suffix)
    references = []
    if previous_root:
        previous_root = Path(previous_root)
        prior = load(previous_root/'plan.json')
        completed = load(previous_root/'completion.json')
        if (completed['status'] != 'complete' or prior['scope'] != 'global' or
                prior['method'] != 'eFAST' or prior['M'] != M or
                prior['n_repeats'] != 5 or prior['n_oxygen'] != 201 or
                prior['phase_seeds'] != list(range(20260923,20260928)) or
                prior['parameters'] != list(ACTIVE) or prior['outputs'] != list(OUTPUTS) or
                prior['oxygen_pct'] != [i/40. for i in range(201)] or
                prior['sif'] != old['sif'] or
                prior['original_plan_sha256'] != sha(initial/'slurm/submission_plan.json')):
            raise ValueError('Previous global batch is incomplete or scientifically different')
        references = resolutions(prior)
        if ns[0] != next_n(references[-1]) or any(b != next_n(a) for a,b in zip(ns,ns[1:])):
            raise ValueError('Refinement must retain consecutive doubling resolutions')
        if sha(previous_root/'inputs/parameter_bounds.tsv') != sha(root/'inputs/parameter_bounds.tsv'):
            raise ValueError('Previous and current bounds differ')
        for name in ['plan.json','completion.json','preparation.json']:
            copy(previous_root/name, 'inputs/previous_batch/'+name)
        copy(previous_root/'validation/method_alignment.json', 'inputs/previous_method_alignment.json')
        for name in ['source_sobol_submission_plan.json','source_sobol_evaluate.R']:
            copy(previous_root/'inputs'/name, 'inputs/'+name)
        for name in ['best_params.tsv','fit_config.rds','parameter_table.csv']:
            copy(previous_root/'inputs/reference_seed25'/name, 'inputs/reference_seed25/'+name)
        for n in references:
            directory = efast_directory(previous_root, prior, n)
            receipt = load(directory/'completion.json')
            if (receipt['status'] != 'complete' or receipt['N'] != n or receipt['M'] != M or
                    receipt['n_repeats'] != 5 or receipt['n_oxygen'] != 201 or
                    receipt['summary_sha256'] != sha(directory/'indices_summary.tsv.gz')):
                raise ValueError('Previous eFAST resolution not certified: '+str(n))
            for path in sorted(directory.iterdir()):
                if path.is_file():
                    copy(path, 'inputs/previous_efast/N%d/%s' % (n,path.name))
        # Reuse the completed Sobol report byte-for-byte; no Sobol job is submitted.
        if load(previous_root/'sobol/completion.json')['status'] != 'complete':
            raise ValueError('Previous Sobol report is incomplete')
        for path in sorted((previous_root/'sobol').rglob('*')):
            if path.is_file():copy(path, 'sobol/'+str(path.relative_to(previous_root/'sobol')))
        save(root/'validation/reference_reuse.json',dict(status='complete',previous_root=str(previous_root),
             reference_resolutions=references,sobol_report_copied=True,sobol_model_calls=0,
             sobol_estimator_calls=0,old_efast_model_calls=0,raw_blocks_copied=False,completed_utc=now()))
    monitor = root/'time_extension'
    for name in ['slurm/logs', 'failures', 'validation', 'time_extension/logs', 'time_extension/bin', 'time_extension/source_scripts']:
        (root/name).mkdir(parents=True, exist_ok=True)
    for name, target in [('time_extension.py', 'monitor.py'), ('run_time_extension.sbatch', 'run_monitor.sbatch'),
                         ('s_add_job_hours_proxy.sh', 'bin/s-add-job-hours')]:
        shutil.copyfile(str(snapshot/name), str(monitor/target))
    (monitor/'bin/s-add-job-hours').chmod(0o755)
    for name in ['submit_hourly_add_hours.sh', 'add_hours_to_running_array_jobs.sh']:
        shutil.copyfile('/home/4482173/'+name, str(monitor/'source_scripts'/name))
    p = dict(scope='global', method='eFAST', revision=revision, created_utc=now(), resolutions=ns,
             previous_root=str(previous_root) if previous_root else None,reference_resolutions=references,
             execution_order='sequential_by_resolution' if previous_root else 'independent_model_arrays',
             N_hard_limit=None, next_resolution_rule='2*N-1; separate immutable subsequent batch', M=M,
             n_repeats=5, phase_seeds=list(range(20260923,20260928)), n_oxygen=201,
             oxygen_pct=[i/40. for i in range(201)], parameters=list(ACTIVE), outputs=list(OUTPUTS),
             source_initial_root=str(initial), source_sobol_root=str(sobol), source_sobol_N=131072,
             original_plan=str(initial/'slurm/submission_plan.json'), original_plan_sha256=sha(initial/'slurm/submission_plan.json'),
             canonical_helpers=str(Path(old['code_snapshot_root'])/'oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol'),
             sif=old['sif'], source_files=source, controller_sha256={name:sha(snapshot/name) for name in FILES},
             qos='xxlarge', initial_time='12:00:00', specified_node=None, array_concurrency_limit=None,
             extension_threshold_hours=6, extension_increment_hours=6, convergence_policy='diagnostic_only',
             repeat_range_tolerance=.05, resolution_delta_tolerance=.02, sample_chunk=1024,
             model_tasks=sum(model_count(n) for n in ns), model_points=sum(ns)*14*201*5,
             neighborhood_jobs_to_submit=0, sobol_model_jobs_to_submit=0,sobol_estimator_jobs_to_submit=0,
             sobol_report_jobs_to_submit=0 if previous_root else 1)
    save(root/'plan.json', p)
    (root/'plan.json').chmod(0o444)
    return p


def sbatch_args(root, stage, dependencies):
    name, role, n, count, memory, _ = stage
    prefix = 'gef_'+hashlib.sha256(str(root).encode()).hexdigest()[:6]+'_'
    args = ['sbatch','--parsable','--hold','--qos=xxlarge','--time=12:00:00',
            '--ntasks=1','--cpus-per-task=1','--mem='+memory,'--job-name='+prefix+name,
            '--comment=global_efast_finalize','--chdir='+str(root),'--export=ALL',
            '--output='+str(root/'slurm/logs'/(name+'_%A_%a.out')),
            '--error='+str(root/'slurm/logs'/(name+'_%A_%a.err'))]
    if count:
        args += ['--array=1-'+str(count)]
    if dependencies:
        args += ['--dependency=afterok:'+':'.join(dependencies)]
    return args+[str(root/'slurm/code_snapshot/run.sbatch'),str(root),role,str(n)]


def submit(repo, root, revision, ns, previous_root=None):
    p=freeze(repo,root,revision,ns,previous_root)
    if (root/'launch_receipt.json').exists():
        print((root/'launch_receipt.json').read_text());return
    ledger=root/'slurm/jobs.json'
    rows=load(ledger) if ledger.exists() else []
    jobs={r['stage']:r['job_id'] for r in rows}
    for stage in stages(p['resolutions'], refinement=bool(p.get('previous_root'))):
        name,role,n,count,memory,deps=stage
        if name in jobs:continue
        args=sbatch_args(root,stage,[jobs[d] for d in deps])
        job_name=[a.split('=',1)[1] for a in args if a.startswith('--job-name=')][0]
        found=set(output(['squeue','-h','-u','4482173','-n',job_name,'-o','%A']).splitlines())
        if len(found)>1:raise ValueError('Duplicate jobs for stage '+name)
        if found:
            job=found.pop()
            details=output(['scontrol','show','job','-o',job])
            if 'WorkDir='+str(root) not in details:raise ValueError('Unrelated same-name job')
        else:
            job=output(args).split(';')[0]
        if not job.isdigit():raise ValueError('Invalid submission receipt')
        rows.append(dict(stage=name,role=role,N=n,job_id=job,array_count=count,dependencies=[jobs[d] for d in deps],
                         sbatch_args=args,submitted_utc=now()))
        save(ledger,rows);jobs[name]=job
        print('SUBMITTED',name,job,flush=True)
    monitor=root/'time_extension'
    targets=dict(target_jobs=[r['job_id'] for r in rows],global_plan_sha256=sha(root/'plan.json'))
    if (monitor/'targets.json').exists() and load(monitor/'targets.json')!=targets:
        raise ValueError('Monitor registry differs')
    save(monitor/'targets.json',targets)
    subprocess.run(['python3',str(monitor/'monitor.py'),'submit'],check=True)
    with (monitor/'submission_jobs.tsv').open() as handle:
        monitors=list(csv.DictReader(handle,delimiter='\t'))
    if {r['target_job_id'] for r in monitors}!=set(targets['target_jobs']):
        raise ValueError('Monitor coverage incomplete')
    for row in rows:
        if output(['squeue','-h','-j',row['job_id'],'-o','%T']):
            subprocess.run(['scontrol','release',row['job_id']],check=True)
    save(root/'launch_receipt.json',dict(status='submitted',scope='global',revision=revision,jobs=rows,monitors=monitors,
        resolutions=p['resolutions'],model_tasks=p['model_tasks'],model_points=p['model_points'],
        n_repeats=5,n_oxygen=201,N_hard_limit=None,convergence_policy='diagnostic_only',submitted_utc=now()))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=Path('/share/lab_crd/taoli/Project/soft_couping_eFAST'))
    parser.add_argument('--root',type=Path)
    parser.add_argument('--revision',required=True)
    parser.add_argument('--resolutions',type=int,nargs='+',default=list(INITIAL_NS))
    parser.add_argument('--previous-root',type=Path,help='Completed immutable global eFAST batch to reuse')
    args=parser.parse_args();root=args.root or args.repo/'oxygen/results/eFAST/global_finalize_20261112'
    root.mkdir(parents=True,exist_ok=True)
    with (root/'submission.lock').open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        submit(args.repo,root,args.revision,args.resolutions,args.previous_root)


if __name__=='__main__':main()
