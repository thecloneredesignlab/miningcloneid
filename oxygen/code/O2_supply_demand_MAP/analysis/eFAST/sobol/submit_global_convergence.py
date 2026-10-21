#!/usr/bin/env python3
"""Idempotently submit only the approved global N65536/N131072 continuation."""
import argparse
import csv
import fcntl
import hashlib
from pathlib import Path
import shutil
import subprocess

from global_convergence import load, save, sha, now, NS, model_count

FILES = ('global_convergence.py', 'global_convergence_plot.py', 'submit_global_convergence.py',
         'test_global_convergence.py', 'run_global_convergence.sbatch', 'retry_global_startup.sh',
         'time_extension.py', 'run_time_extension.sbatch', 's_add_job_hours_proxy.sh')


def output(args):
    return subprocess.check_output(args, universal_newlines=True).strip()


def stages():
    rows = [('prepare', 'prepare', 0, 0, '4G', []), ('smoke', 'smoke', 65536, 0, '4G', ['prepare'])]
    for n in NS:
        key = str(n)
        rows.append(('model_' + key, 'model', n, model_count(n), '4G', ['smoke'] if n == NS[0] else ['report_65536']))
        rows.append(('analyze_' + key, 'analyze', n, 130, '16G', ['model_' + key]))
        rows.append(('report_' + key, 'report', n, 0, '16G', ['analyze_' + key]))
    return rows


def freeze(repo, root, revision):
    planfile = root / 'plan.json'
    if planfile.exists():
        p = load(planfile)
        if p['revision'] != revision:
            raise ValueError('Existing stage belongs to another revision; resume that stage')
        for name, digest in p['controller_sha256'].items():
            if sha(root / 'slurm/code_snapshot' / name) != digest:
                raise ValueError('Frozen controller changed')
        return p
    base = repo / 'oxygen/results/eFAST/sobol'
    oldref, initial, oldglobal = base / 'refinement_20261018', base / 'full_N4096', base / 'global_refinement'
    if load(oldglobal / 'status.json')['status'] != 'complete':
        raise ValueError('Require completed global N32768')
    snapshot = root / 'slurm/code_snapshot'
    snapshot.mkdir(parents=True, exist_ok=True)
    relative = 'oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/'
    for name in FILES:
        data = subprocess.check_output(['git', '-C', str(repo), 'show', revision + ':' + relative + name])
        dest = snapshot / name
        if dest.exists() and dest.read_bytes() != data:
            raise ValueError('Refusing to overwrite a frozen file')
        if not dest.exists():
            dest.write_bytes(data)
            dest.chmod(0o444)
    for d in ['slurm/logs', 'failures', 'time_extension/bin', 'time_extension/logs', 'time_extension/source_scripts', 'inputs']:
        (root / d).mkdir(parents=True, exist_ok=True)
    monitor = root / 'time_extension'
    for name, target in [('time_extension.py', 'monitor.py'), ('run_time_extension.sbatch', 'run_monitor.sbatch'),
                         ('s_add_job_hours_proxy.sh', 'bin/s-add-job-hours')]:
        shutil.copyfile(str(snapshot / name), str(monitor / target))
    (monitor / 'bin/s-add-job-hours').chmod(0o755)
    for name in ['submit_hourly_add_hours.sh', 'add_hours_to_running_array_jobs.sh']:
        shutil.copyfile('/home/4482173/' + name, str(monitor / 'source_scripts' / name))
    source_files = {}
    for path in [oldref / 'slurm/code_snapshot/refinement_run.py', oldglobal / 'status.json',
                 initial / 'global/full/global/completion.json', initial / 'global/full/global/parameter_bounds.tsv',
                 initial / 'global/full/global/indices_summary.npz', initial / 'global/full/global/indices_repeats.npz']:
        source_files[str(path)] = sha(path)
    for n in [8192, 16384, 32768]:
        for rel in ['completion.json', 'cases/global/indices_summary.npz', 'cases/global/indices_summary.tsv.gz']:
            path = oldglobal / ('N%d' % n) / rel
            source_files[str(path)] = sha(path)
    for name in ['figure4_parameter_groups_source.tsv', 'figure4_parameter_group_palette_source.tsv', 'figure4b_o2_classification_source.tsv']:
        src = repo / 'oxygen/results/eFAST' / name
        dest = root / 'inputs' / name
        shutil.copyfile(str(src), str(dest))
        source_files[str(dest)] = sha(dest)
    p = dict(status='frozen', scope='global', revision=revision, created_utc=now(), source_refinement_root=str(oldref),
        source_refinement_plan_sha256=sha(oldref / 'plan.json'), source_initial_root=str(initial), source_global_root=str(oldglobal),
        source_files=source_files, source_N=32768, resolutions=list(NS), n_repeats=5, n_oxygen=201,
        n_parameters=14, n_pairs=91, controller_sha256={name: sha(snapshot / name) for name in FILES},
        qos='xxlarge', initial_time='12:00:00', specified_node=None, array_concurrency_limit=None,
        extension_threshold_hours=6, extension_increment_hours=6, convergence_policy='diagnostic_only',
        original_model_evaluations_to_repeat=0, neighborhood_jobs_to_submit=0)
    save(planfile, p)
    planfile.chmod(0o444)
    return p


def sbatch_args(root, stage, dependencies):
    name, role, n, count, memory, _ = stage
    args = ['sbatch', '--parsable', '--hold', '--qos=xxlarge', '--time=12:00:00', '--ntasks=1', '--cpus-per-task=1',
            '--mem=' + memory, '--job-name=gconv_' + name, '--comment=gconv:' + hashlib.sha256(str(root).encode()).hexdigest()[:16],
            '--chdir=' + str(root), '--export=ALL,GLOBAL_CONVERGENCE_ROOT={},GLOBAL_CONVERGENCE_ROLE={},GLOBAL_CONVERGENCE_N={}'.format(root, role, n or 65536),
            '--output=' + str(root / 'slurm/logs' / (name + '_%A_%a.out')),
            '--error=' + str(root / 'slurm/logs' / (name + '_%A_%a.err'))]
    if count:
        args.append('--array=1-' + str(count))
    if dependencies:
        args.append('--dependency=afterok:' + ':'.join(dependencies))
    return args + [str(root / 'slurm/code_snapshot/run_global_convergence.sbatch')]


def submit(repo, root, revision):
    freeze(repo, root, revision)
    if (root / 'launch_receipt.json').exists():
        print((root / 'launch_receipt.json').read_text())
        return
    ledger = root / 'slurm/jobs.json'
    rows = load(ledger) if ledger.exists() else []
    jobs = {r['stage']: r['job_id'] for r in rows}
    for stage in stages():
        name, role, n, count, memory, deps = stage
        if name in jobs:
            continue
        found = sorted(set(output(['squeue', '-h', '-u', '4482173', '-n', 'gconv_' + name, '-o', '%A']).splitlines()))
        if len(found) > 1:
            raise ValueError('Multiple candidate jobs for ' + name)
        if found:
            details = output(['scontrol', 'show', 'job', '-o', found[0]])
            if 'WorkDir=' + str(root) not in details or 'Command=' + str(root / 'slurm/code_snapshot/run_global_convergence.sbatch') not in details:
                raise ValueError('Unrelated same-name job; refusing adoption')
        args = sbatch_args(root, stage, [jobs[x] for x in deps])
        job = found[0] if found else output(args).split(';')[0]
        if not job.isdigit():
            raise ValueError('Invalid submission receipt')
        rows.append(dict(stage=name, role=role, job_id=job, N=n, array_count=count, memory=memory,
            dependencies=[jobs[x] for x in deps], sbatch_args=args, submitted_utc=now()))
        save(ledger, rows)
        jobs[name] = job
        print('SUBMITTED', name, job, flush=True)
    monitor = root / 'time_extension'
    targets = dict(target_jobs=[r['job_id'] for r in rows], global_plan_sha256=sha(root / 'plan.json'))
    if (monitor / 'targets.json').exists() and load(monitor / 'targets.json') != targets:
        raise ValueError('Monitor registry differs')
    save(monitor / 'targets.json', targets)
    subprocess.run(['python3', str(monitor / 'monitor.py'), 'submit'], check=True)
    with (monitor / 'submission_jobs.tsv').open() as h:
        monitors = list(csv.DictReader(h, delimiter='\t'))
    if {r['target_job_id'] for r in monitors} != set(targets['target_jobs']):
        raise ValueError('Missing monitor receipt')
    for row in rows:
        state = output(['squeue', '-h', '-j', row['job_id'], '-o', '%T'])
        if state and set(state.splitlines()) == {'PENDING'}:
            subprocess.run(['scontrol', 'release', row['job_id']], check=True)
    save(root / 'launch_receipt.json', dict(status='submitted', revision=revision, scope='global', jobs=rows, monitors=monitors,
        model_tasks=3600, analysis_tasks=260, repeats=5, oxygen_points=201, convergence_policy='diagnostic_only',
        extension_threshold_hours=6, extension_increment_hours=6, submitted_utc=now()))


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--repo', type=Path, default=Path('/share/lab_crd/taoli/Project/soft_couping_eFAST'))
    a.add_argument('--revision', required=True)
    args = a.parse_args()
    root = args.repo / 'oxygen/results/eFAST/sobol/global_convergence_20261021'
    root.mkdir(parents=True, exist_ok=True)
    with (root / 'submission.lock').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        submit(args.repo, root, args.revision)


if __name__ == '__main__':
    main()
