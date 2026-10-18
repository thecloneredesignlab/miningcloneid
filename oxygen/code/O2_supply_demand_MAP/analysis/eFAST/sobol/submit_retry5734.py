#!/usr/bin/env python3
"""Precisely restore recovery shard 5734 and its downstream chain (RED Python 3.6+)."""
import argparse
import csv
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
REPO = Path('/share/lab_crd/taoli/Project/soft_couping_eFAST')
SIF = '/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif'
STAGES = [('shard', '5734', '8G'), ('case', '173', '8G'), ('audit', '', '8G'),
          ('aggregate', '', '16G'), ('response', '1-54', '4G'), ('finish', '', '8G')]
OLD_JOBS = ('21129866', '21129867', '21129868')
OLD_MONITORS = ('21138286', '21138287', '21138288')
FILES = ('run_neighborhood_recovery.sbatch', 'retry_container_startup.sh',
         'submit_retry5734.py', 'test_retry5734.py', 'time_extension.py',
         'test_time_extension.py', 'run_time_extension.sbatch', 's_add_job_hours_proxy.sh')


def now():
    return datetime.datetime.utcnow().isoformat() + 'Z'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    temp.write_text(json.dumps(data, indent=2) + '\n'); os.replace(str(temp), str(path))


def output(args):
    return subprocess.check_output(args, universal_newlines=True)


def locations(repo):
    root = repo / 'oxygen/results/eFAST/sobol/full_N4096'
    rec = root / 'recovery_20261018'
    return root, rec, rec / 'retry5734'


def prepare(repo, revision):
    root, rec, retry = locations(repo)
    if (retry / 'preflight.json').exists():
        receipt = json.loads((retry / 'preflight.json').read_text())
        if receipt['controller_revision'] != revision or receipt['plan_sha256'] != sha(rec / 'plan.json'):
            raise ValueError('Existing retry belongs to a different source or plan')
        print('Reusing certified preparation', flush=True)
        return
    rp = json.loads((rec / 'plan.json').read_text())
    task = rp['tasks'][5733]; case = rp['cases'][172]['case']
    if (task['task_id'], task['case_id'], task['repeat'], task['oxygen_indices']) != (
            5734, 'seed283', 2, list(range(160, 176))):
        raise ValueError('Failed shard differs from approved scope')
    if case['case_id'] != 'seed283' or case['case_index'] != 284:
        raise ValueError('Case array index 173 differs')
    missing = [t['task_id'] for t in rp['tasks']
               if not (rec / 'shards' / ('%05d' % t['task_id']) / 'completion.json').is_file()]
    if missing != [5734]:
        raise ValueError('Missing shards changed: ' + str(missing))
    plan = json.loads((root / 'slurm/submission_plan.json').read_text())
    preserved = []; incomplete = []
    for row in plan['cases']:
        receipt = root / row['scope'] / 'full' / row['case_id'] / 'completion.json'
        if receipt.is_file():
            if json.loads(receipt.read_text())['status'] != 'complete':
                raise ValueError('Uncertified existing case')
            preserved.append(dict(path=str(receipt.relative_to(root)), sha256=sha(receipt)))
        else:
            incomplete.append(row['case_id'])
    if incomplete != ['seed283'] or len(preserved) != 500:
        raise ValueError('Expected 499 neighborhoods plus global; missing ' + str(incomplete))
    for row in rp['model_tasks']:
        file = root / 'neighborhood10pct/full' / row['case_id'] / 'runs' / ('R%02d' % row['repeat']) / 'blocks' / row['block'] / 'outputs.npz'
        receipt = Path(str(file) + '.receipt.json')
        rr = json.loads(receipt.read_text())
        if not file.is_file() or rr['status'] != 'complete' or rr['n_points'] != 4096 * 201:
            raise ValueError('A repaired model block is no longer complete')
    snapshot = retry / 'slurm/code_snapshot'; snapshot.mkdir(parents=True, exist_ok=True)
    relative = str(HERE.relative_to(repo))
    for name in FILES:
        data = subprocess.check_output(['git', '-C', str(repo), 'show', revision + ':' + relative + '/' + name])
        dest = snapshot / name
        if dest.exists() and dest.read_bytes() != data:
            raise ValueError('Refusing to replace immutable retry source')
        if not dest.exists(): dest.write_bytes(data)
        dest.chmod(0o444)
    snapshot.chmod(0o555)
    monitor = retry / 'time_extension'
    for sub in ('bin', 'logs', 'source_scripts'):
        (monitor / sub).mkdir(parents=True, exist_ok=True)
    for name, dest in [('time_extension.py', 'monitor.py'),
                       ('run_time_extension.sbatch', 'run_monitor.sbatch'),
                       ('s_add_job_hours_proxy.sh', 'bin/s-add-job-hours')]:
        shutil.copyfile(str(snapshot / name), str(monitor / dest))
    (monitor / 'bin/s-add-job-hours').chmod(0o755)
    for name in ('submit_hourly_add_hours.sh', 'add_hours_to_running_array_jobs.sh'):
        shutil.copyfile(str(rec / 'time_extension/source_scripts' / name), str(monitor / 'source_scripts' / name))
    receipt = dict(prepared_utc=now(), controller_revision=revision,
                   frozen_analysis_revision=rp['recovery_git_commit'],
                   plan_sha256=sha(rec / 'plan.json'), source_revision=rp['source_git_commit'],
                   preserved_cases=preserved, n_preserved_shards=12058,
                   shard_id=5734, case_task_id=173, case_id='seed283',
                   n_preserved_neighborhoods=499, model_blocks_reused=14,
                   code_sha256={name:sha(snapshot / name) for name in FILES},
                   convergence_policy='diagnostic_only')
    save(retry / 'preflight.json', receipt)
    print('Certified 12058 completed shards, 499 neighborhoods, global, and 14 repaired model blocks', flush=True)


def arguments(repo, stage, predecessor=None, hold=False):
    root, rec, retry = locations(repo); role, array, mem = stage
    env = ('ALL,SOBOL_OUT_DIR={},SOBOL_RECOVERY_DIR={},SOBOL_SIF={},'
           'SOBOL_ROLE={},RECOVERY_PLAN_SHA256={}').format(root, rec, SIF, role, sha(rec / 'plan.json'))
    args = ['sbatch', '--parsable', '--qos=xxlarge', '--time=12:00:00', '--ntasks=1',
            '--cpus-per-task=1', '--mem=' + mem, '--job-name=sobol_retry5734_' + role,
            '--chdir=' + str(repo), '--export=' + env,
            '--output=' + str(retry / 'slurm/logs' / (role + '_%A_%a.out')),
            '--error=' + str(retry / 'slurm/logs' / (role + '_%A_%a.err'))]
    if array: args.append('--array=' + array)
    if predecessor: args.append('--dependency=afterok:' + predecessor)
    if hold: args.append('--hold')
    return args + [str(retry / 'slurm/code_snapshot/run_neighborhood_recovery.sbatch')]


def old_pending():
    for jid in OLD_JOBS:
        result = subprocess.run(['squeue', '-h', '-j', jid, '-o', '%T'],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode and 'Invalid job id specified' not in result.stderr:
            raise RuntimeError(result.stderr)
        data = result.stdout.splitlines()
        if any(row.strip() != 'PENDING' for row in data):
            raise ValueError('Old downstream job is no longer pending: ' + jid)


def submit(repo):
    root, rec, retry = locations(repo)
    if (retry / 'launch_receipt.json').exists():
        print('Retry launch already recorded; no duplicate submission.', flush=True)
        return
    receipt = json.loads((retry / 'preflight.json').read_text())
    if sha(rec / 'plan.json') != receipt['plan_sha256']:
        raise ValueError('Recovery plan changed')
    for name, digest in receipt['code_sha256'].items():
        if sha(retry / 'slurm/code_snapshot' / name) != digest:
            raise ValueError('Retry source changed')
    old_pending()
    (retry / 'slurm/logs').mkdir(parents=True, exist_ok=True)
    ledger = retry / 'slurm/submission_jobs.tsv'
    rows = list(csv.DictReader(ledger.open(), delimiter='\t')) if ledger.exists() else []
    if [row['role'] for row in rows] != [s[0] for s in STAGES[:len(rows)]]:
        raise ValueError('Existing retry ledger is not an exact stage prefix')
    previous = rows[-1]['job_id'] if rows else None
    for stage in STAGES[len(rows):]:
        args = arguments(repo, stage, previous, hold=not rows)
        job = output(args).strip().split(';')[0]
        if not job.isdigit(): raise ValueError('Invalid sbatch receipt: ' + job)
        row = dict(role=stage[0], job_id=job, array_spec=stage[1], memory=stage[2],
                   dependency='afterok:' + previous if previous else '',
                   submitted_utc=now())
        with ledger.open('a') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
            if not rows: writer.writeheader()
            writer.writerow(row)
        rows.append(row); previous = job
        print('SUBMITTED', json.dumps(row), flush=True)
    monitor = retry / 'time_extension'; config = monitor / 'targets.json'
    targets = dict(target_jobs=[row['job_id'] for row in rows], recovery=str(rec),
                   recovery_plan_sha256=receipt['plan_sha256'])
    if config.exists() and json.loads(config.read_text()) != targets:
        raise ValueError('Existing monitor registry differs')
    if not config.exists(): save(config, targets)
    config.chmod(0o444)
    subprocess.run(['python3', '-u', str(monitor / 'monitor.py'), 'submit'], check=True)
    monitored = list(csv.DictReader((monitor / 'submission_jobs.tsv').open(), delimiter='\t'))
    if {row['target_job_id'] for row in monitored} != set(targets['target_jobs']):
        raise ValueError('Not all held/pending targets have monitors')
    # Wait for actual monitor startup before releasing the first analysis task.
    ids = [row['monitor_job_id'] for row in monitored]
    for attempt in range(120):
        states = output(['squeue', '-h', '-j', ','.join(ids), '-o', '%i|%T']).splitlines()
        if len(states) == len(ids) and all(line.endswith('|RUNNING') for line in states): break
        if any(line.endswith(('|FAILED', '|CANCELLED', '|TIMEOUT')) for line in states):
            raise ValueError('A new monitor failed before release')
        print('Waiting for monitors to start', attempt + 1, flush=True); time.sleep(5)
    else:
        raise ValueError('Monitors still pending; first shard remains held')
    old_pending()
    subprocess.run(['scancel'] + list(OLD_JOBS), check=True)
    subprocess.run(['scancel'] + list(OLD_MONITORS), check=True)
    subprocess.run(['scontrol', 'release', rows[0]['job_id']], check=True)
    save(retry / 'launch_receipt.json', dict(status='submitted', jobs=rows, monitors=monitored,
         cancelled_old_jobs=list(OLD_JOBS), cancelled_old_monitors=list(OLD_MONITORS),
         source_revision=receipt['controller_revision'], initial_qos='xxlarge',
         initial_time='12:00:00', cpus=1, specified_node=None, array_concurrency_limit=None,
         submitted_utc=now(), convergence_policy='diagnostic_only'))
    print('New recovery chain released; prior blocked chain and its three monitors cancelled.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'submit'))
    parser.add_argument('--repo', type=Path, default=REPO)
    parser.add_argument('--revision')
    args = parser.parse_args(); _, _, retry = locations(args.repo)
    retry.mkdir(parents=True, exist_ok=True)
    with (retry / 'driver.lock').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.mode == 'prepare':
            prepare(args.repo, args.revision or output(['git', '-C', str(args.repo), 'rev-parse', 'HEAD']).strip())
        else: submit(args.repo)


if __name__ == '__main__':
    main()
