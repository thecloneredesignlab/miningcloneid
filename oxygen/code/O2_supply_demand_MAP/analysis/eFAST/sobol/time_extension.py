#!/usr/bin/env python3
"""Supervise the existing RED add-hours scripts for this recovery only (Python 3.6+)."""
import argparse
import contextlib
import datetime
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
TARGETS = ()


def configure():
    global TARGETS
    records = json.loads((ROOT / 'targets.json').read_text())
    targets = records['target_jobs']
    if not targets or len(set(targets)) != len(targets) or any(not isinstance(x, str) or not x.isdigit() for x in targets):
        raise ValueError('Invalid frozen monitor target registry')
    TARGETS = tuple(targets)


WALLTIME = '30-12:00:00'
HOURS = 6
INTERVAL = 300


def utc():
    return datetime.datetime.utcnow().isoformat() + 'Z'


def event(kind, **fields):
    row = dict(utc=utc(), kind=kind, **fields)
    message = json.dumps(row, sort_keys=True)
    print(message, flush=True)
    with (ROOT / 'events.jsonl').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        handle.write(message + '\n')


def command(args, **kwargs):
    return subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, **kwargs)


def seconds(value):
    days, clock = value.split('-', 1) if '-' in value else ('0', value)
    parts = [int(x) for x in clock.split(':')]
    if len(parts) > 3 or not parts:
        raise ValueError(value)
    return int(days) * 86400 + sum(x * 60 ** i for i, x in enumerate(reversed(parts)))


def queue(target, by_name=False):
    args = ['squeue', '-h', '-u', '4482173', '-o', '%i|%T']
    args += ['-n', 'sobol_add6_' + target] if by_name else ['-j', target]
    result = command(args, timeout=60)
    if result.returncode:
        if not by_name and 'Invalid job id specified' in result.stderr:
            return []
        raise RuntimeError(result.stderr or result.stdout)
    return [line.split('|') for line in result.stdout.splitlines() if line.strip()]


@contextlib.contextmanager
def lock(path):
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def sbatch_args(target, predecessor=None):
    args = ['sbatch', '--parsable', '--qos=small', '--ntasks=1',
            '--cpus-per-task=1', '--mem=1G', '--time=' + WALLTIME,
            '--signal=B:USR1@600', '--job-name=sobol_add6_' + target,
            '--chdir=' + str(ROOT),
            '--output=' + str(ROOT / 'logs' / (target + '_%j.out')),
            '--error=' + str(ROOT / 'logs' / (target + '_%j.err')),
            '--export=ALL,THRESHOLD_HOURS=6,CHECK_INTERVAL_SECONDS=300,'
            'DRY_RUN=FALSE,AUTO_CONFIRM=TRUE,WALLTIME=' + WALLTIME]
    if predecessor:
        args.append('--dependency=afterany:' + predecessor)
    return args + [str(ROOT / 'run_monitor.sbatch'), str(ROOT), target]


def submit(target, predecessor=None, dry=False):
    assert target in TARGETS
    with lock(ROOT / 'submission.lock'):
        existing = [r for r in queue(target, by_name=True) if r[0] != predecessor]
        if existing:
            if len(existing) != 1:
                raise RuntimeError('Duplicate monitors for {}: {}'.format(target, existing))
            event('reuse_monitor', target=target, monitor=existing[0][0])
            return existing[0][0]
        if not queue(target):
            event('target_finished_before_submit', target=target)
            return None
        args = sbatch_args(target, predecessor)
        if dry:
            event('dry_submission', target=target, command=args)
            return None
        result = command(args, timeout=60)
        if result.returncode:
            raise RuntimeError(result.stderr or result.stdout)
        job = result.stdout.strip().split(';')[0]
        if not job.isdigit():
            raise RuntimeError('Unrecognized sbatch receipt: ' + result.stdout)
        event('submitted', target=target, monitor=job,
              predecessor=predecessor, command=args)
        ledger = ROOT / 'submission_jobs.tsv'
        empty = not ledger.exists() or ledger.stat().st_size == 0
        with ledger.open('a') as handle:
            if empty:
                handle.write('target_job_id\tmonitor_job_id\tpredecessor\tsubmitted_utc\n')
            handle.write('{}\t{}\t{}\t{}\n'.format(target, job, predecessor or '', utc()))
        return job


def job_details(job):
    result = command(['scontrol', 'show', 'job', '-o', job], timeout=60)
    if result.returncode:
        if 'Invalid job id specified' in result.stderr:
            return None
        raise RuntimeError(result.stderr or result.stdout)
    rows = [dict(part.split('=', 1) for part in row.split() if '=' in part)
            for row in result.stdout.splitlines() if row.strip()]
    # For an array's first task, its internal numeric ID can equal the parent ID;
    # scontrol then returns the whole array. Select that exact internal task only.
    matches = [row for row in rows if row.get('JobId') == job]
    if len(matches) != 1:
        raise RuntimeError('Expected one matching internal job id: ' + job)
    return matches[0]


def extend(job, hours):
    """Called only when the original helper selects a running task; verify readback."""
    target = os.environ['RECOVERY_TIME_TARGET']
    if target not in TARGETS or not job.isdigit() or hours != '6':
        raise ValueError('Extension outside the authorized registered jobs or six-hour increment')
    before = job_details(job)
    if not before:
        event('extension_target_disappeared', target=target, job=job)
        return 0
    if before.get('ArrayJobId', before['JobId']) != target:
        raise ValueError('Internal job does not belong to target ' + target)
    if before.get('JobState') != 'RUNNING':
        event('extension_skipped_state', target=target, job=job, state=before.get('JobState'))
        return 0
    old = seconds(before['TimeLimit'])
    if old - seconds(before['RunTime']) >= HOURS * 3600:
        event('extension_skipped_threshold', target=target, job=job)
        return 0
    result = command(['/app/bin/s-add-job-hours', job, hours])
    print(result.stdout, end='', flush=True)
    print(result.stderr, end='', file=sys.stderr, flush=True)
    after = job_details(job)
    new = seconds(after['TimeLimit']) if after else None
    verified = new is not None and new >= old + HOURS * 3600
    ended = after is None or after.get('JobState') in (
        'COMPLETED', 'CANCELLED', 'FAILED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL')
    event('extension_result', target=target, job=job, returncode=result.returncode,
          before_seconds=old, after_seconds=new, verified=verified,
          target_ended=ended, stdout=result.stdout, stderr=result.stderr)
    return 0 if result.returncode == 0 and (verified or ended) else 1


def stop_child(child):
    if child is not None and child.poll() is None:
        os.killpg(child.pid, signal.SIGTERM)
        try:
            child.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()


def watch(target):
    current = os.environ['SLURM_JOB_ID']
    flags = {'renew': False, 'stop': False}
    signal.signal(signal.SIGUSR1, lambda *_: flags.update(renew=True))
    signal.signal(signal.SIGTERM, lambda *_: flags.update(stop=True))
    signal.signal(signal.SIGINT, lambda *_: flags.update(stop=True))
    lease = (ROOT / ('monitor_' + target + '.lock')).open('a')
    try:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        event('duplicate_monitor_exit', target=target, monitor=current)
        return 0
    env = dict(os.environ)
    env.update(THRESHOLD_HOURS='6', CHECK_INTERVAL_SECONDS='300',
               DRY_RUN='FALSE', AUTO_CONFIRM='TRUE', WALLTIME=WALLTIME,
               MAX_MISSING_CHECKS='3', RECOVERY_TIME_TARGET=target,
               ADD_HOURS_SCRIPT=str(ROOT / 'source_scripts' / 'add_hours_to_running_array_jobs.sh'))
    env['PATH'] = str(ROOT / 'bin') + ':' + env['PATH']
    fallback_end = time.time() + seconds(WALLTIME)
    end = float(env.get('SLURM_JOB_END_TIME') or fallback_end)
    child = None
    event('monitor_started', target=target, monitor=current,
          threshold_hours=HOURS, increment_hours=HOURS, interval_seconds=INTERVAL,
          renewal_end_epoch=end)
    try:
        while not flags['stop']:
            if flags['renew'] or time.time() >= end - 900:
                try:
                    successor = submit(target, predecessor=current)
                except Exception as exc:
                    event('renewal_retry', target=target, monitor=current, error=str(exc))
                    time.sleep(30)
                    continue
                event('handover', target=target, monitor=current, successor=successor)
                return 0
            if child is None or child.poll() is not None:
                if child is not None:
                    try:
                        if not queue(target):
                            event('target_finished', target=target, monitor=current)
                            return 0
                    except Exception as exc:
                        event('queue_check_retry', target=target, error=str(exc))
                        time.sleep(15)
                        continue
                    event('restart_existing_helper', target=target, returncode=child.returncode)
                    time.sleep(15)
                child = subprocess.Popen(
                    ['bash', str(ROOT / 'source_scripts' / 'submit_hourly_add_hours.sh'), target, '6'],
                    env=env, start_new_session=True)
            time.sleep(5)
        return 143
    finally:
        stop_child(child)
        lease.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('submit', 'dry-submit', 'watch', 'extend'))
    parser.add_argument('arguments', nargs='*')
    args = parser.parse_args()
    configure()
    if args.mode in ('submit', 'dry-submit'):
        for target in TARGETS:
            submit(target, dry=args.mode == 'dry-submit')
        return 0
    if args.mode == 'watch' and len(args.arguments) == 1 and args.arguments[0] in TARGETS:
        return watch(args.arguments[0])
    if args.mode == 'extend' and len(args.arguments) == 2:
        return extend(*args.arguments)
    parser.error('Invalid arguments')


if __name__ == '__main__':
    sys.exit(main())
