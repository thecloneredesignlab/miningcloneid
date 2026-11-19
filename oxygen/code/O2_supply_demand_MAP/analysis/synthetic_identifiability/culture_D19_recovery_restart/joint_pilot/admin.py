"""Source, allocation and receipt checks for the finite D19 C/E pilot."""
import csv
import fcntl
import hashlib
import os
from pathlib import Path
import re
import subprocess
import sys
import time

MINUTES = {'C2': 55, 'E2': 55, 'P2': 25, 'V3': 20, 'V4': 20}
ATTEMPTS = {'C2': 800, 'E2': 798, 'P2': 16, 'V3': 4, 'V4': 4}
ORIGINAL_ATTEMPTS = {'T0': 100, 'C': 850, 'E': 850, 'P': 184, 'V1': 8, 'V2': 8}
CAP_SECONDS = 4 * 3600
JOB_FIELDS = ['task', 'job_id', 'reserved_core_seconds']
TERMINAL = {'COMPLETED', 'FAILED', 'TIMEOUT', 'CANCELLED', 'NODE_FAIL', 'OUT_OF_MEMORY', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED'}
WORKER_RDS = {'inputs/public/D19.rds', 'inputs/public/contract.rds',
              'verification/T00/identities/D_C.rds', 'verification/T00/identities/D_E.rds'}


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def table(path, data, fields=None):
    temporary = path.with_name('.pending_' + path.name)
    with temporary.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, delimiter='\t', fieldnames=fields or list(data[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(data)
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def manifest_row(base, path):
    return dict(path=str(path.relative_to(base)), bytes=path.stat().st_size, sha256=sha(path))


def verify_manifest(base, manifest):
    for record in rows(manifest):
        path = base / record['path']
        assert path.stat().st_size == int(record['bytes']) and sha(path) == record['sha256'], str(path)


def verify(archive, pilot):
    release = pilot / 'release_v2'
    assert (archive / 'controller/FINISHED').read_text().strip() == 'DET_OPTIMIZER_UNRESOLVED'
    verify_manifest(release, release / 'manifest.tsv')
    verify_manifest(archive, release / 'archive_inputs.tsv')
    verify_manifest(pilot / 'release', pilot / 'release/manifest.tsv')
    print('FROZEN_ARCHIVE_AND_PILOT_PASS')


def prepare(archive, pilot):
    release = pilot / 'release_v2'
    best = rows(archive / 'reports/best_found_gate_metrics.tsv')
    quality = rows(archive / 'reports/verified_candidates.tsv')
    index = {r['path']: r for r in rows(archive / 'local_collection_manifest.tsv')}
    sources = {'controller/FINISHED', 'reports/selection_lock.tsv', 'reports/verified_candidates.tsv',
               'release/production_bindings.R', 'release/interface.R', 'release/feasibility.R',
               'release/search_helpers.R', 'release/core.R', 'release/policy.py', 'release/common.py',
               'release/seal_phase.py'} | WORKER_RDS
    sources.update(p for p in index if p.startswith(('inputs/production/', 'runtime/')))
    origins = []
    for origin in ('D_C', 'D_E'):
        endpoint = next(r for r in best if r['origin'] == origin)
        record = next(r for r in quality if r['origin'] == origin and r['key'] == endpoint['key'])
        assert record['verified'].upper() == 'TRUE' and float(record['D']) == float(endpoint['D'])
        fields = ['origin', 'key', 'D', 'D_G', 'D_K', 'rms_growth', 'max_growth', 'max_kary_tv', 'state_file']
        origins.append({k: record[k] for k in fields})
        sources.add(record['state_file'])
    table(release / 'origins.tsv', origins)
    table(release / 'archive_inputs.tsv', [index[p] for p in sorted(sources)])
    table(release / 'manifest.tsv', [manifest_row(release, p) for p in sorted(release.iterdir()) if p.name != 'manifest.tsv'])
    verify(archive, pilot)


def audit_reads(archive, pilot, task):
    out = pilot / 'verification' / task
    origins = rows(pilot / 'release/origins.tsv')
    allowed = set(WORKER_RDS)
    arm = task[0] if task in ('C2', 'E2') else None
    if arm:
        allowed &= {p for p in allowed if not p.startswith('verification/') or p.endswith(f'D_{arm}.rds')}
    allowed.update(r['state_file'] for r in origins if not arm or r['origin'] == f'D_{arm}')
    forbidden = []
    archived_rds_reads = 0
    traces = list(out.glob('trace.*'))
    assert traces, 'Scientific worker traces are missing'
    for trace in traces:
        with trace.open(errors='replace') as stream:
            for line in stream:
                if not re.match(r'(open|openat)\(|<\.\.\. (open|openat) resumed>', line):
                    continue
                match = re.search(r'=\s*\d+<(/[^>]+)>\s*$', line)
                if not match or 'O_WRONLY' in line:
                    continue
                p = Path(match[1])
                if archive not in p.parents:
                    if '/oxygen/results/' in str(p) and pilot not in p.parents:
                        forbidden.append(str(p))
                    continue
                relative = str(p.relative_to(archive))
                if relative.startswith(('inputs/production/', 'runtime/', 'release/')):
                    continue
                if relative in allowed:
                    archived_rds_reads += 1
                else:
                    forbidden.append(relative)
    table(out / 'read_audit.tsv', [dict(status='PASS' if not forbidden and archived_rds_reads else 'FAIL',
          allowed_archive_payload_opens=archived_rds_reads, forbidden_reads=len(forbidden),
          findings='|'.join(sorted(set(forbidden))) or 'NONE')])
    assert not forbidden and archived_rds_reads, str(sorted(set(forbidden)))


def scheduler_snapshot(jobs):
    identifiers = {j['job_id'] for j in jobs}
    if not identifiers:
        return []
    assert all(k.isdigit() for k in identifiers), 'An unacknowledged submission requires reconciliation'
    raw = subprocess.check_output(['sacct', '-n', '-P', '-j', ','.join(sorted(identifiers)),
                                  '--format=JobIDRaw,State,ExitCode,AllocCPUS,ElapsedRaw,QOS'], text=True)
    by_id = {r[0]: r for line in raw.splitlines() if len(r := line.split('|')) == 6}
    result = []
    for job in jobs:
        key = job['job_id']
        root = by_id.get(key, [key, 'SUBMITTED', '', '0', '0', ''])
        state = root[1].split()[0].rstrip('+')
        settled = 0
        if state in TERMINAL:
            settled = max(int(r[3] or 0) * int(r[4] or 0) for k, r in by_id.items() if k == key or k.startswith(key + '.'))
        result.append(dict(task=job['task'], job_id=key, state=state, exit_code=root[2], qos=root[5],
                           elapsed_seconds=int(root[4] or 0), settled_core_seconds=settled,
                           live_reserved_core_seconds=0 if state in TERMINAL else int(job['reserved_core_seconds'])))
    return result


def audit_completed(pilot, job, caps):
    task = job['task']
    phase = pilot / 'verification' / task
    assert job['state'] == 'COMPLETED' and job['exit_code'] == '0:0' and job['qos'] == 'xxlarge', job
    assert (phase / 'phase_manifest.tsv').is_file(), f'{task}: phase receipt missing'
    assert all(value == '0' for value in rows(phase / 'phase_exit_codes.tsv')[0].values()), task
    assert all(rows(phase / name)[0]['status'] == 'PASS' for name in ('application_result.tsv', 'policy_gate.tsv', 'read_audit.tsv')), task
    app = rows(phase / 'application_result.tsv')[0]
    events = rows(phase / 'attempts.tsv')
    expected = int(app['forward_attempts'])
    assert int(app['cap']) == caps[task] and expected <= caps[task]
    assert [r['event'] for r in events] == ['START', 'END'] * expected, task
    assert [r['call'] for r in events[::2]] == [str(i) for i in range(1, expected + 1)], task
    assert all(a['key'] == b['key'] and a['call'] == b['call'] for a, b in zip(events[::2], events[1::2])), task
    assert app['N20'] == app['profiles'] == app['tumors'] == '0'
    return expected


def control(archive, pilot):
    marker = pilot / 'accounting/FINISHED'
    assert marker.read_text().strip() == 'DET_OPTIMIZER_UNRESOLVED'
    completed = scheduler_snapshot(rows(pilot / 'accounting/jobs.tsv'))
    prior = {r['task']: audit_completed(pilot, r, ORIGINAL_ATTEMPTS) for r in completed}
    assert set(prior) == set(ORIGINAL_ATTEMPTS)
    assert prior['C'] + ATTEMPTS['C2'] == 850 and prior['E'] + ATTEMPTS['E2'] == 850
    assert sum(prior[t] for t in ('P', 'V1', 'V2')) + sum(ATTEMPTS[t] for t in ('P2', 'V3', 'V4')) <= 200
    assert sum(prior.values()) + sum(ATTEMPTS.values()) <= 2000
    prior_seconds = sum(r['settled_core_seconds'] for r in completed)
    assert prior_seconds + sum(MINUTES.values()) * 60 <= CAP_SECONDS
    release = pilot / 'release_v2'
    acct = pilot / 'accounting_v2'
    ledger = acct / 'jobs.tsv'
    owner = (acct / 'owner.lock').open('a')
    fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (acct / 'FINISHED').exists() and not (acct / 'STOP').exists()
    qos = subprocess.check_output(['sacctmgr', '-nP', 'show', 'qos', 'xxlarge', 'format=Name,MaxWall'], text=True).strip()
    assert qos.startswith('xxlarge|')
    wall = qos.split('|')[1]
    hours, minutes, seconds = map(int, wall.split(':'))
    site_seconds = hours * 3600 + minutes * 60 + seconds
    assert site_seconds >= max(MINUTES.values()) * 60
    table(acct / 'site_qos.tsv', [dict(qos='xxlarge', max_seconds=site_seconds, pilot_cap_core_seconds=CAP_SECONDS)])
    jobs = rows(ledger) if ledger.exists() else []
    audited = {}

    def submit(task):
        assert task not in {r['task'] for r in jobs}
        verify(archive, pilot)
        record = dict(task=task, job_id='SUBMISSION_PENDING', reserved_core_seconds=MINUTES[task] * 60)
        jobs.append(record)
        table(ledger, jobs, JOB_FIELDS)
        output = subprocess.check_output(['sbatch', '--parsable', '--account=tao.li2', '--qos=xxlarge',
            '--ntasks=1', '--cpus-per-task=1', '--hint=nomultithread', '--mem=16G', '--no-requeue',
            '--time=' + str(MINUTES[task]), '--job-name=d19joint_' + task,
            '--output=' + str(pilot / 'logs' / (task + '_%j.txt')), str(release / 'runner.sh'),
            str(archive), str(pilot), task, str(release)], text=True).strip().split(';')[0]
        assert output.isdigit(), 'Ambiguous Slurm submission result'
        record['job_id'] = output
        table(ledger, jobs, JOB_FIELDS)
        print('SUBMITTED', task, output, flush=True)

    try:
        cost = rows(pilot / 'verification/T0/cost.tsv')[0]
        per_attempt = float(cost['seconds_per_evaluation'])
        setup = float(cost['setup_seconds_allowance_per_job'])
        forecasts = {k: setup + per_attempt * n for k, n in ATTEMPTS.items()}
        assert all(0 < s <= MINUTES[t] * 60 for t, s in forecasts.items()), forecasts
        assert prior_seconds + sum(forecasts.values()) <= CAP_SECONDS
        table(acct / 'admission.tsv', [dict(task=t, forecast_core_seconds=s, envelope_core_seconds=MINUTES[t] * 60, status='PASS') for t, s in forecasts.items()])
        for t in ('C2', 'E2'):
            if t not in {r['task'] for r in jobs}:
                submit(t)
        while True:
            status = scheduler_snapshot(jobs)
            table(acct / 'scheduler.tsv', status)
            settled = prior_seconds + sum(r['settled_core_seconds'] for r in status)
            reserved = sum(r['live_reserved_core_seconds'] for r in status)
            future = sum(m * 60 for t, m in MINUTES.items() if t not in {r['task'] for r in status})
            assert settled + reserved + future <= CAP_SECONDS
            table(acct / 'budget.tsv', [dict(cap_core_seconds=CAP_SECONDS, settled_core_seconds=settled,
                  live_reserved_core_seconds=reserved, future_envelope_core_seconds=future)])
            for row in status:
                if row['state'] in TERMINAL and row['task'] not in audited:
                    audited[row['task']] = audit_completed(pilot, row, ATTEMPTS)
            for required, next_task in [(('C2', 'E2'), 'P2'), (('P2',), 'V3'), (('V3',), 'V4')]:
                if all(t in audited for t in required) and next_task not in {r['task'] for r in jobs}:
                    submit(next_task)
            if len(audited) == len(ATTEMPTS):
                counts = dict(prior, **audited)
                caps = dict(ORIGINAL_ATTEMPTS, **ATTEMPTS)
                assert sum(counts.values()) <= 2000 and not reserved and not future
                final = rows(pilot / 'verification/V4/quality.tsv')
                qualified = {r['origin'] for r in final if r['prediction_pass'] == 'TRUE'}
                verdict = 'CONTINUATION_PREDICTION_GATE_MET' if qualified else 'DET_OPTIMIZER_UNRESOLVED'
                table(acct / 'final_attempts.tsv', [dict(task=t, forward_attempts=n, phase_cap=caps[t]) for t, n in counts.items()])
                table(acct / 'final_gate.tsv', [dict(status=verdict, verified_continuation_passes=len(qualified),
                      new_blind_origin_credit=0, jobs=len(counts), active_owned_jobs=0, total_forward_attempts=sum(counts.values()),
                      settled_core_seconds=settled, settled_core_hours=settled / 3600, N20=0, profiles=0, tumors=0)])
                (acct / 'FINISHED').write_text(verdict + '\n')
                print(verdict, flush=True)
                return
            time.sleep(20)
    except Exception as error:
        (acct / 'STOP').write_text(str(error) + '\n')
        raise


mode, source, destination, *extra = sys.argv[1:]
archive = Path(source).resolve()
pilot = Path(destination).resolve()
if mode == 'prepare':
    prepare(archive, pilot)
elif mode == 'verify':
    verify(archive, pilot)
elif mode == 'audit-reads':
    audit_reads(archive, pilot, extra[0])
elif mode == 'control':
    control(archive, pilot)
else:
    raise SystemExit('Unknown pilot administration command')
