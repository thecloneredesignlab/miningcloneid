"""Resume the approved N20 local-stage infrastructure recovery under the original lock."""
import fcntl
from pathlib import Path
import shutil
import subprocess
import sys
import time

root=Path(sys.argv[1]).resolve()
sys.path.append(str(root/'release'))
from common import read, table, verify, sha
from controller import Owner, TERMINAL, RECOVERY_TASK, FAILED_JOB, FAILED_NODE

def resume(owner):
    if owner.terminal.read_text().strip()!='VALIDATION_HOLD':
        raise RuntimeError('The original continuation is not on the approved hold')
    owner.verify_recovery_bundle()
    verify(root)
    owner.refresh()
    failures=[j for j in owner.jobs if (j['state'],j['exit_code'])!=('COMPLETED','0:0')]
    if len(failures)!=1 or not all(j['state'] in TERMINAL for j in owner.jobs):
        raise RuntimeError('Unexpected active or failed job')
    orig,=failures
    if (orig['job_id'],orig['task'],orig['state'],orig['exit_code'],orig['settled_seconds'])!=(
        FAILED_JOB,RECOVERY_TASK,'FAILED','127:0','28'):
        raise RuntimeError('Archived failure differs from the diagnosed allocation')
    if (owner.selection.exists() or owner.reserve.exists() or
        any(j['attempt']!='1' for j in owner.jobs) or
        sum(int(j['settled_seconds']) for j in owner.jobs)!=81253):
        raise RuntimeError('Original scope or charged allocation count changed')
    first=owner.phase(RECOVERY_TASK)
    if (not all(owner.ready(f'N_{letter}_de6') for letter in 'ABCD') or
        (first/'console.log').read_text().strip() != 'FATAL:   user: unknown userid 107865' or
        (first/'evaluations.tsv').exists() or list(first.glob('trace.*'))):
        raise RuntimeError('Expected global checkpoints and zero-forward failure were not reproduced')
    gate,=read(first/'phase_exit_codes.tsv')
    if gate['source_after_exit_code']!='0':
        raise RuntimeError('The failed allocation did not preserve its original immutable source')
    for r in read(first/'phase_manifest.tsv'):
        f=first/r['path']
        if not f.is_file() or f.stat().st_size!=int(r['bytes']) or sha(f)!=r['sha256']:
            raise RuntimeError('Failed phase archive was modified')
    ledger,=read(root/'accounting/final_budget.tsv')
    if ledger['settled_core_seconds']!='81622' or ledger['attempts']!='293605':
        raise RuntimeError('Original terminal usage receipt changed')
    archive=root/'history'/'nss_hold_22065165'
    archive.mkdir(parents=True)
    for f in (owner.path,owner.control/'controller.log',owner.control/'controller_pid.txt',
              root/'REPORT.md',root/'source_manifest.tsv',owner.accounting/'candidate_attempts.tsv',
              owner.accounting/'usage_ledger.tsv',owner.accounting/'slurm_snapshot.txt'):
        shutil.copy2(f,archive/f.name)
    failed_phase=first.with_name('N_A_local__attempt1_22065165')
    first.rename(failed_phase)
    for f in (owner.hold,owner.terminal,owner.accounting/'final_budget.tsv'):
        f.rename(archive/f.name)
    table(owner.accounting/'infrastructure_retry.tsv',[dict(task=RECOVERY_TASK,prior_job=FAILED_JOB,
        prior_node=FAILED_NODE,prior_charged_forwards=0,prior_core_seconds=28,
        reason='APPTAINER_HOST_UID_LOOKUP_FAILURE_BEFORE_R',
        archived_phase=failed_phase.relative_to(root).as_posix(),
        archived_hold=archive.relative_to(root).as_posix(),
        start_from='N_A_de6',approved_infrastructure_retries=1,core_hour_cap='NONE')])
    owner.submit(RECOVERY_TASK,'infrastructure_retry',2)

def main():
    with (root/'controller/owner.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        owner=Owner(root)
        resume(owner)
        while not owner.terminal.exists():
            try:
                owner.tick()
            except subprocess.CalledProcessError as exc:
                if any(r['job_id']=='submission_pending' for r in owner.jobs):
                    owner.stop('SUBMISSION_RECONCILIATION_HOLD',str(exc))
                    break
                print('SCHEDULER_QUERY_RETRY',str(exc),flush=True)
            except (OSError,ValueError,RuntimeError,KeyError) as exc:
                if owner.hold.exists():
                    raise
                owner.stop('VALIDATION_HOLD',str(exc))
            if not owner.terminal.exists():
                time.sleep(8)

if __name__=='__main__':
    main()
