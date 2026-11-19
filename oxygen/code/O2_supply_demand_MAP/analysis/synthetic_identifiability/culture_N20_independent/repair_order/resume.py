"""Reopen an isolated pilot after the diagnosed bijective observation-order mismatch."""
import fcntl
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(sys.argv[1]).resolve()
PATCH = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'release'))
from common import read, sha, table, verify
from controller import Owner, TERMINAL

FAILED = '22058909'
FIXED = ('n20_model.R','assess.R')

def reopen(owner):
    if owner.terminal.read_text().strip()!='VALIDATION_HOLD':
        raise RuntimeError('Unexpected terminal state')
    verify(ROOT)
    owner.refresh()
    failed = owner.job('E')
    if len(owner.jobs)!=3 or not all(j['state'] in TERMINAL for j in owner.jobs) or failed['job_id']!=FAILED:
        raise RuntimeError('Unexpected startup ledger or live allocation')
    pred = owner.phase('E')
    message = (pred/'console.log').read_text()
    if failed['state']!='FAILED' or 'identical(names(pred$growth), names(c$case$data$growth))' not in message:
        raise RuntimeError('Failure differs from the diagnosed order mismatch')
    if (owner.selection.exists() or owner.reserve.exists() or
        any(j['attempt']!='1' for j in owner.jobs) or not owner.ready('T00') or not owner.ready('G')):
        raise RuntimeError('Unexpected prior state')
    for row in read(pred/'phase_manifest.tsv'):
        item=pred/row['path']
        if item.stat().st_size!=int(row['bytes']) or sha(item)!=row['sha256']:
            raise RuntimeError('Failed attempt provenance changed')
    forward = sum(r['event']=='START' for r in read(pred/'forward_attempts.tsv'))
    if forward != 1 or any(r['event']=='END' for r in read(pred/'forward_attempts.tsv')):
        raise RuntimeError('Unexpected original attempt count')
    files = {r['path']:r for r in read(PATCH/'manifest.tsv')}
    for name,row in files.items():
        p=PATCH/name
        if sha(p)!=row['sha256'] or p.stat().st_size!=int(row['bytes']):
            raise RuntimeError('Sealed patch source mismatch')
    history = ROOT/'history/order_hold_22058909'
    history.mkdir(parents=True)
    for file in (owner.path, ROOT/'source_manifest.tsv', owner.control/'controller.log',
                 owner.control/'controller_pid.txt', ROOT/'REPORT.md',
                 owner.accounting/'candidate_attempts.tsv', owner.accounting/'slurm_snapshot.txt'):
        shutil.copy2(file,history/file.name)
    original = read(ROOT/'source_manifest.tsv')
    for name in FIXED:
        p=ROOT/'release'/name
        shutil.copy2(p,history/name)
        tmp=p.with_name(p.name+'.new')
        shutil.copy2(PATCH/name,tmp)
        tmp.replace(p)
    for r in original:
        item=Path(r['path'])
        if item.parent==Path('release') and item.name in FIXED:
            f=ROOT/item
            r.update(bytes=f.stat().st_size,sha256=sha(f))
    table(ROOT/'source_manifest.tsv',original)
    verify(ROOT)
    prior = pred.with_name('E__attempt1_'+FAILED)
    pred.rename(prior)
    for f in (owner.hold,owner.terminal,owner.accounting/'final_budget.tsv'):
        f.rename(history/f.name)
    table(owner.accounting/'technical_retry.tsv',[dict(task='E', previous_job=FAILED,
        reason='BIJECTIVE_SAMPLE_ID_ORDER_DIFFERENCE',previous_forward_attempts=forward,
        original_phase=str(prior.relative_to(ROOT)),old_manifest_sha=sha(history/'source_manifest.tsv'),
        new_manifest_sha=sha(ROOT/'source_manifest.tsv'),
        scope='Align production predictions and archived assessment reference by existing unique IDs')])
    owner.submit('E','retry',2)

def main():
    with (ROOT/'controller/owner.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        owner=Owner(ROOT)
        reopen(owner)
        while not owner.terminal.exists():
            try:owner.tick()
            except subprocess.CalledProcessError as exc:
                if any(r['job_id']=='submission_pending' for r in owner.jobs):
                    owner.stop('SUBMISSION_RECONCILIATION_HOLD',str(exc))
                    break
                print('SCHEDULER_QUERY_RETRY',str(exc),flush=True)
            except (ValueError,OSError,RuntimeError,KeyError) as exc:
                if owner.hold.exists():raise
                owner.stop('VALIDATION_HOLD',str(exc))
            if not owner.terminal.exists():time.sleep(8)

if __name__=='__main__':
    main()
