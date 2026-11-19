"""Resume the sealed D19 study after the diagnosed Slurm-ID propagation failure."""
import fcntl
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(sys.argv[1]).resolve()
PATCH = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'release'))
from common import read, sha, table, verify
import controller

PHASE = 'D_B_replay1b'
ORIGINAL = '22050987'
FIRST = '22050947'
PREDECESSOR = 'D_B_replay1a'
HISTORY = 'replay_env_hold_22050987'

def verify_patch():
    for row in read(PATCH/'manifest.tsv'):
        p = PATCH/row['path']
        if sha(p) != row['sha256'] or p.stat().st_size != int(row['bytes']):
            raise RuntimeError('Recovery release mismatch: '+str(p))

class Owner(controller.Owner):
    def submit(self, task, purpose, attempt=1):
        if attempt == 1 and self.job(task):
            return
        if attempt == 2 and (task != PHASE or purpose != 'retry' or
            self.job(task)['job_id'] != ORIGINAL or any(j['purpose']=='retry' for j in self.jobs)):
            raise RuntimeError('The recovery authorizes one retry of the diagnosed phase')
        verify(self.root)
        verify_patch()
        left = controller.ENVELOPES[purpose] - self.used(purpose)
        seconds = min(self.site, left)//60*60
        if seconds < 60:
            raise RuntimeError('No scheduler allocation remains for '+task)
        job = dict(task=task, attempt=str(attempt), purpose=purpose, job_id='submission_pending',
                   state='RESERVED', exit_code='', allocated_cpus='', elapsed_seconds='',
                   settled_seconds='0', reserved_seconds=str(seconds))
        self.jobs.append(job)
        self.save()
        self.budget()
        args = ['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1',
                '--cpus-per-task=1','--hint=nomultithread','--mem=16G','--no-requeue',
                '--time='+str(seconds//60),'--job-name=d19blind_'+task,
                '--output='+str(self.root/'logs'/(task+'_%j.txt')),
                str(PATCH/'launch.sh'),str(self.root),task]
        raw = subprocess.check_output(args,text=True).strip().split(';')[0]
        if not raw.isdigit():
            raise RuntimeError('Uncertain sbatch identifier. Retained reservation for reconciliation.')
        job.update(job_id=raw,state='SUBMITTED')
        self.save()
        print('SUBMITTED',task,raw,'attempt',attempt,flush=True)

def open_recovery(owner):
    history = ROOT/'history'/HISTORY
    if history.exists():
        raise RuntimeError('Recovery history already exists')
    if owner.terminal.read_text().strip()!='VALIDATION_HOLD':
        raise RuntimeError('Unexpected terminal marker')
    hold, = read(owner.hold)
    if PHASE not in hold['detail'] or ORIGINAL not in hold['detail']:
        raise RuntimeError('Unexpected halt reason')
    verify(ROOT)
    verify_patch()
    owner.refresh()
    if any(j['state'] not in controller.TERMINAL for j in owner.jobs):
        raise RuntimeError('An existing allocation remains live')
    retry = owner.job(PHASE)
    first = owner.job(PREDECESSOR)
    if (retry['job_id'],retry['state'],retry['exit_code'],first['job_id']) != (ORIGINAL,'FAILED','1:0',FIRST):
        raise RuntimeError('The failed or prerequisite allocation changed')
    if any(j['purpose']=='retry' for j in owner.jobs) or owner.reserve.exists() or owner.selection.exists():
        raise RuntimeError('Retry, reserve or result state changed')
    for job in owner.jobs:
        if job['state']=='COMPLETED' and not owner.ready(job['task']):
            raise RuntimeError('Previous successful phase is unsealed: '+job['task'])
    if not owner.ready(PREDECESSOR):
        raise RuntimeError('Predecessor must pass independent source and policy checks')
    failed = owner.phase(PHASE)
    checks = {r['check']:controller.true(r['passed']) for r in read(failed/'checks.tsv')}
    if checks != {'replay_1_score':True,'replay_1_native_parameters':True,
                   'replay_1_independent_components':True,'two_fresh_processes_1':False}:
        raise RuntimeError('Failure differs from diagnosed metadata defect')
    if read(failed/'policy_gate.tsv')[0]['status']!='PASS':
        raise RuntimeError('Failed run had a scientific read policy error')
    for row in read(failed/'phase_manifest.tsv'):
        p = failed/row['path']
        if p.stat().st_size!=int(row['bytes']) or sha(p)!=row['sha256']:
            raise RuntimeError('Failed allocation receipt changed: '+str(p))
    counts = controller.phase_buckets(failed/'evaluations.tsv')
    if counts['verify']!=1 or sum(counts.values())!=1:
        raise RuntimeError('Failed replay did not consume the single audited call')
    attempted = read(ROOT/'accounting/candidate_attempts.tsv')
    charged_verify, = [int(r['attempts']) for r in attempted if r['origin']=='D_B' and r['bucket']=='verify']
    # Three candidate keys per replay at most, then the second 150-direction poll and its two replays.
    if charged_verify+3+150+6 > controller.ATTEMPT_CAPS['verify']:
        raise RuntimeError('Original B verification cap cannot fit the complete remaining work')
    provenance, = read(PATCH/'diagnosis.tsv')
    for phase, key in ((PREDECESSOR,'first_point_sha256'),(PHASE,'failed_point_sha256')):
        path = owner.phase(phase)/('point_'+provenance['key']+'.rds')
        if sha(path)!=provenance[key]:
            raise RuntimeError('Diagnosed replay payload changed: '+phase)
    if sum(int(j['settled_seconds']) for j in owner.jobs) > controller.CAP:
        raise RuntimeError('Previously spent compute exceeds cap')

    history.mkdir(parents=True)
    for path in (owner.path, ROOT/'accounting/candidate_attempts.tsv', ROOT/'accounting/slurm_snapshot.txt',
                 ROOT/'controller/controller.log', ROOT/'controller/controller_pid.txt', ROOT/'REPORT.md'):
        shutil.copy2(path,history/path.name)
    prior = failed.with_name(PHASE+'__attempt1_'+ORIGINAL)
    failed.rename(prior)
    for path in (owner.hold,owner.terminal,ROOT/'accounting/final_budget.tsv'):
        path.rename(history/path.name)
    table(ROOT/'accounting/technical_retry.tsv', [dict(
        task=PHASE, previous_job=ORIGINAL, prior_first_job=FIRST, reason='APPTAINER_CLEANENV_REMOVED_SLURM_JOB_ID',
        failed_scientific_calls=1, charged_verify_at_hold=charged_verify,
        failed_phase=str(prior.relative_to(ROOT)), prior_history=str(history.relative_to(ROOT)),
        first_point_sha256=provenance['first_point_sha256'], failed_point_sha256=provenance['failed_point_sha256'],
        resume_from=PREDECESSOR, original_sources_preserved=True)])
    owner.submit(PHASE,'retry',2)

def main():
    with (ROOT/'controller/owner.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        owner = Owner(ROOT)
        open_recovery(owner)
        while not owner.terminal.exists():
            try:
                owner.tick()
            except subprocess.CalledProcessError as exc:
                if any(j['job_id']=='submission_pending' for j in owner.jobs):
                    owner.stop('SUBMISSION_RECONCILIATION_HOLD',str(exc))
                    break
                print('SCHEDULER_QUERY_RETRY',str(exc),flush=True)
            except (ValueError,OSError,RuntimeError,KeyError) as exc:
                if owner.hold.exists():
                    raise
                owner.stop('VALIDATION_HOLD',str(exc))
            if not owner.terminal.exists():
                time.sleep(8)

if __name__=='__main__':
    main()
