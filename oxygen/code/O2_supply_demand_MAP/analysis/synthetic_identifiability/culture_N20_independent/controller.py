"""Finite N20 release controller. It handles metadata and never reads parameter values."""
import fcntl
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys
import time

from common import read, sha, table, verify

ORIGINS = tuple('N_'+x for x in 'ABCDEF')
PHASES = ('import', 'bench', 'de1b', 'de2', 'de3', 'de4', 'de5', 'de6', 'local',
          'pattern', 'verify1', 'replay1a', 'replay1b',
          'verify2', 'replay2a', 'replay2b')
ENVELOPES = {**dict.fromkeys(ORIGINS, 100800), 'engineering':14400, 'assessment':14400, 'retry':28800}
TERMINAL = {'COMPLETED','FAILED','TIMEOUT','CANCELLED','NODE_FAIL','OUT_OF_MEMORY',
            'PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED','SPECIAL_EXIT'}
INFRASTRUCTURE_FAILURES = {'NODE_FAIL', 'BOOT_FAIL', 'PREEMPTED'}
FIELDS = ('task','attempt','purpose','job_id','state','exit_code','allocated_cpus',
          'elapsed_seconds','settled_seconds','reserved_seconds')
CAP = 662400
ATTEMPT_CAPS = dict(global_=73400, local=6000, pattern=2000, verify=1800, handoff=800)
ATTEMPT_CAPS['global'] = ATTEMPT_CAPS.pop('global_')

def true(value):
    return str(value).lower() == 'true'

def site_seconds(raw):
    days, sep, clock = raw.strip().rstrip('|').partition('-')
    if not sep:
        clock, days = days, '0'
    h, m, s = map(int, clock.split(':'))
    return int(days) * 86400 + h * 3600 + m * 60 + s

def benchmark(rows):
    if len(rows) != 400:
        raise ValueError('A cost gate requires exactly 400 complete callbacks')
    all_times = [float(r['seconds']) for r in rows]
    feasible = sorted(float(r['seconds']) for r in rows if true(r['feasible']))
    if len(feasible) < 200 or any(not math.isfinite(t) or t < 0 for t in all_times):
        raise ValueError('Invalid complete callback timings or imported feasibility')
    return max(1.0, 1.25 * statistics.mean(all_times), feasible[math.ceil(.95 * len(feasible))-1])

def phase_buckets(path):
    counts = dict.fromkeys(ATTEMPT_CAPS, 0)
    if path.is_file():
        for row in read(path):
            if row['event'] == 'START':
                counts[row['bucket']] += 1
    return counts

class Owner:
    def __init__(self, root):
        self.root = root.resolve()
        self.accounting = self.root / 'accounting'
        self.reports = self.root / 'reports'
        self.control = self.root / 'controller'
        self.path = self.accounting / 'jobs.tsv'
        self.jobs = read(self.path) if self.path.exists() else []
        self.admission = self.accounting / 'cost_admission.tsv'
        self.reserve = self.accounting / 'reserve_release.tsv'
        self.selection = self.reports / 'selection_lock.tsv'
        self.hold = self.accounting / 'stop_record.tsv'
        self.terminal = self.control / 'TERMINAL'
        self.cache = {}
        assert sum(ENVELOPES.values()) == CAP and sum(ATTEMPT_CAPS.values()) == 84000
        raw = subprocess.check_output(['sacctmgr','-nP','show','qos','xxlarge','format=MaxWall'], text=True)
        self.site = site_seconds(raw)
        table(self.accounting / 'site_limits.tsv', [dict(qos='xxlarge', max_seconds=self.site)])
        self.save()

    def save(self):
        table(self.path, self.jobs, FIELDS)

    def phase(self, task):
        return self.root / 'verification' / task

    def job(self, task):
        return next((j for j in reversed(self.jobs) if j['task'] == task), None)

    def used(self, purpose):
        return sum(int(j['settled_seconds']) + int(j['reserved_seconds']) for j in self.jobs if j['purpose'] == purpose)

    def budget(self):
        rows = []
        for purpose, cap in ENVELOPES.items():
            booked = self.used(purpose)
            if booked > cap:
                raise RuntimeError('Budget envelope exceeded: ' + purpose)
            rows.append(dict(purpose=purpose, cap_core_seconds=cap, booked_core_seconds=booked,
                             protected_remaining_core_seconds=cap-booked))
        table(self.accounting / 'budget_ledger.tsv', rows)

    def refresh(self):
        if not self.jobs:
            return
        ids = [j['job_id'] for j in self.jobs]
        if not all(s.isdigit() for s in ids):
            raise RuntimeError('Unreconciled submission reservation. No duplicate submission is permitted.')
        raw = subprocess.check_output(['sacct','-nP','-j',','.join(ids),
              '--format=JobIDRaw,State%40,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,QOS'],text=True)
        (self.accounting / 'slurm_snapshot.txt').write_text(raw)
        rows = {v[0]:v for line in raw.splitlines() if len(v:=line.split('|')) >= 7}
        for j in self.jobs:
            jid = j['job_id']
            if jid not in rows:
                continue
            r = rows[jid]
            if r[6] not in ('', 'xxlarge'):
                raise RuntimeError('Unexpected scheduler QOS: ' + jid)
            state = r[1].split()[0]
            j.update(state=state, exit_code=r[2], allocated_cpus=r[3], elapsed_seconds=r[4])
            if state in TERMINAL:
                actual = max(int(a[3] or 0)*int(a[4] or 0) for k,a in rows.items()
                             if k == jid or k.startswith(jid+'.'))
                j.update(settled_seconds=str(actual), reserved_seconds='0')
        self.save()
        self.budget()

    def ready(self, task):
        j = self.job(task)
        if not j or j['state'] != 'COMPLETED':
            return False
        phase = self.phase(task)
        manifest = phase / 'phase_manifest.tsv'
        if not manifest.is_file():
            return False
        h = sha(manifest)
        if j['exit_code'] != '0:0':
            raise RuntimeError('Nonzero Slurm exit: ' + task)
        if self.cache.get(task) == h:
            return True
        for r in read(manifest):
            path = phase / r['path']
            if not path.is_file() or path.stat().st_size != int(r['bytes']) or sha(path) != r['sha256']:
                raise RuntimeError('Phase integrity mismatch: ' + str(path))
        for filename in ('application_result.tsv','gate.tsv','policy_gate.tsv'):
            if not all(r['status']=='PASS' and ('job' not in r or r['job']==j['job_id']) for r in read(phase/filename)):
                raise RuntimeError('Failed phase gate: ' + task + ' ' + filename)
        if any(int(x) for x in read(phase/'phase_exit_codes.tsv')[0].values()):
            raise RuntimeError('Nonzero application, policy or source exit: ' + task)
        deny = phase/'read_denials.tsv'
        if task.startswith('N_') and deny.exists() and deny.read_text().strip():
            raise RuntimeError('Unexpected blind read denial: ' + task)
        self.cache[task] = h
        return True

    def submit(self, task, purpose, attempt=1):
        if attempt == 1 and self.job(task):
            return
        verify(self.root)
        left = ENVELOPES[purpose] - self.used(purpose)
        seconds = min(self.site, left) // 60 * 60
        if seconds < 60:
            raise RuntimeError('No scheduler allocation remains for ' + task)
        job = dict(task=task, attempt=str(attempt), purpose=purpose, job_id='submission_pending',
                   state='RESERVED', exit_code='', allocated_cpus='', elapsed_seconds='',
                   settled_seconds='0', reserved_seconds=str(seconds))
        self.jobs.append(job)
        self.save()
        self.budget()
        args = ['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1',
                '--cpus-per-task=1','--hint=nomultithread','--mem=16G','--no-requeue',
                '--time='+str(seconds//60),'--job-name=n20blind_'+task,
                '--output='+str(self.root/'logs'/(task+'_%j.txt')),
                str(self.root/'release/runner.sh'),str(self.root),task]
        raw = subprocess.check_output(args, text=True).strip().split(';')[0]
        if not raw.isdigit():
            raise RuntimeError('Uncertain sbatch identifier. Retained reservation for reconciliation.')
        job.update(job_id=raw, state='SUBMITTED')
        self.save()
        print('SUBMITTED',task,raw,'attempt',attempt,flush=True)

    def handle_failures(self):
        for task in dict.fromkeys(j['task'] for j in self.jobs):
            j = self.job(task)
            if j['state'] not in TERMINAL or (j['state']=='COMPLETED' and j['exit_code']=='0:0'):
                continue
            phase = self.phase(task)
            no_calls = (not sum(phase_buckets(phase/'evaluations.tsv').values()) and
                        not (phase/'forward_attempts.tsv').exists())
            used_retry = any(r['purpose']=='retry' for r in self.jobs)
            if j['state'] in INFRASTRUCTURE_FAILURES and not used_retry and no_calls:
                if task.startswith('N_'):
                    index = PHASES.index(task[4:])
                    if index and not self.ready(task[:4]+PHASES[index-1]):
                        raise RuntimeError('Retry has no verified parent checkpoint')
                    if task[4:].startswith('de') and task[4:]!='de1b':
                        if not self.admission.exists() or float(read(self.admission)[0]['seconds_per_complete_callback'])*12200 + 120 > ENVELOPES['retry']:
                            raise RuntimeError('Technical retry envelope cannot fit a full DE block')
                if phase.exists():
                    phase.rename(self.phase(f"{task}__attempt1_{j['job_id']}"))
                table(self.accounting/'technical_retry.tsv',[dict(task=task, previous_job=j['job_id'],
                    reason=j['state'], first_attempt_scientific_calls=0, restart_from='verified phase parent')])
                self.submit(task,'retry',2)
                return True
            raise RuntimeError(f"Phase {task} job {j['job_id']} ended {j['state']} {j['exit_code']}.")
        return False

    def next_phase(self, task, purpose):
        if not self.job(task):
            self.submit(task,purpose)
            return False
        return self.ready(task)

    def chain(self, origin):
        for stage in PHASES:
            if not self.next_phase(origin+'_'+stage,origin):
                return False
        return True

    def cost_gate(self):
        if self.admission.exists():
            return read(self.admission)[0]['status']=='PASS'
        rows = read(self.phase('N_A_bench')/'callback_timings.tsv')
        c = benchmark(rows)
        forecast = math.ceil(84000*c + 3600)
        bench_job = self.job('N_A_bench')
        rss = [r[5] for line in (self.accounting/'slurm_snapshot.txt').read_text().splitlines()
               if len(r:=line.split('|')) >= 7 and r[0].split('.')[0]==bench_job['job_id'] and r[5]]
        if not rss:
            raise RuntimeError('Benchmark peak memory accounting is missing')
        elapsed = int(bench_job['settled_seconds'])
        timing_sum = sum(float(r['seconds']) for r in rows)
        setup = max(0, elapsed - timing_sum)
        if setup > 200:
            forecast += math.ceil(18*(setup-200))
        status = 'PASS' if forecast <= ENVELOPES['N_A'] and c*12200+setup < self.site else 'RESOURCE_ADMISSION_HOLD'
        table(self.admission,[dict(status=status, seconds_per_complete_callback=c,
              forecast_core_seconds_per_origin=forecast, benchmark_job_id=bench_job['job_id'],
              benchmark_allocated_core_seconds=elapsed, callbacks=len(rows),
              feasible_callbacks=sum(true(r['feasible']) for r in rows),
              measured_overhead_seconds=setup, peak_memory=';'.join(rss))])
        print('COST_GATE',status,forecast,flush=True)
        return status=='PASS'

    def classify(self, origins):
        rows, counts = [], []
        for origin in origins:
            phase = self.phase(origin+'_replay2b')
            quality = read(phase/'quality.tsv')
            if not quality or any(r['origin'] != origin or not true(r['verified']) for r in quality):
                raise RuntimeError('Foreign or unverified final quality record')
            rows.extend(quality)
            b, = read(phase/'best_identity.tsv')
            verified = min(quality,key=lambda r:float(r['Q']))
            counts.append(dict(origin=origin,best_found_key=b['key'],best_found_Q=b['Q'],
                best_verified_key=verified['key'],best_verified_Q=verified['Q']))
        raw = min(counts,key=lambda r:float(r['best_found_Q']))
        selected = min(rows,key=lambda r:float(r['Q']))
        exact = next((r for r in rows if r['key']==raw['best_found_key'] and r['origin']==raw['origin']),None)
        agreeing = sum(float(r['best_verified_Q'])-float(selected['Q']) <= 1e-4 for r in counts)
        local = bool(exact and true(exact['local_pass']))
        qualified = agreeing >= 2 and local and abs(float(raw['best_found_Q'])-float(selected['Q']))<=1e-8
        return rows,counts,raw,selected,agreeing,local,qualified

    def lock_selection(self, origins):
        quality,counts,raw,selected,agreeing,local,qualified = self.classify(origins)
        status = 'N20_BOUNDED_OPTIMIZER_QUALIFIED' if qualified else 'NOISY_ESTIMATE_PROVISIONAL_OPTIMIZER_UNRESOLVED'
        table(self.reports/'origin_verdicts.tsv',counts)
        table(self.reports/'optimizer_quality.tsv',quality)
        table(self.selection,[dict(status=status,independent_origins_within_1e_4=agreeing,
            best_found_origin=raw['origin'],best_found_key=raw['best_found_key'],best_found_Q=raw['best_found_Q'],
            best_found_local_pass=local,selected_origin=selected['origin'],selected_key=selected['key'],
            selected_Q=selected['Q'],state_file=selected['state_file'],profiles=0,tumor=0,repeated_noisy=0)])
        (self.control/'SELECTION_LOCKED').write_text(status+'\n')
        print('SELECTION_LOCKED',status,flush=True)

    def tally(self):
        attempts = []
        for origin in ORIGINS:
            totals = dict.fromkeys(ATTEMPT_CAPS, 0)
            for phase in self.root.joinpath('verification').glob(origin+'_*'):
                for bucket, n in phase_buckets(phase/'evaluations.tsv').items():
                    totals[bucket] += n
            for key, n in totals.items():
                if n > ATTEMPT_CAPS[key]:
                    raise RuntimeError('Origin attempt cap exceeded: ' + origin + '/' + key)
                attempts.append(dict(origin=origin,bucket=key,attempts=n,cap=ATTEMPT_CAPS[key]))
        external = 0
        for name in self.root.joinpath('verification').glob('E*'):
            path = name/'forward_attempts.tsv'
            if path.is_file():
                external += sum(r['event']=='START' for r in read(path))
        if external > 400:
            raise RuntimeError('External engineering attempt cap exceeded')
        attempts.append(dict(origin='engineering',bucket='external',attempts=external,cap=400))
        table(self.accounting/'candidate_attempts.tsv',attempts)
        return sum(r['attempts'] for r in attempts)

    def close(self, status, detail):
        if any(j['state'] not in TERMINAL for j in self.jobs):
            return
        calls = self.tally()
        actual = sum(int(j['settled_seconds']) for j in self.jobs)
        if actual > CAP:
            raise RuntimeError('Settlement exceeds approved cap')
        table(self.accounting/'final_budget.tsv',[dict(status=status, settled_core_seconds=actual,
          settled_core_hours=actual/3600, cap_core_seconds=CAP, live_reserved_core_seconds=0,
          attempts=calls, max_attempts=504400, technical_retries=sum(j['purpose']=='retry' for j in self.jobs),
          scientific_reserve_rounds=int(self.reserve.exists()), owned_jobs=len(self.jobs), active_owned_jobs=0)])
        lines = ['# N20 single-dataset pilot', '', f'**{status}**', '', detail, '',
                 f'Settled allocated use was {actual/3600:.6f} of 184 core-hours across {len(self.jobs)} xxlarge jobs.',
                 f'All charged forwards totaled {calls:,} of 504,400.', '']
        if self.selection.exists():
            selection = read(self.selection)[0]
            lines.extend([f"Independent origins within 1e-4 Q: {selection['independent_origins_within_1e_4']}.",
                    f"Lowest found Q: {selection['best_found_Q']} ({selection['best_found_origin']}).",
                    f"Local certificate at the exact lowest point: {selection['best_found_local_pass']}.",
                    f"Selected independently verified Q: {selection['selected_Q']} ({selection['selected_origin']}).",
                    'The post-selection assessor describes growth, karyotype, parameter, kernel and response recovery.', ''])
        lines.extend(['The pilot uses one frozen noisy growth-plus-karyotype dataset and fits all 20 coordinates.',
                      'The result concerns numerical agreement and descriptive recovery for this dataset.',
                      'Repeated synthetic recovery, profiles and tumor fitting require separate releases.', ''])
        (self.root/'REPORT.md').write_text('\n'.join(lines))
        self.terminal.write_text(status+'\n')
        print('TERMINAL',status,actual,calls,flush=True)

    def tick(self):
        self.refresh()
        if self.hold.exists():
            row = read(self.hold)[0]
            self.close(row['status'],row['detail'])
            return
        if self.handle_failures():
            return
        if self.selection.exists():
            if not self.next_phase('P','assessment'):
                return
            record = read(self.selection)[0]
            self.close(record['status'],'Post-selection assessment and independent source, role and numerical gates passed.')
            return
        for task in ('T00','G','E'):
            if not self.next_phase(task,'engineering'):
                return
        for stage in ('import','bench'):
            if not self.next_phase('N_A_'+stage,'D_A'):
                return
        if not self.cost_gate():
            self.stop('RESOURCE_ADMISSION_HOLD','The measured complete-work forecast exceeds a registered allowance.')
            return
        complete = [self.chain(origin) for origin in ORIGINS[:4]]
        if not all(complete):
            return
        if self.classify(ORIGINS[:4])[-1]:
            self.lock_selection(ORIGINS[:4])
            return
        if not self.reserve.exists():
            table(self.reserve,[dict(round=1,families='N_E,N_F',reason='Primary score agreement or exact local quality unresolved')])
        if all([self.chain(origin) for origin in ORIGINS[4:]]):
            self.lock_selection(ORIGINS)

    def stop(self, status, detail):
        table(self.hold,[dict(status=status,detail=detail)])
        print(status,detail,flush=True)
        self.close(status,detail)

def main():
    root = Path(sys.argv[1])
    with (root/'controller/owner.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        owner = Owner(root)
        while not owner.terminal.exists():
            try:
                owner.tick()
            except subprocess.CalledProcessError as exc:
                if any(j['job_id']=='submission_pending' for j in owner.jobs):
                    owner.stop('SUBMISSION_RECONCILIATION_HOLD',str(exc))
                    break
                print('SCHEDULER_QUERY_RETRY',str(exc),flush=True)
            except (ValueError, OSError, RuntimeError, KeyError) as exc:
                if owner.hold.exists():
                    raise
                owner.stop('VALIDATION_HOLD',str(exc))
            if not owner.terminal.exists():
                time.sleep(8)

if __name__=='__main__':
    main()
