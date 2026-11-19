"""Finite N20 release controller. It handles metadata and never reads parameter values."""
import fcntl
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

from common import read, sha, table, verify

RECOVERY_TASK = "N_A_local"
FAILED_JOB = "22065165"
FAILED_NODE = "hpctpa3pc0027"

ORIGINS = tuple('N_'+x for x in 'ABCDEF')
PHASES = ('import', 'bench', 'de1b', 'de2', 'de3', 'de4', 'de5', 'de6', 'local',
          'pattern', 'verify1', 'replay1a', 'replay1b',
          'verify2', 'replay2a', 'replay2b')
TERMINAL = {'COMPLETED','FAILED','TIMEOUT','CANCELLED','NODE_FAIL','OUT_OF_MEMORY',
            'PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED','SPECIAL_EXIT'}
FIELDS = ('task','attempt','purpose','job_id','state','exit_code','allocated_cpus',
          'elapsed_seconds','settled_seconds','reserved_seconds')
HISTORICAL_CORE_SECONDS = 369
HISTORICAL_ENGINEERING_FORWARDS = 5
HANDOFFS = ('N_A_import', 'N_A_bench')
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
        assert sum(ATTEMPT_CAPS.values()) == 84000
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

    def accounting_snapshot(self):
        actual=sum(int(j['settled_seconds']) for j in self.jobs)
        requested=sum(int(j['reserved_seconds']) for j in self.jobs)
        table(self.accounting/'usage_ledger.tsv',[dict(
            historical_core_seconds=HISTORICAL_CORE_SECONDS,stage_core_seconds=actual,
            cumulative_core_seconds=HISTORICAL_CORE_SECONDS+actual,
            live_requested_core_seconds=requested,core_hour_cap='NONE')])

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
        self.accounting_snapshot()

    def ready(self, task):
        if task in HANDOFFS:
            return True
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

    def verify_recovery_bundle(self):
        recovery=self.root/'recovery_nss_22065165'
        for row in read(recovery/'manifest.tsv'):
            path=recovery/row['path']
            if path.stat().st_size!=int(row['bytes']) or sha(path)!=row['sha256']:
                raise RuntimeError('Infrastructure recovery manifest changed: '+str(path))

    def node_recovery_ready(self):
        j=self.job(RECOVERY_TASK)
        if j is None or j['attempt']!='2':
            return False
        path=self.phase(RECOVERY_TASK)/'node_preflight.tsv'
        if not path.is_file():
            if j['state']=='COMPLETED':
                raise RuntimeError('Completed retry has no verified node preflight')
            return False
        receipt,=read(path)
        if (receipt['status']!='PASS' or receipt['job_id']!=j['job_id'] or
            receipt['uid']!='107865' or receipt['host_user']!='4482173' or
            receipt['container_uid']!='107865' or receipt['node']==FAILED_NODE):
            raise RuntimeError('Node recovery preflight mismatch')
        return True

    def submit(self, task, purpose, attempt=1):
        if attempt == 1 and self.job(task):
            return
        if attempt != 1 and (attempt!=2 or task!=RECOVERY_TASK or purpose!='infrastructure_retry'
            or not self.job(task) or self.job(task)['job_id']!=FAILED_JOB
            or any(r['attempt']=='2' for r in self.jobs)):
            raise RuntimeError('Only the separately approved infrastructure retry may be submitted')
        verify(self.root)
        self.verify_recovery_bundle()
        # Slurm requires a reservation. Use the configured site limit.
        seconds = self.site // 60 * 60
        if seconds < 60:
            raise RuntimeError('Cannot determine a valid Slurm reservation')
        job = dict(task=task, attempt=str(attempt), purpose=purpose, job_id='submission_pending',
                   state='RESERVED', exit_code='', allocated_cpus='', elapsed_seconds='',
                   settled_seconds='0', reserved_seconds=str(seconds))
        self.jobs.append(job)
        self.save()
        self.accounting_snapshot()
        args = ['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1',
                '--cpus-per-task=1','--hint=nomultithread','--mem=16G','--no-requeue',
                '--time='+str(seconds//60),'--job-name=n20blind_'+task,
                '--output='+str(self.root/'logs'/(task+'_%j.txt'))]
        if attempt==2:
            args.append('--exclude='+FAILED_NODE)
        args.extend([str(self.root/'recovery_nss_22065165/launcher.sh'),str(self.root),task,str(attempt)])
        raw = subprocess.check_output(args, text=True).strip().split(';')[0]
        if not raw.isdigit():
            raise RuntimeError('Uncertain sbatch identifier. Retained reservation for reconciliation.')
        job.update(job_id=raw, state='SUBMITTED')
        self.save()
        print('SUBMITTED',task,raw,'attempt',attempt,flush=True)

    def handle_failures(self):
        for j in self.jobs:
            if self.job(j['task']) is not j:
                continue
            if j['state'] in TERMINAL and (j['state'],j['exit_code'])!=('COMPLETED','0:0'):
                raise RuntimeError(f"Phase {j['task']} job {j['job_id']} ended {j['state']} {j['exit_code']}.")
        return False

    def next_phase(self, task, purpose):
        if task in HANDOFFS:
            return True
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
        saved, = read(self.root/'continuity/benchmark_N.tsv')
        rows = read(self.phase('N_A_bench')/'callback_timings.tsv')
        c = benchmark(rows)
        if (saved['status']!='PASS' or saved['job_id']!='22059685' or int(saved['callbacks'])!=400 or
            c!=float(saved['admission_seconds_per_callback']) or
            phase_buckets(self.phase('N_A_bench')/'evaluations.tsv')['global']!=400):
            raise RuntimeError('A400 handoff did not reproduce the admitted checkpoint evidence')
        overhead=float(saved['overhead_seconds'])
        block=12200*c+overhead
        forecast=math.ceil(84000*c+3600+18*max(0,overhead-200))
        status='PASS' if block<self.site else 'SCHEDULER_COMPATIBILITY_HOLD'
        if not self.admission.exists():
            table(self.admission,[dict(status=status, benchmark_job_id=saved['job_id'], callbacks=400,
                forecast_core_seconds_per_origin=forecast, estimated_largest_block_seconds=block,
                site_max_seconds=self.site, core_hour_cap='NONE')])
            print('VERIFIED_A400_HANDOFF',status,forecast,flush=True)
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
        external = HISTORICAL_ENGINEERING_FORWARDS
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
        new = sum(int(j['settled_seconds']) for j in self.jobs)
        actual = HISTORICAL_CORE_SECONDS + new
        if calls>504400:
            raise RuntimeError('Registered forward-attempt allowance exceeded')
        table(self.accounting/'final_budget.tsv',[dict(status=status, settled_core_seconds=actual,
          settled_core_hours=actual/3600, historical_core_seconds=HISTORICAL_CORE_SECONDS,
          continuation_core_seconds=new, core_hour_cap='NONE', live_reserved_core_seconds=0,
          attempts=calls, max_attempts=504400,
          technical_retries=1+sum(j['purpose']=='infrastructure_retry' for j in self.jobs),
          infrastructure_retries=sum(j['purpose']=='infrastructure_retry' for j in self.jobs),
          supplemental_engineering_allocations=1,
          scientific_reserve_rounds=int(self.reserve.exists()), owned_jobs=len(self.jobs), active_owned_jobs=0)])
        lines = ['# N20 single-dataset pilot', '', f'**{status}**', '', detail, '',
                 f'Cumulative settled use was {actual/3600:.6f} core-hours, including the original and A400 receipts.',
                 f'The continuation used {new/3600:.6f} core-hours across {len(self.jobs)} xxlarge allocations.',
                 f'All charged forwards totaled {calls:,} of 504,400, including the archived engineering and A400 calls.', '']
        retry = next((j for j in self.jobs if j['purpose']=='infrastructure_retry'),None)
        if retry:
            lines.extend([f"Original infrastructure failure: {FAILED_JOB} on {FAILED_NODE}. "
                f"The separate retry was allocation {retry['job_id']} ({retry['state']}).",
                'The original 28 core-seconds and the retry allocation remain in cumulative accounting.', ''])
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
        if not self.node_recovery_ready():
            return
        if self.selection.exists():
            if not self.next_phase('P','assessment'):
                return
            record = read(self.selection)[0]
            self.close(record['status'],'Post-selection assessment and independent source, role and numerical gates passed.')
            return
        if not self.admission.exists():
            verify(self.root)
            if not self.cost_gate():
                self.stop('SCHEDULER_COMPATIBILITY_HOLD','Forecasted block duration exceeds the Slurm site limit.')
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
