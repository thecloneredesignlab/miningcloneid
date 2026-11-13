#!/usr/bin/env python3
"""Budgeted culture-only review DAG. No application concurrency throttle."""
import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
from datetime import datetime, timezone

SLURM = '/cm/shared/apps/slurm/current/bin/'
REL = 'oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability'
TERMINAL = ('accepted', 'failed', 'rejected', 'screen_failed')
P1 = ('pm_00', 'bb_00', 'pm_05', 'pm_06', 'bb_07', 'bb_08', 'pm_09')


def now():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.pending')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def commitment(reservations, kind=None):
    return sum(r.get('actual_core_hours', r['maximum_core_hours'])
               for r in reservations.values() if kind is None or r['kind'] == kind)


def recovered(values, local_improvement_required=False):
    if not values:
        return dict(pass_gate=False, reason='no_verified_candidate')
    best = min(values, key=lambda v: v['objective'])
    fresh = {v['task_id'] for v in values if v.get('pure_independent') and
             v.get('mode') == 'global' and v['objective'] - best['objective'] <= .001}
    local = [v for v in values if v.get('tested_directional_stability') and
             v['objective'] - best['objective'] <= .0001 and
             (not local_improvement_required or v.get('local_improvement', float('inf')) < .0001)]
    passed = len(fresh) >= 2 and bool(local)
    return dict(pass_gate=passed, reason='independent_and_local_checks_pass' if passed else
                'independent_recovery_unresolved' if len(fresh) < 2 else 'local_stability_unresolved',
                best_task=best['task_id'], objective=best['objective'],
                fresh_populations_near_minimum=len(fresh), stable_local_checks=len(local))


class Review:
    def __init__(self, project, study):
        self.project = Path(project).resolve(); self.study = Path(study).resolve()
        self.config = read(self.study / 'study.json'); self.state = read(self.study / 'state.json')
        self.parent = self.project / self.config['parent_relative_path']
        self.old = read(self.parent / 'state.json')
        self.wrapper = self.study / 'implementation' / REL / 'run_culture_review.sh'
        self.script = Path(__file__).resolve()
        self.started = time.time()

    def save(self):
        self.state['updated_utc'] = now()
        self.state['nonfit_reserved_core_hours'] = commitment(self.state['reservations'], 'nonfit')
        write(self.study / 'state.json', self.state)

    def event(self, kind, **details):
        self.state['events'].append(dict(kind=kind, utc=now(), **details))

    def verify(self):
        cfg = self.config
        if not cfg['authorized'] or cfg['qos'] != 'xxlarge' or cfg['core_hour_cap'] != 10000 or cfg['tumor_fitting']:
            raise RuntimeError('Wrong study authorization')
        approved = read(self.study / 'manifests/approved_inputs.json')
        if sha(self.study / 'study.json') != approved['study_sha256'] or sha(self.parent / 'state.json') != cfg['parent_state_sha256']:
            raise RuntimeError('Approved study or immutable parent changed')
        if sha(self.study / 'manifests/code_manifest.tsv') != approved['code_manifest_sha256']:
            raise RuntimeError('Approved implementation manifest changed')
        with (self.study / 'manifests/code_manifest.tsv').open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                if sha(self.study / 'implementation' / row['path']) != row['sha256']:
                    raise RuntimeError('Frozen implementation changed: ' + row['path'])
        registry = read(self.study / 'manifests/orchestrator_versions.json')['active']
        if sha(self.script) != registry['sha256']:
            raise RuntimeError('Use the registered controller')

    def submit(self, key, args, maximum_core_hours, kind):
        folder = self.study / 'jobs' / key
        receipt = folder / 'submission.json'; intent = folder / 'submission_intent.json'
        if receipt.exists():
            job = read(receipt)['job_id']
            self.state['reservations'][key]['job_id'] = job
            return job
        if intent.exists():
            raise RuntimeError('Ambiguous submission requires scheduler audit: ' + key)
        if commitment(self.state['reservations']) + maximum_core_hours > 10000 + 1e-8:
            raise RuntimeError('Total allocation ceiling reached')
        if kind == 'nonfit' and commitment(self.state['reservations'], kind) + maximum_core_hours > 120 + 1e-8:
            raise RuntimeError('Screening/validation/controller budget reached')
        self.state['reservations'][key] = dict(kind=kind, maximum_core_hours=maximum_core_hours,
                                              receipt=str(receipt.relative_to(self.study)))
        self.save()
        command = [SLURM + 'sbatch', '--parsable', '--partition=red', '--qos=xxlarge', '--ntasks=1', '--no-requeue'] + args
        write(intent, dict(arguments=command, utc=now(), maximum_core_hours=maximum_core_hours, kind=kind))
        p = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=60)
        if p.returncode:
            write(folder / 'submission_error.json', dict(stdout=p.stdout, stderr=p.stderr, returncode=p.returncode))
            raise RuntimeError('Scheduler rejected submission: ' + p.stderr)
        job = p.stdout.strip().split(';')[0]
        if not job.isdigit():
            raise RuntimeError('Ambiguous scheduler response')
        write(receipt, dict(job_id=job, arguments=command, maximum_core_hours=maximum_core_hours, kind=kind, utc=now()))
        self.state['reservations'][key]['job_id'] = job
        self.event('submitted', key=key, job_id=job); self.save()
        return job

    def controller(self, dependency, tag):
        key = 'controller_' + self.script.parent.name + '_' + tag
        if key in self.state['reservations']:
            return
        self.submit(key, ['--cpus-per-task=1', '--mem=8G', '--time=00:10:00',
            '--dependency=afterany:' + dependency, '--job-name=cr1110_reduce',
            '--output=' + str(self.study / 'logs' / (key + '_%j.out')),
            '--error=' + str(self.study / 'logs' / (key + '_%j.err')),
            '--wrap=python3 ' + str(self.script) + ' advance ' + str(self.project) + ' ' + str(self.study)], 1/6, 'nonfit')

    def accounting(self):
        ids = [v['job_id'] for v in self.state['reservations'].values() if v.get('job_id')]
        if not ids:
            return {}
        p = subprocess.run([SLURM + 'sacct', '-j', ','.join(ids), '-P',
            '--format=JobID,State,ExitCode,ElapsedRaw,AllocCPUS,MaxRSS,TotalCPU'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=60)
        if p.returncode:
            raise RuntimeError('Scheduler accounting unavailable')
        import io
        rows = list(csv.DictReader(io.StringIO(p.stdout), delimiter='|'))
        parents = {r['JobID']: r for r in rows if r['JobID'] in ids}
        for r in self.state['reservations'].values():
            a = parents.get(r.get('job_id'))
            if a and a['State'].split()[0] not in ('PENDING', 'RUNNING', 'COMPLETING', 'CONFIGURING', 'SUSPENDED'):
                r['actual_core_hours'] = int(a['ElapsedRaw']) * int(a['AllocCPUS']) / 3600
                r['scheduler_state'] = a['State']
        write(self.study / 'manifests/latest_scheduler_accounting.json', dict(utc=now(), records=rows))
        return parents

    def old_values(self, pid):
        ids = ([k for k, t in self.old['tasks'].items() if t['stage'] == 'baseline']
               if pid == 'base_original' else self.old['points'][pid]['task_ids'])
        values = []
        for tid in ids:
            if self.old['tasks'][tid]['status'] != 'accepted':
                continue
            r = read(self.parent / 'fits' / tid / 'validated_result.json')
            values.append(dict(task_id='parent_' + tid, objective=r['objective'], pure_independent=False,
                parent_independent_anchor=any(label.startswith('independent_') for label in self.old['tasks'][tid].get('starts_labels', [])),
                source_fit=str((self.parent / 'fits' / tid / 'accepted_fit.rds').relative_to(self.project)),
                params={k: float(v) for k, v in read(self.parent / 'fits' / tid / 'artifact_validation_v05.json')['numeric_params_exact'].items()}))
        return values

    def values(self, pid):
        point = self.state['points'][pid]
        values = self.old_values(point['parent_point']) if point['domain'] == 'original' and point.get('include_parent', True) else []
        for tid in point['task_ids']:
            if self.state['tasks'][tid]['status'] != 'accepted':
                continue
            result = read(self.study / 'fits' / tid / 'validated_result.json')
            result['source_fit'] = str((self.study / 'fits' / tid / 'accepted_fit.rds').relative_to(self.project))
            result['params'] = {k: float(v) for k, v in result['numeric_params_exact'].items()}
            values.append(result)
        return values

    def add_point(self, pid, domain='original', parent_point=None, priority='P1'):
        if pid in self.state['points']:
            return
        parent_point = parent_point or pid
        original = self.old['points'].get(parent_point)
        self.state['points'][pid] = dict(point_id=pid, parent_point=parent_point, domain=domain,
            fixed={original['target']: original['value']} if original else {}, priority=priority,
            task_ids=[], rounds=0)

    def create(self, pid, mode, source=None, pure=False, reserve=False, slot=None):
        if slot:
            existing = next((t['task_id'] for t in self.state['tasks'].values() if t['point_id'] == pid and t.get('slot') == slot), None)
            if existing:
                return existing
        count = mode + '_attempts'
        cap = self.config[mode + '_attempt_cap']
        reserve_count = 'reserve_' + count
        if self.state[count] >= cap or (reserve and self.state[reserve_count] >= self.config['reserve_' + mode + '_cap']):
            return None
        tid = 'review_{:03d}'.format(self.state['next_ordinal'])
        point = self.state['points'][pid]
        task = dict(task_id=tid, point_id=pid, domain=point['domain'], fixed=point['fixed'], mode=mode,
            priority=point['priority'], pure_independent=pure,
            start_class='fresh_independent_population' if pure else 'transferred_accepted_vector',
            source_fit=source['source_fit'] if source else None,
            source_fit_sha256=sha(self.project / source['source_fit']) if source else None,
            source_objective=source['objective'] if source else None,
            seed=self.config['master_seed'] + self.state['next_ordinal'] * 100003,
            min_generations=300, max_generations=1000, reserve=reserve, status='planned', slot=slot)
        screening_registry = self.study / 'manifests/screening_versions.json'
        if mode == 'global' and screening_registry.exists():
            screen = read(screening_registry)['active']
            validation = read(self.study / screen['validation_path'])
            if not validation['pass'] or validation['script_sha256'] != screen['script_sha256']:
                raise RuntimeError('Initialization-only proposal is not validated')
            for field in ('script', 'wrapper'):
                if sha(self.study / screen[field + '_path']) != screen[field + '_sha256']:
                    raise RuntimeError('Initialization source changed')
            task.update(screening_script=screen['script_path'],screening_script_sha256=screen['script_sha256'],
                screening_wrapper=screen['wrapper_path'],screening_wrapper_sha256=screen['wrapper_sha256'],
                screening_validation_sha256=sha(self.study / screen['validation_path']))
        write(self.study / 'tasks' / tid / 'task.json', {k:v for k,v in task.items() if k != 'status'})
        task['task_sha256'] = sha(self.study / 'tasks' / tid / 'task.json')
        self.state['tasks'][tid] = task; point['task_ids'].append(tid)
        self.state[count] += 1; self.state['next_ordinal'] += 1
        if reserve:
            self.state[reserve_count] += 1
        self.save(); self.launch(task)
        return tid

    def launch(self, task):
        tid = task['task_id']
        mode = 'screen' if task['mode'] == 'global' and task['status'] == 'planned' else task['mode']
        if mode == 'screen':
            cores, memory, wall, charge, kind = 1, '8G', '01:00:00', 1, 'nonfit'
        elif mode == 'global':
            cores, memory, wall, charge, kind = 8, '16G', '12:00:00', 96, 'global'
        else:
            cores, memory, wall, charge, kind = 1, '8G', '06:00:00', 6, 'local'
        key = tid + '_' + mode
        wrapper = self.wrapper
        if mode == 'screen' and task.get('screening_wrapper'):
            wrapper = self.study / task['screening_wrapper']
            if sha(wrapper) != task['screening_wrapper_sha256']:
                raise RuntimeError('Screening wrapper changed')
        job = self.submit(key, ['--cpus-per-task=' + str(cores), '--mem=' + memory, '--time=' + wall,
            '--job-name=cr1110_' + tid, '--output=' + str(self.study / 'logs' / (key + '_%j.out')),
            '--error=' + str(self.study / 'logs' / (key + '_%j.err')),
            str(wrapper), str(self.project), str(self.study), mode, tid], charge, kind)
        task['job_id'] = job; task['execution_mode'] = mode; task['status'] = mode + '_submitted'
        self.save(); self.controller(job, key)

    def refresh(self, accounting):
        more_validation = False
        for tid, task in list(self.state['tasks'].items()):
            if task['status'] in TERMINAL:
                continue
            if sha(self.study / 'tasks' / tid / 'task.json') != task['task_sha256']:
                raise RuntimeError('Task specification changed')
            if task.get('source_fit') and sha(self.project / task['source_fit']) != task['source_fit_sha256']:
                raise RuntimeError('Source vector payload changed')
            if task.get('screening_script') and sha(self.study / task['screening_script']) != task['screening_script_sha256']:
                raise RuntimeError('Screening proposal source changed')
            if task['status'] in ('planned', 'ready'):
                if not self.state.get('scientific_hold'):
                    self.launch(task)
                continue
            a = accounting.get(task.get('job_id'))
            if not a or a['State'].split()[0] in ('RUNNING', 'PENDING', 'COMPLETING', 'CONFIGURING'):
                continue
            task['scheduler'] = a
            if a['State'] != 'COMPLETED' or a['ExitCode'] != '0:0':
                task['status'] = 'screen_failed' if task['execution_mode'] == 'screen' else 'failed'
                task['reason'] = 'scheduler_' + a['State']; self.event('task_failed', task_id=tid, reason=task['reason'])
                if a['State'] == 'FAILED':
                    log = self.study / 'logs' / (tid + '_' + task['execution_mode'] + '_' + task['job_id'] + '.err')
                    expected_screen_failure = task['execution_mode'] == 'screen' and log.exists() and 'Independent feasible-bank screening unresolved' in log.read_text(errors='replace')
                    if not expected_screen_failure:
                        self.state['scientific_hold'] = 'Unexpected execution failure; inspect task ' + tid
                continue
            if task['execution_mode'] == 'screen':
                result = read(self.study / 'populations' / tid / 'screening_result.json')
                if not result['pass'] or result['feasible_members'] < 8:
                    task['status'] = 'screen_failed'; continue
                if sha(self.study / 'populations' / tid / 'population.rds') != result['population_sha256']:
                    raise RuntimeError('Screened population changed')
                if task.get('screening_script') and (result['script_sha256'] != task['screening_script_sha256'] or
                        sha(self.study / 'populations' / tid / 'member_origins.tsv') != result['member_origins_sha256']):
                    raise RuntimeError('Initialization provenance changed')
                task['status'] = 'ready'; self.save()
                if not self.state.get('scientific_hold'):
                    self.launch(task)
                continue
            if not os.getenv('SLURM_JOB_ID') or time.time() - self.started > 430:
                more_validation = True; continue
            out = self.study / 'fits' / tid
            if not (out / 'candidate_result.json').exists():
                task['status'] = 'rejected'; task['reason'] = 'missing_complete_candidate'; continue
            fatal = False
            for log in (out / 'evaluation_logs').glob('*.tsv'):
                with log.open() as stream:
                    for row in csv.DictReader(stream, delimiter='\t'):
                        if any(token in (row.get('error') or '') for token in
                               ('checkpoint_write_failed', 'Stale candidate', 'Immutable culture configuration', 'Wrong free-vector coordinates')):
                            fatal = True; break
                if fatal:
                    break
            if fatal:
                task['status'] = 'rejected'; task['reason'] = 'fatal_bookkeeping_error_in_evaluations'
                self.state['scientific_hold'] = 'Bookkeeping integrity failure; inspect task ' + tid
                continue
            run = subprocess.run(['bash', str(self.wrapper), str(self.project), str(self.study), 'check', tid],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=150)
            write(out / 'external_check_process.json', dict(returncode=run.returncode, stdout=run.stdout, stderr=run.stderr))
            if run.returncode:
                task['status'] = 'rejected'; task['reason'] = 'fresh_process_check_failed';
                self.state['scientific_hold'] = 'Fresh-process validation failed; dependent release stopped'; continue
            guard = read(out / 'artifact_validation.json'); result = read(out / 'candidate_result.json')
            if not guard['pass'] or guard['candidate_result_sha256'] != sha(out / 'candidate_result.json'):
                raise RuntimeError('Independent guard changed')
            for name, digest in result['files_sha256'].items():
                if sha(out / name) != digest:
                    raise RuntimeError('Accepted payload changed')
            result.update(accepted=True, artifact_validation_sha256=sha(out / 'artifact_validation.json'))
            if task['mode'] == 'local':
                result['local_improvement'] = max(0., result.get('local_improvement', task['source_objective'] - result['objective']))
            write(out / 'validated_result.json', result)
            task['status'] = 'accepted'; self.event('accepted', task_id=tid, objective=result['objective'])
            self.save()
        return more_validation

    def initial(self, pid, baseline=False, domain=False):
        point = self.state['points'][pid]
        if point.get('initial_registered'):
            return
        # The parent endpoints are verified feasible sources; new-domain transfer
        # vectors are rescored by the new optimizer and its external checker.
        sources = sorted(self.old_values(point['parent_point']), key=lambda x:x['objective'])
        primary = sources[0]
        alternative = next((x for x in sources[1:] if x.get('parent_independent_anchor') and x.get('params') != primary.get('params')),
                           next((x for x in sources[1:] if x.get('params') != primary.get('params')), sources[-1]))
        for index in range(2):
            self.create(pid, 'global', pure=True, slot='initial_independent_' + str(index))
        if baseline or domain:
            self.create(pid, 'global', source=primary, slot='initial_warm_1')
        if baseline and not domain:
            self.create(pid, 'global', source=alternative, slot='initial_warm_2')
        self.create(pid, 'local', source=primary, slot='initial_local_1')
        if not domain:
            self.create(pid, 'local', source=alternative, slot='initial_local_2')
        point['initial_registered'] = True; self.save()

    def closed(self, pid):
        return all(self.state['tasks'][tid]['status'] in TERMINAL for tid in self.state['points'][pid]['task_ids'])

    def escalate(self, pid, baseline=False, source=None):
        point = self.state['points'][pid]
        if point['rounds'] >= 2 or not self.closed(pid):
            return False
        values = self.values(pid)
        if not values:
            return False
        source = source or min(values, key=lambda x:x['objective'])
        gate = recovered(values, baseline)
        tasks = []
        fresh = gate.get('fresh_populations_near_minimum', 0)
        if source['objective'] < gate['objective'] - .001:
            fresh = 0
        for index in range(max(0, 2 - fresh)):
            tasks.append(self.create(pid, 'global', pure=True, reserve=True,
                slot='reserve_global_' + str(point['rounds']) + '_' + str(index)))
        if not gate.get('stable_local_checks') or source['objective'] < gate['objective'] - 1e-8:
            tasks.append(self.create(pid, 'local', source=source, reserve=True, slot='reserve_local_' + str(point['rounds'])))
        if any(tasks):
            point['rounds'] += 1; self.event('escalated', point_id=pid, round=point['rounds']); self.save(); return True
        return False

    def advance(self):
        self.verify(); self.state.pop('controller_error', None)
        accounting = self.accounting()
        validation = self.study / 'tests/validation_status.json'
        if not validation.exists() or not read(validation).get('integration_pass'):
            self.state['status'] = 'validation_pending'; self.save(); self.report(); return
        policy_path = self.study / 'tests/controller_policy_status.json'
        if not policy_path.exists() or not read(policy_path).get('pass') or read(policy_path).get('controller_sha256') != sha(self.script):
            self.state['status'] = 'controller_validation_pending'; self.save(); self.report(); return
        if read(validation)['code_manifest_sha256'] != sha(self.study / 'manifests/code_manifest.tsv'):
            raise RuntimeError('Validation fingerprints are stale')
        more = self.refresh(accounting)
        if self.state.get('scientific_hold'):
            self.state['status'] = 'validation_hold'; self.save(); self.report(); return
        self.add_point('base_original')
        for pid in P1:
            self.add_point(pid)
        self.save()
        self.initial('base_original', baseline=True)
        for pid in P1:
            self.initial(pid)
        # Repair any domain's unrestricted reference from a better restricted
        # complete vector; keep every negative delta until the repair succeeds.
        domains = sorted({p['domain'] for p in self.state['points'].values()})
        for domain in domains:
            bid = 'base_' + domain
            bvalues = self.values(bid)
            if not bvalues:
                continue
            bbest = min(bvalues, key=lambda x:x['objective'])
            constrained = [v for pid,p in self.state['points'].items() if p['domain'] == domain and p['fixed'] for v in self.values(pid)]
            better = min(constrained, key=lambda x:x['objective']) if constrained else None
            if better and better['objective'] < bbest['objective'] - 1e-8:
                if self.closed(bid):
                    if not self.escalate(bid, True, better):
                        self.state['scientific_hold'] = 'Unrestricted repair unresolved within reserve'; self.state['status'] = 'baseline_repair_unresolved'
                        self.save(); self.report(); return
                elif not self.state['points'][bid].get('early_repair_started'):
                    # The complete restricted vector is already available;
                    # fresh baseline globals need not finish before this local
                    # unrestricted repair starts. Do not duplicate globals.
                    point=self.state['points'][bid]
                    if point['rounds'] < 2 and self.create(bid,'local',source=better,reserve=True,slot='early_baseline_repair'):
                        point['early_repair_started']=True;point['rounds']+=1
                        self.event('early_baseline_repair',point_id=bid,source_task=better['task_id']);self.save()
        self.recover_failed_screening()
        b_gate = recovered(self.values('base_original'), True)
        self.state['baseline_gate'] = b_gate
        if self.closed('base_original') and not b_gate['pass_gate']:
            if not self.escalate('base_original', True):
                self.state['status'] = 'baseline_optimizer_unresolved'; self.save(); self.report(); return
        p1_closed = all(self.closed(pid) for pid in P1) and self.closed('base_original')
        if p1_closed and b_gate['pass_gate']:
            for pid in P1:
                if not recovered(self.values(pid))['pass_gate']:
                    self.escalate(pid)
            self.replacement_anchors()
            # P2 is released after the P1 review; remaining P1 flags stay visible.
            for pid in self.old['points']:
                if pid not in P1:
                    self.add_point(pid, priority='P2'); self.initial(pid)
            if all(recovered(self.values(pid))['pass_gate'] for pid in P1):
                self.release_domains()
            elif all(self.closed(pid) for pid in P1):
                self.state['domain_release'] = 'Deferred: high-priority independent/local recovery remains unresolved'
        if all(pid in self.state['points'] and self.closed(pid) for pid in self.old['points']):
            for pid in self.old['points']:
                if not recovered(self.values(pid))['pass_gate']:
                    self.escalate(pid)
        for pid,point in list(self.state['points'].items()):
            if point.get('replacement_anchor') and self.closed(pid) and not recovered(self.values(pid))['pass_gate']:
                self.escalate(pid)
            if point['domain'] != 'original' and self.closed(pid) and not recovered(self.values(pid), not point['fixed'])['pass_gate']:
                self.escalate(pid, not point['fixed'])
        active = any(t['status'] not in TERMINAL for t in self.state['tasks'].values())
        self.state['status'] = 'review_running' if active else 'review_fitting_finished'
        # Recheck accounting after potentially long external validations. A
        # finished task whose original trigger timed out waiting for the lock
        # must still get a fresh reconciliation job.
        latest = self.accounting()
        more = more or any(t['status'] not in TERMINAL and latest.get(t.get('job_id'),{}).get('State') == 'COMPLETED'
                           for t in self.state['tasks'].values())
        if more:
            self.controller(os.getenv('SLURM_JOB_ID', next(iter(accounting))), 'continuation_' + str(len(self.state['events'])))
        self.save(); self.report()

    def recover_failed_screening(self):
        # An independent population screen needs no local endpoint or finished
        # baseline global, once its replacement initializer has been validated.
        if not (self.study/'manifests/screening_versions.json').exists():
            return
        for pid in P1:
            point=self.state['points'][pid]
            if point['rounds'] or point.get('screen_recovery_registered'):
                continue
            failed=[self.state['tasks'][tid]for tid in point['task_ids']
                    if self.state['tasks'][tid]['status']=='screen_failed']
            if not failed:
                continue
            created=[]
            for index in range(min(2,len(failed))):
                created.append(self.create(pid,'global',pure=True,reserve=True,slot='screen_recovery_'+str(index)))
            if any(created):
                point['rounds']+=1;point['screen_recovery_registered']=True
                self.event('screen_recovery_submitted',point_id=pid,independent_populations=sum(bool(x)for x in created));self.save()

    def replacement_anchors(self):
        import math
        best = min(self.values('base_original'), key=lambda x:x['objective'])
        for target, upper, scale in [('p_misseg', .4, .015), ('buffer_beta', 10., 1.)]:
            pid = 'replacement_' + target
            if pid in self.state['points']:
                continue
            value = best['params'][target]
            grid = [p['value'] for p in self.old['points'].values() if p['target'] == target]
            distance = min(abs(math.log1p(value/scale)-math.log1p(v/scale)) for v in grid) / math.log1p(upper/scale)
            if distance <= .001:
                continue
            if self.state['reserve_global_attempts'] + 3 > 12 or self.state['reserve_local_attempts'] + 1 > 20:
                self.state.setdefault('unfunded_replacement_anchors', {})[target] = value
                continue
            self.state['points'][pid] = dict(point_id=pid, parent_point='base_original', domain='original',
                fixed={target:value}, priority='P1', task_ids=[], rounds=0, replacement_anchor=True,
                include_parent=False, initial_registered=True)
            self.save()
            for index in range(2):
                self.create(pid, 'global', pure=True, reserve=True, slot='anchor_independent_' + str(index))
            self.create(pid, 'global', source=best, reserve=True, slot='anchor_transferred')
            self.create(pid, 'local', source=best, reserve=True, slot='anchor_local')
            self.event('replacement_anchor_registered', point_id=pid, value=value, transformed_range_distance=distance)
            self.save()

    def release_domains(self):
        triggers = {}
        for domain, pids, parameter, upper in [('D_base', ('pm_00','bb_00'), 'p_mis_base', .005),
                                             ('D_sat', ('pm_08','pm_09'), 'k_o_mis', .02)]:
            if any(pid not in self.state['points'] or not self.closed(pid) for pid in pids):
                continue
            trigger = any(abs(min(self.values(pid), key=lambda x:x['objective'])['params'][parameter] - upper) <= 1e-8 for pid in pids)
            triggers[domain] = trigger
            if not trigger:
                continue
            bid = 'base_' + domain
            self.add_point(bid, domain, 'base_original', 'P3'); self.initial(bid, baseline=True, domain=True)
            for parent in pids:
                pid = domain + '_' + parent
                self.add_point(pid, domain, parent, 'P3'); self.initial(pid, domain=True)
        self.state['domain_triggers'] = triggers

    def report(self):
        points = []
        for pid,p in self.state['points'].items():
            values = self.values(pid); gate = recovered(values, not p['fixed'])
            points.append(dict(point_id=pid, domain=p['domain'], priority=p['priority'],
                attempts=len(p['task_ids']), active=sum(self.state['tasks'][tid]['status'] not in TERMINAL for tid in p['task_ids']), **gate))
        write(self.study / 'reports/progress.json', dict(status=self.state['status'], utc=now(),
            baseline_gate=self.state.get('baseline_gate'), task_status_counts={status:sum(t['status']==status for t in self.state['tasks'].values()) for status in sorted({t['status']for t in self.state['tasks'].values()})},
            global_attempts=self.state['global_attempts'], local_attempts=self.state['local_attempts'],
            total_allocated_or_committed_core_hours=commitment(self.state['reservations']),
            nonfit_allocated_or_committed_core_hours=commitment(self.state['reservations'],'nonfit'),
            scientific_hold=self.state.get('scientific_hold'), points=points))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['advance','status']); parser.add_argument('project'); parser.add_argument('study')
    args = parser.parse_args(); study = Path(args.study).resolve()
    if args.mode == 'status':
        print(json.dumps(read(study/'reports/progress.json'),indent=2)); return
    with (study/'controller.lock').open('a') as lock:
        # Short, bounded lock waits avoid charging many idle metadata allocations.
        deadline = time.time() + 45
        while True:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB); break
            except BlockingIOError:
                if time.time() >= deadline:
                    print('Another controller is advancing the same study; this trigger exits.'); return
                time.sleep(1)
        review = Review(args.project,args.study)
        try:
            review.advance()
        except Exception as error:
            review.state['controller_error'] = repr(error); review.save(); review.report(); raise


if __name__ == '__main__':
    main()
