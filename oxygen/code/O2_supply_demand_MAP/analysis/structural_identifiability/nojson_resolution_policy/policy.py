#!/usr/bin/env python3
"""Resolution-aware task policy over an unchanged frozen numerical release."""
from pathlib import Path
import argparse
import contextlib
import math
import os
import shutil
import subprocess
import sys
import time
import uuid

sys.dont_write_bytecode = True
POLICY_ID = 'culture_resolution_aware_v1'
INTEGRITY = 1e-8
RESOLUTION = .001
REPRESENTATIVE = 1e-4
FINE_SHAPE = .002
W = L = None


def initialize(study, release=None):
    global W, L, RELEASE, BASE_CANDIDATES, OLD_GATES, OLD_EMIT, OLD_SUBMIT
    study = Path(study).resolve()
    sys.path.insert(0, str(study / 'implementation'))
    import workflow
    W = workflow
    L = W.L
    if Path(W.__file__).resolve().parent != study / 'implementation':
        raise ValueError('Wrong frozen workflow imported')
    RELEASE = Path(release or __file__).resolve()
    if RELEASE.is_file():
        RELEASE = RELEASE.parent
    BASE_CANDIDATES = W.candidates
    OLD_GATES, OLD_EMIT, OLD_SUBMIT = W.gates, W.emit, L.submit


def verify_release(study):
    L.verify_code(study)
    for row in L.read_table(RELEASE / 'source_manifest.tsv'):
        if L.sha(RELEASE / row['path']) != row['sha256']:
            raise ValueError('Policy release fingerprint mismatch')
    c = L.read_kv(RELEASE / 'controls.tsv')
    assert c == dict(policy_id=POLICY_ID, integrity=INTEGRITY, resolution=RESOLUTION,
                    representative=REPRESENTATIVE, fine_shape=FINE_SHAPE,
                    total_cap=4000, nonfit_cap=120, evidence_cap=6,
                    global_attempt_cap=36, local_attempt_cap=70, rounds=2,
                    qos='xxlarge', language='English')


def flag(value):
    if isinstance(value, bool):
        return value
    return {'True': True, 'TRUE': True, 'true': True,
            'False': False, 'FALSE': False, 'false': False}[value]


def repair_trigger(anchor_q, candidate_q):
    if not all(math.isfinite(q) for q in (anchor_q, candidate_q)):
        raise ValueError('Nonfinite reference decision')
    return anchor_q - candidate_q > RESOLUTION


def endpoint_evidence(directory, result):
    """Pre-screen records; the vector and scores still require fresh R verification."""
    if result['local_stable'] or result['budget_limited']:
        return None
    sweeps = L.read_table(directory / 'local_sweeps.tsv')
    if not sweeps or not all(flag(r['complete']) for r in sweeps):
        return None
    for h in (.001, .0001):
        fine = [r for r in sweeps if float(r['mesh']) == h]
        if not fine:
            return None
        r = fine[-1]
        if float(r['maximum_tested_improvement']) > 1e-4 or float(r['improvement']) >= 1e-5:
            return None
    method = L.read_kv(directory / 'method_sequence.tsv')
    if method['total_calls'] > 20000:
        return None
    q = float(sweeps[-1]['after'])
    probes = directory / 'boundary_probes.tsv'
    if probes.exists():
        for r in L.read_table(probes):
            if abs(float(r['endpoint']) - q) > INTEGRITY:
                raise ValueError('Wrong endpoint in inward evidence')
            if q - float(r['probe']) > 1e-5:
                return None
    gap = q - result['objective']
    if gap < -INTEGRITY:
        raise ValueError('Endpoint is below the accepted minimum')
    return q if gap <= REPRESENTATIVE else None


def compatible_evidence_manifests(study):
    allowed={L.sha(RELEASE/'source_manifest.tsv')}
    path=RELEASE/'evidence_compatibility.tsv'
    if not path.exists():return allowed
    for r in L.read_table(path):
        old=study/r['release_path']
        if L.sha(old/'source_manifest.tsv')!=r['sha256']:
            raise ValueError('Changed compatible evidence release')
        for item in L.read_table(old/'source_manifest.tsv'):
            if L.sha(old/item['path'])!=item['sha256']:
                raise ValueError('Changed compatible release payload')
        for name in ['evidence.R','evidence.sh','controls.tsv']:
            if L.sha(old/name)!=L.sha(RELEASE/name):
                raise ValueError('Numerical evidence compatibility changed')
        allowed.add(r['sha256'])
    return allowed


def evidence_results(data, study):
    rows = []
    allowed=compatible_evidence_manifests(study)
    for j in data['jobs']:
        if j['mode'] not in {'policy_validate', 'policy_evidence'} or not L.job_ok(j):
            continue
        folder = study / 'policy/evidence' / j['task']
        if not (folder / 'COMPLETE').exists():
            raise ValueError('Missing policy evidence completion')
        receipt = L.read_kv(folder / 'validation.tsv')
        if (not receipt['pass'] or receipt['policy_sha256'] not in allowed
                or receipt['requests_sha256'] != L.sha(study/'policy/requests'/(j['task']+'.tsv'))
                or receipt['payload_manifest_sha256'] != L.sha(folder/'payload_manifest.tsv')):
            raise ValueError('Stale policy evidence release')
        for item in L.read_table(folder / 'payload_manifest.tsv'):
            if L.sha(folder / item['path']) != item['sha256']:
                raise ValueError('Changed policy evidence payload')
        for r in L.read_table(folder / 'results.tsv'):
            assert r['policy_id'] == POLICY_ID
            source = study.parents[3] / r['source_fit']
            if L.sha(source) != r['source_sha256']:
                raise ValueError('Changed policy evidence source')
            r = dict(r)
            for k in ['objective', 'source_objective', 'endpoint_gap']:
                r[k] = float(r[k])
            r['local_evidence'] = flag(r['local_evidence'])
            r['batch'] = j['task']
            r['evidence_job_id'] = j['job_id']
            r['evidence_directory'] = str(folder.relative_to(study))
            rows.append(r)
    return rows


def candidates(data, study, project):
    cs = BASE_CANDIDATES(data, study, project)
    for r in evidence_results(data, study):
        if r['action'] == 'free_map':
            original = next(z for z in cs if z['task_id'] == r['task_id'])
            cs.append(dict(task_id=r['task_id'] + '_free_verified', point='base',
                           origin=original['origin'], family=original['family'],
                           mode='mapped', phase='policy_free_mapping',
                           objective=r['objective'], local_stable=False,
                           alternative=False, start_objective=r['source_objective'],
                           source_fit=r['verified_fit']))
    return sorted(cs, key=lambda r: (r['objective'], r['task_id']))


def witnesses(data, study, cs):
    """Never transfer an endpoint flag to a different checkpoint vector."""
    reps = []
    for r in cs:
        if r['mode'] in {'local', 'import'} and r['local_stable']:
            reps.append(dict(r, witness_objective=r['objective'],
                             witness_fit=r['source_fit'], evidence='original_accepted_local_flag'))
    index = {r['task_id']: r for r in cs}
    for ev in evidence_results(data, study):
        if ev['action'] == 'endpoint' and ev['local_evidence']:
            r = index[ev['task_id']]
            assert r['mode'] == 'local'
            reps.append(dict(r, witness_objective=ev['objective'], witness_fit=ev['verified_fit'],
                             evidence='fresh_verified_endpoint', evidence_job_id=ev['evidence_job_id']))
    return reps


def root_witnesses(data, study, reps, minimum, point):
    roots = {}
    for r in sorted(reps, key=lambda x: (x['witness_objective'], x['task_id'])):
        if r['point'] != point or r['family'] not in {'A', 'B'}:
            continue
        if r['witness_objective'] - minimum > RESOLUTION:
            continue
        if point != 'base' and L.accepted(data, study, r['origin'] + '_global') is None:
            continue
        if r['origin'] not in roots:
            roots[r['origin']] = r
    return list(roots.values())


def anchor(data, study):
    records = [r for r in data['events'] if r['event'] == 'resolution_anchor_installed']
    if not records:
        raise ValueError('No completed reference anchor')
    path = study / records[-1]['detail']
    receipt = L.read_kv(path / 'anchor.tsv')
    for r in L.read_table(path / 'manifest.tsv'):
        if L.sha(path / r['path']) != r['sha256']:
            raise ValueError('Changed reference anchor evidence')
    return receipt


def has_families(roots):
    return len(roots) >= 2 and {r['family'] for r in roots} == {'A', 'B'}


def gates(data, study, project):
    cs = candidates(data, study, project)
    reps = witnesses(data, study, cs)
    b = W.best(cs, 'base'); q = b['objective']; a = anchor(data, study)
    near = root_witnesses(data, study, reps, q, 'base')
    local = sorted([r for r in reps if r['point'] == 'base' and r['witness_objective'] - q <= REPRESENTATIVE],
                   key=lambda r:(r['witness_objective'],r['task_id']))
    repeat = [r for r in local if r['mode'] == 'local' and r.get('alternative', False)
              and r['start_objective'] - r['witness_objective'] < 1e-4]
    bg = dict(objective=q, task_id=b['task_id'], families=';'.join(sorted({r['family'] for r in near})),
              recovered_origins=len(near), repeat_pass=bool(repeat),
              representative_pass=bool(local), policy_id=POLICY_ID,
              anchor_objective=a['objective'], anchor_task_id=a['task_id'],
              cumulative_gain=a['objective'] - q, repair_trigger=repair_trigger(a['objective'], q),
              lowest_candidate_local_stable=b['local_stable'])
    bg['pass'] = has_families(near) and bool(local) and bool(repeat)
    bg['representative_task'] = local[0]['task_id'] if local else 'none'
    bg['representative_objective'] = local[0]['witness_objective'] if local else None
    bg['representative_gap'] = local[0]['witness_objective'] - q if local else None
    rows = []
    for c in L.read_table(study / 'cases.tsv'):
        point = c['point']; best = W.best(cs, point); minimum = best['objective']
        roots = root_witnesses(data, study, reps, minimum, point)
        raw = [r for r in cs if r['point'] == point and r['mode'] == 'local'
               and r['family'] in {'A', 'B'} and r['objective'] - minimum <= RESOLUTION
               and L.accepted(data, study, r['origin'] + '_global') is not None]
        local = sorted([r for r in reps if r['point'] == point
                        and r['witness_objective'] - minimum <= REPRESENTATIVE],
                       key=lambda x: (x['witness_objective'], x['task_id']))
        repeat = [r for r in local if r['mode'] == 'local' and r.get('alternative', False)
                  and (r['phase'] == 'repeat' or r['phase'].startswith('reserve_'))
                  and r['start_objective'] - r['witness_objective'] < 1e-4]
        reference = bg['pass'] and not bg['repair_trigger'] and minimum >= q - INTEGRITY
        rows.append(dict(point=point, target=c['target'], value=float(c['value']), objective=minimum,
                         delta=minimum-q, best_task=best['task_id'], independent_pass=has_families(roots),
                         recovered_origins=len(roots), families=';'.join(sorted({r['family'] for r in roots})),
                         raw_near_families=';'.join(sorted({r['family'] for r in raw})),
                         local_pass=bool(local), repeat_pass=bool(repeat), reference_pass=reference,
                         pass_gate=has_families(roots) and bool(local) and bool(repeat) and reference,
                         policy_id=POLICY_ID, delta_vs_anchor=minimum-a['objective'],
                         reference_candidate_gap=a['objective']-q,
                         lowest_candidate_local_stable=best['local_stable'],
                         representative_task=local[0]['task_id'] if local else 'none',
                         representative_objective=local[0]['witness_objective'] if local else None,
                         representative_gap=local[0]['witness_objective']-minimum if local else None,
                         representative_source_fit=local[0]['witness_fit'] if local else 'none',
                         shape_below_resolution=minimum-q <= FINE_SHAPE))
    return bg, rows, cs


@contextlib.contextmanager
def legacy_candidates():
    current = W.candidates
    W.candidates = BASE_CANDIDATES
    try:
        yield
    finally:
        W.candidates = current


def legacy_gates(data, study, project):
    with legacy_candidates():
        return OLD_GATES(data, study, project)


def save_anchor(state, data, project, bg, cs, proof_generation, reason, proof_data=None):
    assert bg['pass']
    proof_data = data if proof_data is None else proof_data
    reps = witnesses(proof_data, state.study, cs)
    roots = root_witnesses(proof_data, state.study, reps, bg['objective'], 'base')
    assert has_families(roots)
    repeats = [r for r in reps if r['point']=='base' and r['mode']=='local'
               and r.get('alternative',False) and r['witness_objective']-bg['objective'] <= REPRESENTATIVE
               and r['start_objective']-r['witness_objective'] < 1e-4]
    assert repeats
    best = W.best(cs, 'base')
    folder = state.study / 'policy/anchors' / ('anchor_' + uuid.uuid4().hex)
    folder.mkdir(parents=True)
    L.write_kv(folder / 'anchor.tsv', dict(policy_id=POLICY_ID, objective=bg['objective'],
               task_id=best['task_id'], source_fit=best['source_fit'],
               source_sha256=L.sha(project / best['source_fit']), utc=L.now(),
               proof_generation=proof_generation, reason=reason,
               repeat_task=repeats[0]['task_id'], recovery_band=RESOLUTION))
    fields = ['task_id','origin','family','witness_objective','witness_fit','evidence']
    L.write_table(folder / 'witnesses.tsv', [{k:r[k] for k in fields} for r in roots], fields)
    L.write_table(folder / 'manifest.tsv', [dict(path=p.name,sha256=L.sha(p)) for p in sorted(folder.iterdir())], ['path','sha256'])
    state.event(data, 'resolution_anchor_installed', str(folder.relative_to(state.study)))
    state.save(data)


def repair_baseline(state, data, project):
    bg, points, cs = gates(data, state.study, project)
    deficit = not bg['pass'] or bg['repair_trigger']
    n = W.round_number(state.study)
    if not deficit:
        # Complete an early intervention, without repeatedly moving a quiet anchor.
        if n is not None and (state.study/'reserve'/str(n)/'decision.tsv').exists():
            decision=L.read_kv(state.study/'reserve'/str(n)/'decision.tsv')
            if decision['reason']=='baseline_priority':
                active=[L.read_kv(state.study/'tasks'/(r['task_id']+'.tsv')) for r in data['tasks']
                        if r['phase']=='reserve_'+str(n) and r['point']=='base']
                complete=[W.chain(state,data,project,t) if t['mode']=='global' else W.ready(state,data,project,t) for t in active]
                if not all(complete):
                    return None
                save_anchor(state,data,project,bg,cs,(state.study/'state/CURRENT').read_text().strip(),'completed_intervention')
                L.atomic_text(state.study/'reserve'/str(n)/'CLOSED',L.now()+' baseline repaired; early intervention consumes round\n')
                state.event(data,'resolution_round_closed',str(n));state.save(data)
        return True
    if n is None:
        return False
    phase='reserve_'+str(n);W.open_round(state,data,n,'baseline_priority')
    active=[L.read_kv(state.study/'tasks'/(r['task_id']+'.tsv')) for r in data['tasks'] if r['phase']==phase and r['point']=='base']
    pending=False
    for t in active:
        done=W.chain(state,data,project,t) if t['mode']=='global' else W.ready(state,data,project,t)
        pending=not done or pending
    if pending:
        return None
    bg,points,cs=gates(data,state.study,project)
    if bg['pass']:
        save_anchor(state,data,project,bg,cs,(state.study/'state/CURRENT').read_text().strip(),'completed_reference_checkpoint')
        L.atomic_text(state.study/'reserve'/str(n)/'CLOSED',L.now()+' baseline repaired; early intervention consumes round\n')
        state.event(data,'resolution_round_closed',str(n));state.save(data)
        return True
    if W.counts(data,phase)[1]>=8:
        L.atomic_text(state.study/'reserve'/str(n)/'CLOSED',L.now()+' baseline unresolved at round cap\n')
        return None
    if not bg['repeat_pass'] or not bg['representative_pass']:
        r=W.best(cs,'base');tid=W.next_local_id(data,phase)
        t=W.task(state.study,tid,'local',phase,'base',r['origin'],r['family'],r['source_fit'],True)
        W.add(state,data,t);W.ready(state,data,project,t);return None
    missing={'A','B'}-set(bg['families'].split(';'))
    # Missing local evidence has priority over a new global origin. A completed
    # global plus an unstable near-minimum local may use remaining shared locals.
    eligible=[r for r in cs if r['point']=='base' and r['mode']=='local'
              and r['family'] in missing and not r['local_stable']
              and r['objective']-bg['objective']<=RESOLUTION]
    if eligible:
        r=min(eligible,key=lambda x:(x['objective'],x['task_id']))
        tid=W.next_local_id(data,phase)
        t=W.task(state.study,tid,'local',phase,'base',r['origin'],r['family'],r['source_fit'],True)
        state.event(data,'resolution_missing_family_local_repair',r['task_id']+' '+tid)
        W.add(state,data,t);W.ready(state,data,project,t);return None
    launched=False
    for family in sorted(missing):
        tid='base_R'+str(n)+'_'+family+'_global'
        if any(r['task_id']==tid for r in data['tasks']):
            continue
        if W.counts(data,phase)[0]>=4:
            break
        t=W.task(state.study,tid,'global',phase,'base','base_R'+str(n)+'_'+family,family)
        W.add(state,data,t);W.chain(state,data,project,t);launched=True
    if not launched:
        L.atomic_text(state.study/'reserve'/str(n)/'CLOSED',L.now()+' baseline independent recovery unresolved\n')
        state.event(data,'resolution_round_closed_unresolved',str(n));state.save(data)
    return None


def evidence_commitment(data):
    return sum(float(j['actual_ch']) if L.is_closed(j) and j['actual_ch'] else float(j['reserved_ch'])
               for j in data['jobs'] if j['phase'] in {'policy_validation','policy_evidence'})


def policy_submit(state,data,project,tid,mode,role,phase):
    if mode not in {'controller','policy_validate','policy_evidence'}:
        existed=L.get_job(data,tid,mode) is not None
        result=OLD_SUBMIT(state,data,project,tid,mode,role,phase)
        if not existed:
            state.event(data,'resolution_dispatch',result['token']+' '+POLICY_ID+' '+L.sha(RELEASE/'source_manifest.tsv'))
            state.save(data)
        return result
    old=L.get_job(data,tid,mode)
    if old:
        return old
    cpus,seconds,mem=L.resources(role)
    request=L.budget_check(data,role,cpus,seconds)
    if mode in {'policy_validate','policy_evidence'} and evidence_commitment(data)+request>6+1e-9:
        raise ValueError('Policy evidence allocation ceiling reached')
    token='cpr_'+uuid.uuid4().hex[:18]
    row=dict(slot=uuid.uuid4().hex,task=tid,mode=mode,role=role,phase=phase,job_id='',token=token,
             state='INTENT',exit_code='',cpus=str(cpus),wall_seconds=str(seconds),reserved_ch=format(request,'.17g'),actual_ch='',submitted_utc=L.now())
    data['jobs'].append(row);state.event(data,'submit_intent',token);state.save(data)
    command=[str(L.SLURM/'sbatch'),'--parsable','--no-requeue','--partition=red','--qos=xxlarge',
             '--cpus-per-task='+str(cpus),'--mem='+str(mem)+'G',
             '--time=%02d:%02d:00'%(seconds//3600,(seconds%3600)//60),'--job-name='+token,
             '--output='+str(state.study/'logs')+'/'+token+'_%j.out',
             '--error='+str(state.study/'logs')+'/'+token+'_%j.err',
             str(RELEASE/'batch.sh'),str(project),str(state.study),str(RELEASE),mode,tid]
    L.write_table(state.study/'submissions'/(token+'_arguments.tsv'),[dict(order=i,argument=a) for i,a in enumerate(command)],['order','argument'])
    try:
        p=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,timeout=55)
        L.atomic_text(state.study/'submissions'/(token+'_stdout.txt'),p.stdout)
        L.atomic_text(state.study/'submissions'/(token+'_stderr.txt'),p.stderr)
        text=p.stdout.strip().split(';')[0]
        if p.returncode or not text.isdigit():
            raise ValueError('Ambiguous policy submission; reconcile token')
        row['job_id']=text;row['state']='SUBMITTED'
        state.event(data,'submission_receipt',token+' '+text)
        state.event(data,'resolution_dispatch',token+' '+POLICY_ID+' '+L.sha(RELEASE/'source_manifest.tsv'))
        state.save(data)
    except BaseException:
        row['state']='AMBIGUOUS';state.save(data);raise
    return row


REQUEST_FIELDS=['action','task_id','source_fit','source_sha256','source_objective','source_result_sha256']


def request_batch(state,data,project,items,validation=False):
    tid=('policy_validation_' if validation else 'policy_evidence_')+uuid.uuid4().hex[:12]
    rows=[]
    for action,r in items:
        assert L.accepted(data,state.study,r['task_id']) is not None
        rows.append(dict(action=action,task_id=r['task_id'],source_fit=r['source_fit'],
                         source_sha256=L.sha(project/r['source_fit']),source_objective=format(r['objective'],'.17g'),
                         source_result_sha256=L.sha(state.study/'fits'/r['task_id']/'result.tsv')))
    L.write_table(state.study/'policy/requests'/(tid+'.tsv'),rows,REQUEST_FIELDS)
    state.event(data,'policy_evidence_request',tid+' '+L.sha(state.study/'policy/requests'/(tid+'.tsv'))+' '+L.sha(RELEASE/'source_manifest.tsv'))
    state.save(data)
    policy_submit(state,data,project,tid,'policy_validate' if validation else 'policy_evidence','check',
                  'policy_validation' if validation else 'policy_evidence')
    return tid


def prepare_evidence(state,data,project):
    for j in data['jobs']:
        if j['mode'] not in {'policy_validate','policy_evidence'}:
            continue
        if not L.is_closed(j):
            return False
        if not L.job_ok(j):
            raise ValueError('Policy evidence job failed; no numerical retry')
    existing={(r['action'],r['task_id']) for r in evidence_results(data,state.study)}
    cs=candidates(data,state.study,project);free=W.best(cs,'base')['objective'];items=[]
    for r in cs:
        if r['mode']=='local' and not r['local_stable'] and ('endpoint',r['task_id']) not in existing:
            if endpoint_evidence(state.study/'fits'/r['task_id'],r) is not None:
                items.append(('endpoint',r))
    for c in L.read_table(state.study/'cases.tsv'):
        r=W.best(cs,c['point'])
        if r['objective']<free-INTEGRITY and ('free_map',r['task_id']) not in existing:
            items.append(('free_map',r))
    if items:
        request_batch(state,data,project,items);return False
    return True


def emit(state,data,project,status):
    result=OLD_EMIT(state,data,project,status)
    old,points,_=legacy_gates(data,state.study,project)
    L.write_kv(state.study/'policy/tables/original_baseline_gate_current.tsv',old)
    L.write_table(state.study/'policy/tables/original_point_quality_current.tsv',points,list(points[0]))
    bg,ps,cs=result
    L.write_kv(state.study/'policy/tables/decision.tsv',dict(policy_id=POLICY_ID,utc=L.now(),
               anchor_task_id=bg['anchor_task_id'],anchor_objective=bg['anchor_objective'],
               lowest_free_task=bg['task_id'],lowest_free_objective=bg['objective'],
               cumulative_gain=bg['cumulative_gain'],repair_trigger=bg['repair_trigger'],
               independent_recovery=bg['recovered_origins']>=2 and bg['families']=='A;B',
               representative_pass=bg['representative_pass'],repeat_pass=bg['repeat_pass'],
               policy_evidence_commitment=evidence_commitment(data)))
    reps=witnesses(data,state.study,cs)
    fields=['task_id','point','origin','family','witness_objective','witness_fit','evidence']
    L.write_table(state.study/'policy/tables/checked_representatives.tsv',[{k:r[k] for k in fields} for r in reps],fields)
    return result


def install():
    W.candidates=candidates;W.gates=gates;W.repair_baseline=repair_baseline;W.emit=emit;L.submit=policy_submit


def assert_owner(data,study):
    active=[r for r in data['events'] if r['event']=='resolution_policy_activated']
    if not active or active[-1]['detail']!=str(RELEASE.relative_to(study)):
        raise ValueError('Inactive policy dispatcher')
    fingerprints=[r['detail'] for r in data['events'] if r['event']=='resolution_policy_fingerprint']
    if not fingerprints or fingerprints[-1]!=L.sha(RELEASE/'source_manifest.tsv'):
        raise ValueError('Active policy manifest changed')
    live=[j for j in data['jobs'] if j['role']=='controller' and not L.is_closed(j)]
    own=os.environ.get('SLURM_JOB_ID')
    if len(live)==1 and live[0]['job_id']==own:
        return True
    renew=[r for r in data['events'] if r['event']=='resolution_controller_renewal']
    if len(live)==2 and any(r['detail'].split()[-1]==own for r in renew):
        return False
    raise ValueError('Controller ownership is not unique')


def controller(project,study):
    install();state=L.State(study);start=time.monotonic()
    while True:
        with state.locked():
            data=state.load()
            try:
                L.refresh(state,data);verify_release(study)
                if not assert_owner(data,study):
                    status='waiting_for_predecessor_controller_exit'
                else:
                    # Accept closed unchanged payloads even while reference dispatch is held.
                    for row in list(data['tasks']):
                        if row['mode'] not in {'global','local'} or row['phase']=='fixture':
                            continue
                        j=L.get_job(data,row['task_id'],row['mode'])
                        if j is not None and L.is_closed(j) and L.job_ok(j):
                            try:
                                W.ready(state,data,project,L.read_kv(study/'tasks'/(row['task_id']+'.tsv')))
                            except ValueError as exc:
                                if not (study/'tables/resource_stop.tsv').exists() or not W.resource_limit(exc):raise
                    stopped=(study/'tables/resource_stop.tsv').exists()
                    status=W.drive(state,data,project) if stopped or prepare_evidence(state,data,project) else 'policy_evidence'
                    if status.startswith('complete_'):
                        drained=[]
                        for row in list(data['tasks']):
                            if row['mode'] in {'global','local'} and row['phase']!='fixture' and L.get_job(data,row['task_id'],row['mode']):
                                try:
                                    drained.append(W.ready(state,data,project,L.read_kv(study/'tasks'/(row['task_id']+'.tsv'))))
                                except ValueError as exc:
                                    if not (study/'tables/resource_stop.tsv').exists() or not W.resource_limit(exc):raise
                                    drained.append(L.is_closed(L.get_job(data,row['task_id'],row['mode'])))
                        if not all(drained):status='closing_submitted_fits_and_checks'
                        elif data['reservations']:
                            state.event(data,'unlaunched_dependent_slots_released_at_stop',';'.join(r['task'] for r in data['reservations']))
                            data['reservations']=[];state.save(data)
                    if status.startswith('complete_') and not (study/'tables/resource_stop.tsv').exists():
                        snapshot=study/'manifests/diagnostic_sources.tsv'
                        if not snapshot.exists():
                            cs=candidates(data,study,project)
                            fields=['point','origin','family','source_fit']
                            L.write_table(snapshot,[{k:W.best(cs,p)[k] for k in fields} for p in ['base']+[r['point'] for r in L.read_table(study/'cases.tsv')]],fields)
                        done=[]
                        for r in L.read_table(snapshot):
                            tid='diagnostic_'+r['point'];t=W.task(study,tid,'diagnostic','diagnostic',r['point'],r['origin'],r['family'],r['source_fit']);W.add(state,data,t)
                            j=L.get_job(data,tid,'diagnostic')
                            if j is None:W.submit(state,data,project,tid,'diagnostic','check','diagnostic');done.append(False)
                            else:done.append(L.is_closed(j))
                        if not all(done):status='final_event_diagnostics'
                    if status.startswith('complete_') and any(not L.is_closed(j) for j in data['jobs'] if j['role']!='controller'):
                        status='closing_submitted_work'
                emit(state,data,project,status)
                if status.startswith('complete_'):
                    L.write_kv(study/'tables/final_status.tsv',dict(status=status,policy_id=POLICY_ID,utc=L.now(),static_delivery_pending=True))
                    L.atomic_text(study/'CLOSURE.md',(study/'STATUS.md').read_text()+'\nScientific dispatch closed. Final accounting, static presentation and verified delivery remain required.\n')
                    state.event(data,'resolution_scientific_dispatch_closed',status);state.save(data)
                    return
                if time.monotonic()-start>42300:
                    child=policy_submit(state,data,project,'controller_policy_'+uuid.uuid4().hex[:10],'controller','controller','control')
                    state.event(data,'resolution_controller_renewal',os.environ['SLURM_JOB_ID']+' '+child['job_id']);state.save(data)
                    return
            except BaseException as exc:
                if W.resource_limit(exc):
                    state.event(data,'resource_stop',str(exc));state.save(data)
                    L.write_kv(study/'tables/resource_stop.tsv',dict(utc=L.now(),reason=str(exc)))
                    emit(state,data,project,'resource_stop_closing')
                else:
                    state.event(data,'dependent_submissions_held',str(exc));state.save(data)
                    L.atomic_text(study/'BLOCKED.txt',L.now()+' '+str(exc)+'\n')
                    emit(state,data,project,'blocked');raise
        time.sleep(30)


def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['controller','audit']);p.add_argument('project',type=Path);p.add_argument('study',type=Path);a=p.parse_args()
    initialize(a.study);verify_release(a.study)
    if a.mode=='controller':controller(a.project,a.study)
    else:
        with L.State(a.study).locked():L.audit(a.study)


if __name__=='__main__':main()
