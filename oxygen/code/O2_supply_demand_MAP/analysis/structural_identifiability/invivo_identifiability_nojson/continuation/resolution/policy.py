"""Resolution01 atomic reservation policy. Called inside the owned state lock."""
from pathlib import Path
import math
import workflow as w

VERSION = 'resolution01'
CAPS = {'baseline': 220*3600, 'profile': 360*3600}
CLOSURE = 4*3600
LIMITS = {'existing_local':3, 'baseline_global':1, 'baseline_polish':1,
          'baseline_repair':2, 'profile_global':2, 'profile_polish':1,
          'profile_ascending':2, 'profile_descending':2, 'profile_repair':1}
FIELDS = ['task','phase','slot_class','slot_key','source_task','decision_fit',
          'reservation_seconds','polish_escrow','check_escrow','consumed_parent',
          'policy_version','recorded_utc']

class PhaseBudgetStop(ValueError):
    """A complete allocation and its dependent evidence no longer fit."""


def gain_compare(a, b):
    # Only a floating subtraction tie guard, not the 1e-8 integrity tolerance.
    scale=max(1., abs(float(a)), abs(float(b)))
    delta=float(a)-float(b)
    return 0 if abs(delta)<=32*2.220446049250313e-16*scale else (1 if delta>0 else -1)


def local_decision(gain, probe, complete, passed=False):
    if passed:return 'accepted'
    if not complete:return 'unresolved_check'
    if gain_compare(gain,.001)>=0 or gain_compare(probe,.002)>0:return 'continue'
    return 'resolution_stop'


def batch_eligible(before_uniform, after_uniform, before_best, after_best):
    return gain_compare(max(before_uniform-after_uniform,before_best-after_best),.002)>0


def completed_batch_eligible(s, policy):
    """Batch two is earned by batch one's frozen completion, never later work."""
    s=Path(s);path=s/'tables/resolution_batch_outcomes.tsv'
    if not path.exists():raise ValueError('Missing frozen complete first-batch outcome')
    matches=[r for r in w.read_table(path) if r['policy']==policy and int(r['batch'])==1]
    if len(matches)!=1:raise ValueError('Ambiguous frozen first-batch outcome')
    r=matches[0]
    values=[float(r[key]) for key in ('before_uniform','before_feasible','after_uniform','after_feasible','protected_gain')]
    if not all(math.isfinite(value) for value in values):raise ValueError('Non-finite frozen protected batch score')
    if r['eligible'].lower() not in ('true','false'):raise ValueError('Invalid frozen eligibility flag')
    before=next(b for b in w.read_table(s/'tables/resolution_batches.tsv') if b['policy']==policy and int(b['batch'])==1)
    for key in ('before_uniform','before_feasible'):
        if gain_compare(float(r[key]),float(before[key]))!=0:raise ValueError('Frozen protected batch reference changed')
    expected=batch_eligible(float(r['before_uniform']),float(r['after_uniform']),float(r['before_feasible']),float(r['after_feasible']))
    gain=max(float(r['before_uniform'])-float(r['after_uniform']),float(r['before_feasible'])-float(r['after_feasible']))
    if gain_compare(float(r['protected_gain']),gain)!=0:raise ValueError('Inconsistent frozen protected gain')
    if (r['eligible'].lower()=='true')!=expected:raise ValueError('Inconsistent frozen first-batch eligibility')
    w.verify_manifest(s,s/r['evidence_manifest'])
    collection=s/'collections'/r['collection']
    if w.sha(collection/'manifest.tsv')!=r['collection_manifest_sha256']:raise ValueError('Completed batch snapshot changed')
    w.verify_manifest(collection,collection/'manifest.tsv')
    return expected


def repair_continuation_allowed(s, rows, previous):
    if not previous:return True
    latest=previous[-1];job=next(j for j in rows if j['task']==latest['task'])
    if job['state']!='COMPLETED' or job['exit_code']!='0:0':return False
    evidence=[r for r in w.read_table(Path(s)/'tables/resolution_local_decisions.tsv') if r['task']==job['task']]
    if not evidence:return False
    r=evidence[-1];complete=r['complete'].lower()=='true'
    return complete and r['decision']=='continue' and local_decision(float(r['window_gain']),float(r['full_check_gain']),complete)=='continue'


def assignments(s):
    p=Path(s)/'tables/resolution_assignments.tsv'
    return w.read_table(p) if p.exists() else []


def positions(rows, ledger):
    by={j['task']:j for j in rows};out={}
    consumed={a['consumed_parent'] for a in ledger if a['consumed_parent']}
    for phase in CAPS:
        actual=live=escrow=closure_used=0
        for a in ledger:
            if a['phase']!=phase:continue
            j=by.get(a['task'])
            if j and w.closed(j):actual+=int(j['actual_core_seconds'])
            else:live+=int(a['reservation_seconds'])
            if a['slot_class']=='closure':closure_used+=int(j['actual_core_seconds']) if j and w.closed(j) else int(a['reservation_seconds'])
            if a['task']+'/polish' not in consumed:escrow+=int(a['polish_escrow'])
            if a['task']+'/check' not in consumed:escrow+=int(a['check_escrow'])
        closure=max(0,CLOSURE-closure_used)
        out[phase]=dict(actual=actual,live=live,escrow=escrow,closure=closure,
                        cap=CAPS[phase],committed=actual+live+escrow+closure)
    return out


def check_capacity(rows, ledger, addition):
    proposed=ledger+[addition]
    pos=positions(rows,proposed)[addition['phase']]
    if pos['committed']>pos['cap']:
        raise PhaseBudgetStop('Phase commitment including dependent work/closure exceeds cap')
    cls=addition['slot_class']
    if cls in LIMITS:
        count=sum(a['phase']==addition['phase'] and a['slot_class']==cls and
                  a['slot_key']==addition['slot_key'] for a in proposed)
        if count>LIMITS[cls]:raise ValueError('Non-resetting role/root allowance exhausted')
    return pos


def reserve(s, rows, task, spec):
    s=Path(s);marker=w.kv(s/'tables/resolution_policy.tsv')
    if marker['policy_version']!=VERSION:raise ValueError('Wrong active policy')
    if w.sha(s/marker['manifest'])!=marker['manifest_sha256']:raise ValueError('Migration policy manifest changed')
    w.verify_manifest(s,s/marker['manifest'])
    phase=spec.get('phase');cls=spec.get('slot_class');key=spec.get('slot_key')
    if phase not in CAPS or not cls or not key:raise ValueError('Explicit phase/class/key required for EVERY allocation')
    if cls not in LIMITS and cls not in ('validation','check','audit','closure','reference'):
        raise ValueError('Unknown allocation class')
    ledger=assignments(s)
    if any(a['task']==task for a in ledger):raise ValueError('Existing reservation; reconcile, never duplicate')
    if any(a['task'] not in {j['task'] for j in rows} for a in ledger):
        raise ValueError('Unattached reservation requires reconciliation')
    expected_phase='profile' if cls.startswith('profile_') else 'baseline'
    if cls in LIMITS and phase!=expected_phase:raise ValueError('Wrong fitting phase')
    if phase=='profile' and not (s/'tables/profile_epoch.tsv').exists():
        raise ValueError('Profile phase not released')
    fit=cls in LIMITS
    stop_path=s/'tables/resolution_budget_stops.tsv'
    if fit and stop_path.exists() and any(a.get('phase','baseline')==phase for a in w.read_table(stop_path)):
        raise PhaseBudgetStop('Fitting held after a recorded phase budget stop')
    if cls=='reference':
        if spec['cores']!=1 or spec['max_calls']!=1 or spec['wall_seconds']>900:
            raise ValueError('Reference reconstruction permits only one point, not hidden fitting')
        purpose=spec.get('reference_purpose')
        if any(a['slot_class']=='reference' and a['slot_key']==key and
               any(j['task']==a['task'] and j['state']=='COMPLETED' and j['exit_code']=='0:0' for j in rows) for a in ledger):
            raise ValueError('This reference is already verified; no reset by epoch/name')
        if purpose=='target_replay':
            if int(spec.get('density_index',0))!=0:raise ValueError('Target replay must reconstruct the unrestricted map')
            if spec['root_family']!='informed':raise ValueError('A reconstructed target point is not an archive/uniform witness')
            gate=w.kv(s/'tables/resolution_review_gate.tsv');w.verify_manifest(s,s/gate['evidence_manifest'])
            queue=w.read_table(s/'tables/resolution_target_replays.tsv')
            if not any(r['policy']==spec['policy'] and r['selection']==spec['selection_name'] and r['root_id']==spec['root_id'] for r in queue):
                raise ValueError('No verified lower target-policy/free-map candidate')
        elif purpose=='profile_initial':
            if spec['root_family']!='profile':raise ValueError('A fixed-density reference belongs to its profile branch')
            w.gate(s,'baseline_gate.tsv')
            if any(not w.closed(j) and any(a['task']==j['task'] and a['phase']=='baseline' for a in ledger) for j in rows):
                raise ValueError('Baseline phase still active')
            epoch=w.kv(s/'tables/profile_epoch.tsv')['epoch']
            chosen=s/'tables'/(spec['selection_name']+'.rds')
            if w.sha(chosen)!=w.sha(s/'profile_epochs'/epoch/('baseline_'+spec['policy']+'.rds')):
                raise ValueError('Profile starts require the frozen free reference, not another branch winner')
        elif purpose=='profile_counterexample':
            if spec['root_family']!='informed':raise ValueError('A counterexample cannot manufacture an independent free witness')
            if w.kv(s/'tables/profile_gate.tsv')['free_baseline_invalidated']!='true':
                raise ValueError('No material constrained/free counterexample')
        else:raise ValueError('Undeclared reference purpose')
    if fit:
        native=w.gate(s,'resolution_native_identity_gate.tsv')
        if native.get('schema')!='2':raise ValueError('Original-parent native identity schema2 is required')
        if cls in ('existing_local','baseline_global','baseline_polish','baseline_repair','profile_repair') and int(spec.get('density_index',0))!=0:
            raise ValueError('Free-baseline fitting cannot retain a constrained density target')
        origins=[j for j in rows if j['root_id']==spec['root_id'] and j['kind']!='engineering']
        if origins and any(j['policy']!=spec['policy'] or j['root_family']!=spec['root_family'] for j in origins):
            raise ValueError('A root cannot change its score policy or lineage family')
        allowed=(8,7200,200000) if cls=='baseline_global' else (8,3600,100000) if cls=='profile_global' else (1,3600,10000) if cls in ('profile_ascending','profile_descending') else (1,7200,20000)
        if spec['cores']!=allowed[0] or spec['wall_seconds']>allowed[1] or spec['max_calls']>allowed[2]:
            raise ValueError('Per-window resource/callback contract exceeded')
        if any(not w.closed(j) and any(a['task']==j['task'] and a['phase']==phase and
                a['slot_class']==cls and a['slot_key']==key for a in ledger) for j in rows):
            raise ValueError('The same fitting component already has an active window')
        if cls=='existing_local':
            roots=w.read_table(s/'policies/resolution01/bridge.tsv')
            if key!=spec['root_id'] or key not in {r['root_id'] for r in roots}:
                raise ValueError('Unknown existing root allowance')
            original=next(r for r in roots if r['root_id']==key)
            if original['root_family']!=spec['root_family']:raise ValueError('Inherited root family cannot be relabelled')
            gate=w.kv(s/'tables/resolution_review_gate.tsv')
            if gate['pass']!='true':raise ValueError('No completed resolution review')
            w.verify_manifest(s,s/gate['evidence_manifest'])
            d=next((r for r in w.read_table(s/'tables/resolution_local_decisions.tsv')
                    if r['task']==spec.get('decision_fit')),None)
            if not d or d['decision']!='continue' or d['root_id']!=key:
                raise ValueError('No current checked whole-window continuation decision')
            if any(a['decision_fit']==spec['decision_fit'] and a['slot_class']==cls for a in ledger):
                raise ValueError('This endpoint already issued its continuation')
        if cls=='baseline_global':
            w.gate(s,'resolution_global_entry_gate.tsv')
            batch=int(spec.get('batch',0));branch=spec.get('branch','')
            if batch not in (1,2) or branch not in ('hybrid','uniform'):
                raise ValueError('Finite batch identity required')
            if key!='{}:{}:{}'.format(spec['policy'],batch,branch):raise ValueError('Wrong global slot key')
            if spec['root_family']!='uniform' or spec['init_mode']!=branch:raise ValueError('Wrong independent initialization')
            if branch=='uniform' and any(j['root_id']==spec['root_id'] for j in rows):
                raise ValueError('A new full-box uniform branch must have a genuinely new root identity')
            if branch=='hybrid' and not origins:raise ValueError('Hybrid restart requires its existing own uniform history')
            w.verify_manifest(s,s/'tables/resolution_dispatch_evidence.tsv')
            decision=next(r for r in w.read_table(s/'tables/strategy_review.tsv') if r['policy']==spec['policy'])
            if decision['decision']!=('initial_global_batch' if batch==1 else 'next_global_batch') or int(decision['batch'])!=batch:
                raise ValueError('No verified finite batch authorization')
            if batch==2 and not completed_batch_eligible(s,spec['policy']):
                raise ValueError('The frozen complete first batch did not earn a second batch')
        if cls in ('baseline_polish','profile_polish'):
            parent=next((a for a in ledger if a['task']==spec.get('polish_parent',spec['source_task'])),None)
            if not parent or parent['slot_class'] not in ('baseline_global','profile_global') or parent['slot_key']!=key:
                raise ValueError('A polish must have its actual global parent')
            native_parent=next(j for j in rows if j['task']==parent['task'])
            if native_parent['root_id']!=spec['root_id'] or native_parent['policy']!=spec['policy'] or native_parent['root_family']!=spec['root_family']:
                raise ValueError('A polish cannot copy another root or family')
        if cls in ('baseline_repair','profile_repair') and key!=spec['policy']:
            raise ValueError('Repair allowance belongs to policy, not renamed component')
        if cls in ('baseline_repair','profile_repair') and not origins:
            raise ValueError('Repair requires an existing genuine root or verified informed reconstruction')
        if cls=='baseline_repair':
            w.verify_manifest(s,s/'tables/resolution_dispatch_evidence.tsv')
            d=next(r for r in w.read_table(s/'tables/strategy_review.tsv') if r['policy']==spec['policy'])
            if d['decision']!='finite_baseline_repair':raise ValueError('No verified finite repair condition')
            previous=[a for a in ledger if a['slot_class']==cls and a['slot_key']==key]
            if not repair_continuation_allowed(s,rows,previous):
                raise ValueError('No complete whole-window/probe improvement justifying another repair')
        if cls=='profile_repair' and w.kv(s/'tables/profile_gate.tsv')['free_baseline_invalidated']!='true':
            raise ValueError('No material profile counterexample')
        if cls.startswith('profile_') and cls!='profile_repair':
            index=int(spec.get('density_index',0))
            if index not in range(1,8) or key!='{}:{}'.format(spec['policy'],index):
                raise ValueError('Profile quota belongs to the original policy/density case, not epoch or root rename')
            pg=w.kv(s/'tables/profile_gate.tsv')
            if pg['free_baseline_invalidated']=='true':raise ValueError('Free baseline veto holds new profile fits')
            cases=w.read_table(s/'tables/profile_anchor_gates.tsv')
            case=next(r for r in cases if r['policy']==spec['policy'] and int(r['index'])==index)
            if case['pass'].lower()=='true':raise ValueError('Accepted profile case has no additional fitting entitlement')
            if cls=='profile_global':
                w.gate(s,'resolution_global_entry_gate.tsv')
                previous=[a for a in ledger if a['slot_class']==cls and a['slot_key']==key]
                if previous:
                    first=next(j for j in rows if j['task']==previous[0]['task'])
                    if first['root_id']!=spec['root_id']:raise ValueError('A profile population continuation cannot change its branch root')
                    proof=w.kv(s/'tables/resolution_profile_review_gate.tsv');w.verify_manifest(s,s/proof['evidence_manifest'])
                    choices=w.read_table(s/'tables/resolution_profile_global_decisions.tsv')
                    decision=next((r for r in choices if r['policy']==spec['policy'] and int(r['index'])==index),None)
                    if not decision or decision['eligible'].lower()!='true':raise ValueError('No meaningful own-branch global progress authorization')
    consumed='';polish=check=0
    if cls=='reference' and spec.get('kind')=='baseline':check=900
    if fit:
        check=900
        if cls=='baseline_global':polish=8100
        if cls=='profile_global' and not any(a['slot_class']=='profile_polish' and a['slot_key']==key for a in ledger):
            # The single polish is escrowed only once across two global windows.
            if not any(a['slot_class']=='profile_global' and a['slot_key']==key for a in ledger):polish=8100
    source=spec.get('source_task','')
    if cls in ('baseline_polish','profile_polish'):
        parent=next(a for a in ledger if a['task']==spec.get('polish_parent',source))
        # A second profile global can use the case's original polish escrow.
        if not int(parent['polish_escrow']):parent=next(a for a in ledger if a['slot_key']==key and int(a['polish_escrow']))
        consumed=parent['task']+'/polish'
    elif cls=='check':
        parent=next((a for a in ledger if a['task']==source),None)
        if parent and int(parent['check_escrow']):consumed=source+'/check'
    if consumed and any(a['consumed_parent']==consumed for a in ledger):
        raise ValueError('Dependent reservation already consumed')
    addition=dict(task=task,phase=phase,slot_class=cls,slot_key=key,source_task=source,
        decision_fit=spec.get('decision_fit',''),reservation_seconds=str(spec['cores']*spec['wall_seconds']),
        polish_escrow=str(polish),check_escrow=str(check),consumed_parent=consumed,
        policy_version=VERSION,recorded_utc=w.utc())
    check_capacity(rows,ledger,addition)
    w.table(s/'tables/resolution_assignments.tsv',ledger+[addition],FIELDS)


def export_accounting(s, rows):
    data=positions(rows,assignments(s))
    w.table(Path(s)/'tables/resolution_phase_accounting.tsv',
            [dict(phase=p,**values) for p,values in data.items()],
            ['phase','actual','live','escrow','closure','cap','committed'])
    return data
