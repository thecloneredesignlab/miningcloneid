#!/usr/bin/env python3
"""Read-only numerical evidence audit and typed closure proof; no submissions."""
import argparse
from pathlib import Path
import workflow as w

IMAGE_SHA='dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372'

def main():
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);a=p.parse_args();s=a.study.resolve()
    rows=w.read_table(s/'tables/jobs.tsv');budget=w.read_table(s/'manifests/budget.tsv')
    if not rows or not all(w.closed(j) for j in rows): raise ValueError('Closure requires every owned allocation closed')
    if not all(j['state']=='COMPLETED' and j['exit_code']=='0:0' for j in rows): raise ValueError('An allocation did not complete successfully')
    if len({j['job_id'] for j in rows})!=len(rows) or len({j['task'] for j in rows})!=len(rows): raise ValueError('Duplicate job/task')
    state=w.State(s).load()
    if state!=rows: raise ValueError('Collected state differs from the exported coherent ledger')
    raw={r.split('|')[0]:r.split('|') for r in (s/'tables/raw_accounting.txt').read_text().splitlines() if r}
    w.verify_manifest(s,s/'manifests/input_manifest.tsv');prior=w.prior_debit(s)
    for release in {j['release'] for j in rows}: w.verify_manifest(s/'releases'/release,s/'releases'/release/'manifest.tsv')
    total=0;fit_count=0;improvement_count=0
    for j in rows:
        parent=raw[j['job_id']];actual=int(parent[3])*int(parent[4]);total+=actual
        if actual!=int(j['actual_core_seconds']) or actual>int(j['reserved_core_seconds']): raise ValueError('Parent accounting mismatch')
        runtime=s/'runtime'/j['task']
        if (runtime/'image_sha256.txt').read_text().strip().split()[0]!=IMAGE_SHA: raise ValueError('Image mismatch')
        if (runtime/'launch_exit_code.txt').read_text().strip()!='0': raise ValueError('Launch mismatch')
        receipt=(s/'logs'/(j['task']+'_receipt.txt')).read_text().strip().splitlines()[0].split(';')[0]
        if receipt!=j['job_id']: raise ValueError('Receipt ownership mismatch')
        command=(s/'logs'/(j['task']+'_submission.txt')).read_text().splitlines()
        for expected in ('--partition=red','--qos=xxlarge','--no-requeue','--cpus-per-task='+j['cores'],'--time='+w.slurm_time(int(j['wall_seconds']))):
            if expected not in command: raise ValueError('Submission contract mismatch: '+expected)
        dispatch=s/'dispatch'/j['task']
        if w.sha(dispatch/'manifest.tsv')!=j['dispatch_sha256']: raise ValueError('Dispatch proof mismatch')
        w.verify_manifest(dispatch,dispatch/'manifest.tsv')
        scope='fits' if j['mode'] in ('global','local') else 'tests'
        root=s/scope/j['task'];w.verify_payload(root)
        if scope=='fits': fit_count+=1
        if scope=='fits' or j['mode'] in ('directions','profile_fresh'):
            improvement_count+=1
            evidence={r['key']:r['value'] for r in w.read_table(s/'tables/record_audits'/(j['task']+'_summary.tsv'))}
            if evidence.get('pass')!='true' or evidence['payload_sha256']!=w.sha(root/'payload_manifest.tsv'): raise ValueError('Improvement audit missing or stale')
    for b in budget:
        group=[j for j in rows if j['role']==b['role']]
        if len(group)>int(b['attempts']) or sum(map(w.charge,group))>int(b['cap_core_seconds']): raise ValueError('Role cap exceeded')
    if sum(int(b['cap_core_seconds']) for b in budget)!=446400 or 192497+total>648000 or prior+total>691200: raise ValueError('Formal/lifetime aggregate cap mismatch')
    phases={b['role']:b['phase'] for b in budget}
    for phase,cap in (('repair',129600),('profile',316800)):
        if sum(w.charge(j) for j in rows if phases[j['role']]==phase)>cap: raise ValueError('Phase cap exceeded')
    for release in {j['release'] for j in rows}:
        if release=='repair02': w.verify_manifest(s,s/'tables/release_bridge_repair02.tsv')
    for name in ('engineering_gate.tsv','engineering_gate_repair02.tsv','baseline_gate.tsv'):
        proof={r['key']:r['value'] for r in w.read_table(s/'tables'/name)}
        if name.startswith('engineering') and proof.get('pass')!='true': raise ValueError('Engineering gate no longer passes')
        if proof['input_manifest_sha256']!=w.sha(s/'manifests/input_manifest.tsv'): raise ValueError('Stale gate input')
        if proof['release_manifest_sha256']!=w.sha(s/'releases'/proof['release']/'manifest.tsv'): raise ValueError('Stale gate release')
        if proof.get('evidence_manifest'): w.verify_manifest(s,s/proof['evidence_manifest'])
        if name=='baseline_gate.tsv' and any(j['role'].startswith(('density_','reserve_density_')) for j in rows):
            # A later constrained counterexample can legitimately invalidate a
            # previously passing gate. Dispatch copies retain its release proof.
            if proof.get('pass')!='true' and not (s/'tables/profile_counterexamples.rds').exists(): raise ValueError('Profile without a valid baseline or counterexample history')
    reserve=[j for j in rows if j['role'] in ('reserve_global','reserve_local','reserve_fresh','reserve_directions')]
    if reserve:
        activated={r['task']:r for r in w.read_table(s/'tables/reserve_activation.tsv')}
        decisions={r['task']:r for r in w.read_table(s/'tables/reserve_decision.tsv')}
        for j in reserve:
            if j['task'] not in activated: raise ValueError('Unrecorded reserve activation')
            context={r['key']:r['value'] for r in w.read_table(s/'dispatch'/j['task']/'context.tsv')}
            if j['mode'] in ('global','local'):
                choice=decisions[j['task']]
                if choice['selection_sha256']!=w.sha(s/'dispatch'/j['task']/'selection.rds'): raise ValueError('Reserve selection changed after decision')
                if any(choice[a]!=context[b] for a,b in (('policy','policy'),('root_id','root_id'),('family','root_family'))): raise ValueError('Reserve decision/dispatch ancestry differs')
        starts=w.read_table(s/'tables/local_start_identity.tsv')
        if {r['task'] for r in starts}!={j['task'] for j in rows if j['mode']=='local'} or any(r['pass']!='TRUE' for r in starts): raise ValueError('First local evaluation proof missing')
    w.audit_formats(s)
    w.table(s/'tables/closure_audit.tsv',[
        dict(key='pass',value='true'),dict(key='owned_closed_allocations',value=str(len(rows))),
        dict(key='audited_fit_payloads',value=str(fit_count)),dict(key='audited_improvement_payloads',value=str(improvement_count)),
        dict(key='profile_allocations',value=str(sum(j['role'].startswith(('density_','reserve_density_')) for j in rows))),
        dict(key='new_actual_core_seconds',value=str(total)),dict(key='formal_actual_core_seconds',value=str(192497+total)),
        dict(key='prior_actual_core_seconds',value=str(prior)),dict(key='lifetime_actual_core_seconds',value=str(prior+total)),
        dict(key='repair_cap_core_seconds',value='129600'),dict(key='profile_cap_core_seconds',value='316800'),dict(key='formal_cap_core_seconds',value='648000'),dict(key='lifetime_cap_core_seconds',value='691200'),
        dict(key='image_sha256',value=IMAGE_SHA),dict(key='model_evaluations',value='0'),dict(key='created_utc',value=w.utc())],['key','value'])
    print('Closure audit passed:',len(rows),'closed allocations;',fit_count,'fit payloads;',total,'new and',prior+total,'lifetime core-seconds.')

if __name__=='__main__': main()
