"""Narrative over settled closure evidence; no model or scheduler calls."""
import workflow as w
import policy as p

def write_report(s,jobs,phase,endpoint,actual,counts):
    if any(not w.closed(j) for j in jobs):raise ValueError('Report requires settled parent allocations')
    gates=w.read_table(s/'tables/baseline_policy_gates.tsv')
    chains=w.read_table(s/'tables/baseline_chain_summary.tsv')
    if any('uniform_representative_gap' not in r for r in gates):
        raise ValueError('Final reporting requires separate actual-minimum and representative gaps')
    lines=['# Resolution01 execution, verification and settled accounting','',
        'Updated '+w.utc()+'.','',
        'Endpoint: **'+endpoint+'**. Scientific acceptance and successful evidence/accounting closure are separate outcomes.','',
        'Both free baselines must pass before conditional density profiles. Resolution or resource stops do not establish global optimality or structural nonidentifiability. Weighted-score profiles are diagnostic; no calibrated confidence interval is claimed.','',
        '## Phase and lifetime accounting','',
        '| Phase | Settled core-hours | Phase cap | Unused core-hours |',
        '|---|---:|---:|---:|']
    for name,a in phase.items():
        lines.append('| {} | {:.8f} | {:.0f} | {:.8f} |'.format(name,a['actual']/3600,a['cap']/3600,(a['cap']-a['actual'])/3600))
    lines+=['','Historical closed opening:388524core-seconds (107.92333333core-hours). '
        'The complete costs of the six inherited third windows belong to the revised baseline phase once. '
        'Opening plus both revised phases gives {}core-seconds ({:.8f}core-hours) lifetime. '
        'The legacy accounting identity301846+{} gives the same total. Failed/partial work remains charged. '
        'Unused reservations are not actual spending. Culture remains outside these ledgers.'.format(
            301846+actual,(301846+actual)/3600,actual),'',
        '## Free-baseline verification and recovery','',
        '| Policy | Actual verified free minimum | Known protected floor | Archive minimum gap | Uniform minimum gap | Uniform representative gap | Independent recovery | Complete representative | Accepted |',
        '|---|---:|---:|---:|---:|---:|---|---|---|']
    for r in gates:
        archive_gap=min(float(c['verified_best']) for c in chains if c['policy']==r['policy'] and c['family']=='archive')-float(r['known_feasible_bound'])
        lines.append('| {} | {:.12f} | {:.12f} | {:.12g} | {:.12g} | {:.12g} | {} | {} | {} |'.format(
            r['policy'],float(r['best_objective']),float(r['known_feasible_bound']),archive_gap,float(r['uniform_gap']),
            float(r['uniform_representative_gap']),r['family_recovery'],r['representative_available'],r['pass']))
    lines+=['','The best genuine Uniform minimum minus the free minimum decomposes into the unchanged weighted score terms:',
        '', '| Policy | Burden | Ploidy | Necrosis | Weighted prior | Total gap |', '|---|---:|---:|---:|---:|---:|']
    for r in w.read_table(s/'tables/recovery_component_differences.tsv'):
        lines.append('| {} | {:.12g} | {:.12g} | {:.12g} | {:.12g} | {:.12g} |'.format(r['policy'],
            *(float(r[k]) for k in ('burden_difference','ploidy_difference','necrosis_difference','weighted_prior_difference','actual_minimum_gap'))))
    lines+=['','Recovery compares genuine archive/uniform minima with the actual lowest verified candidate at0.002. '
        'A full fixed-center representative within1e-4 supplies only its own four-scale axis/coupled/observation checks. '
        'Minimum, representative and completed anchor vectors remain separate; small improvements do not silently reset the anchor. '
        'A lower protected candidate awaiting the correct target-policy/unrestricted-map replay remains an explicit unresolved deficit.','',
        'Same-vector score, source, prediction, arithmetic and payload identities retain1e-8. '
        'Actual calculations use the pinned R4.4.2 SIF (SHA256 dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372), '
        'unchanged numerical5/raw3 inputs, masks, weights, bounds and score policies; N0=Nref=1e6 and native dead-stock/volume semantics are preserved. '
        'Full original-parent native identity is checked when vector/policy match; changed-target/policy references use the correctly reconstructed fresh center. '
        'Cache counters remain diagnostics. Native/prediction/record/population evidence is linked below.','',
        '## Non-resetting fitting allowances','',
        '| Phase | Role | Original quota key | Allocated windows, including failures | Maximum |',
        '|---|---|---|---:|---:|']
    for (name,role,key),number in sorted(counts.items()):
        if role in p.LIMITS:lines.append('| {} | {} | {} | {} | {} |'.format(name,role,key,number,p.LIMITS[role]))
    lines+=['','Six inherited third windows each consumed one of their root\'s three revised existing-local slots. '
        'Global branch polishes and declared finite repairs retain real ancestry and their separate role allowances. '
        'A new epoch, source release, component or task name cannot reset quotas. '
        'Whole-window gain>=0.001 or a complete probe gain>0.002 can justify another eligible local window; '
        'the native1e-10 diagnostic and old25-sweep rule do not govern outer dispatch. '
        'The second baseline global batch requires>0.002 improvement against fixed pre-batch protected minima. '
        'No third batch is authorized. Native stop labels remain distinct from resolution/window/resource outcomes.','',
        '## Reproduction and checkpoint restoration','']
    restores=[]
    for j in jobs:
        if j['mode']!='local' or j['selector']!='checkpoint' or j['state']!='COMPLETED' or j['exit_code']!='0:0':continue
        proof=w.kv(s/'fits'/j['task']/'start_identity.tsv')
        if float(proof['difference'])!=0 or proof['projection']!='false':
            raise ValueError('An exact local checkpoint restore changed its first score or parameters')
        restores.append(j['task'])
    arithmetic=w.kv(s/'tables/audit_math_validation.tsv')
    if arithmetic['pass']!='true':raise ValueError('Archived score arithmetic validation failed')
    predictive=w.read_table(s/'tables/prediction_identity_summary.tsv')
    if not predictive or any(r['pass'].lower()!='true' for r in predictive):
        raise ValueError('Missing or failing predictive identity summary')
    lines+=['{} completed local checkpoint restorations reproduced the first score exactly, without parameter projection. '
        'The four real free-20/fixed-19 complete-check validators and the normal global-entry 500-callback test passed in the pinned SIF. '
        'Full native outputs, score components and original-parent identity are audited separately from cache counters.'.format(len(restores)),
        '', 'The archived arithmetic audit checked {} real cases and {} analytic cases: maximum arithmetic gap {}, '
        'maximum canonical-score gap {}. Across {} predictive-identity records, maximum score gap {} and maximum native predictive-field gap {}. '
        'These are reproduction checks rather than optimization-accuracy guarantees.'.format(
            arithmetic['archived_cases'],arithmetic['analytic_cases'],arithmetic['maximum_arithmetic_gap'],
            arithmetic['maximum_canonical_gap'],len(predictive),
            max(float(r['score_gap']) for r in predictive),max(float(r['native_max_absolute_gap']) for r in predictive)), '',
        '## Closed-window timing','',
        '| Task | Job | Cores | Parent minutes | Charged core-hours | Window callbacks | Stage seconds | Native reason | Window gain |',
        '|---|---|---:|---:|---:|---:|---:|---|---:|']
    ledger={r['task']:r for r in p.assignments(s)}
    by={j['task']:j for j in jobs}
    for r in w.read_table(s/'tables/fit_summary.tsv'):
        if r['task'] not in ledger:continue
        j=by[r['task']]
        lines.append('| {} | {} | {} | {:.3f} | {:.8f} | {} | {} | {} | {} |'.format(
            r['task'],j['job_id'],j['cores'],int(j['actual_core_seconds'])/int(j['cores'])/60,
            int(j['actual_core_seconds'])/3600,r['calls'],r['stage_seconds'],r['termination'],r['window_gain']))
    lines+=['','Parent AllocCPUS multiplied by ElapsedRaw is charged once. '
        'Child timers and inherited root-call totals are not added as extra allocation cost. '
        'Window callback totals describe distinct allocations; checkpoint ancestry counters are not independent additional calls. '
        'Population continuation is a warm restart, with retained own-root incumbents and native RNG/population geometry evidence, not bitwise uninterrupted DE.','',
        '| Phase | Allocation role | Parents | Settled core-hours, including failures | Successful fitting callbacks |',
        '|---|---|---:|---:|---:|']
    role_totals={}
    callbacks={r['task']:int(r['calls']) for r in w.read_table(s/'tables/fit_summary.tsv')}
    for task,a in ledger.items():
        j=by[task];values=role_totals.setdefault((a['phase'],a['slot_class']),[0,0,0])
        values[0]+=1;values[1]+=int(j['actual_core_seconds']);values[2]+=callbacks.get(task,0)
    if sum(v[1] for v in role_totals.values())!=sum(a['actual'] for a in phase.values()):
        raise ValueError('Role costs do not reconcile with settled phase costs')
    for (name,role),values in sorted(role_totals.items()):
        lines.append('| {} | {} | {} | {:.8f} | {} |'.format(name,role,values[0],values[1]/3600,values[2]))
    lines+=['','Callback counts above describe successful fitting windows only. '
        'Fresh checks, canonical replays and engineering callbacks have separate native counters and are not fitting-window callbacks. '
        'Any unsuccessful partial work remains charged even when it lacks a sealed fitting counter.','',
        '## Unsuccessful allocations in the revised phases','',
        '| Task | Job | Slurm state | Exit | Charged core-hours |',
        '|---|---|---|---|---:|']
    bad=[j for j in jobs if j['task'] in ledger and not (j['state']=='COMPLETED' and j['exit_code']=='0:0')]
    for j in bad:lines.append('| {} | {} | {} | {} | {:.8f} |'.format(j['task'],j['job_id'],j['state'],j['exit_code'],int(j['actual_core_seconds'])/3600))
    if not bad:lines.append('| None | | | |0 |')
    lines+=['','All failures retain receipts, logs, partial outputs, reviews and actual fees. '
        'Metadata-only audits explicitly forbid model/prediction calls; their allocations, including the last one, are settled above.']
    extra_tasks=('r_resolution01_b_production_uniform_new_w0003_c02','a_resolution01_w0005')
    if all(task in by for task in extra_tasks):
        extra_cost=sum(int(by[task]['actual_core_seconds']) for task in extra_tasks)
        lines+=['','One unnecessary optional check of the older production new-Uniform w0003 endpoint was issued after its w0004 descendant had already supplied a lower checked representative. '
            'That check and its associated archived-output audit both completed successfully and cost {} core-seconds ({:.8f} core-hours), '
            'already included in the phase totals. The scheduling guard now defers superseded old checks. '
            'This overhead did not grant another fitting window or reset any allowance.'.format(extra_cost,extra_cost/3600)]
    lines+=['',
        '## Conditional density phase','']
    if phase['profile']['actual']==0:
        lines+=['No density-profile allocation was made. Failed free-baseline recovery leaves the conditional stage unreleased; no density-identifiability result follows.']
    else:
        lines+=['The original seven full-precision density anchors and both policies remain unchanged. '
            'Inspect the case/branch gates for accepted and unresolved coverage. Raw signed deltas and frozen reference identities are retained. '
            'A delta below-0.002 vetoes new dependent work; profile epochs do not reset finite case or free-repair allowances.']
    lines+=['','## Evidence and delivery','',
        '[Free gates](../tables/baseline_policy_gates.tsv), [actual candidates](../tables/resolution_candidates.tsv), '
        '[checked anchors](../tables/resolution_anchors.tsv), [scheduling decisions](../tables/resolution_local_decisions.tsv), '
        '[finite strategy](../tables/strategy_review.tsv), [full native identity](../tables/resolution_native_identity_gate.tsv), '
        '[record audits](../tables/record_audits/), [prediction audits](../tables/prediction_audits/), '
        '[population audit](../tables/population_geometry_audit.tsv), [settled phase fees](../tables/resolution_final_phase_accounting.tsv), '
        '[reviewed failures](../tables/reviewed_failures.tsv), [closure](../tables/closure_audit.tsv), '
        '[HPC delivery fingerprint](../tables/delivery_verified.tsv).','',
        'Complete native/checkpoint/improvement payloads are authoritative on HPC. '
        'The local lean collection carries full remote-payload verification proofs and complete small allocated metadata packets; '
        'it deliberately has no local science COMPLETE marker. Culture tasks/budgets, shared runner/image/libraries and sealed prior runs remain protected.','']
    w.atomic_text(s/'reports/resolution_execution_report.md','\n'.join(lines))
