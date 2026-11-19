"""Read-only scientific evidence collection and metadata closeout, without CAS."""
from pathlib import Path
import shutil
import sys
from protocol import read,table,sha,verify_archive,validate_roles,PARENTS,CAP

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve()
if '--verify-only' in sys.argv[3:]:
    rows=verify_archive(study)
    actual={str(p.relative_to(study)) for p in study.rglob('*') if p.is_file() and p.name not in ('COMPLETE','collection_manifest.tsv')}
    assert actual=={r['path'] for r in rows}
    print('READ_ONLY_ARCHIVE_VERIFICATION_PASSED',len(rows));raise SystemExit(0)
assert not (study/'COMPLETE').exists(),'Closed archives are immutable'
checks=[]
def check(name,condition=True):
    assert condition,name;checks.append(dict(check=name,passed=True))

parent_files=0
for name,digest in PARENTS.items():
    parent_files+=len(verify_archive(project/'oxygen/results/structural_identifiability'/name,digest))
check('all_four_parent_archives_unchanged')
for r in read(study/'inputs/parent_verification.tsv'):
    if r['namespace'] not in PARENTS:
        parent_files+=len(verify_archive(project/'oxygen/results/structural_identifiability'/r['namespace'],r['manifest_sha256']))
check('failed_preparation_archive_preserved')
for r in read(study/'inputs/frozen_source_verification.tsv'):check('frozen_production:'+r['path'],sha(project/r['path'])==r['sha256'])
for r in read(study/'inputs/protected_worktree_hashes.tsv'):check('protected:'+r['path'],sha(project/r['path'])==r['sha256'])
for r in read(study/'source_manifest.tsv'):check('admitted_source:'+r['path'],sha(study/r['path'])==r['sha256'])
jobs=read(study/'accounting/jobs.tsv')
check('all_owned_allocations_terminal',bool(jobs) and all(int(j['reserved_core_seconds'])==0 and j['state'].split()[0] in ('COMPLETED','FAILED','TIMEOUT','NODE_FAIL','OUT_OF_MEMORY','CANCELLED','PREEMPTED','BOOT_FAIL','DEADLINE') for j in jobs))
check('owned_queue_empty',not (study/'accounting/owned_queue.txt').read_text().strip())
used=sum(int(j['actual_core_seconds']) for j in jobs)
check('independent_budget',used<=CAP)
auth={r['key']:r['value'] for r in read(study/'inputs/authorization.tsv')}
retry_used=int(auth['technical_retries_used'])+len(jobs)-1
check('retry_limit',retry_used<=1)
result=read(study/'verification/application_result.tsv')[0]
controls=read(study/'verification/control_results.tsv')
exact=read(study/'verification/exact_checks.tsv')
policy=read(study/'verification/policy_gate.tsv')[0]
passed=result['status']=='PASS_EXACT_INTERFACE_CONTROLS' and policy['status']=='PASS' and all(j['state']=='COMPLETED' and j['exit_code']=='0:0' for j in jobs) and {r['case'] for r in controls if r['status']=='PASS'}=={'T'+str(i).zfill(2) for i in range(10)} and all(r['residual']=='0' for r in exact)
status='interface_pass_project_cost_unverified' if passed else 'HOLD_interface_validation'
roles=read(study/'inputs/parameter_roles.tsv');validate_roles(roles)
mapped=[]
for row in roles:
    r=dict(row)
    if r['parameter'].startswith('init_'):validation='exact_initial_chain_and_sharing_controls_passed' if passed else 'held'
    elif r['parameter']=='sigma_kary':validation='exact_normalized_observer_controls_passed' if passed else 'held'
    elif r['parameter']=='sigma_growth':validation='positive_likelihood_SD_role_verified_outside_mean_classification' if passed else 'held'
    else:validation='static_original_mechanism_mapping_preserved'
    r.update(interface_evidence=validation,project_structural_result='not_analyzed',Stage_B_released=False)
    mapped.append(r)
table(study/'reports/parameter_roles.tsv',mapped)

# This bridge is source inspection only. No R function or simulator is invoked.
source_paths={Path(r['original_path']).name:study/r['frozen_copy'] for r in read(study/'inputs/production_sources.tsv')}
lineage=source_paths['o2_supply_demand_map_invitro_lineage_simulation_utils.R'].read_text()
objective=source_paths['o2_supply_demand_map_invitro_objective_utils.R'].read_text()
source_checks={
 'finite_initial_normalization':'ivt_gaussian_initial_state <- function' in lineage and 'w / s' in lineage,
 'initial_override_and_inheritance':'init_state_override' in lineage and 'reseeded_state = chosen_frac * reseed_live_cells' in lineage,
 'target_count_selected_day':'order(abs(live_cells[candidate_idx] - target_live_cells_use), candidate_idx)' in lineage,
 'karyotype_column_normalization':'kernel_colsum <- colSums(kernel_mat)' in objective and 'kernel_mat <- sweep(kernel_mat, 2, kernel_colsum' in objective,
 'flow_post_mixture_grid_normalization':'norm_const <- sum(density_vals) * step' in objective,
 'growth_standard_deviation':'sd = sigma_use' in objective and 'ivt_growth_loglik_df' in objective,
 'flow_configured_not_fitted':'sigma_flow_ploidy <- if (is.null(flow_kernel_sd_ploidy))' in objective,
 'flow_default_width':'2 * step' in objective and '1 / 22' in objective,
 'complete_executable_parameter_names':all('"'+r['parameter']+'"' in objective for r in roles)}
for k,v in source_checks.items():check('static_bridge:'+k,v)
table(study/'reports/static_source_contract.tsv',[dict(check=k,passed=v,evidence_route='static_source_contract') for k,v in source_checks.items()])
tool_map=[]
for r in read(study/'inputs/model_tool_map.tsv'):
    q=dict(r);q.update(StrikePy='native_initial_manifold_support_unverified; exact_standalone_interface_available' if passed else 'interface_hold',project_structural_result='not_analyzed',original_parameter_fixed_or_dropped=False)
    tool_map.append(q)
table(study/'reports/model_tool_map.tsv',tool_map)
readiness=[]
for r in read(study/'inputs/project_case_readiness.tsv'):
    q=dict(r);q['compute_released']='False';q['project_structural_result']='not_analyzed'
    if q['tool']=='StrikePy':
        q['method_gate']='mixed_prior_controls_plus_exact_interface; unchanged_native_exponent_path_incomplete' if passed else 'interface_hold'
        q['protocol_map']='exact_initial_observer_fixture_pass; original_joint_model_and_branches_unverified' if q['model_form']=='C' and passed else q['protocol_map']
    readiness.append(q)
table(study/'reports/project_case_readiness.tsv',readiness)
acceptance=[dict(gate='A3_scoped_interface_controls',passed=passed,status=status),
            dict(gate='A4_static_readiness_audit',passed=True,status='completed_scoped_assessment'),
            dict(gate='native_initial_manifold_support',passed=False,status='not_verified'),
            dict(gate='production_reduction_tail_and_discrete_mapping',passed=False,status='existing_holds_preserved'),
            dict(gate='Stage_B',passed=False,status='not_released'),dict(gate='Step_5',passed=False,status='not_released')]
table(study/'reports/acceptance.tsv',acceptance)
all_results=controls+[dict(case='T10',status='PASS_STATIC_REVIEW',evidence_route='static_source_contract',detail='Original sources and all roles mapped; project solver and bridge remain unverified')]
table(study/'reports/control_results.tsv',all_results)
table(study/'accounting/budget_summary.tsv',[dict(key=k,value=v) for k,v in dict(cap_core_hours=8,cap_core_seconds=CAP,new_allocated_core_seconds=used,used_core_hours=used/3600,live_reserved_core_seconds=0,technical_retries_used=retry_used,technical_retry_max=1,HPC_repeated_allocations=len(jobs)-1,prior_package_charge_carried=0).items()])
table(study/'reports/final_status.tsv',[dict(key=k,value=v) for k,v in dict(delivery_complete=True,interface_status=status,original_mechanisms_retained=14,original_nuisances_retained=6,original_culture_parameters_classified=0,native_initial_manifold_support='not_verified',new_native_calls=0,new_randomized_rank_calls=0,Stage_B=False,Step_5=False,fits=0,synthetic_culture_datasets=0,production_replays=0,author_contacts=0).items()])
(study/'reports/custom_module_decision.md').write_text(
 '# Scoped adapter decision\n\nA finite initial-substitution/chain-rule and normalized-observer adapter was tested with independent analytical references using existing SymPy operations. No general-purpose structural-identifiability solver, native rank patch or discrete-time engine was built.\n\nThe evidence route is an exact standalone interface certificate. Native StrikePy parameterized-initial support remains unverified, and its historical unknown-exponent automated return remains incomplete. A compliant MATLAB route was not reopened.\n\nAny actual culture structural case requires separate Stage B review and authorization, preserving all14+6 roles, finite support, branch assumptions and the continuous/discrete distinction.\n')
text='''# Stage A initial-distribution and observation-interface validation

Date: 2026-11-16, America/New_York.
Outcome: STATUS.

## What was checked

The registered small controls tested finite normalized Gaussian initial weights,
unknown mean/width chain rules, state-space output jets before initial
substitution, normalized karyotype kernel-width derivatives, finite-grid flow
normalization, shared cohort parameters, variance-only nuisance roles, exact-zero
activity expressions, and invalid/stale interface rejection. Independent
references use exponential-family moment formulas, closed-form toy jets,
quotient derivatives, explicit inverses and output-preserving families.

The three-point initial-family minor is 2*p_minus*p_zero*p_plus/s^5, positive on
the declared domain. A toy common-rate model has a nonzero constrained witness.
Live-total output alone preserves arbitrary initial mean/width changes. The
unbounded Gaussian negative control preserves outputs when initial and observer
variances change by opposite amounts. This last symmetry is not attributed to
the finite production kernels. Flow is a finite-grid density normalized after
mixing; karyotype kernels normalize each column. Their edge normalization
counterexample is nonzero exactly.

## Parameter and evidence scope

All fourteen original mechanisms and six nuisances remain in the inventory.
Four initial-family nuisances retain cohort sharing, sigma_kary remains an
unknown applicable kernel SD, and sigma_growth remains a positive likelihood SD
outside deterministic mean-output classification. No new fitted flow width was
introduced. No original culture parameter was fixed, estimated or classified.

This delivery is an exact standalone interface certificate plus a static source
bridge, not a claim that the unchanged native toolbox supports this initial
manifold or completes every nonrational problem. Original C04 native timeout and
all earlier attempts remain preserved. Only the previously admitted private
Python/SymPy path was used; no native structural rank solver was called.

## Remaining gates

Full culture equations, nonlinear branch domains, retained accounts, finite
upper/lower tails and full-state symbolic cost still need project review.
Failed class closures are not repaired by an interface pass. Continuous ODE jets
do not establish DT=0.05 updates, target-count selected end days, inherited
passage state or sparse sampling. Flow processing and physical sampling/live-dead
documentation remain experimental limitations. Stage B and Step 5 remain
not_released. No synthetic culture data, propagation, replay, fit, profile or
tumor analysis was run.

## Resources and runtime policy

The independent ceiling is8 allocated core-hours, with at most one technical
retry. There is no discretionary case/process/derivative timer and no native
maxLietime path. Scheduler reservation is derived from the aggregate budget,
explicit retry/settlement reserves and inspected QoS constraints. No application
concurrency cap or automatic requeue was introduced. Only this package's job IDs
were managed. Prior budget ledgers were not reset or spent.

The scientific worker and children have file/process traces, a preload access
guard and imported-module inventory. The recorded policy gate must pass for an
interface PASS. Environment module bootstrap is outside this worker trace;
its frozen module recipe is inherited, not a universal whole-platform audit.
Authored outputs are English and no Java/JavaScript/JSON artifacts are used.
No author contact, old controller change, production change, commit or push
occurred. Raw exact checks, witnesses and collection records are retained.
'''.replace('STATUS',status)
text+='\nNew allocation charge: '+str(used)+' core-seconds ('+format(used/3600,'.9f')+' core-hours). Technical retries used: '+str(retry_used)+'. The retained pre-submission failure concerned a missing metadata destination directory; its correction consumed the one technical retry and no HPC allocation. HPC repeated allocations: '+str(len(jobs)-1)+'. No live owned reservation remains.\n'
(study/'REPORT.md').write_text(text)
for p in study.rglob('*'):
    if p.is_file():check('artifact_format:'+str(p.relative_to(study)),p.suffix.lower() not in ('.json','.json5','.jsonl','.js','.java','.jar','.ipynb'))
check('English_report',not any('\u4e00'<=c<='\u9fff' for c in text))
table(study/'reports/delivery_verification.tsv',checks)
table(study/'reports/verification_counts.tsv',[dict(parent_files_verified=parent_files,frozen_source_inputs_verified=len(read(study/'inputs/frozen_source_verification.tsv')),protected_worktree_files_verified=len(read(study/'inputs/protected_worktree_hashes.tsv')),static_preparation_checks=len(read(study/'reports/static_tests.tsv')),exact_worker_checks=len(exact),delivery_checks=len(checks))])
manifest=[]
for p in sorted(study.rglob('*')):
    if p.is_file() and p.name not in ('COMPLETE','collection_manifest.tsv'):manifest.append(dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)))
table(study/'collection_manifest.tsv',manifest)
(study/'COMPLETE').write_text(sha(study/'collection_manifest.tsv')+'\n'+str(len(manifest))+' manifest-listed files\n')
print('DELIVERY_CLOSED',status,len(exact),'exact checks;',used,'new core-seconds;',len(manifest),'files')
