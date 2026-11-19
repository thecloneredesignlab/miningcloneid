"""Collect bounded Step5 evidence and design decisions; metadata only."""
import ast
from collections import Counter
from pathlib import Path
import shutil
import sys
from common import *


def sealed(folder, manifest):
    rows = read(folder / manifest)
    for row in rows:
        p = folder / row['path']
        assert p.is_file() and sha(p) == row['sha256'], str(p)
        if 'bytes' in row:
            assert p.stat().st_size == int(row['bytes']), str(p)
    return len(rows)


def original_numeric_parameters(rows):
    wanted = set(MECHANISMS + NUISANCES)
    selected = [r for r in rows if r['parameter'] in wanted]
    assert len(selected) == 20 and {r['parameter'] for r in selected} == wanted
    return {r['parameter']: float(r['value']) for r in selected}


def main():
    project = Path(sys.argv[1]).resolve()
    study = Path(sys.argv[2]).resolve()
    assert not (study / 'CLOSED').exists(), 'Archive is already closed'
    verify(study)
    sealed(study, 'validation_source_manifest.tsv')
    roles = read(study / 'inputs/parameter_roles.tsv')
    roles_valid(roles)
    out = study / 'reports'
    jobs = read(study / 'accounting/jobs.tsv')
    expected = {'M01', 'O01', 'R01', 'R02', 'O02', 'R03', 'G01'}
    assert {r['task'] for r in jobs} == expected and len(jobs) == 7
    phase_rows = []
    for job in jobs:
        d = study / 'verification' / job['task']
        app = read(d / 'application_result.tsv')[0]['status']
        policy = read(d / 'policy_gate.tsv')[0]['status']
        gate = read(d / 'gate.tsv')[0]['status']
        assert job['state'] == 'COMPLETED' and job['exit_code'] == '0:0'
        assert app == policy == 'PASS' and gate.startswith('PASS')
        n = sealed(d, 'phase_manifest.tsv')
        phase_rows.append(dict(task=job['task'], job_id=job['job_id'],
            slurm_state=job['state'], exit_code=job['exit_code'],
            application_status=app, policy_status=policy, scientific_gate=gate,
            verified_files=n, actual_core_seconds=job['actual_core_seconds']))
    table(out / 'phase_acceptance.tsv', phase_rows)
    spent = sum(int(j['actual_core_seconds']) for j in jobs)
    assert spent <= CAP and not any(int(j['reserved_core_seconds']) for j in jobs)
    retry = int(kv(study / 'accounting/retry_ledger.tsv')['used'])
    assert retry <= 1
    old_checks = []
    for name, marker, digest in (
        (ORIGIN, 'COMPLETE', ORIGIN_SHA),
        ('20261116_culture_step04_stageB_B6_closeout_01', 'CLOSED',
         '4a2f843ad0249211acc74690332254abad5c2788da04d7ee0bb69e2dcf9f19d9')):
        archive = project / 'oxygen/results/structural_identifiability' / name
        assert sha(archive / 'collection_manifest.tsv') == digest
        assert (archive / marker).read_text().splitlines()[0] == digest
        n = sealed(archive, 'collection_manifest.tsv')
        old_checks.append(dict(archive=name, files=n, manifest_sha256=digest,
                              status='PASS_UNCHANGED'))
    table(out / 'prior_archive_preservation.tsv', old_checks)
    changes = []
    for r in read(study / 'preservation_baseline.tsv'):
        p = project / r['project_path']
        after = sha(p) if p.exists() else 'MISSING'
        changes.append(dict(**r, current_sha256=after,
                            status='UNCHANGED' if after == r['sha256'] else 'EXTERNAL_CHANGE_REQUIRES_REVIEW'))
    table(out / 'unrelated_work_preservation.tsv', changes)
    assert all(r['status'] == 'UNCHANGED' for r in changes), 'Review concurrent external changes; never overwrite'

    ledger = []
    for task in ('R01', 'R02', 'R03'):
        ledger.extend(dict(task=task, **r) for r in read(study / 'verification' / task / 'forward_ledger.tsv'))
    assert len(ledger) <= 74 and len({r['eval_id'] for r in ledger}) == len(ledger)
    assert all(r['status'] in ('returned', 'protocol_infeasible') for r in ledger)
    points = read(study / 'verification/M01/bridge_points.tsv')
    attempted = {r['eval_id'] for r in ledger if r['task'] == 'R03'}
    assert attempted == {r['eval_id'] for r in points if r['status'] == 'registered'}
    table(out / 'forward_evaluation_registry.tsv', ledger)
    table(out / 'registered_point_dispositions.tsv', [dict(**r,
          disposition=next((x['status'] for x in ledger if x['eval_id'] == r['eval_id']), r['status'])) for r in points])
    observer = read(study / 'verification/O02/observer_call_ledger.tsv')
    contrasts = read(study / 'verification/O01/contrast_attempts.tsv')
    certs = read(study / 'verification/O01/observer_extension_certificates.tsv')
    assert len(observer) <= 27 and len(contrasts) <= 6
    assert len(certs) == 3 and all(r['status'] == 'CERTIFIED_LOCAL_IDEAL_EXTENSION' for r in certs)
    assert all(int(r['new_full_determinant_evaluations']) == 0 for r in certs)
    table(out / 'observer_extension_manifest.tsv', certs)
    unresolved = read(out / 'case_status_crosswalk.tsv')
    assert len(unresolved) == 12 and all(r['historical_direct_status'] == 'UNRESOLVED_REGISTERED_SCOPE' for r in unresolved)
    table(out / 'case_inference_addendum.tsv', [dict(**r,
        inference=('Local ideal18 by L-output inclusion' if r['case'] == 'B02' else
                   'Local ideal14 by known-initial restriction and output inclusion' if r['case'] == 'B06' else
                   'Local ideal19 by new block observer extension' + (' and flow-output inclusion' if r['case'] == 'B05' else '')),
        old_order2_result='UNCHANGED; single-O nomination remains unresolved',
        actual_sparse_discrete_identification='NOT_PROVED') for r in unresolved])
    for name, task in [('floating_observer_audit.tsv', 'O02'),
                       ('variance_role_audit.tsv', 'O02'),
                       ('local_observation_geometry.tsv', 'G01'),
                       ('selection_branch_audit.tsv', 'R03'),
                       ('anchor_replay_checks.tsv', 'R02'),
                       ('cost_measurement.tsv', 'R01')]:
        shutil.copy2(study / 'verification' / task / name, out / name)
    shutil.copy2(study / 'verification/O01/observer_extension_proof.md', out / 'observer_extension_proof.md')
    floating = read(out / 'floating_observer_audit.tsv')
    params = original_numeric_parameters(read(study / 'verification/M01/anchor_parameters.tsv'))
    seed_rows = [r for r in floating if abs(float(r['width']) - params['sigma_kary']) < 1e-14]
    assert len(seed_rows) == 3 and all(r['all_reference_offdiagonal_zero'] == 'TRUE' and float(r['identity_max_gap']) == 0 for r in seed_rows)
    geometry = read(out / 'local_observation_geometry.tsv')
    assert len(geometry) == 19 and {r['parameter'] for r in geometry} == set(MECHANISMS + NUISANCES) - {'sigma_growth'}
    geom = {r['parameter']: r for r in geometry}
    assert geom['sigma_kary']['status'] == 'NUMERICALLY_FLAT_AT_REGISTERED_PROBES'
    table(out / 'branch_switch_summary.tsv', [dict(
        flag=flag, affected_coordinates=sum(flag in r['observed_branch_flags'] for r in geometry),
        interpretation=('score-floor branch; does not itself imply a mean-output discontinuity'
                        if flag == 'score_floor_sha_switch' else
                        'protocol/pruning/floating-support branch flag; not a nonidentifiability proof'))
        for flag in ('selection_sha_switch', 'pruning_sha_switch',
                     'initial_zero_sha_switch', 'score_floor_sha_switch',
                     'clamp_branch_not_certified')])
    columns = read(study / 'verification/G01/column_scaling.tsv')
    sigma_column = next(r for r in columns if r['parameter'] == 'sigma_kary')
    table(out / 'sigma_local_derivative_diagnostic.tsv', [dict(
        parameter='sigma_kary', normalized_coordinate_secant_norm=sigma_column['unscaled_column_norm'],
        probes=4, status=geom['sigma_kary']['status'],
        scope='two registered local scales on saved latent bundle; not a global derivative certificate or a derivative scan at all9 widths')])
    table(out / 'sigma_kary_status.tsv', [dict(regime=r, status=s, implication=i) for r, s, i in (
        ('ideal continuous L plus initial karyotype; three interior witnesses', 'CERTIFIED_LOCAL_19',
         'sigma_kary joins14 mechanisms and4 unknown initial coordinates; no global or actual-schedule claim'),
        ('seed144 binary64 production kernel', 'IDENTITY_KERNEL_AT_ANCHOR',
         'all offdiagonal kernel entries underflow at the anchor; raw width recovery cannot be required here'),
        ('four registered width perturbation passes on saved latent bundle', geom['sigma_kary']['status'],
         'retain a free diagnostic coordinate; report kernel/distribution recovery and calibration need'))])

    # Original20-coordinate identities are retained, with current dispositions separate from old columns.
    current_roles = []
    for r in roles:
        name = r['parameter']
        status = 'variance-law target; outside deterministic mean map' if name == 'sigma_growth' else geom[name]['status']
        current_roles.append(dict(parameter=name, scientific_role=r['scientific_role'],
            units=r['units'], lower=r['lower'], upper=r['upper'], transform=r['transform'],
            seed144_value=params[name], future_primary_role='unknown; not silently fixed',
            ideal_structural_evidence='conditional local19 observer extension' if name != 'sigma_growth' else 'positive Gaussian variance law under stated observation assumption',
            production_bridge_status=status, full_sparse_D_structural_status='NOT_PROVED',
            source_file=r['source_file'], source_line=r['source_line'], source_sha256=r['source_sha256']))
    table(out / 'parameter_roles.tsv', current_roles)
    estimands = []
    for r in current_roles:
        n = r['parameter']
        estimands.append(dict(target=n, type='original_coordinate', inventory_role=r['scientific_role'],
            primary_future_status=r['future_primary_role'], evidence=r['production_bridge_status'],
            error_report='natural and normalized transformed error; absolute error at zero',
            release='HELD_FOR_SEPARATE_DESIGN_AND_EXECUTION_APPROVAL'))
    for target, definition, warning in (
        ('effective_missegregation', 'p_base+p_misseg*mu(N,O)/(k_o_mis+mu(N,O))', 'track induction, basal risk and death-response compensation'),
        ('post_error_survival', 'production S(N,m), including smax,beta,n_exp and error size', 'report the full retained production function, not only its ploidy prefactor'),
        ('growth_death_responses', 'lambda(N,O) and mu(N,O)', 'retain ploidy_related death and original unknown exponents'),
        ('observation_kernel', 'column-normalized finite integer K_sigma', 'raw sigma may be resolution limited in binary64'),
        ('predicted_distributions', 'complete karyotype and flow mass plus selected-endpoint growth', 'absolute population and normalized composition are separate'),
        ('single_O_combinations', 'alpha_o2*H(O;C,n) and mu_hp*H(O;C,n)', 'single-constant-O quotient only; never remove original coordinates in the multi-O protocol')):
        estimands.append(dict(target=target, type='function_or_combination', inventory_role='derived; no new fitted parameter',
            primary_future_status=definition, evidence=warning, error_report='function/distribution error at declared conditions',
            release='REGISTERED_FOR_REVIEW; NOT_FITTED'))
    table(out / 'recovery_estimands.tsv', estimands)
    configs = [
        dict(configuration='N0_known6', free_mechanisms=14, free_initials=0, free_errors=0, total_free=14, purpose='conditional engineering control'),
        dict(configuration='Ninit', free_mechanisms=14, free_initials=4, free_errors=0, total_free=18, purpose='initial-composition confounding comparison'),
        dict(configuration='Nerror', free_mechanisms=14, free_initials=0, free_errors=2, total_free=16, purpose='error comparison; noisy probability-law study only'),
        dict(configuration='Nall', free_mechanisms=14, free_initials=4, free_errors=2, total_free=20, purpose='intended unknown-nuisance comparison; unresolved coordinates reported honestly')]
    table(out / 'nuisance_configurations.tsv', [dict(**r, execution='NOT_STARTED',
        initial_policy='root-specific Gaussian composition shared through descendants; inherited propagation; no passage reinitialization',
        abundance_policy='counts remain conditioning; no invented seventh nuisance',
        sigma_growth_policy='positive variance target in noisy design; fixed in noiseless engineering controls') for r in configs])
    designs = [
        ('production_full_mask', 'growth+karyotype+processed flow; death off', 'REVIEW_COMPLETE_DESIGN_HELD', 'actual sparse discrete identification and physical flow/count contracts remain open'),
        ('conditional_growth_karyotype', 'observed histories/counts fixed; endpoint growth and integer-cell sampling assumption', 'CANDIDATE_REQUIRES_EXPLICIT_SCOPE_DECISION', 'new observation mask, initial-interface and noise assumptions need separate review and approval'),
        ('physical_experiment_faithful', 'shared raw counts; OLS, stopping and reseeding generated consistently', 'HELD_MISSING_EXPERIMENTAL_CONTRACT', 'passage SOP, count-error dependence and physical sample mapping'),
        ('finite_noise_flow', 'event-level or separately registered approximate flow law', 'HELD_MISSING_EXPERIMENTAL_CONTRACT', 'FCS, gating, retained G0/G1 counts and512/500/density processing'),
        ('live_dead_addition', 'assay-to-live/dead/buffer observation map', 'HELD_HYPOTHETICAL', 'classification, detached material, fractional timing and numerical-loss buffer')]
    table(out / 'observation_design_cases.tsv', [dict(design=d, operator=o, status=s, blocker=b, generated_observations=0, fits=0) for d,o,s,b in designs])
    missing = read(study / 'inputs/contracts/unresolved_inputs.tsv')
    table(out / 'missing_information.tsv', [dict(**r, current_refresh='No new release found in bounded local filename scan; not a global absence claim', step5_status='OPEN_RETAINED') for r in missing])
    table(out / 'estimand_decision_log.tsv', [dict(topic=t, decision=d, reason=r) for t,d,r in (
        ('14+6 inventory', 'retain all20 original identities', 'weak or flat local behavior does not authorize fixing'),
        ('12 prior direct cases', 'retain raw unresolved statuses; add scoped output-map implications', 'new proof route is not a successful old order2 solve'),
        ('sigma_kary', 'free diagnostic plus kernel/distribution targets', 'ideal local identifiability coexists with binary64 identity behavior'),
        ('sigma_growth', 'noise-law variance target only', 'free SD on noiseless means tends to its lower bound and is not valid recovery'),
        ('initial composition', 'four unknown shared root coordinates', 'fresh Gaussian at each passage changes the protocol'),
        ('abundance', 'keep observed counts as declared conditioning', 'physical count uncertainty needs its own model; not an implicit extra parameter'),
        ('numeric geometry', 'diagnostic secants only', 'branch changes and scale disagreement prevent structural-rank interpretation'),
        ('future truth', 'seed144 remains disclosed candidate, not released generating truth', 'do not choose a favorable point or inflate sigma to manufacture recovery'),
        ('future fitting', 'four blind origins plus at most one reserve round of two, proposal only', 'separate cost-based authorization; no current allocation released'),
        ('confidence intervals', 'no1.92/3.84 or routine likelihood-ratio thresholds', 'weighted modality means require separate statistical calibration'))])
    table(out / 'source_switches.tsv', read(study / 'inputs/fixed_settings.tsv'))
    table(out / 'boundary_activity.tsv', read(study / 'inputs/contracts/parameter_activity_rules.tsv') + [dict(
        condition='O=0', inactive_parameters='O2_crit;n_O', existing_map_tested='False',
        policy='Scoped Step4 algebraic forward-dependence argument; not a new optimizer map test', penalty_contribution='None in production culture objective')])
    for name in ('growth_observation_contract.md', 'flow_observation_contract.md', 'conditioning_graph.md'):
        shutil.copy2(study / 'inputs/contracts' / name, out / name)
    (out / 'conditioning_map.md').write_text('''# Production conditioning and observation map

Seed144 uses v2,133 live states22:154, DT0.05 day, ploidy_related death,
fixed oxygen per segment, original pruning/drop support and logistic crowding.
The two root compositions are Gaussian families; descendant fractions come
from the selected parent live state. External initial/target counts remain
conditioning inputs. Infeasible inocula are failures, not silently repaired.
Retained stocks reset at each segment. The selected daily endpoint and the
parent boundary are recorded together; target counts can move both decisions.

Growth score predictions are selected-endpoint log ratios. Physical stored
growth was OLS on raw count records:78/114 differences, maximum0.05657 in the
accepted audit. These operators are not interchangeable. Removing growth from
a score would still leave count information in protocol conditioning.

Karyotype uses a column-normalized integer Gaussian kernel. Flow uses mixing
then normalization with a fixed derived bandwidth; it adds no fitted nuisance.
Root Gaussian tails can underflow; transition pruning and score-floor masks
are retained in branch receipts. A1e-12 score floor is not a normalized sampling
law. Continuous fixed-mask certificates do not establish the full sparse D map.
Exact same-input repeats pass; source sufficient invariant-region bounds address
clamps in real arithmetic, not literal floating-bit equivalence.

The new adapter loads unchanged frozen function definitions and compiles the
byte-identical C++ directly into a private cache. It bypasses only the old timed
compile-lock bootstrap. Production-state, score, selection, source and interface
checks document equivalence; the biological propagation is unchanged.
''')
    reference_count = len(read(study / 'verification/O02/production_reference_checks.tsv'))
    selection_count = len(read(study / 'verification/R01/selection.tsv'))
    table(out / 'full_model_bridge_checks.tsv', [dict(check=k, status='PASS', detail=v) for k,v in (
        ('source_and_SIF_identity', SIF_SHA), ('original20_parameter_map', 'all active; roundtrip and injection checks'),
        ('original_production_saved_state', str(reference_count)+' numeric field comparisons and'+str(selection_count)+' selected segments'),
        ('fresh_process_repeat', 'states, distributions, score and selected indices'),
        ('finite_scope', str(len(ledger))+' full bundles; no fit or data generation'),
        ('piecewise_diagnostics', 'branch and two-scale diagnostics retained; no positive full-D structural claim'))])
    table(out / 'structural_review_gate.tsv', [dict(scope=s, status=g, reason=r) for s,g,r in (
        ('source_and_engineering', 'PASS', 'frozen production reproduction and bounded interfaces verified'),
        ('new_ideal_observer_extension', 'PASS_SCOPED_LOCAL_19', 'three original interior witnesses; ideal continuous L plus initial karyotype'),
        ('actual_sparse_discrete_identification', 'NOT_PROVED', 'finite secants and output inclusion do not close this gap'),
        ('synthetic_design_review', 'REVIEW_COMPLETE_DESIGN_HELD', 'sampling/conditioning and unresolved sparse-map assumptions remain explicit'),
        ('Step6_generation_and_fitting', 'NOT_AUTHORIZED_OR_STARTED', 'requires separate scientific scope and execution confirmation'))])
    (out / 'structural_informed_synthetic_design.md').write_text('''# Structural-informed synthetic design review

Review completed; design remains held. All14 mechanisms and6 nuisances remain.
The ideal19-coordinate extension is local at three interior witnesses, using
continuous absolute live totals across three oxygens and an initial karyotype
moment. Its assumptions differ from the actual sparse selected-endpoint outputs.
The12 old direct order2 attempts remain unresolved; output inclusion/restriction
and the new observer factorization are separate evidence routes.

Retain individual-coordinate and function/combination recovery side by side.
At single constant positive oxygen, alpha_o2*H and mu_hp*H are identifiable-map
candidates while oxygen-response individuals have an exact ambiguity. The real
multi-O study must retain every original coordinate. Exact-zero inactive
coordinates require absolute-error/activity reporting, not generic relative
error. Functions include p(N,O), S(N,m), growth/death curves, K_sigma and complete
predicted distributions. Numerical weak vectors are hypotheses, not proved
compensation or global combinations.

Initial composition has four shared root Gaussian coordinates, with inherited
states at passages. sigma_kary stays a free diagnostic coordinate: its ideal
nonzero derivative does not prevent a production binary64 identity kernel at
seed144. Report raw width and kernel/distribution recovery separately; independent
assay calibration would provide information absent at that resolution. Do not
inflate widths or replace the integer kernel. sigma_growth is an unknown positive
variance only in a later noisy study; do not assess its recovery with noiseless
mean data. N0_known6, Ninit, Nerror and Nall are registered comparisons, not run.
Nall is the intended unknown-nuisance comparison and cannot be replaced silently
by a known-nuisance success. Counts/inocula remain separate conditioning inputs.

A growth+karyotype production-compatible conditional candidate is prepared for
review: fixed declared histories/counts, endpoint-mean growth with explicitly
assumed passage noise, and integer karyotype sampling. This is an estimator
study, not verified wet-lab generation. It requires a new mask/design decision;
it is not automatically released by this report. Keep the remaining modality
averaging/normalization when registering any mask change. Full flow-noise and
live/dead designs remain held pending assay contracts. A processed200-point
flow curve is not200 independent events. Model numerical-loss buffer is not
automatically biological death. Physical-study growth must generate shared raw
counts once and derive OLS, stopping and reseeding consistently.

Seed144 remains a provenance-based candidate with disclosed boundary and assay
resolution issues, not a released truth. Do not choose truths/noise for favorable
recovery or demand another real-data global optimum. Future blind fitting must
exclude generating active values from fit bundles, inject all free nuisance
coordinates from independent streams, select by independently verified objective
only and report boundary, optimization and predictive errors separately.

Proposed later pilot: generator/interface controls, one noisy case, then a cost
benchmark before replicated studies. Four blind origins and at most one reserve
round of two are proposals requiring separate approval. Evaluation/generation
ceilings must follow a measured benchmark; no arbitrary research time cutoff.
Agreement of two origins is optimizer evidence, not a global proof. A noisy fit
can score below generating truth, so old noiseless Q-Qtruth gates do not apply.
Register natural/transformed and zero-boundary errors, function/distribution
errors, variance recovery and optimization quality independently. No1.92/3.84
threshold, standard likelihood-ratio p value or confidence interval follows from
the weighted modality-mean score. No generation, fit, profile or tumor release
is made here.
''')
    (out / 'future_execution_requirements.md').write_text('''# Requirements before any Step6 execution

1. Review the scoped continuous-to-discrete gap and choose an explicit scientific
   conditional design or supply missing experimental contracts; do not declare
   the full sparse production map structurally proved.
2. Confirm observation mask, initial-root interface and conditioning/noise laws.
   Resolve physical growth/SOP/shared-count contracts for faithful experiments;
   flow and live/dead remain held without their own interfaces/assumptions.
3. Freeze provenance-selected truth,20 roles, activity/boundaries, all four nuisance
   configurations, function estimands and seed/stream/blinding manifests.
4. Implement and validate a new generator and unknown-nuisance evaluator after
   separate authorization. Reuse biological functions, not old D0 fixed-error
   objective, mutable counters, controllers or deadline policy.
5. Obtain separate generation/fitting budget and confirm a measured single-case
   benchmark before any replicated batch. Current unused allowance transfers
   neither budget nor authority. No real-data optimum search is prerequisite.
''')
    ref = read(study / 'verification/O02/production_reference_checks.tsv')
    maxgap = max(float(r['max_scaled_gap']) for r in ref)
    counts = Counter(r['status'] for r in geometry)
    sv = read(study / 'verification/G01/singular_values.tsv')
    ranks = {r['scaling']: (r['nominal_rank_1e6'], r['nominal_rank_1e8'], r['nominal_rank_1e10']) for r in sv}
    budget = [dict(key=k, value=v) for k,v in dict(
        ceiling_core_seconds=CAP, settled_core_seconds=spent,
        remaining_core_seconds=CAP-spent, technical_retries_used=retry,
        active_or_queued_own_allocations=0, complete_forward_bundles=len(ledger),
        outside_domain_points=sum(r['status']=='SKIP_OUTSIDE_DOMAIN' for r in points),
        standalone_observer_calls=len(observer), observer_contrasts=len(contrasts),
        new_full_determinants=0, synthetic_observations=0, optimizer_calls=0,
        profiles=0, tumor_analyses=0).items()]
    table(study / 'accounting/final_resource_summary.tsv', budget)
    (study / 'REPORT.md').write_text(f'''# Step5 bounded execution report

Collection: COMPLETE. Source/engineering: PASS. Design:
REVIEW_COMPLETE_DESIGN_HELD. Step6 generation/fitting was not started.

## Evidence and results

- Preserved14 mechanism and6 nuisance identities;120 case-coordinate map rows.
  Rehashed prior760-file Stage B and426-file B6 packages unchanged. Retained all
  12 unresolved direct order2 runs and their old spend/failure/retry records.
- All three prespecified interior witnesses passed the new ideal observer block
  extension: local rank19 for14 mechanisms,4 unknown root coordinates and
  sigma_kary, using ideal continuous L plus initial karyotype. Used{len(contrasts)}
  of6 allowed contrasts; no new full determinant or state-case solve. This is
  local fixed-support/mask evidence, not global or actual sparse-D identification.
- Verified pinned production SIF/source, ploidy_related death, original20-active
  mapping and independent-score agreement. Original seed144 saved-state comparison
  maximum scaled difference{maxgap:.6g}; fresh-process repeat discrepancies were
  zero. Production score3.85564658889903 was reproduced, not optimized.
- Executed{len(ledger)} of74 allowed full two-background/history bundles;
  {sum(r['status']=='SKIP_OUTSIDE_DOMAIN' for r in points)} outside-domain points
  were skipped, never clipped/replaced. Anchor propagation0.333 seconds excludes
  compilation. Allocated CPU and compilation are counted separately below.
- Completed{len(observer)} standalone observer calls and4 observer-only saved-
  latent perturbation passes. At full-precision seed144 sigma_kary, all kernel
  offdiagonals underflow and the binary64 kernel is identity. Registered local
  sigma probes are numerically flat. Wider registered controls show observer
  response; no width was substituted as a new truth.
- Production geometry has19 mean columns; sigma_growth is outside that map.
  Status counts: {dict(counts)}. Nominal secant ranks at relative thresholds
  1e-6/1e-8/1e-10: {ranks}. These are scaling/branch-sensitive diagnostics, not
  structural proofs, confidence intervals or proof of parameter compensation.
- Sigma_growth positive-Gaussian score wiring passed literal residual controls.
  It is a future noise-law variance target, not recoverable from exact means by
  treating a lower-bound SD optimum as success.

## Design decision

All20 original coordinates remain visible and intended unknown in Nall. Known-
nuisance configurations are controls. Report individual targets and meaningful
functions/distributions separately, especially sigma_kary at seed144. Initial
Gaussian families are shared at roots; actual descendant states and counts
remain propagated/conditioned under the audited production interface.

The full physical/noisy study remains held: passage/count-error contracts,
physical sample mapping, raw flow/gating/retained counts and live/dead operator
remain open. The bounded material scan did not supply them. A narrower conditional
growth+karyotype candidate is documented for an explicit scientific decision;
it is not automatically released. The full sparse/discrete positive structural
gap is still unproved. No demand for a new real-data global optimum is introduced.

## Resources, policy and preservation

Seven own allocations completed with Slurm0:0, internal application PASS and
dependency/write-trace policy PASS. Conservatively settled{spent} allocated
core-seconds ({spent/3600:.6f} core-hours), versus12 core-hours authorized.
Technical retries:{retry}/1. No own queued/running allocation remains. MaxRSS
in the raw Slurm snapshot is scheduler-sampled, not a guaranteed peak. The actual
scientific workers/children were guarded and traced; launchers are outside that
coverage. Denied terminal probes stayed denied and are explicitly classified.
Site/budget-derived reservations are not arbitrary research deadlines. No
application concurrency cap, requeue or optimizer was used. Old530/50 core-second
Stage B/B6 spend and old retry records remain historical, not reset or charged
to this independent package.

Files are English and use Python/R/C++/shell with TSV/CSV/RDS/Markdown evidence;
no Java, JavaScript or JSON workflow/file was introduced. New scientific math
ran only in xxlarge: production checks inside the pinned SIF; outward-arithmetic
observer/geometry work in the separately recorded host runtime. No old controller,
task or archive was changed. No synthetic observations, random draws, fitting,
profiles, tumor analysis or author contact occurred. Unused allowance authorizes
no extra work. Collection completion and scientific design HOLD are separate.

Read reports/structural_review_gate.tsv, recovery_estimands.tsv,
parameter_roles.tsv, case_inference_addendum.tsv, nuisance_configurations.tsv,
observation_design_cases.tsv, missing_information.tsv and
structural_informed_synthetic_design.md for scoped next decisions.
''')
    table(out / 'collection_gate.tsv', [dict(status='PASS_COLLECTION_COMPLETE',
        phases=7, full_forward_bundles=len(ledger), settled_core_seconds=spent,
        retries_used=retry, design_status='REVIEW_COMPLETE_DESIGN_HELD',
        generation_and_fitting='NOT_STARTED')])
    code = Path(__file__).parent
    # Check final metadata sources without importing any scientific module.
    tests = []
    for p in sorted(code.glob('*.py')):
        tree = ast.parse(p.read_text())
        bad = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import): bad.extend(a.name.split('.')[0] for a in node.names)
            if isinstance(node, ast.ImportFrom): bad.append((node.module or '').split('.')[0])
        assert not set(bad) & {'json', '_json', 'jpype', 'jnius'}
        tests.append(dict(check='final_metadata_syntax_imports_'+p.name, passed=True))
    table(out / 'collection_static_tests.tsv', tests)
    table(out / 'final_implementation_manifest.tsv', [dict(project_path=str(p.relative_to(project)), bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(code.iterdir()) if p.is_file()])
    forbidden = [p for p in study.rglob('*') if p.is_file() and (p.suffix.lower() in ('.json','.jsonl','.js','.mjs','.cjs','.java','.class','.jar') or '__pycache__' in p.parts)]
    assert not forbidden, forbidden
    if retry:
        retry_record = read(study / 'accounting/technical_retry.tsv')
        assert len(retry_record) == 1 and retry_record[0]['scientific_reruns'] == '0'
        retry_record[0]['status'] = 'PASS_COLLECTION_RETRY; NO_SCIENTIFIC_RERUN'
        table(study / 'accounting/technical_retry.tsv', retry_record)
    table(study / 'collection_manifest.tsv', [dict(path=str(p.relative_to(study)), bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(study.rglob('*')) if p.is_file() and p.name not in ('collection_manifest.tsv','CLOSED')])
    (study / 'CLOSED').write_text(sha(study / 'collection_manifest.tsv')+'\nPASS_COLLECTION_COMPLETE\nREVIEW_COMPLETE_DESIGN_HELD\n')
    print('STEP5_CLOSED', spent, 'allocated core-seconds;', len(ledger), 'forward bundles; DESIGN_HELD')


if __name__ == '__main__':
    main()
