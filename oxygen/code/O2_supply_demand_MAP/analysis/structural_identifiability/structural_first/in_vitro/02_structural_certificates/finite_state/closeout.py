"""Metadata-only stopped-package collection; never imports a mathematical runtime."""
import ast
from datetime import datetime
from pathlib import Path
import sys
from zoneinfo import ZoneInfo
from protocol import read, table, sha, kv, roles_valid, verify_archive, CAP, MECHANISMS, INITIALS, NUISANCES

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve()
assert not (study/'COMPLETE').exists(),'Do not rewrite a closed package'
jobs=read(study/'accounting/jobs.tsv');roles=read(study/'inputs/parameter_roles.tsv');roles_valid(roles)
terminal={'COMPLETED','FAILED','TIMEOUT','CANCELLED','NODE_FAIL','OUT_OF_MEMORY','PREEMPTED','BOOT_FAIL','DEADLINE'}
assert len(jobs)==23 and all(j['job_id'].isdigit() and j['state'].split()[0] in terminal for j in jobs)
assert not (study/'accounting/owned_queue.txt').read_text().strip()
assert not any(int(j['reserved_core_seconds']) for j in jobs)
used=sum(int(j['actual_core_seconds']) for j in jobs);assert used<=CAP
retry=kv(study/'accounting/technical_retry.tsv');assert retry['used']=='1' and retry['maximum']=='1'
failed={j['task']:j for j in jobs if j['state']!='COMPLETED' or j['exit_code']!='0:0'}
assert set(failed)=={'B2','B6'}
assert failed['B2']['exit_code']=='1:0' and failed['B6']['exit_code']=='77:0'
assert 'ModuleNotFoundError' in (study/'verification/B6/console.log').read_text()
assert 'from algebra import decode, independent_minor_check' in (study/'verification/B6/console.log').read_text()
assert read(study/'verification/B6/policy_gate.tsv')[0]['status']=='FAIL'
assert not (study/'verification/B6/boundary_gate.tsv').exists()
assert not (study/'verification/B6/application_result.tsv').exists()
table(study/'accounting/stop_record.tsv',[dict(status='STOPPED_NO_TECHNICAL_RETRY_REMAINING',
    triggering_task='B6',job_id=failed['B6']['job_id'],reason='Import before private runtime initialization; missing final module inventory',
    additional_allocations_authorized=False,remaining_technical_retries=0)])

preserved=[]
for name in ('frozen_source_verification.tsv','protected_worktree_hashes.tsv'):
    for r in read(study/'inputs'/name):
        p=project/r['path'];actual=sha(p) if p.is_file() else 'MISSING'
        preserved.append(dict(group=name,path=r['path'],original_sha256=r['sha256'],current_sha256=actual,unchanged=actual==r['sha256']))
table(study/'reports/preservation_verification.tsv',preserved)
assert all(r['unchanged'] for r in preserved),'Preservation discrepancy requires explicit review'
old_count=verify_archive(project/'oxygen/results/structural_identifiability/20261116_culture_step04_initial_observer_02',
    '9693b628dc0afa3f7e02a62023588298620f33f4b6ab9a46b01d2802c710f309')
assert old_count==69
source_checks=[]
for m in sorted(study.glob('source_manifest_r*.tsv')):
    for row in read(m):
        p=study/row['path'];assert p.stat().st_size==int(row['bytes']) and sha(p)==row['sha256']
    source_checks.append(dict(manifest=m.name,sha256=sha(m),files=len(read(m)),status='PASS'))
table(study/'reports/release_manifest_verification.tsv',source_checks)
table(study/'implementation_manifest.tsv',[dict(project_path=str(p.relative_to(project)),bytes=p.stat().st_size,sha256=sha(p))
    for p in sorted(Path(__file__).parent.iterdir()) if p.is_file()])

forbidden={'.json','.json5','.jsonl','.js','.mjs','.cjs','.java','.jar','.ipynb'}
assert not [p for p in study.rglob('*') if p.is_file() and p.suffix.lower() in forbidden]
assert not [p for p in Path(__file__).parent.iterdir() if p.is_file() and p.suffix.lower() in forbidden]
for p in Path(__file__).parent.glob('*.py'):ast.parse(p.read_text())

workers=[];results=[];attempts=[];nominations=[]
for j in jobs:
    task=j['task'];directory=study/'verification'/task
    phase='B2_native' if task=='B2' else task
    policy=read(study/'verification'/phase/'policy_gate.tsv')[0]
    app='MISSING'
    if (directory/'application_result.tsv').exists():app=read(directory/'application_result.tsv')[0]['status']
    if task=='B0':
        assert read(directory/'semantic_gate.tsv')[0]['status']=='PASS'
        app='PASS_SEMANTIC_GATE_NO_GENERIC_APPLICATION_FILE'
    workers.append(dict(task=task,job_id=j['job_id'],slurm_state=j['state'],slurm_exit=j['exit_code'],
        application_status=app,policy_status=policy['status'],policy_scope=phase,
        successful_trace_events=policy['successful_events'],allocated_core_seconds=j['actual_core_seconds']))
    if task not in ('B2','B6'):
        assert j['state']=='COMPLETED' and j['exit_code']=='0:0' and policy['status']=='PASS'
        assert app=='PASS' or task=='B0' and app=='PASS_SEMANTIC_GATE_NO_GENERIC_APPLICATION_FILE'
    if not task.startswith(('B01_','B02_','B03_','B04_','B05_','B06_')):continue
    r=read(directory/'case_result.tsv')[0];results.append(r)
    assert r['native_tool_initial_support']=='not_verified' and r['Step_5']=='False' and r['fits']=='0'
    a=read(directory/'certificate_attempts.tsv')
    calls=read(directory/'determinant_call_ledger.tsv')
    assert len(calls)==int(r['determinant_evaluations'])==len(a)<=8
    assert all(x['status']=='returned' for x in calls)
    for record in a:attempts.append(dict(task=task,**record))
    registry={int(x['index']):x for x in read(directory/'row_registry.tsv')}
    columns=read(directory/'column_registry.tsv')
    for c in read(directory/'candidate_registry.tsv'):
        ids=[int(v) for v in c['rows'].split(',')];oxygens=sorted({registry[i]['oxygen'] for i in ids})
        fullcols=c['columns'].split(',');nomination_only_single_O=len(oxygens)==1
        oxygen_shape_selected={'O2_crit','n_O'}.issubset(fullcols)
        original_d=len(columns)
        singular=nomination_only_single_O and oxygen_shape_selected and len(ids)>original_d-2
        nominations.append(dict(task=task,candidate=c['candidate'],size=len(ids),selected_oxygens=','.join(oxygens),
            all_mean_unknowns=original_d,exact_single_O_upper_bound=original_d-2 if nomination_only_single_O else '',
            candidate_singular_by_known_family=singular,
            note='Midpoint pivots are nomination only; rounded full rank cannot override exact invariance'))
    expected='CERTIFIED_RANK_16_WITH_EXACT_SYMMETRY' if task.startswith('B01_') else 'CERTIFIED_FULL_LOCAL' if task.startswith('B04_') else 'UNRESOLVED_REGISTERED_SCOPE'
    assert r['status']==expected
assert len(results)==18 and len(attempts)==102 and len(nominations)==33
assert all(r['candidate_singular_by_known_family'] for r in nominations if r['task'].startswith(('B02_','B03_','B05_','B06_')))
table(study/'reports/worker_outcomes.tsv',workers)
table(study/'reports/structural_results.tsv',sorted(results,key=lambda r:(r['case'],r['witness'])))
table(study/'reports/certificate_attempts.tsv',attempts)
table(study/'reports/nomination_diagnostics.tsv',nominations)
policy_records=[]
for p in sorted((study/'verification').glob('*/policy_gate.tsv')):
    policy_records.append(dict(phase=p.parent.name,**read(p)[0]))
assert len(policy_records)==25 and sum(r['status']=='PASS' for r in policy_records)==24
assert [r['phase'] for r in policy_records if r['status']!='PASS']==['B6']
table(study/'reports/runtime_policy_summary.tsv',policy_records)
assert len(read(study/'verification/B2_native/kernel_calls.tsv'))==9
assert read(study/'verification/B2_review/operator_gate.tsv')[0]['status']=='PASS'

statuses=[]
oxygen_family={'O2_crit','n_O','alpha_o2','mu_hp'}
for case in ('B01','B02','B03','B04','B05','B06'):
    for name in MECHANISMS+NUISANCES:
        if name=='sigma_growth':
            status='OUTSIDE_DETERMINISTIC_MEAN';evidence='Variance-only coordinate; no noise law identifiability analysis';global_scope='not_analyzed'
        elif name=='sigma_kary':
            status='UNRESOLVED_OBSERVER_COORDINATE' if case in ('B03','B05') else 'ABSENT_FROM_DECLARED_OUTPUT'
            evidence='No valid 19-column direct certificate' if case in ('B03','B05') else 'Retained inventory; assay width not an unknown mean coordinate in this case'
            global_scope='not_analyzed'
        elif case=='B06' and name in INITIALS:
            status='KNOWN_IN_DIAGNOSTIC';evidence='Original nuisance retained but conditioned as known in B06';global_scope='not_analyzed'
        elif case=='B01' and name in oxygen_family:
            status='NONIDENTIFIABLE_EXPLICIT_SINGLE_O_FAMILY';evidence='B1 exact response invariance plus B2 event/operator correspondence; original parameters vary';global_scope='continuous_local_family_nonuniqueness'
        elif case=='B01':
            status='LOCAL_QUOTIENT_IDENTIFIABLE';evidence='B01 rank16 and exact quotient (alpha*H,mu_hp*H,other mechanisms,initial coordinates)';global_scope='only ideal full-initial Gaussian inverse if an initial coordinate'
        elif case=='B04':
            status='JOINT_LOCAL_IDENTIFIABLE';evidence='Three direct rank18 order4 certificates for 14 mechanisms plus 4 unknown initial coordinates';global_scope='not_claimed'
        else:
            status='LOCAL_BY_OBSERVED_L_SUBMAP'
            evidence='Logical implication of B04: same model/inputs/initial convention and L is an observed function in this case; not a direct order2 certificate'
            global_scope='only ideal full-initial Gaussian inverse for B02 initial coordinates; otherwise not_claimed'
        statuses.append(dict(case=case,parameter=name,status=status,evidence=evidence,global_scope=global_scope,
            witnesses='W1,W2,W3',domain='conditional continuous 133-state finite support; separate strict-mask neighborhoods',
            actual_sparse_culture_inference='not_established',Step_5='not_released'))
assert len(statuses)==120
table(study/'reports/project_parameter_status.tsv',statuses)

# Slurm child-step rounding is accounted by max(parent,steps), not summed twice.
raw=(study/'accounting/slurm_snapshot.txt').read_text().splitlines();resources=[]
for j in jobs:
    entries=[line.split('|') for line in raw if line.split('|')[0]==j['job_id'] or line.split('|')[0].startswith(j['job_id']+'.')]
    rss=[r[5] for r in entries if r[5]]
    resources.append(dict(task=j['task'],job_id=j['job_id'],requested_cpus=1,allocated_cpus=j['allocated_cpus'],
        elapsed_seconds=j['elapsed_seconds'],conservative_core_seconds=j['actual_core_seconds'],MaxRSS_reports=';'.join(rss),qos='xxlarge'))
table(study/'accounting/resource_summary.tsv',resources)
table(study/'accounting/budget_final.tsv',[dict(cap_core_seconds=CAP,used_core_seconds=used,used_core_hours=used/3600,
    unspent_core_hours=(CAP-used)/3600,live_reservation_core_seconds=0,technical_retries_used=1,
    technical_retry_max=1,failed_allocations=2,owned_jobs=23,owned_active_jobs=0)])
acceptance=[
 dict(gate='source_SIF_actual_death_mode',status='PASS',evidence='B0; actual ploidy_related; SIF digest verified'),
 dict(gate='arithmetic_response_Taylor_controls',status='PASS',evidence='B1:124 exact/interface controls'),
 dict(gate='operator_source_derivative_bridge',status='PASS_AFTER_ONE_TECHNICAL_REVIEW',evidence='B2_native9 returned kernels; B2_review; original B2 failure retained'),
 dict(gate='single_O_local_quotient',status='PASS_CONDITIONAL',evidence='B01 W1-W3 rank16 with exact two-dimensional family'),
 dict(gate='multi_O_live_mean_local_rank',status='PASS_CONDITIONAL',evidence='B04 W1-W3 rank18; order4'),
 dict(gate='direct_B02_B03_B05_B06_minor_tests',status='UNRESOLVED_FINITE_SCOPE',evidence='24 nominated minors remain singular under one-O family; 12 case/witness runs exhausted scope'),
 dict(gate='complete_boundary_and_discrete_audit',status='NOT_EXECUTED_TECHNICAL_STOP',evidence='B6 failed before runtime boot; no boundary_gate or application_result'),
 dict(gate='saved_certificate_additional_closeout_review',status='NOT_EXECUTED_TECHNICAL_STOP',evidence='Per-case independent certificate checks passed; added B6 saved-matrix rereview did not execute'),
 dict(gate='runtime_policy',status='PARTIAL_HOLD',evidence='24 scientific phases have PASS; B6 lacks completed module inventory and has FAIL; no compliant success claimed for B6'),
 dict(gate='metadata_preservation_collection',status='PASS',evidence='Old69-file archive, all frozen manifests and protected worktree files verified; queue empty'),
 dict(gate='Stage_B_delivery',status='CLOSED_WITH_PARTIAL_RESULTS_AND_HOLD',evidence='Bounded stop is honored; no technical retry remains'),
 dict(gate='Step_5',status='NOT_RELEASED',evidence='No mapping/design release, fitting or new synthetic data generation')]
table(study/'acceptance.tsv',acceptance)
table(study/'package_status.tsv',[dict(status='CLOSED_WITH_PARTIAL_RESULTS_AND_HOLD',scientific_all_gates_pass=False,
    collection_complete=True,Step_5=False,new_fits=0,new_trajectories=0,new_synthetic_datasets=0)])

now=datetime.now(ZoneInfo('America/New_York')).isoformat(timespec='seconds')
report=f'''# Culture Step 4 Stage B: finite analysis and stopped closeout

Collected at {now}. The package is **closed with partial structural results and a HOLD**. Collection completion does not mean every scientific gate passed. No Step 5 release is made.

## Results and their exact scope

All 18 registered case/witness jobs returned with application PASS, policy PASS and Slurm COMPLETED/0:0. Six produced positive certificates; twelve exhausted their registered candidate/precision scope. The question is the conditional continuous finite-support model, with 133 live states N=22:154, separate positive oxygen experiments, fixed drop support and strict pruning masks. It is not the actual sparse, inherited-passage experimental protocol.

| Case | Output / initial convention | Direct result at W1,W2,W3 |
|---|---|---|
| B01 | One O=0.2%; ideal full live state; four shared initial unknowns | Rank16/18; matching exact two-dimensional invariance |
| B02 | O=0.2%,1%,20.5%; ideal full live state; four initial unknowns | Registered order2 minors unresolved |
| B03 | Same three O; L and complete normalized karyotype; four initial unknowns and sigma_kary | Registered order2 minors unresolved |
| B04 | Same three O; ideal continuous absolute L only; four initial unknowns | Rank18/18 through order4; joint local identifiability |
| B05 | Conditional B03 plus fixed-interface normalized flow | Registered order2 minors unresolved |
| B06 | Conditional B02 with known initial distributions | Registered order2 minors unresolved |

For a single positive O, let H=H(O;C,n), choose nearby positive C',n' and write H'=H(O;C',n'). Set alpha'=alpha*H/H' and mu_hp'=mu_hp*H/H'. All other coordinates are unchanged. Lambda, mu, p, survival, masks, operator coefficients and retained inflows are identical. Thus O2_crit, n_O, alpha_o2 and mu_hp are individually nonunique under this experiment. Effective alpha*H and mu_hp*H, the other ten mechanisms and four initial coordinates form a locally identifiable 16-dimensional quotient near each witnessed branch. A positive rank16 minor and the exact upper argument are both required; the four original coordinates were retained, not fixed out of the inventory.

The three B04 certificates establish joint local identifiability of all14 mechanisms and the four shared Gaussian initial coordinates under the declared ideal continuous absolute-live-total experiments. These are local statements near three distinct interior witnesses with their own masks. There is no proof over the entire production box or of global mechanism uniqueness. Orders0:4 use Taylor coefficients with the factorial convention preserved; the informative rows, scaled matrices, finite determinant balls and independently reconstructed event-route Neumann checks are saved in each case directory.

B02/B03/B05/B06 direct certificate attempts must remain marked unresolved. Their deterministic midpoint nomination selected only one O in both registered candidates at every witness. Exact single-O invariance bounds such candidate ranks two below the original mean dimension, so a nominally full midpoint pivot is not a valid structural rank. Increasing precision cannot cure that selection. No new candidates, orders or witnesses were added after exhaustion. This is a known nomination limitation, not evidence of nonidentifiability of the full multi-input maps. See `reports/nomination_diagnostics.tsv`.

At the continuous output-map level, B04 also has a simple consequence: B02/B03/B05 observe L directly or as sum(x); their same18 mechanism/initial coordinates are therefore locally determined by that observed submap, independently of an assay-width coordinate. B06 conditions on known initial values and observes an output containing L, so the mechanism restriction inherits this local result. These are logical output-inclusion/restriction consequences of B04, **not successful direct order2 certificates** and not a newly computed higher-order B02/B03/B05/B06 test. Sigma_kary remains unresolved in B03/B05. No claim of a full19-dimensional certificate is made.

## Initial distributions and observation nuisance coordinates

The original14+6 inventory is preserved in `inputs/parameter_roles.tsv`; case-specific conclusions appear in the120 rows of `reports/project_parameter_status.tsv`.

- Four initial Gaussian mean/SD coordinates remain unknown and shared across the ideal independent oxygen experiments in B01-B05; there is one family per2N/4N starting background. B06 is explicitly a known-initial diagnostic. No fresh Gaussian is imposed at an actual inherited passage.
- For an ideal true full initial vector, three adjacent positive weights yield s^2=-1/log(x_minus*x_plus/x_center^2) and mean=center+s^2*log(x_plus/x_minus)/2. Normalization cancels. The separate complete inversion establishes initial-family global uniqueness only for those true full-initial outputs; it does not apply to live-only, blurred assays or sparse measurements.
- Sigma_kary is an unknown shared observer coordinate in B03/B05 and unresolved by this bounded direct analysis. It is absent from the other declared mean outputs, while retained in the inventory.
- Sigma_growth is the growth likelihood SD. It does not enter the deterministic mean map; its variance-law identifiability was not studied and it was not treated as an identified/fitted mean parameter.
- Flow has fixed audited smoothing and normalization; no new free flow-width nuisance was invented. Actual physical calibration and sampling assumptions remain held.

## Production source and numerical validation

The actual seed144 cfg and run parameters use **ploidy_related**, verified within the hash-checked production SIF. The new death law is mu_N=mu_hp*H*(N/44)^gamma_mu. Older inherited diploid_NULL descriptions were incorrect for that production configuration; the old closed Stage A archive was not edited. See `inputs/semantic_correction.md` and `verification/B0/`.

SIF SHA256: `dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372`.

B1 passed124 arithmetic, response-derivative, small-system Taylor, Gaussian, quotient and collector controls. The method uses pinned SymPy1.12, python-flint0.9.0/Arb and a fixed-model exporter, not a new general symbolic SI solver. The native initial-manifold support status of earlier general tools remains **not verified**. No prior native control was rerun.

B2 saved exactly9 instantaneous production G kernels (3 witnesses by3 oxygens), not trajectories. All133-state operators, rates and buffer ledgers passed componentwise1e-12-scaled comparisons in the accepted saved-result review. Event derivatives for all14 mechanisms were independently checked; biology, probability-pruning, lower-drop and upper-drop contributions are separate. Retained buffer is not pure biological death. Numeric agreement alone is not a symbolic equivalence proof; source event identities, response derivatives and the fixed recurrence are separate obligations and evidence routes.

The original B2 allocation failed in a generic symbolic validator after all9 kernels returned. Its one authorized technical retry reviewed the saved kernels, with no native repeats, using an equivalent explicit nonnegative event representation; equations, masks and tolerance stayed unchanged. Original failure and review are both retained.

## Stop event and incomplete gates

B6 job21971224 failed before runtime boot: `boundary_audit.py` imported `algebra`, which imports flint, before `boot()` installed the private runtime path. This was an implementation import-order error. Slurm records FAILED/77:0 because the wrapper's policy audit also detected a missing completed Python module inventory; the application itself exited1. There is no boundary result, application PASS, saved-matrix additional review or complete runtime-policy acceptance for this phase. No forbidden dependency was reported in its trace, but incomplete coverage is not policy PASS.

The single technical retry had already been consumed by B2. The confirmed stopping rule therefore prohibited another computation or allocation. B6 was not repaired and rerun, and the failed release is immutable. Per-case primary and independently reconstructed certificates already returned successfully; the additional B6 saved-enclosure rereview remains unexecuted. No claim is made that the B6 failure disproves them.

Three response-level boundary controls had already passed in B1: p_misseg=0 removes the k_o_mis derivative, buffer_beta=0 removes the buffer_n_exp derivative, and O=0 gives H=1 for positive shape parameters. The **complete compositional boundary and production discrete/passage audit is not closed**. The planned constant-O real-arithmetic induction must not be reported as an executed B6 proof. There is no positive identifiability certificate for DT=0.05, clamps, target-day selection and inherited passages, nor any literal floating-point equivalence claim.

Production target-count selection, first-index ties, inoculum fallback/failure and inherited fractions remain recorded in frozen source. Observed final counts used for conditioning must not be mistaken for held-out measurements. Observed OLS growth versus scored endpoint log-ratio, physical flow processing and eligible live/dead denominators remain experimental/design issues. This package does not resolve them.

## Resources, preservation and next boundary

All23 owned allocations used xxlarge and one CPU/16GiB. Settled conservative use, including both failed jobs, is **{used} core-seconds = {used/3600:.6f} core-hours**, below64 core-hours. Active/queued commitments are zero. No application concurrency cap, GNU timeout, subprocess timeout or arbitrary research deadline was used. The unavoidable xxlarge site maximum was12hours; scheduler reservations were derived from that limit and the remaining aggregate authorization with settlement headroom. The budget is a ceiling, not permission to spend the unspent balance after the retry stop.

There were102 registered determinant evaluations across33 nominated candidates and18 cases (27 candidates actually evaluated), below the36-candidate/144-evaluation ceilings; exactly9 native instantaneous kernels and zero forward trajectories, new synthetic culture datasets or fits. MaxRSS and parent/step rounding are retained in `accounting/resource_summary.tsv`. The old69-file archive, all six immutable release manifests,76 frozen production-source records and147 pre-existing dirty/untracked worktree files were rehashed unchanged. No old controller or unrelated task was changed.

All scientific result-bearing phases except failed B6 have successful file/process traces and module inventories. The result tree has no Java/JavaScript/JSON-format files. B6 is explicitly a coverage HOLD. English artifacts and raw failure records are retained. `COMPLETE` seals the collected archive only; inspect `package_status.tsv` and `acceptance.tsv` before using it as scientific evidence.

Step5 remains unreleased. A future separately confirmed action would first move the algebra import after guarded runtime initialization, add an import-order static control, and complete only the missing B6 boundary/discrete and saved-certificate reviews under an explicit new retry/resource allowance. Repairing direct candidate nomination or identifying sigma_kary is a distinct proposal; it must not silently reuse the exhausted candidate/retry scope. No practical recovery, profiles or tumor fitting may follow automatically from the present conditional ranks.
'''
(study/'REPORT.md').write_text(report)
doc=project/'oxygen/results/analysis_documentation/frozen_inputs/in_vitro/structural_certificates/culture_structural_step04_stageB_finite_execution_20261116.md'
doc.write_text(f'''# Culture Step 4 Stage B execution record

Date: {now}. Status: **closed with partial results and HOLD**; Step5 not released.

The64-core-hour bounded package reconciled the production SIF and actual ploidy_related death mode, passed124 prerequisite controls, and accepted9 saved instantaneous native kernels after its one authorized technical review. It then completed all18 predeclared case/witness runs. B01 obtained rank16 with an exact single-O parameter family at all three witnesses. B04 obtained full local rank18 for14 mechanisms and4 unknown initial coordinates under three ideal constant-oxygen absolute-live-total experiments. These continuous finite-support findings do not establish practical identifiability or global production/passage identifiability.

B02/B03/B05/B06 direct order2 tests exhausted their registered minors; both candidates sampled one oxygen and therefore could not break the exact oxygen-response family. Sigma_kary remains unresolved; sigma_growth is outside the deterministic mean map. All14 mechanisms and6 nuisance rows are preserved.

The B6 closeout worker imported flint before private runtime initialization and failed. Its missing final module inventory also prevents policy PASS. The one retry was already used, so no further computation was submitted. Complete boundary/discrete composition and added saved-certificate review remain unexecuted. Existing per-case independent certificate checks are retained; neither failure nor a missing review is called nonidentifiability.

Conservative use: {used} core-seconds ({used/3600:.6f} core-hours),23 terminal owned jobs, no active commitments; xxlarge, no application concurrency cap or discretionary timer. Both failures remain in the ledger. No fits, new culture datasets, trajectories, Step5 or tumor work ran.

Detailed English report: [REPORT.md](../oxygen/results/structural_identifiability/{study.name}/REPORT.md).
See `acceptance.tsv`, `reports/project_parameter_status.tsv`, `reports/nomination_diagnostics.tsv`, `reports/worker_outcomes.tsv`, `accounting/budget_final.tsv` and the raw verification directories. The archive completion marker seals collection, not a scientific all-gates PASS. Pre-existing147 files and frozen production inputs were preserved.
''')
manifest=[]
for p in sorted(study.rglob('*')):
    if p.is_file() and p.name not in ('collection_manifest.tsv','COMPLETE'):
        manifest.append(dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)))
table(study/'collection_manifest.tsv',manifest)
(study/'COMPLETE').write_text(sha(study/'collection_manifest.tsv')+'\nCLOSED_WITH_PARTIAL_RESULTS_AND_HOLD\nCOLLECTION_ONLY_NOT_SCIENTIFIC_ALL_GATES_PASS\n')
print('STOPPED_PACKAGE_COLLECTED',len(manifest),'files;',used,'core-seconds; Step5 held')
