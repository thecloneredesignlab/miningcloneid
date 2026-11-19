"""Independent standard-library collection; never regenerates scientific rows."""
from datetime import datetime
from pathlib import Path
import sys
from zoneinfo import ZoneInfo
from common import *

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();code=Path(__file__).parent
assert not (study/'CLOSED').exists()
verify(study)
for r in read(study/'retry_source_manifest.tsv'):
    p=study/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256']
for r in read(study/'reports/failed_P03_manifest.tsv'):
    p=study/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256']
old=project/'oxygen/results/structural_identifiability'/ORIGIN
assert sha(old/'collection_manifest.tsv')==ORIGIN_SHA==(old/'COMPLETE').read_text().splitlines()[0]
for r in read(old/'collection_manifest.tsv'):
    p=old/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256']
preservation=[]
for r in read(study/'preservation_baseline.tsv'):
    p=project/r['project_path'];current=sha(p)
    preservation.append(dict(**r,current_sha256=current,unchanged=current==r['sha256']))
assert all(r['unchanged'] for r in preservation)
table(study/'reports/preservation_final.tsv',preservation)
jobs=read(study/'accounting/jobs.tsv')
assert {j['task'] for j in jobs}=={'P01','P02','P03','P03_R1'} and len(jobs)==4
assert not (study/'accounting/owned_queue.txt').read_text().strip()
outcomes=[]
for j in jobs:
    folder=study/'verification'/j['task'];exits=read(folder/'phase_exit_codes.tsv')[0]
    app=read(folder/'application_result.tsv')[0]['status'];policy=read(folder/'policy_gate.tsv')[0]['status']
    if j['task']=='P03':
        assert j['state']=='FAILED' and j['exit_code']=='77:0' and app=='PASS' and policy=='FAIL'
        assert exits==dict(application_exit_code='0',policy_exit_code='77',source_after_exit_code='0')
    else:
        assert accepted(j['state'],j['exit_code'],app,policy,'PASS')
        assert all(v=='0' for v in exits.values())
    assert int(j['allocated_cpus'])==1 and not int(j['reserved_core_seconds'])
    outcomes.append(dict(task=j['task'],job_id=j['job_id'],slurm_state=j['state'],slurm_exit=j['exit_code'],
        application_exit=exits['application_exit_code'],policy_exit=exits['policy_exit_code'],source_after_exit=exits['source_after_exit_code'],application=app,policy=policy))
table(study/'reports/worker_outcomes.tsv',outcomes)
inventories=[]
excluded={'json','_json','simplejson','orjson','ujson','jpype','jnius','execjs','py_mini_racer','symbtools','strike_goldd','rJava','V8','jsonlite','rjson','RJSONIO'}
for pattern,field in (('imported_modules.tsv','module'),('loaded_namespaces.tsv','namespace')):
    for p in sorted((study/'verification').rglob(pattern)):
        rows=read(p)
        assert rows and not any(r[field].split('.')[0] in excluded for r in rows)
        inventories.append(dict(path=str(p.relative_to(study)),sha256=sha(p),entries=len(rows),status='PASS'))
for p in (study/'verification').rglob('blocked*'):
    if p.is_file() and p.stat().st_size:
        assert p.parent.name in ('P03','P03_R1') and all(line=='WRITE_OUTSIDE_PHASE\topen\t/dev/tty' for line in p.read_text().splitlines()),'Unexpected blocked scientific activity: '+str(p)
assert len(inventories)==5  # Original SIF attempt, accepted retry, P01/P02 and expected-failure child.
table(study/'reports/complete_inventory_review.tsv',inventories)
P01=study/'verification/P01';P02=study/'verification/P02';P03=study/'verification/P03_R1'
assert read(P01/'runtime_gate.tsv')[0]['status']=='PASS'
for p in ('saved_certificate_gate.tsv','boundary_gate.tsv'):assert read(P02/p)[0]['status']=='PASS'
assert read(P03/'passage_gate.tsv')[0]['status']=='PASS'
certs=read(P02/'saved_certificate_inventory.tsv');matrices=read(P02/'saved_certificate_review.tsv');boundaries=read(P02/'boundary_activity.tsv');fixtures=read(P03/'passage_contract_tests.tsv')
assert len(certs)==len({(r['case'],r['witness']) for r in certs})==6
assert len(matrices)==len({(r['case'],r['witness'],r['route']) for r in matrices})==12
assert len(boundaries)==3 and len(fixtures)==9 and all(r['passed']=='TRUE' for r in fixtures)
roles=read(study/'inputs/parameter_roles.tsv');roles_valid(roles)
assert (P03/'image_sha256.txt').read_text().strip()==SIF_SHA
for p in study.rglob('*'):
    if p.is_file():assert p.suffix.lower() not in {'.json','.jsonl','.json5','.js','.mjs','.cjs','.java','.jar','.ipynb'},str(p)
used=sum(int(j['actual_core_seconds']) for j in jobs);retry=int(kv(study/'accounting/retry_ledger.tsv')['used'])
assert used<=CAP and retry<=1
table(study/'accounting/budget_final.tsv',[dict(independent_cap_core_seconds=CAP,used_core_seconds=used,used_core_hours=used/3600,
    unspent_core_hours=(CAP-used)/3600,live_reservations=0,technical_retries_used=retry,technical_retry_max=1,failed_allocations=1,
    original_closed_package_used_core_seconds=530,original_used_retries=1,combined_historical_and_new_core_seconds=530+used)])
raw=(study/'accounting/slurm_snapshot.txt').read_text().splitlines();resources=[]
for j in jobs:
    rows=[line.split('|') for line in raw if line.split('|')[0]==j['job_id'] or line.split('|')[0].startswith(j['job_id']+'.')]
    resources.append(dict(task=j['task'],job_id=j['job_id'],qos='xxlarge',allocated_cpus=j['allocated_cpus'],elapsed_seconds=j['elapsed_seconds'],
        conservative_core_seconds=j['actual_core_seconds'],reported_MaxRSS=';'.join(r[5] for r in rows if r[5])))
table(study/'accounting/resource_summary.tsv',resources)
table(study/'B6_gate.tsv',[dict(status='PASS_BOUNDED_B6_SUPPLEMENT',saved_certificates=6,saved_matrices=12,boundary_arguments=3,passage_fixtures=9,
    new_project_cases=0,new_project_determinants=0,new_G_kernels=0,new_trajectories=0,new_culture_datasets=0,fits=0,Step_5=False)])
table(study/'scientific_acceptance.tsv',[
    dict(gate='guarded_bootstrap_and_regression_tests',status='PASS',evidence='P01 plus static_tests.tsv; original import-before-boot rejected'),
    dict(gate='six_saved_certificates_twelve_matrices',status='PASS',evidence='P02 fixed256-bit exact-decoding/Neumann rereview; no determinant/operator regeneration'),
    dict(gate='three_boundary_dependence_arguments',status='PASS_SCOPED',evidence='exact scalar identities and source-reviewed complete dependency composition'),
    dict(gate='constant_O_discrete_passage_induction',status='PASS_UNDER_DECLARED_ASSUMPTIONS',evidence='P02 discrete_induction.md plus accepted P03_R1 nine unchanged-function contracts'),
    dict(gate='SIF_identity',status='PASS',evidence='P03_R1 production image digest; no native kernel/semantic rerun'),
    dict(gate='scientific_runtime_policy',status='PASS_ACCEPTED_PHASES',evidence='P01/P02/P03_R1 accepted; original P03 FAILED policy retained; failed terminal probes denied, not opened'),
    dict(gate='old_evidence_budget_and_tasks_preserved',status='PASS',evidence='all760 sealed files and preservation baseline rehashed'),
    dict(gate='new_B6_delivery',status='PASS_BOUNDED_SUPPLEMENT',evidence='B6_gate.tsv; all owned jobs terminal; independent accounting'),
    dict(gate='Step_5',status='NOT_RELEASED',evidence='separate review/confirmation still required')])
table(study/'origin_status_crosswalk.tsv',[
    dict(item='original_package',historical='CLOSED_WITH_PARTIAL_RESULTS_AND_HOLD',new_status='UNCHANGED',meaning='archive/failed jobs/old used retry remain immutable'),
    dict(item='B6_import_order_and_closeout',historical='failed before mathematical work',new_status='PASS_BOUNDED_SUPPLEMENT',meaning='new independent evidence; not a rewrite of job21971224'),
    dict(item='B02_B03_B05_B06_direct_cases',historical='12 unresolved finite-scope runs',new_status='UNCHANGED_UNRESOLVED',meaning='no new candidates or reruns'),
    dict(item='native_initial_manifold_support',historical='not_verified',new_status='UNCHANGED',meaning='saved fixed-model certificates do not validate general native tool support'),
    dict(item='sigma_kary',historical='unresolved in B03/B05',new_status='UNCHANGED_UNRESOLVED',meaning='not part of B6 repair'),
    dict(item='actual_growth_flow_live_dead_contracts',historical='held experimental/design issues',new_status='UNCHANGED',meaning='pure selectors do not repair experimental metadata'),
    dict(item='Step5',historical='not released',new_status='NOT_RELEASED',meaning='no fitting or design release')])
table(study/'parameter_status_addendum.tsv',[dict(parameter=r['parameter'],scientific_role=r['scientific_role'],
    existing_structural_status='retained; see original case-specific status table',
    boundary_activity='inactive k_o_mis at p_misseg=0' if r['parameter']=='k_o_mis' else 'inactive buffer_n_exp at buffer_beta=0' if r['parameter']=='buffer_n_exp' else 'inactive at zero O only' if r['parameter'] in ('O2_crit','n_O') else 'not reclassified at boundary',
    nuisance_note='outside deterministic mean; variance law not analyzed' if r['parameter']=='sigma_growth' else 'unresolved B03/B05 observer width' if r['parameter']=='sigma_kary' else 'original case roles unchanged',Step_5='not_released') for r in roles])
table(study/'collection_status.tsv',[dict(status='COMPLETE',scientific_scope='bounded B6 supplement only',all_required_B6_gates_pass=True,Step_5=False)])
table(study/'implementation_manifest.tsv',[dict(project_path=str(p.relative_to(project)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(code.iterdir()) if p.is_file()])
now=datetime.now(ZoneInfo('America/New_York')).isoformat(timespec='seconds')
report=f'''# B6 repair and closeout validation supplement

Collected at {now}. Status: **PASS for this bounded B6 supplement; Step5 not released**.

## What passed

The import-order repair installs guarded private-runtime initialization before
the transitive algebra/flint import. A static regression rejects the old entry;
actual-entry stub tests run locally without CAS imports. The xxlarge preflight
passed real private-runtime, dyadic-decoding, nonsingular/singular Neumann,
Python/native read-only write-guard and expected-failure receipt controls.
The intentional missing-dependency child is a passing negative fixture, not a
failed project allocation. Its exit1 and completed inventory are retained.

All6 saved certificates were reused: B01 W1/W2/W3 rank16 and B04 W1/W2/W3
rank18. All12 saved primary/reference matrix enclosures passed an exact-rational
midpoint preconditioner plus outward Neumann bound strictly below1 at256 bits.
The6 saved determinant enclosures exclude zero. No project determinant was
recomputed, no derivative row/operator was regenerated and no structural case
was rerun. Prior source/derivative validity evidence is reused separately from
this new nonsingularity check. See verification/P02/saved_certificate_review.tsv.

Three source-reviewed forward dependence arguments passed, with exact scalar
identities: p_misseg=0 removes k_o_mis while retaining basal missegregation;
buffer_beta=0 removes buffer_n_exp while retaining smax^m and m=0 survival1;
zero O with positive oxygen-shape coordinates removes O2_crit/n_O for that
zero-input experiment only. Priors/penalties and positive boundary ranks are
not claimed. See boundary_proofs.md, boundary_proof_obligations.tsv and the
complete source_dependency_map.tsv under verification/P02.

Nine literal observation/passage fixtures passed inside the digest-verified
production SIF. Only unchanged allowlisted selection/extraction functions and
their pure weighted-mean helper were evaluated; no full bootstrap, simulator,
kernel compilation or trajectory ran. Tests cover positive-day selection,
nearest-target/first-index ties, missing target, eligible-inoculum fallback,
protocol infeasibility, fraction-preserving reseeding, invalid shapes/inoculum
and fallback tie/state bookkeeping. See verification/P03_R1/passage_contract_tests.tsv.

The first P03 application passed those9 fixtures but its policy gate failed:
five shell terminal probes were denied, and R could not access the intended
scratch directory and created/deleted an empty fallback directory under /tmp.
That FAILED/77:0 allocation and its application0/policy77 distinction remain
unchanged. The one technical retry repaired only the harness: an explicitly
writable phase-owned temporary bind, an early tempdir assertion, and a narrow
classification of denied /dev/tty probes as nonpersistent queries. The native
guard still rejects them. Other blocked accesses and any successful out-of-scope
file/directory/terminal mutation still fail, with negative tests for each.
The retry passed the same9 fixtures, so there were9 unique interface fixtures
evaluated twice; no successful structural case or saved matrix was rerun.

## Meaning and limits of the discrete argument

For the same one constant O with feedback disabled, the known family
C,n -> C',n', alpha'=alpha*H/H', mu_hp'=mu_hp*H/H' preserves rates and event
arguments. Source-reviewed induction preserves the real-arithmetic DT=0.05
indexed live/death/buffer updates and clamps. Equal daily arrays yield equal
target selection, ties, inoculum fallback/failure, fractions and known reseeding.
The production wrapper resets retained compartments to zero at each new segment;
it does not carry them as inherited stocks. Parent live overrides preserve the
selected fractions rather than imposing a fresh Gaussian.

This proof requires identical external inputs, valid original branches and
fresh/correct same-O memoization without hash collision or alias from a previous
distinct oxygen in the same bin. These restrictions are explicit. It proves
neither literal floating-point bit equality nor invariance of multi-O histories,
and it is not a positive discrete identifiability certificate. Source-token
matches serve as provenance anchors; the proof combines exact response algebra,
source-reviewed event/update composition and independent pure interface tests.

The old continuous rank findings retain their conditional local finite-support
scope. They do not prove global uniqueness or practical recovery from sparse
data. Sigma_kary remains unresolved; sigma_growth remains outside the
deterministic mean map. All14 mechanisms and6 nuisances and their original
roles are preserved. The12 unresolved direct cases were not modified or rerun.
Observed OLS growth versus endpoint log-ratio, physical flow calibration and
eligible live/dead denominators remain separate experimental/design questions.
Observed counts used in passage selection remain conditioning information.

## Environment, resources and preservation

Python/Arb rereview used the previously validated host runtime; the pure R
contract tests used the production SIF with SHA256 {SIF_SHA}.
Do not describe the host mathematics as running inside the SIF. The3 accepted
scientific phases P01/P02/P03_R1 passed policy gates; original P03 remains FAIL.
File/process traces and completed inventories include the expected-failure child;
container/scheduler launchers are excluded from that scientific coverage.
Inputs were guarded read-only and rehashed before/after each phase. Outputs are
English, with no Java/JavaScript/JSON-format files or successful excluded runtime
access. Denied terminal probes are preserved as such, with no successful terminal
access allowed. Native import/write guards and SIF read-only mounts are explicit
controls. The accepted SIF tempdir is inside its phase-owned scratch root.

New independent use: **{used} core-seconds = {used/3600:.6f} core-hours** of the
4-core-hour ceiling. Technical retries used: **{retry}/1**. All4 allocations
are terminal (3 COMPLETED,1 FAILED); no new-package job remains active. All used xxlarge, one CPU,
16GiB, single-thread math and no automatic requeue. There was no application
concurrency cap or discretionary timeout; reservations followed the site's
12-hour hard limit and remaining aggregate authority with retry/settlement
headroom. Actual use includes conservative parent/step rounding. MaxRSS is a
scheduler report, not a guarantee that short-lived memory peaks were sampled.

The old530-core-second spend, used retry, FAILED job21971224 and historical
HOLD remain unchanged. Its760-file archive and all preservation-baseline files
were rehashed. New authorization does not reopen its remaining budget. The old
collection erratum (30 evaluated candidate identities,102 precision evaluations)
is retained. No old controller/task, production fit, paper or figure was changed.

This addendum closes the missing B6 obligations only. It does not turn the old
package into an all-questions PASS and does not release Step5. No fits, profiles,
new synthetic culture datasets, trajectories or tumor analysis were started.
'''
(study/'REPORT.md').write_text(report)
doc=project/'oxygen/results/analysis_documentation/frozen_inputs/in_vitro/structural_certificates/culture_structural_step04_B6_closeout_execution_20261116.md'
doc.write_text(f'''# B6 repair and closeout execution

Date: {now}. The bounded B6 supplement passed; Step5 remains NOT RELEASED.

The import-order regression and guarded bootstrap passed,6 saved certificates
and12 matrices passed fixed256-bit rereview,3 source-reviewed boundary arguments
passed, and9 pure production selection/reset fixtures passed in the verified SIF.
No project structural case, determinant, kernel or trajectory was rerun. The
discrete result is a same-constant-O real-arithmetic induction with explicit
cache/branch/input assumptions, not positive discrete identifiability or literal
floating-point equality. All14 mechanisms and6 nuisances are retained.

New use: {used} core-seconds ({used/3600:.6f} core-hours); retries {retry}/1;
xxlarge,4 terminal owned jobs (3 completed,1 failed retained), no live commitments, no application concurrency
cap or discretionary timer. Old530-second spend, failed jobs, retry and HOLD
are preserved. Sigma_kary,12 unresolved direct cases and experimental contracts
remain unchanged. No Step5 or fitting work was released.

Detailed report: [REPORT.md](../oxygen/results/structural_identifiability/{study.name}/REPORT.md).
See B6_gate.tsv, scientific_acceptance.tsv, origin_status_crosswalk.tsv,
verification/P02/saved_certificate_review.tsv and verification/P03_R1/passage_contract_tests.tsv.
''')
manifest=[dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(study.rglob('*')) if p.is_file() and p.name not in ('collection_manifest.tsv','CLOSED')]
table(study/'collection_manifest.tsv',manifest)
(study/'CLOSED').write_text(sha(study/'collection_manifest.tsv')+'\nPASS_BOUNDED_B6_SUPPLEMENT\nSTEP5_NOT_RELEASED\n')
print('B6_SUPPLEMENT_COLLECTED',len(manifest),'files;',used,'core-seconds; retry',retry,'; Step5 held')
