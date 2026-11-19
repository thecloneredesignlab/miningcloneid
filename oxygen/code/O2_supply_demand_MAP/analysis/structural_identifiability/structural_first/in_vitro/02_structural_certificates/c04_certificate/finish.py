"""Metadata closeout of the accepted saved-matrix certificate route."""
from pathlib import Path
import shutil
import sys
from protocol import read,table,sha

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve()
if '--verify-only' in sys.argv[3:]:
    rows=read(study/'collection_manifest.tsv')
    for r in rows:
        p=study/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256']
    assert {r['path'] for r in rows}=={str(p.relative_to(study)) for p in study.rglob('*') if p.is_file() and p.name not in ('collection_manifest.tsv','COMPLETE')}
    assert (study/'COMPLETE').read_text().splitlines()[0]==sha(study/'collection_manifest.tsv')
    print('READ_ONLY_ARCHIVE_VERIFICATION_PASSED',len(rows));raise SystemExit(0)
assert not (study/'COMPLETE').exists()
checks=[]
def check(name,condition):
    assert condition,name;checks.append(dict(check=name,passed=True))
old=project/'oxygen/results/structural_identifiability/20261116_culture_step04_saved_review_01'
for r in read(study/'inputs/parent_collection_manifest.tsv'):check('parent:'+r['path'],sha(old/r['path'])==r['sha256'])
for r in read(study/'source_manifest.tsv'):check('frozen:'+r['path'],sha(study/r['path'])==r['sha256'])
for r in read(study/'inputs/protected_worktree_hashes.tsv'):check('protected:'+r['path'],sha(project/r['path'])==r['sha256'])
cert=read(study/'verification/certificate.tsv')[0]
check('exact_certificate_passed',cert['status']=='PASS_EXACT_SAVED_MATRIX_CERTIFICATE' and cert['rank']=='3')
check('no_new_native_or_random_rank',cert['new_native_calls']==cert['new_randomized_rank_calls']=='0')
check('all_nine_rows_exact',len(read(study/'verification/row_checks.tsv'))==9 and all(r['residual']=='0' for r in read(study/'verification/row_checks.tsv')))
check('both_recovery_identities_exact',len(read(study/'verification/recovery_identities.tsv'))==2 and all(r['residual']=='0' for r in read(study/'verification/recovery_identities.tsv')))
check('runtime_policy_passed',read(study/'verification/policy_gate.tsv')[0]['status']=='PASS')
check('no_solver_loaded',not any(r['module'].split('.')[0] in ('strike_goldd','native_strike','symbtools','json','_json','jpype','jnius') for r in read(study/'verification/imported_modules.tsv')))
check('original_native_incomplete_preserved',cert['native_backend_status']=='original_timeout_and_incomplete_return_preserved' and (study/'inputs/original_application_status.tsv').read_text().strip()=='application_exit_code\t124')
roles=read(study/'inputs/parameter_roles.tsv')
check('all_twenty_roles_preserved',len(roles)==20 and sum(r['scientific_role']=='mechanism' for r in roles)==14 and all(r['structural_result']=='not_analyzed' for r in roles))
job=read(study/'accounting/jobs.tsv')[0]
check('scheduler_success',job['state']=='COMPLETED' and job['exit_code']=='0:0' and int(job['reserved_core_seconds'])==0)
check('owned_queue_empty',not (study/'accounting/owned_queue_final.txt').read_text().strip())
used=int(job['actual_core_seconds']);previous=313;total=previous+used
check('cumulative_budget_preserved',total<=7200)
table(study/'accounting/cumulative_budget.tsv',[dict(key=k,value=v) for k,v in dict(cap_total_core_hours=2,previous_core_seconds=previous,new_core_seconds=used,cumulative_core_seconds=total,cumulative_core_hours=total/3600,reserved_core_seconds=0,new_native_calls=0,new_randomized_rank_calls=0,new_technical_retries=0,remaining_technical_retry_allowance=1).items()])
table(study/'reports/control_status.tsv',[
    dict(case='C02',mathematical_control_status='PASS',evidence='accepted earlier native result and exact checks',native_backend_status='completed'),
    dict(case='C03',mathematical_control_status='PASS',evidence='accepted saved native result and exact checks',native_backend_status='native completed; historical collector failure retained'),
    dict(case='C04',mathematical_control_status='PASS_EXACT_CERTIFICATE',evidence='saved native matrix equals independent jet Jacobian; nonzero determinant and explicit identities',native_backend_status='original automated native return still incomplete after timeout'),
    dict(case='C05',mathematical_control_status='PASS',evidence='accepted native result and exact checks',native_backend_status='completed')])
table(study/'reports/final_status.tsv',[dict(key=k,value=v) for k,v in dict(delivery_complete=True,C04_exact_local_certificate=True,C04_native_backend_completed=False,original_culture_parameters_classified=0,mechanisms_retained=14,nuisances_retained=6,new_native_calls=0,new_randomized_rank_calls=0,Stage_B=False,Step_5=False,synthetic_datasets=0,recovery_fits=0,profiles=0,tumor_fits=0,author_contacts=0,new_general_solver=False).items()])
text='''# C04 reprocessing by exact saved-matrix certificate

Date: 2026-11-16, America/New_York.
Outcome: PASS_EXACT_SAVED_MATRIX_CERTIFICATE for the C04 toy ODE.
Original automated native return: still incomplete; original timeout retained.

## Authority and scope

The user requested reprocessing the timed-out case without discretionary time
cutoffs and explicitly accepted an exact certificate as sufficient, while
retaining the incomplete native record. No new C04 native solve or randomized
rank calculation was performed. No other accepted control was repeated.

## Exact evidence

The archived order-2 native matrix has three rows/columns for (x1,a,b). Its input
model remains x_dot=-a*x1^b and y=x1, with a,b,x1 positive and b unknown real.
Every matrix entry agrees exactly with the independently specified output-jet
Jacobian for [x1,-a*x1^b,a^2*b*x1^(2*b-1)]. All nine residuals are zero.

Its determinant equals -a^2*x1^(3*b-1) exactly. With a>0 and x1>0, a^2 is positive
and x1 raised to any real exponent is positive. The determinant is therefore
negative and nonzero on the registered domain. The analytic augmented
observability-identifiability rank condition supplies a local certificate for
both a and b, with initial x1 directly observed through y=x1.

The explicit identities b=x1*x_ddot/x_dot^2 and a=-x_dot/x1^b also reduce exactly
to their declared parameters; x_dot is nonzero on this domain. These are
closed-form control checks using existing SymPy operations, not a new general
structural-analysis engine. No global-identifiability claim is needed here.

## Interpretation and historical record

The original native job exported this matrix but did not complete its automated
rank/return path. Its application exit124, call marker, source/input identities
and closed-package hashes are preserved. The new certificate does not rewrite
that job as a native success, establish why it was slow, or validate the original
randomized-rank arithmetic for unknown real exponents. The producer/core column
ordering and constant-parameter dynamics were audited as data, not reexecuted.

All four toy mathematical control expectations now have accepted evidence using
their recorded routes. This is mixed native/certificate evidence; it is not a
claim that the unchanged native toolbox completed all four benchmarks. It does
not classify the fourteen original culture mechanisms or six nuisance parameters,
validate their Gaussian initial/observer manifold or release Stage B/Step5.

## Timing, resource and policy

There is no 300-second process cutoff, per-derivative timer or native invocation
in this certificate worker. Historical maxLietime=60 is retained only as input
provenance and is not executed. The original two core-hour aggregate budget is
continued, not reset:313 core-seconds had already been charged. The one-CPU
scheduler reservation was derived from the full remaining6887 core-seconds,
rounded down to114 scheduler minutes with47 seconds left for settlement. This
is enforcement of the confirmed aggregate resource budget, not an independently
chosen case cutoff. Site xxlarge maximum12 hours was audited separately.

No forbidden Java/JavaScript/JSON file/process/module activity was observed in
the admitted path; policy, module inventory and trace checks passed. Artifacts
and authored reports are English. No author contact, new dataset, production
replay, fit, profile, tumor analysis, production-model change, old controller
change, commit or push occurred.

'''
text+='New allocated use: '+str(used)+' core-seconds. Cumulative use: '+str(total)+' core-seconds ('+format(total/3600,'.9f')+' core-hours), below2. No new technical retry, owned live task or reservation remains.\n'
(study/'REPORT.md').write_text(text)
for p in study.rglob('*'):
    if p.is_file():check('allowed_format:'+str(p.relative_to(study)),p.suffix.lower() not in ('.json','.json5','.jsonl','.js','.java','.jar','.ipynb'))
check('English_authored_report',not any('\u4e00'<=c<='\u9fff' for c in text))
table(study/'reports/delivery_verification.tsv',checks)
release=study/'releases/r02_metadata';release.mkdir()
for p in Path(__file__).parent.glob('*.py'):shutil.copy2(p,release/p.name)
manifest=[]
for p in sorted(study.rglob('*')):
    if p.is_file() and p.name not in ('collection_manifest.tsv','COMPLETE'):manifest.append(dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)))
table(study/'collection_manifest.tsv',manifest)
(study/'COMPLETE').write_text(sha(study/'collection_manifest.tsv')+'\n'+str(len(manifest))+' manifest-listed files\n')
print('DELIVERY_CLOSED',len(checks),'checks;',len(manifest),'files;',used,'new and',total,'cumulative core-seconds')
