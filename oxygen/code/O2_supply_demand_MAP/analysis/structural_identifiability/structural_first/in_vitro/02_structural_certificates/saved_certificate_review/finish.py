"""Finite saved-review closeout; standard-library metadata/hash checks only."""
from pathlib import Path
import datetime
import shutil
import sys
from base_contract import read,table,sha,policy_findings

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve()
if '--verify-only' in sys.argv[3:]:
    rows=read(study/'collection_manifest.tsv')
    for r in rows:
        p=study/r['path'];assert p.is_file() and sha(p)==r['sha256'] and p.stat().st_size==int(r['bytes'])
    expected={r['path'] for r in rows}
    actual={str(p.relative_to(study)) for p in study.rglob('*') if p.is_file() and p.name not in ('collection_manifest.tsv','COMPLETE')}
    assert expected==actual
    assert sha(study/'collection_manifest.tsv')==(study/'COMPLETE').read_text().splitlines()[0]
    print('READ_ONLY_ARCHIVE_VERIFICATION_PASSED',len(rows));raise SystemExit(0)
assert not (study/'COMPLETE').exists()
checks=[]
def check(name,value):
    assert value,name
    checks.append(dict(check=name,passed=True))
old=project/'oxygen/results/structural_identifiability/20261116_culture_step04_four_controls_01'
ancestor=project/'oxygen/results/structural_identifiability/20261116_culture_step04_public_tools_01'
for row in read(study/'inputs/parent_collection_manifest.tsv'):check('parent:'+row['path'],sha(old/row['path'])==row['sha256'])
for row in read(study/'inputs/ancestor_manifest.tsv'):check('ancestor:'+row['path'],sha(ancestor/row['path'])==row['sha256'])
for row in read(ancestor/'reports/frozen_source_verification.tsv'):check('production:'+row['path'],sha(project/row['path'])==row['sha256'])
for row in read(study/'source_manifest.tsv'):check('frozen:'+row['path'],sha(study/row['path'])==row['sha256'])
for row in read(study/'inputs/protected_worktree_hashes.tsv'):check('protected:'+row['path'],sha(project/row['path'])==row['sha256'])
roles=read(study/'inputs/parameter_roles.tsv')
check('all_twenty_roles',len(roles)==20 and len({r['parameter'] for r in roles})==20)
check('fourteen_mechanisms',sum(r['scientific_role']=='mechanism' for r in roles)==14)
check('six_nuisances',sum(r['scientific_role']!='mechanism' for r in roles)==6)
check('no_original_parameter_classification',all(r['structural_result']=='not_analyzed' for r in roles))
jobs=read(study/'accounting/jobs.tsv')
check('finite_job_count',len(jobs)<=4)
check('all_terminal',all(r['state'] in ('COMPLETED','FAILED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','CANCELLED') for r in jobs))
check('zero_reservations',all(int(r['reserved_core_seconds'])==0 for r in jobs))
check('empty_owned_queue',not (study/'accounting/owned_queue_final.txt').read_text().strip())
actual=sum(int(r['actual_core_seconds']) for r in jobs);retries=sum(r['retry']=='True' for r in jobs)
check('two_hour_cap',actual<=7200);check('one_retry_cap',retries<=1)
review=study/'review/C03';rr=read(review/'result.tsv')[0]
check('saved_C03_accepted',rr['status']=='PASS' and rr['new_native_C03_calls']=='0')
check('C03_no_new_native_marker',not (study/'controls/C03').exists())
check('C03_no_new_randomized_rank',read(review/'call_ledger.tsv')[0]['randomized_rank_calls']=='0')
check('C03_no_backend_loaded',not any(r['module'].split('.')[0] in ('symbtools','strike_goldd','native_strike','functions') for r in read(review/'imported_modules.tsv')))
check('C03_policy_pass',not policy_findings(review) and read(review/'policy_gate.tsv')[0]['status']=='PASS')
for folder in (review,review/'intermediate_order2'):
    check(str(folder.relative_to(study))+':exact_rows',all(r['residual']=='0' for r in read(folder/'row_checks.tsv')))
summary=[dict(case='C02',status='PASS_REUSED_PARENT',rank=2,derivative_order=1,new_native_calls=0,origin='accepted parent; not rerun'),dict(case='C03',status='PASS_SAVED_REVIEW',rank=2,derivative_order=3,new_native_calls=0,origin='saved original job 21966795; original FAILED collector receipt retained')]
errata=[]
for case in ('C04','C05'):
    folder=study/'controls'/case
    call=read(folder/'call_started.tsv')[0];check(case+':one_native_call_marker',call['native_calls']=='1' and call['case']==case)
    check(case+':core_unchanged',sha(folder/'native_app/strike_goldd.py')==sha(study/'upstream/StrikePy/StrikePy/strike_goldd.py'))
    for name in ('elim_and_recalc.py','rationalize.py'):check(case+':helper_unchanged:'+name,sha(folder/'native_app/functions'/name)==sha(study/'upstream/StrikePy/StrikePy/functions'/name))
    state=next(r for r in jobs if r['task']=='native_'+case)
    if (folder/'result.tsv').exists():
        r=read(folder/'result.tsv')[0]
        check(case+':accepted_receipt',r['status']=='PASS' and r['run_id']==study.name and state['state']=='COMPLETED' and state['exit_code']=='0:0')
        check(case+':policy_pass',not policy_findings(folder) and read(folder/'policy_gate.tsv')[0]['status']=='PASS')
        check(case+':all_rows_exact',all(v['residual']=='0' for v in read(folder/'row_checks.tsv')))
        completion=read(folder/'native_completion.tsv')[0]
        check(case+':durable_native_return',completion['native_returned']=='True' and completion['order_provenance']=='persisted_native_return_frame')
        for item in read(folder/'native_artifact_inventory.tsv'):check(case+':native_artifact:'+item['path'],sha(folder/item['path'])==item['sha256'])
        summary.append(dict(case=case,status='PASS_NEW_NATIVE',rank=r['rank'],derivative_order=r['derivative_order'],new_native_calls=1,origin='immutable call marker and native completion receipt'))
        if 'native_calls_this_process' in r:
            errata.append(dict(case=case,raw_field='native_calls_this_process',raw_value=r['native_calls_this_process'],correct_interpretation='zero native calls during collection subroutine only',authoritative_new_native_calls=1,evidence='call_started.tsv plus native_completion.tsv; raw file preserved'))
    else:
        code=(folder/'application_status.tsv').read_text().strip().split('\t')[1]
        summary.append(dict(case=case,status='UNRESOLVED_TIMEOUT' if code=='124' else 'UNRESOLVED_APPLICATION_OR_VALIDATION',rank='',derivative_order='',new_native_calls=1,origin='app_exit_'+code+'; no complete accepted receipt'))
        check(case+':no_false_completion',not (folder/'native_completion.tsv').exists() or read(folder/'native_completion.tsv')[0]['native_returned']!='True' or not (folder/'result.tsv').exists())
table(study/'reports/control_results.tsv',summary)
table(study/'reports/call_count_field_errata.tsv',errata,['case','raw_field','raw_value','correct_interpretation','authoritative_new_native_calls','evidence'])
table(study/'accounting/call_ledger.tsv',[dict(case=r['case'],new_native_calls=r['new_native_calls'],evidence=r['origin']) for r in summary])
allpass=all(r['status'].startswith('PASS') for r in summary)
status='METHOD_CONTROLS_PASSED_INITIAL_OBSERVER_GATE_HELD' if allpass else 'HOLD_incomplete_nonrational_method_control'
table(study/'accounting/budget_summary.tsv',[dict(key=k,value=v) for k,v in dict(cap_core_hours=2,conservative_actual_core_seconds=actual,conservative_actual_core_hours=actual/3600,reserved_core_seconds=0,technical_retries_used=retries,technical_retries_remaining=1-retries,task_allocations=len(jobs),new_C03_native_calls=0,new_C04_native_calls=1,new_C05_native_calls=1).items()])
table(study/'reports/final_status.tsv',[dict(key=k,value=v) for k,v in dict(delivery_complete=True,scientific_status=status,source_mechanisms=14,source_nuisances=6,new_C03_native_calls=0,new_C04_native_calls=1,new_C05_native_calls=1,project_structural_cases=0,nuisance_interface_cases=0,synthetic_datasets=0,recovery_fits=0,profiles=0,tumor_fits=0,author_contacts=0,Stage_B=False,Step_5=False,general_solver_created=False,closed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()).items()])
lines=['# Saved C03 review and C04/C05 continuation report','', 'Date: 2026-11-16, America/New_York.', 'Status: '+status+'.','',
'## Authority and scope','',
'The user authorized collector repair/tests, saved C03 verification without a native repeat, and C04/C05 only after acceptance. Independent ceiling: two allocated core-hours, xxlarge, no application concurrency cap and at most one technical retry. All authored outputs are English. Java/JavaScript/JSON, author contact, Stage B, Step 5, synthetic recovery, profiles and tumor fitting were excluded.','',
'## Outcomes','', '| Case | Outcome | Augmented rank | New native calls |','|---|---|---:|---:|']
for r in summary:lines.append('| '+str(r['case'])+' | '+r['status']+' | '+str(r['rank'])+' | '+str(r['new_native_calls'])+' |')
lines+=['', 'C02 is reused from the accepted parent. C03 is a separately verified saved native result; the original Slurm FAILED/collector-exit record is unchanged. C03 final order three was reconstructed from the saved log, one-state/two-parameter input, pinned stopping branch and final writer. Both order-2 and order-3 matrices match closed-form references exactly. Identical a/b columns, a nonzero -x1 minor and the admissible a+b invariance establish the expected rank-two control reference without a new randomized rank or native solver call.','',
'The saved-review process rejected solver/backend imports by policy. Its inventory contains no StrikePy/symbtools backend module. C04/C05 interfaces were checked without solving before their admission. The shared collector selects the unique matrix using the authoritative final derivative order, retaining intermediate files.','',
'C04 tests unknown symbolic exponent b in -a*x1^b; its expected reference is rank three. C05 tests -exp(a+b)*x1; its expected reference is rank two and individual a,b ambiguity. These are toy control expectations; a missing/timeout result is unresolved and is never called nonidentifiable on that basis.','']
for r in summary[2:]:
    if not r['status'].startswith('PASS'):
        folder=study/'controls'/r['case']
        lines += [r['case']+' did not produce a complete validated result. Application status: '+(folder/'application_status.tsv').read_text().strip()+'. The admitted timeout/resource stop is enforced; no native repeat, seed change, exponent fixing, extra derivatives or algorithm patch was used. An exported intermediate matrix is retained as partial evidence, not accepted as a completed native classification. The precise internal bottleneck has not been established.','']
lines+=['## Collector and runtime verification','',
'Fifty-two local standard-library tests passed before submission. They cover actual saved C03 order reconstruction, multi-file and lowered-final-order selection, incomplete/corrupt summaries and matrices, symbol/column/domain errors, stale receipts, no scientific side effects on import, forbidden C03 native admission, duplicate markers and durable completion metadata surviving a mocked post-return collector failure.','',
'Native call-start markers are exclusive. Native frame/metadata, artifact hashes and rank/substitution records are persisted before collection assertions. The accepted native results have exact row checks and analytic rank witnesses. The mathematical core/helpers remain byte-identical to the pinned source. No new general structural solver was built.','',
'An auxiliary raw result field native_calls_this_process was mislabeled zero by the collection subroutine in the executed worker. It means zero additional native calls during collection, not zero calls in the whole worker. Original files are preserved. The canonical call ledger uses immutable call-start markers and completion receipts: C03 zero, C04 one, C05 one. The source label is corrected for future use; no scientific result or execution count was changed.','',
'Accepted paths passed import/process/file trace policy checks with no observed forbidden runtime activity. Partial/failed paths are not promoted to compliance or scientific completion merely because a log exists. No installer, MATLAB probe or author communication was used.','',
'## Accounting and protection','',
'Conservative actual charge: '+str(actual)+' core-seconds ('+format(actual/3600,'.9f')+' core-hours), below two. Technical retries used: '+str(retries)+'. No live reservation or owned queued/running task remains. Failed/time-limited work is included and step charges are not double-counted with parents. Scheduler, application, native return, collection and validation states are recorded separately.','',
'Both earlier closed packages, all76 frozen production/input references and unrelated tracked changes retain their hashes. All14 mechanisms and six nuisance roles remain not_analyzed for the original culture model. Four Gaussian initial coordinates retain sharing/normalization; sigma_kary is the unknown assay kernel width; sigma_growth is variance-only and outside deterministic mean-output identification. No production model/SIF or old controller/task was changed; no commit/push was requested.','',
'## Remaining gate','',
'Method controls do not validate the production Gaussian initial/observer manifold, the full133-state problem, boundary/reduction adequacy or the discrete propagation/passage protocol. Stage B and Step 5 remain held. No original culture parameter or combination is classified by these toy controls. A further separately reviewed decision is required for remaining unresolved control work or initial/observer mapping.']
(study/'REPORT.md').write_text('\n'.join(lines)+'\n')
for p in study.rglob('*'):
    if p.is_file():check('allowed_format:'+str(p.relative_to(study)),p.suffix.lower() not in ('.json','.jsonl','.json5','.js','.java','.jar','.ipynb'))
for p in [study/'REPORT.md',*list((study/'reports').glob('*.md'))]:check('English_authored:'+p.name,not any('\u4e00'<=c<='\u9fff' for c in p.read_text()))
table(study/'reports/delivery_verification.tsv',checks)
release=study/'releases/r02_metadata';release.mkdir()
for p in Path(__file__).parent.glob('*.py'):shutil.copy2(p,release/p.name)
manifest=[]
for p in sorted(study.rglob('*')):
    if p.is_file() and p.name not in ('collection_manifest.tsv','COMPLETE'):manifest.append(dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)))
table(study/'collection_manifest.tsv',manifest)
(study/'COMPLETE').write_text(sha(study/'collection_manifest.tsv')+'\n'+str(len(manifest))+' manifest-listed files\n')
print('DELIVERY_CLOSED',status,actual,'core-seconds;',len(manifest),'files;',len(checks),'checks')
