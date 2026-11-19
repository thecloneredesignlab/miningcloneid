"""English evidence closeout and independent metadata checks; no scientific imports."""
from pathlib import Path
import datetime
import shutil
import subprocess
import sys
from contract import CASES,read,table,sha,receipt_valid,policy_findings

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve()
if '--verify-only' in sys.argv[3:]:
    for row in read(study/'collection_manifest.tsv'):
        p=study/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==int(row['bytes'])
    assert (study/'COMPLETE').read_text().splitlines()[0]==sha(study/'collection_manifest.tsv')
    print('READ_ONLY_ARCHIVE_VERIFICATION_PASSED');raise SystemExit(0)
assert not (study/'COMPLETE').exists()
old=project/'oxygen/results/structural_identifiability/20261116_culture_step04_public_tools_01'
checks=[]
def check(name,value):
    assert value,name
    checks.append(dict(check=name,passed=True))
for row in read(study/'inputs/parent_collection_manifest.tsv'):
    check('closed_parent:'+row['path'],sha(old/row['path'])==row['sha256'])
check('closed_parent_COMPLETE_unchanged',(old/'COMPLETE').read_text().splitlines()[0]=='6ea489717a56b25bb1932635b91f8ff8de91a71036bfe0184ea32203088d00ef')
for row in read(old/'reports/frozen_source_verification.tsv'):
    check('production_reference:'+row['path'],sha(project/row['path'])==row['sha256'])
for row in read(study/'inputs/protected_worktree_hashes.tsv'):
    check('protected_worktree:'+row['path'],sha(project/row['path'])==row['sha256'])
for filename in ('source_manifest.tsv','releases/r02/manifest.tsv'):
    for row in read(study/filename):check('frozen_source:'+row['path'],sha(study/row['path'])==row['sha256'])
roles=read(study/'inputs/parameter_roles.tsv')
check('all_twenty_roles_retained',len(roles)==20 and len({r['parameter'] for r in roles})==20)
check('fourteen_mechanisms_retained',sum(r['scientific_role']=='mechanism' for r in roles)==14)
check('six_nuisances_retained',sum(r['scientific_role']!='mechanism' for r in roles)==6)
check('original_parameters_not_analyzed',all(r['structural_result']=='not_analyzed' for r in roles))
jobs=read(study/'accounting/jobs.tsv')
check('no_live_reservation',all(int(r['reserved_core_seconds'])==0 for r in jobs))
check('finite_task_count',len(jobs)<=6)
check('one_technical_retry',sum(r['retry']=='True' for r in jobs)<=1)
check('terminal_jobs',all(r['state'] in ('COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL') for r in jobs))
actual=sum(int(r['actual_core_seconds']) for r in jobs)
check('three_core_hour_cap',actual<=10800)
summaries=[]
for case in CASES:
    folder=study/'controls'/case
    if (folder/'native_app/strike_goldd.py').exists():
        check(case+':native_core_preserved',sha(folder/'native_app/strike_goldd.py')==sha(study/'upstream/StrikePy/StrikePy/strike_goldd.py'))
        for name in ('elim_and_recalc.py','rationalize.py'):
            check(case+':native_helper_preserved:'+name,sha(folder/'native_app/functions'/name)==sha(study/'upstream/StrikePy/StrikePy/functions'/name))
    if (folder/'result.tsv').exists():
        row=read(folder/'result.tsv')[0]
        good=receipt_valid(row,case,sha(study/'fixtures'/case/'model.py'),study.name)
        good=good and not policy_findings(folder) and read(folder/'policy_gate.tsv')[0]['status']=='PASS'
        j=next(r for r in jobs if r['task']==case)
        good=good and j['state']=='COMPLETED' and j['exit_code']=='0:0'
        check(case+':complete_validated_receipt',good)
        names=read(folder/'native_named_results.tsv')[0]
        check(case+':all_symbolic_rows_exact',all(r['residual']=='0' for r in read(folder/'row_checks.tsv')))
        check(case+':native_core_unchanged',sha(folder/'native_app/strike_goldd.py')==sha(study/'upstream/StrikePy/StrikePy/strike_goldd.py'))
        for name in ('elim_and_recalc.py','rationalize.py'):
            check(case+':native_helper_unchanged:'+name,sha(folder/'native_app/functions'/name)==sha(study/'upstream/StrikePy/StrikePy/functions'/name))
        summaries.append(dict(case=case,status='PASS',identifiable=names['identifiable'],nonidentifiable=names['nonidentifiable'],rank=names['rank'],derivative_order=names['derivative_order'],exact_witness=names['exact_witness'],scope='toy_control_only'))
    elif case=='C03' and (study/'reports/C03_saved_native_evidence.tsv').exists():
        check('C03_saved_trace_policy',not policy_findings(folder))
        summaries.append(dict(case=case,status='NATIVE_COMPLETED_VALIDATION_INCOMPLETE',identifiable='',nonidentifiable='a,b',rank=2,derivative_order=3,exact_witness='',scope='native_expected_disposition_only; exact_row_validation_incomplete'))
    else:summaries.append(dict(case=case,status='NOT_RUN_GATE_HELD',identifiable='',nonidentifiable='',rank='',derivative_order='',exact_witness='',scope='no_structural_conclusion'))
passed=all(r['status']=='PASS' for r in summaries)
table(study/'reports/control_results.tsv',summaries)
for case in CASES:
    pre=study/'preflight_retry'/case
    check(case+':preflight_passed',read(pre/'result.tsv')[0]['status']=='PASS' and not policy_findings(pre))
table(study/'accounting/budget_summary.tsv',[dict(key=k,value=v) for k,v in dict(cap_core_hours=3,conservative_actual_core_seconds=actual,conservative_actual_core_hours=actual/3600,reserved_core_seconds=0,technical_retries_used=sum(r['retry']=='True' for r in jobs),technical_retries_remaining=1-sum(r['retry']=='True' for r in jobs),task_allocations=len(jobs),native_method_calls=sum(r['task'] in CASES for r in jobs)).items()])
scientific='METHOD_CONTROLS_PASSED_INITIAL_OBSERVER_GATE_HELD' if passed else 'HOLD_incomplete_method_controls'
table(study/'reports/final_status.tsv',[dict(key=k,value=v) for k,v in dict(delivery_complete=True,scientific_status=scientific,method_controls_passed=sum(r['status']=='PASS' for r in summaries),mechanisms_retained=14,nuisances_retained=6,project_structural_cases=0,nuisance_interface_controls=0,production_replays=0,synthetic_datasets=0,recovery_fits=0,profiles=0,tumor_fits=0,author_contacts=0,new_general_solver=False,Stage_B_released=False,Step_5_released=False,closed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()).items()])
lines=['# Four-control finite rerun report','', 'Date: 2026-11-16, America/New_York.', 'Status: '+scientific+'.','',
'The user authorized only C02-C05, preceding shape/interface checks, xxlarge with no application concurrency cap, an independent three allocated core-hour ceiling and at most one technical retry. All authored outputs are English. No author was contacted. Stage B, Step 5, synthetic recovery, profiles and tumor fitting were excluded.','',
'## Results','', '| Control | Status | Individually identifiable | Individually nonidentifiable | Augmented rank | Derivatives |','|---|---|---|---|---|---|']
for r in summaries:lines.append('| '+ ' | '.join(str(r[k]) or 'none' for k in ('case','status','identifiable','nonidentifiable','rank','derivative_order'))+' |')
lines += ['', 'C03/C05 have a predeclared analytic a+b invariance. This is the control reference, not an executed C05 result. C04/C05 were not admitted. The saved C03 native answer matches its expected classification, but its complete validation gate did not pass. All statements are scoped to the toy controls; none classifies an original culture mechanism or nuisance.','',
'## Input, interface and evidence validation','',
'Fifty-four initial standard-library static checks passed before submission. The first preflight stopped on an indentation error in the adapter-extracted packing span, before any upstream method call. The one technical retry corrected only extraction padding/dedentation. The replacement span has the same AST as the pinned native branch. Thirteen additional extraction/generated-source checks and shell syntax validation preceded the retry. The failed attempt, original inputs and charges were retained; no retry remains.','',
'The retried preflight passed all four shapes, positive symbol domains, parameter roles, options resolution, zero parameter dynamics, initial output Jacobians, helper signature binding and column-selection primitives. The previously invalid nested dynamics were rejected. A two-state C01 packing regression passed without rerunning its solver.','',
'C02 and C03 were admitted sequentially. C02 completed all checks: its raw matrix was read through a bounded AST vocabulary and matched the independent closed-form Jacobian exactly, with the nonzero -x determinant witness. Native core/helpers remained byte-identical. Profiling recorded native named results and randomized-rank seeds/substitutions without changing the algorithm.','',
'C03 completed its native call, reported a,b individually nonidentifiable, and recorded four generic-rank results of two. The collector then incorrectly required exactly one matrix file, although native derivative extension legitimately produced both order-2 and order-3 files. It stopped before independent exact-row validation or a complete receipt. This is our collector error, not a scientific disagreement or native solver failure. The new selector uses the authoritative final native order and preserves intermediate files. Five static filename-selection checks passed; the corrected candidate is prepared but not executed. C03 saved output and failed validation status are retained.','',
'The one technical retry was already spent on preflight extraction. No second retry, verification allocation or C03 native repeat was submitted. C04/C05 remained blocked by the incomplete C03 receipt. Any next continuation requires separate review and should validate the saved C03 output before considering further native cases. Randomized native ranks are not promoted to global project proofs.','',
'No Java/JVM, JavaScript or JSON runtime/file findings were observed in the admitted private Python paths. The native guard, import hook, process/file traces and loaded-module inventories were reviewed; no forbidden dependency was relaxed and no MATLAB probe or installer was run. The guard alone is not claimed as a universal sandbox.','',
'The native derivative-timeout fallback has a source-level helper-arity inconsistency. It was documented, not patched, and not used to infer a result. Timeouts would have remained unresolved.','',
'## Resources and protected work','',
'Conservative allocated use: '+str(actual)+' core-seconds ('+format(actual/3600,'.9f')+' core-hours), below 3. All failed attempts and the retry are included; parent and step charges are not summed twice. No live reservation remains. Slurm states, application exits, validation receipts and collection status are separate in the evidence package. See accounting/jobs.tsv and slurm_snapshot.txt.','',
'The old closed package, all 76 referenced frozen production/input files and unrelated tracked changes retained their hashes. Old controllers/tasks, production code/SIF, parameter bounds and scientific roles were not modified. No commit or push was requested.','',
'All fourteen mechanism and six nuisance rows remain not_analyzed for the original culture problem. Four Gaussian initial coordinates retain cohort sharing and normalization; sigma_kary remains an unknown observation-kernel width, and sigma_growth remains a variance-only nuisance outside deterministic mean-output identification.','',
'## Remaining release boundary','',
'A complete four-control pass validates only this bounded method/control interface. The production Gaussian initial/observer mapping, 133-state scalability, reduction/tail qualifications, discrete propagation and passage protocol remain unvalidated or held as previously recorded. No N01-N06 or full project structural case was executed. No general-purpose structural solver was created. A separately reviewed initial/observer mapping proposal is the next decision; Stage B and Step 5 are not released automatically.']
(study/'REPORT.md').write_text('\n'.join(lines)+'\n')
for p in study.rglob('*'):
    if p.is_file():check('allowed_format:'+str(p.relative_to(study)),p.suffix.lower() not in ('.json','.jsonl','.json5','.js','.java','.jar','.ipynb'))
for p in [study/'REPORT.md',*list((study/'reports').glob('*.md'))]:check('English_authored:'+p.name,not any('\u4e00'<=c<='\u9fff' for c in p.read_text()))
table(study/'reports/delivery_verification.tsv',checks)
release=study/'releases/r03_metadata';release.mkdir(exist_ok=True)
for p in Path(__file__).parent.glob('*.py'):shutil.copy2(p,release/p.name)
manifest=[]
for p in sorted(study.rglob('*')):
    if p.is_file() and p.name not in ('collection_manifest.tsv','COMPLETE'):manifest.append(dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)))
table(study/'collection_manifest.tsv',manifest)
(study/'COMPLETE').write_text(sha(study/'collection_manifest.tsv')+'\n'+str(len(manifest))+' manifest-listed files\n')
print('DELIVERY_CLOSED',len(checks),'checks;',len(manifest),'files;',sum(r['status']=='PASS' for r in summaries),'of four controls passed;',actual,'core-seconds')
