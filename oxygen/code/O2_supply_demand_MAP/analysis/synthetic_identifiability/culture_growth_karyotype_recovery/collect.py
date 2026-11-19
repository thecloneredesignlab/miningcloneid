"""Close a finite screening stop. No simulation, RNG draws, fitting or dispatch."""
from pathlib import Path
import math
import sys
from common import *
s=Path(sys.argv[1]).resolve();project=Path(sys.argv[2]).resolve();verify(s)
assert (s/'controller/STOP').exists()
jobs=read(s/'accounting/jobs.tsv');assert len(jobs)==5
assert all(j['state'] in ('COMPLETED','FAILED') and int(j['reserved_core_seconds'])==0 for j in jobs)
assert [j['task'] for j in jobs]==['I01','I02','I02_R1','I03','D_A_screen']
assert next(j for j in jobs if j['task']=='I02')['exit_code']=='1:0'
for task in ('I01','I02_R1','I03','D_A_screen'):
 j=next(j for j in jobs if j['task']==task);assert j['state']=='COMPLETED' and j['exit_code']=='0:0'
 p=s/'verification'/task
 for r in read(p/'phase_manifest.tsv'):assert (p/r['path']).stat().st_size==int(r['bytes']) and sha(p/r['path'])==r['sha256']
 for name in ('application_result.tsv','gate.tsv','policy_gate.tsv'):assert read(p/name)[0]['status']=='PASS'
 for r in read(p/'checks.tsv') if (p/'checks.tsv').exists() else []:assert r['passed']=='TRUE'
p=s/'verification/I02'
for r in read(p/'phase_manifest.tsv'):assert sha(p/r['path'])==r['sha256']
assert read(p/'policy_gate.tsv')[0]['status']=='PASS'
assert any(r['check']=='representative_nontruth_scores' and r['passed']=='FALSE' for r in read(p/'checks.tsv'))
external=0;raw=[];totals=[]
for task in ('I02','I02_R1','D_A_screen'):
 records=read(s/'verification'/task/'evaluations.tsv');starts={int(r['call']):r for r in records if r['event']=='START'};ends={int(r['call']):r for r in records if r['event']=='END'}
 assert len(starts)==len(ends) and set(starts)==set(ends)==set(range(1,len(starts)+1))
 for i,r in starts.items():
  e=ends[i];assert r['candidate_sha']==e['candidate_sha'] and all(r[k]==e[k] for k in r if k.startswith('u'))
  assert all(0<=float(r[k])<=1 for k in r if k.startswith('u'))
  if task=='D_A_screen':raw.append(dict(draw=i,candidate_sha=r['candidate_sha'],feasible=e['feasible'],score=e['score'],**{k:v for k,v in r.items() if k.startswith('u')}))
 feasible=sum(r['feasible']=='TRUE' for r in ends.values());totals.append(dict(task=task,attempts=len(starts),feasible=feasible,infeasible=len(starts)-feasible))
 if task!='D_A_screen':external+=len(starts)
 else:assert len(starts)==1000 and feasible==0
assert external==134 and external<=512
assert int(read(s/'verification/I02_R1/engineering_receipt.tsv')[0]['external_attempts'])==external
assert float(read(s/'verification/I02_R1/engineering_receipt.tsv')[0]['positive_D'])<=1e-12
assert kv(s/'accounting/retry_ledger.tsv')['used']=='1'
assert all(r['rounds_used']=='0' for r in read(s/'accounting/reserve_ledger.tsv'))
assert sha(s/'inputs/public/synthetic_dataset.rds')=='cfe0803a09baf31c5c3244a752ca899395185301f57e939a148ba38e30ac9205'
assert all(not(p/'read_denials.tsv').exists() or (p/'read_denials.tsv').stat().st_size==0 for p in [s/'verification/D_A_screen'])
assert len(set((s/'verification/I03/read_denials.tsv').read_text().splitlines()))==4
old=project/'oxygen/results/synthetic_identifiability/20261116_culture_growth_karyotype_conditional_01'
assert sha(old/'collection_manifest.tsv')=='d8d9e0407a8260f0510ec8e09f0c1d94978c2c2ffaf9a82999ce1e04e5ddbea0'
for r in read(old/'collection_manifest.tsv'):assert sha(old/r['path'])==r['sha256']
preserved=[]
for r in read(s/'preservation_baseline.tsv'):
 f=project/r['project_path'];preserved.append(dict(project_path=r['project_path'],baseline_sha=r['sha256'],current_sha=sha(f) if f.exists() else 'MISSING',unchanged=f.exists() and sha(f)==r['sha256']))
table(s/'reports/preservation_audit.tsv',preserved)
assert all(r['unchanged'] for r in preserved),'Unrelated drift needs explicit review before closure'
for p in s.rglob('*'):
 if p.is_file():assert p.suffix.lower() not in ('.json','.jsonl','.js','.java','.jar')
table(s/'reports/attempt_reconciliation.tsv',totals)
table(s/'reports/raw_screen_bank_from_journal.tsv',raw)
actual=sum(int(j['actual_core_seconds']) for j in jobs);assert actual==182 and actual<CAP
for name,status,detail in [('adapter_gate.tsv','PASS_AFTER_ONE_TECHNICAL_REPAIR','Source, production SIF, mapping, scoring, exact-zero activity, optimizer APIs and native read isolation passed'),('D_gate.tsv','SCREENING_UNRESOLVED','D_A found0 feasible candidates in1000 registered attempts; no population or biological optimizer call'),('N_quality_gate.tsv','NOT_RELEASED','D19 recovery gate not passed; no N20 job or observation regeneration')]:table(s/'reports'/name,[dict(status=status,detail=detail)])
table(s/'reports/final_status.tsv',[dict(status='SCREENING_UNRESOLVED',closed=True,actual_core_seconds=actual,actual_core_hours=actual/3600,maximum_core_hours=400,technical_retries_used=1,scientific_reserve_rounds_used=0,D_origins_screened=1,D_optimizer_origins=0,N_origins=0,screen_attempts=1000,screen_feasible=0,external_forward_attempts=external,profiles=0,tumor=0)])
(s/'REPORT.md').write_text('''# D19/N20 conditional recovery: finite screening stop

Final status: **SCREENING_UNRESOLVED; N20 NOT_RELEASED**.

## What passed

The production SIF/source identities, all14 mechanism and6 nuisance roles, D19/N20 mapping, transform endpoints, exact-zero activity, repeated protocol/state replay, independent scoring, native optimizer APIs, population/RNG fixtures and blinded-read isolation passed. D19 positive-control discrepancy was1.66105023660953e-30. The4 private positive-control perturbations validate wiring only; they are not independent blind recovery origins. Eighty-two local static checks passed. The original closed generator and frozen noisy dataset remained unchanged.

## Finite failure and meaning

D_A's registered blind mixture completed1000 charged screening attempts. All1000 were recognized production-protocol infeasible;0 of the required8 feasible anchors were found. No NP200 population was formed, no DE/Subplex biological recovery optimization began, and the400-global-callback benchmark was not admissible. Under the confirmed R2 insufficient-screening stop rule, remaining D origins/reserves and every N20 origin stayed unreleased. There are no recovered parameter estimates or multi-origin agreement results to interpret.

Slurm21992565 completed0:0 and its application/runtime/source receipts passed: the program completed its finite screen correctly. Its scientific screening target failed. This is not evidence that all admissible parameters are infeasible or that the model is non-identifiable. The positive control and its small engineering perturbations were feasible. The accepted broad blind proposal did not reach that feasible region at its registered resolution. The swallowed recognized protocol conditions were logged as infeasible, but their detailed passage/reason fields were not retained; no unsupported biological attribution is made.

## One technical repair, fully retained

I01 job21992476 passed. I02 job21992484 failed1:0 after its128 independent random engineering proposals were all protocol-infeasible, leaving no nontruth vector for the score-comparison fixture. Its source/mapping/positive/zero checks had passed. The sole authorized technical retry21992545 corrected that engineering fixture to4 private small perturbations; none was injected into blind screening. Original129 attempts and35 core-seconds remained charged. I03 job21992558 then passed4 explicit denial canaries, including the retry namespace. No second retry was used or released.

## Resources and artifacts

Five scientific allocations, all xxlarge/one CPU, settled182 allocated core-seconds (0.0505555556 core-hours) of the400-hour independent ceiling. No live allocation remains. There were134 cumulative external forward attempts and1000 origin screening attempts. Scientific reserve rounds used:0. Scheduler reservations respected the site's12-hour cap and remaining resource budget; no arbitrary application time cutoff was used. All artifacts are English; runtime/file checks found no forbidden Java/JavaScript/JSON use. Existing tasks/controllers and all preservation-baseline files remained unchanged.

`reports/attempt_reconciliation.tsv` reconciles every START/END pair. `reports/raw_screen_bank_from_journal.tsv` preserves all19-dimensional raw proposals from their durable exact-coordinate journals. Since screening failed, its RDS population bank was never constructed; the per-draw mixture-component labels were not persisted by that failed path. Original RNG bank/state and exact raw coordinates are retained. `accounting/jobs.tsv` separates scheduler states/exits from application and scientific gates. The original failed phase and repaired release are preserved separately.

## Next decision

Do not spend unused authority automatically. First design a focused audit of the exact infeasibility conditions/passages, then propose a truth-independent feasible-origin construction with its own validation and finite budget. Do not change parameter bounds, conditioning counts, first-eight rule or blind proposal bank retrospectively. No profiles, tumor, new noisy dataset or noisy fit has been started. Structural and experimental-contract gaps remain unchanged; sigma_kary flatness cannot be assessed by a recovery run that never reached optimization.
''')
(s/'CLOSED').write_text('SCREENING_UNRESOLVED; finite insufficient-screening stop; N20 NOT_RELEASED\n')
(s/'reports/validation_final.tsv').write_text('check\tstatus\nsource_inputs_runtime\tPASS\nphase_receipts\tPASS\nSTART_END_and_caps\tPASS\nblind_screening_target\tFAIL_0_OF_8\nD19_recovery\tNOT_RUN\nN20\tNOT_RELEASED\npreservation\tPASS\nold_generator\tPASS\nfrozen_dataset\tPASS\nforbidden_artifacts\tPASS\nno_live_allocations\tPASS\n')
table(s/'collection_manifest.tsv',[dict(path=str(p.relative_to(s)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(s.rglob('*')) if p.is_file() and p.name!='collection_manifest.tsv'])
print('CLOSED_SCREENING_UNRESOLVED',actual,'core_seconds',len(read(s/'collection_manifest.tsv')),'files')
