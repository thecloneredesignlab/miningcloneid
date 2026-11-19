"""Assemble an explicitly partial local collection backed by a verified HPC archive."""
from pathlib import Path
from collections import defaultdict
import csv,math,sys
from common import *
s=Path(sys.argv[1]).resolve();project=Path(sys.argv[2]).resolve();verify(s)
for name in ('control_manifest.tsv','closure_manifest.tsv','audit_manifest.tsv'):
 for r in read(s/name):assert sha(s/r['path'])==r['sha256']
assert read(s/'verification/C1/audit_gate.tsv')[0]['status']=='PASS'
assert read(s/'reports/bank_release_gate.tsv')[0]['status']=='BANKS_READY_FOR_D19_REVIEW'
assert kv(s/'accounting/retry_ledger.tsv')['used']=='1'
jobs=read(s/'accounting/jobs.tsv');assert len(jobs)==42 and all(int(j['reserved_core_seconds'])==0 for j in jobs)
assert sum(j['state']=='FAILED' for j in jobs)==1
assert next(j for j in jobs if j['state']=='FAILED')['task']=='P0'
assert all(j['allocated_cpus']=='1' for j in jobs)
actual=sum(int(j['actual_core_seconds']) for j in jobs);assert actual<=CAP
# C1 verified every original scientific phase receipt and audited actual opens.
location=[]
for r in read(s/'reports/authoritative_hpc_inventory.tsv'):
 p=s/r['path'];local=p.is_file()
 if local:assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],'Local copied bytes differ '+r['path']
 assert p.suffix.lower() not in ('.json','.jsonl','.js','.java','.jar')
 location.append(dict(path=r['path'],bytes=r['bytes'],sha256=r['sha256'],location='HPC_AND_LOCAL' if local else 'HPC_ONLY'))
table(s/'reports/payload_location_manifest.tsv',location)
banks=read(s/'verification/C0/bank_manifest.tsv');assert len(banks)==4 and all(int(r['members'])==200 and int(r['construction_attempts'])<=12000 for r in banks)
metrics=[];attempts=0
for b in banks:
 origin=b['origin'];assert (s/'verification'/b['bank_path']/'bank.rds').is_file()
 primary=[j for j in jobs if j['task'].startswith(origin+'_') and not j['task'].endswith(('check1','check2','robust'))]
 starts=ends=0
 for j in primary:
  records=read(s/'verification'/j['task']/'evaluations.tsv');a={int(r['call']):r for r in records if r['event']=='START'};z={int(r['call']):r for r in records if r['event']=='END'}
  assert set(a)==set(z)==set(range(1,len(a)+1));starts+=len(a);ends+=len(z)
  assert all(a[i]['key']==z[i]['key'] and all(a[i]['u'+str(k)]==z[i]['u'+str(k)] for k in range(1,20)) for i in a)
 assert starts==ends==int(b['construction_attempts']);attempts+=starts
 validation=0;robust=0
 for stage in ('check1','check2','robust'):
  p=s/'verification'/(origin+'_'+stage);assert kv(p/'result.tsv')['status']=='PASS'
  assert all(r['passed']=='TRUE' for r in read(p/'checks.tsv'))
  validation+=sum(r['event']=='START' for r in read(p/'evaluations.tsv'))
  if stage=='robust':
   probes=read(p/'bank_robustness.tsv');assert len(probes)<=288;robust=sum(r['feasible']=='TRUE' for r in probes)
  else:assert len(read(p/'bank_replay_checks.tsv'))==200
 assert validation<=688
 metrics.append(dict(origin=origin,construction_attempts=starts,bank_members=200,fresh_member_replays=400,validation_attempts=validation,robust_feasible=robust,robust_probes=len(probes),bank_path=b['bank_path'],status='PASS'))
assert attempts==25715
values=read(s/'verification/C0/bank_parameters_natural_and_unit.tsv');assert len(values)==16000
coords=defaultdict(dict);names=set(MECHANISMS+NUISANCES)
for r in values:coords[(r['origin'],int(r['member']))][r['parameter']]=r
assert len(coords)==800 and all(set(v)==names for v in coords.values())
physical=MECHANISMS+NUISANCES[:4]
for origin in ('F_A','F_B','F_C','F_D'):
 u=[[float(coords[(origin,i)][n]['unit_coordinate']) for n in physical] for i in range(1,201)]
 for i in range(200):
  for j in range(i+1,200):
   distance=max(abs(x-y) for x,y in zip(u[i],u[j]));assert distance>=1e-6
   if i<8 and j<8:assert distance>=.02
 assert all(coords[(origin,i)]['sigma_growth']['status']=='neutral_placeholder_not_estimated' for i in range(1,201))
supply=read(s/'verification/C0/bank_supply_margins.tsv');assert len(supply)==96000
minimum=min(float(r['supply_ratio']) for r in supply)
source_gap=0.;daily_count=0
for i in range(1,5):
 for r in read(s/'verification'/('P2_'+str(i))/'daily_supply.tsv'):
  a=float(r['reported_live']);b=float(r['state_total']);source_gap=max(source_gap,abs(a-b)/max(1,abs(a),abs(b)));daily_count+=1
assert source_gap<=1e-10
old=project/'oxygen/results/synthetic_identifiability/20261116_culture_growth_karyotype_recovery_01'
assert sha(old/'collection_manifest.tsv')=='16d38d3b5536d627926843f679620fb601d4ff508b1ca7f413601a012d58f170'
preservation=[]
for r in read(s/'preservation_baseline.tsv'):
 p=project/r['path'];preservation.append(dict(path=r['path'],baseline_sha=r['sha256'],current_sha=sha(p) if p.exists() else 'MISSING',unchanged=p.exists() and sha(p)==r['sha256']))
table(s/'reports/preservation_audit.tsv',preservation);assert all(r['unchanged'] for r in preservation),'Unrelated drift requires review'
table(s/'reports/constructor_summary.tsv',metrics)
table(s/'reports/count_state_consistency.tsv',[dict(records=daily_count,max_relative_gap=source_gap,status='PASS',scope='Recorded historical parent trajectories; not external experimental unit validation')])
table(s/'reports/bank_margin_summary.tsv',[dict(nominal_boundary_records=len(supply),minimum_supply_ratio=minimum,below_one=sum(float(r['supply_ratio'])<1 for r in supply),scope='Selected original parent boundaries, not experimental uncertainty')])
table(s/'reports/final_status.tsv',[dict(status='BANKS_READY_FOR_D19_REVIEW',validated_banks=4,members=800,fresh_replays=1600,robust_probes=sum(r['robust_probes'] for r in metrics),robust_feasible=sum(r['robust_feasible'] for r in metrics),construction_attempts=attempts,diagnostic_attempts=1064,test_forward_attempts=10,test_toy_attempts=3,validation_attempts=sum(r['validation_attempts'] for r in metrics),full_candidate_attempts=attempts+1064+10+sum(r['validation_attempts'] for r in metrics),maximum_full_candidate_attempts=51944,allocated_core_seconds=actual,allocated_core_hours=actual/3600,maximum_core_hours=32,technical_retries=1,scientific_reserves=0,live_allocations=0,D19_fit=0,N20=0,profiles=0,tumor=0,local_collection='verified_subset_HPC_full_scientific_payload_authority')])
table(s/'accounting/budget_ledger.tsv',[dict(cap_core_seconds=CAP,settled_core_seconds=actual,live_reserved_core_seconds=0,technical_retry_hold=0,collection_hold=0,validation_hold=0,unused_closed_core_seconds=CAP-actual,status='CLOSED_NO_AUTOMATIC_COMPUTE_RELEASE')])
(s/'reports/BANKS_README.md').write_text('''# Frozen initializer banks, not fitted estimates

Four bank RDS objects are in verification/F_A_de4/bank.rds through F_D_de4/bank.rds. Each stores200 selected complete19-dimensional normalized vectors, full20-role natural values, complete construction predictions and origin/RNG/contract provenance. The first8 are the registered coordinate-diverse anchors. Unit distances are in18 declared physical coordinates; this is not proof of functionally distinct mechanisms or different optimum basins. Members within a family are correlated search outputs, not posterior/prior samples.

The natural-and-unit table is verification/C0/bank_parameters_natural_and_unit.tsv. The original bank_members.tsv columns named after parameters are normalized unit coordinates, not natural values. Use the explicit C0 table or native RDS for interpretation.

sigma_growth is a public neutral placeholder used only in unscored physical calls. It is NOT the generating/control SD and must NOT be reused as the D19 reference scale. A separately authorized D19 recovery must inject its original positive known control SD after bank locking. sigma_kary remains free over its complete original bounds and is not constrained by passage feasibility. Its old local binary64-flat observer evidence remains unchanged.

These are only D19 restart-review banks. No D19 data fitting or N20 was performed. Do not silently repurpose these families for N20, import D19 fits/known scales into N20, change full bounds, add truth-centered starts or turn ranges into confidence intervals. All bank members passed two independent unmodified production replays; all1152 registered small feasibility probes passed. This is finite local evidence, not a whole-neighborhood/global certificate.

The initial runtime guard protected the frozen private exports and foreign families. Closure completed alias/historical/fit-object barrier tests using a stricter template. That template was not retroactively deployed in old workers. All97860 successful payload opens in the32 actual construction/validation worker phases were independently audited against the stronger read exclusions: no forbidden read occurred. Future D19 dispatch needs its own reviewed role-specific guard. The real original fit-object path was absent at closure; its absence is not called a successful access-denial test. Existing canonical aliases/historical data and a non-scientific forbidden-filename fixture were actually denied.

Full original traces, compiler caches, member states and continuation checkpoints remain authoritative on HPC. Local storage contains a checksum-verified deliverable subset. Consult payload_location_manifest.tsv and authoritative_hpc_inventory.tsv; do not treat missing local raw artifacts as missing HPC results. No model was rerun for collection/storage repair.
''')
(s/'REPORT.md').write_text('''# D19 feasibility diagnosis and independent start banks

Final status: **BANKS_READY_FOR_D19_REVIEW**. No D19 response fitting, N20, profiles, tumor analysis or new observations were started.

The1000 old global failures reproduced.2N failed974/1000;4N failed1000/1000. The old sequential Boolean concealed26 candidates that could complete2N but still failed4N. Structured source evidence identifies insufficient eligible positive-day parent supply for the requested child inoculum. Later boundaries are NOT_REACHED, not additional demonstrated failures. This does not prove a specific mechanism caused the shortfall or establish non-identifiability.

The unchanged production contract has62/60 segments and61/59 required parent boundaries for2N/4N. Eight missing recorded initial counts use the original cfg fallback;4 missing targets use the next-array recorded initial count and4 use the last positive day. Raw/effective values and sources are frozen. No count, horizon, threshold, inheritance, retained-stock reset or response function was changed. Reported live totals and state sums agreed over172893 parent/day records, with maximum relative discrepancy1.0043749839979664e-14; this checks internal bookkeeping, not physical SOP or external units.

The observation-free, balanced progress/shortfall merit found four new fully feasible NP200 populations without reading truth, expected/noisy responses or old endpoints. Construction used6414/6975/6413/5913 attempts for F_A/F_B/F_C/F_D, totaling25715. Each bank retains all14 mechanisms,4 unknown root coordinates and a free sigma_kary coordinate. The sixth nuisance, sigma_growth, is retained as a neutral nonphysical placeholder in the20-role object and is not estimated or disclosed from truth.

All800 members passed two fresh-process unmodified production replays (1600 total), including exact protocol decisions and registered state/probability tolerances. All1152 registered local physical-coordinate probes stayed feasible. The minimum nominal selected-boundary supply ratio was1.00043724473532 across96000 member/boundary records. These probes/ranges establish neither global optimization nor parameter identification. No parameter truth-error or D19 recovery objective was evaluated.

One technical retry was used: the initial exporter incorrectly rejected recorded missing counts despite the original cfg fallback. Initial P0's12 seconds and repaired P0_R1's23 seconds remain charged; the failure preceded propagation. Original and repaired releases are both retained. No scientific reserve or replacement family was used. The only native DE interruptions were explicit bank-target endings, not claims of convergence or data-fit optimization.

The archive/read audit verified3703 authoritative HPC scientific payload files, all scientific phase receipts and the old closed collection. Actual worker/child read traces contained97860 successful payload opens and0 prohibited reads. Initial direct canaries preceded construction; closure additionally completed canonical-alias/historical/fit-object barrier completeness tests with a stronger template. That later template was not retroactively claimed as the runtime guard. Experimental growth/count SOP and flow/live-dead gaps, sigma_kary observer flatness and unresolved structural cases remain unchanged.

Local full transfer hit ENOSPC. Only this study's HPC-backed local raw/duplicate copies were removed; original HPC outputs and all unrelated work remained protected. A filtered local collection retains four bank RDS objects, natural/unit parameter and supply tables, diagnostics, verification tables and checksum inventories. The complete raw scientific archive remains on HPC. C0 export and C1 archive/read audit were metadata-only, with zero model/optimizer calls; no numerical work was repeated for storage repair.

Use reports/final_status.tsv, constructor_summary.tsv, BANKS_README.md and payload_location_manifest.tsv for exact status, resource accounting and handoff. The independent32-core-hour package closes here. Any D19 recovery restart needs a separately reviewed explicit budget, known-SD injection, guard/counter adaptation and unchanged prediction gates. N20 remains held even though the initializer gate passed.
''')
with (s/'REPORT.md').open('a') as f:f.write('\nFinal allocated use: '+str(actual)+' core-seconds ('+format(actual/3600,'.9f')+' core-hours) of32; live allocations0; unused authority is closed.\n')
table(s/'reports/final_validation.tsv',[dict(check=n,status='PASS') for n in ('source_and_code_layers','remote_scientific_payload_and_phase_hashes','old_collection_unchanged','local_deliverable_subset_hashes','four_independent_rng_families','four_banks_200_members','two_replays_per_member','1152_feasibility_probes','actual_worker_no_forbidden_reads','late_alias_completeness_canaries','attempt_and_resource_caps','local_preservation','no_forbidden_artifact_extensions','no_D19_fit_N20_profiles_or_tumor')])
(s/'CLOSED').write_text('BANKS_READY_FOR_D19_REVIEW; no automatic D19 recovery or N20 release; HPC full scientific payload authority, local verified subset\n')
print('LOCAL_SUBSET_AND_BANKS_VALIDATED',actual,'core_seconds;',sum(r['location']=='HPC_ONLY' for r in location),'HPC-only payload files')
