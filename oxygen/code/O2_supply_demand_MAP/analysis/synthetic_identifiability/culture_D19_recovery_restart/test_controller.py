"""Metadata-only safety fixtures. No model, cases, RDS data or real dispatch."""
from pathlib import Path
import sys,tempfile
from common import table,read
from controller import Owner,FIELDS,cost

def run():
 assert cost([dict(seconds=.4,feasible='TRUE') for _ in range(8)])==1
 with tempfile.TemporaryDirectory(prefix='d19_budget_fixture_') as d:
  s=Path(d);(s/'controller').mkdir();table(s/'accounting/site_limits.tsv',[dict(max_seconds=43200)]);table(s/'accounting/jobs.tsv',[],FIELDS);table(s/'accounting/retry_ledger.tsv',[dict(key='used',value=0)]);o=Owner(s);o.ready=lambda t:False
  assert o.future()[0]==720000
  o.jobs=[dict(task='T00',actual_core_seconds=100,reserved_core_seconds=0)];o.ready=lambda t:t in {'T01','T02D','T03F'};table(s/'accounting/cost_admission.tsv',[dict(forecast_primary_seconds=87600)])
  for l in 'ABCD':o.jobs.append(dict(task='D_'+l+'_bench',actual_core_seconds=0,reserved_core_seconds=43200))
  assert o.used()+o.future()[0]<720000
  (s/'controller/SELECTION_LOCKED').write_text('FIXTURE\n');assert o.future()[0]==28800+28800+14400
  o.jobs.append(dict(task='C0',actual_core_seconds=0,reserved_core_seconds=28800));assert o.future()[0]==28800+14400
  o.jobs.append(dict(task='C1',actual_core_seconds=0,reserved_core_seconds=14400));assert o.future()[0]==28800
  print('FIXTURE_PASS complete_scope_cost nonoverlapping_reservations protected_retry postselection_pool_transfer')
  o.jobs=[];o.ready=lambda t:t in {'D_A_replay2b','D_B_replay2b'}
  for l,val,lp in [('A',1e-5,False),('B',2e-5,True)]:
   p=s/'verification'/('D_'+l+'_replay2b');table(p/'best_identity.tsv',[dict(origin='D_'+l,key='key'+l,D=val)]);table(p/'quality.tsv',[dict(origin='D_'+l,key='key'+l,D=val,verified=True,prediction_pass=True,local_pass=lp)])
  o.close_selection('AB');r=read(s/'reports/selection_lock.tsv')[0];assert r['status']=='D19_PREDICTION_RECOVERED_LOCAL_UNRESOLVED' and not (s/'controller/FINISHED').exists()
  print('FIXTURE_PASS exact_lowest_point_not_hidden truth_assessment_after_lock no_premature_finished')
  assert cost([dict(seconds=.3,feasible='TRUE') for _ in range(7)]+[dict(seconds=.4,feasible='TRUE')])==1
  print('FIXTURE_PASS short_final_reserve_block_uses_entire_family_cost_history')
 with tempfile.TemporaryDirectory(prefix='d19_close_fixture_') as d:
  s=Path(d);(s/'controller').mkdir();table(s/'accounting/site_limits.tsv',[dict(max_seconds=43200)]);jobs=[]
  for t in ['D_A_replay2b','D_B_replay2b','T01','C0','C1']:jobs.append(dict(job_id=str(len(jobs)+1),task=t,state='COMPLETED',exit_code='0:0',allocated_cpus=1,elapsed_seconds=10,actual_core_seconds=10,reserved_core_seconds=0))
  table(s/'accounting/jobs.tsv',jobs,FIELDS);table(s/'accounting/retry_ledger.tsv',[dict(key='used',value=0)]);table(s/'reports/selection_lock.tsv',[dict(status='D19_PREDICTION_RECOVERED_LOCAL_UNRESOLVED',best_found_origin='D_A',best_found_key='keyA',best_found_D=1e-5)])
  rows=[]
  for l in 'AB':
   r=dict(origin='D_'+l,key='key'+l,D=1e-5,verified=True,prediction_pass=True,local_pass=False);rows.append(r);table(s/'verification'/('D_'+l+'_replay2b')/'quality.tsv',[r]);table(s/'verification'/('D_'+l+'_replay2b')/'evaluations.tsv',[dict(event='START'),dict(event='END')])
  table(s/'reports/verified_candidates.tsv',rows);table(s/'verification/C0/assessment_gate.tsv',[dict(status='PASS')]);table(s/'verification/C1/audit_gate.tsv',[dict(status='PASS')]);table(s/'verification/T01/attempt_receipt.tsv',[dict(full_forward_attempts=10)])
  Owner(s).finish_report();assert (s/'controller/FINISHED').exists();assert read(s/'accounting/final_budget.tsv')[0]['settled_core_seconds']=='50';assert sum(int(r['full_bundle_attempts']) for r in read(s/'accounting/candidate_attempts.tsv'))==12;assert read(s/'accounting/budget_ledger.tsv')[0]['uncommitted_protected_future_core_seconds']=='0'
  print('FIXTURE_PASS English_report attempt_counts settled_allocations closed_unused_authority')
 print('ALL_METADATA_FIXTURES_PASS; real_jobs=0; physical_forwards=0')
if __name__=='__main__':run()
