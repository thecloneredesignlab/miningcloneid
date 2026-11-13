"""Close owned results with integer allocation and immutable payload evidence.

No aggregate spending ceiling. Run only after the accepted scientific endpoint
or the predeclared, evidence-based baseline strategy stagnation endpoint.
"""
import argparse,re,os
from pathlib import Path
import workflow as w

def main():
 p=argparse.ArgumentParser();p.add_argument('study',type=Path);a=p.parse_args();s=a.study.resolve();jobs=w.read_table(s/'tables/jobs.tsv')
 if not jobs or any(not w.closed(j) or not j['job_id'] for j in jobs):raise ValueError('An owned allocation/receipt remains active or uncertain')
 if len({j['job_id'] for j in jobs})!=len(jobs):raise ValueError('Duplicate allocated parent identity')
 w.verify_manifest(s,s/'manifests/input_manifest.tsv');w.verify_manifest(s,s/'tables/cache_source_manifest.tsv');w.audit_formats(s/'implementation')
 prior_costs={};prior_ids=set()
 for phase in ('pilot','formal','repair'):
  d=s/'prior'/phase;old=w.read_table(d/'jobs.tsv');account={r.split('|')[0]:r.split('|') for r in (d/'raw_accounting.txt').read_text().splitlines() if r};cost=0
  for j in old:
   if j['job_id'] in prior_ids:raise ValueError('Historical allocation counted twice')
   prior_ids.add(j['job_id']);r=account[j['job_id']];cost+=int(r[3])*int(r[4])
   if not w.closed(j) or r[1]!=j['state'] or r[2]!=j['exit_code']:raise ValueError('Historical accounting snapshot is not closed/coherent')
  prior_costs[phase]=cost
 if sum(prior_costs.values())!=301846 or prior_costs['formal']+prior_costs['repair']!=271355:raise ValueError('Historical integer opening ledger changed')
 raw={r.split('|')[0]:r.split('|') for r in (s/'tables/raw_accounting.txt').read_text().splitlines() if r}
 reviewed={r['task']:r for r in w.read_table(s/'tables/reviewed_failures.tsv')};actual=0;failed=0;success=0
 geometry=w.read_table(s/'tables/population_geometry_audit.tsv') if (s/'tables/population_geometry_audit.tsv').exists() else []
 for j in jobs:
  record=raw[j['job_id']];cost=int(record[3])*int(record[4]);actual+=cost
  if cost!=int(j['actual_core_seconds']) or int(record[3])!=int(j['cores']) or cost>int(j['reserved_core_seconds']):raise ValueError('Allocation accounting mismatch: '+j['task'])
  if record[1]!=j['state'] or record[2]!=j['exit_code']:raise ValueError('Slurm terminal identity mismatch')
  dispatch=s/'dispatch'/j['task'];w.verify_manifest(dispatch,dispatch/'manifest.tsv')
  if w.sha(dispatch/'manifest.tsv')!=j['dispatch_sha256'] or w.sha(dispatch/'spec.tsv')!=j['spec_sha256']:raise ValueError('Dispatch fingerprint mismatch')
  release=s/'releases'/j['release'];w.verify_manifest(release,release/'manifest.tsv')
  receipt=(s/'logs'/(j['task']+'_receipt.txt')).read_text().splitlines()[0].split(';')[0]
  if receipt!=j['job_id']:raise ValueError('Submission receipt mismatch')
  runtime=s/'runtime'/j['task'];image=runtime/'image_sha256.txt'
  if not image.exists() or image.read_text().strip()!='dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372':raise ValueError('Missing pinned-image launch proof')
  if j['state']=='COMPLETED' and j['exit_code']=='0:0':
   success+=1;root=s/('fits' if j['mode'] in ('local','global') else 'tests')/j['task'];w.verify_payload(root)
   if (runtime/'launch_exit_code.txt').read_text().strip()!='0':raise ValueError('Successful parent has inconsistent launcher exit')
   result=w.kv(root/'result.tsv')
   if result['pass']!='true':raise ValueError('Internal numerical result did not pass')
   if j['mode'] in ('local','global','directions','profile_fresh') and j['kind']!='engineering':
    audit=w.kv(s/'tables/record_audits'/(j['task']+'_summary.tsv'))
    if audit['pass']!='true' or audit['payload_sha256']!=w.sha(root/'payload_manifest.tsv'):raise ValueError('Missing full native-record audit')
   if j['mode']=='global' and j['kind']!='engineering':
    proof=[r for r in geometry if r['task']==j['task']];populations=list(root.rglob('actual_initial_population.rds'))
    if len(proof)!=len(populations) or not proof:raise ValueError('Missing archived population geometry/RNG audit')
    for row in proof:
     if row['pass'].lower()!='true' or row['population_sha256']!=w.sha(root/row['block']/'actual_initial_population.rds'):raise ValueError('Population audit fingerprint mismatch')
   if j['mode'] in ('fresh','directions','profile_fresh','replay'):
    proof=w.kv(s/'tables/prediction_audits'/(j['task']+'_summary.tsv'))
    if proof['pass']!='true' or proof['payload_sha256']!=w.sha(root/'payload_manifest.tsv'):raise ValueError('Missing full native predictive-output identity audit')
  else:
   failed+=1
   if j['task'] not in reviewed or int(reviewed[j['task']]['charged_core_seconds'])!=cost:raise ValueError('Unreviewed failed/partial allocation')
 baseline=w.kv(s/'tables/baseline_gate.tsv');strategy=w.read_table(s/'tables/strategy_review.tsv')
 if baseline['pass']=='true':
  w.gate(s,'baseline_gate.tsv');profile=w.gate(s,'profile_gate.tsv');endpoint='accepted_density_profiles'
 else:
  states={r['policy']:r['decision'] for r in strategy}
  if set(states)!={'production','data_only'} or not any(v=='strategy_stagnation' for v in states.values()) or any(v not in ('accepted','strategy_stagnation') for v in states.values()):raise ValueError('Authorized numerical work remains; no closure at submission, a resource window or a historical budget')
  if any(j['kind']=='profile' for j in jobs):
   pg=w.kv(s/'tables/profile_gate.tsv')
   if pg['free_baseline_invalidated']!='true':raise ValueError('Profiles ran without a valid accepted/invalidated-epoch record')
  endpoint='baseline_strategy_stagnation'
 # Scientific/billing status remain separate, including unsuccessful allocations.
 from final_accounting import finalize
 if finalize(s)!=actual:raise ValueError('Final billing total differs from closure audit')
 w.typed(s/'tables/closure_audit.tsv',{'pass':True,'endpoint':endpoint,'owned_parents':len(jobs),'successful_parents':success,'reviewed_unsuccessful_parents':failed,'new_actual_core_seconds':actual,'lifetime_actual_core_seconds':301846+actual,'historical_formal_plus_new_core_seconds':271355+actual,'active_parents':0,'aggregate_core_hour_ceiling':'none','model_evaluations':0})
 files=set()
 for directory in ('manifests','reference','source','prior','releases','metadata_releases','implementation','collections','dispatch','tables','reports','plots','profile_epochs','logs'):
  d=s/directory
  if d.exists():files.update(f for f in d.rglob('*') if f.is_file() and f.name!='.DS_Store' and '__pycache__' not in f.parts and f.name not in ('delivery_manifest.tsv','delivery_verified.tsv'))
 # Include required hidden frozen inputs, even though Finder sidecars elsewhere
 # are not scientific delivery content.
 files.update(s/r['path'] for r in w.read_table(s/'manifests/input_manifest.tsv'))
 files.update(s/r['path'] for r in w.read_table(s/'tables/cache_source_manifest.tsv'))
 for j in jobs:
  runtime=s/'runtime'/j['task']
  files.update(f for f in runtime.glob('*.txt') if f.is_file())
  root=s/('fits' if j['mode'] in ('local','global') else 'tests')/j['task']
  if (root/'COMPLETE').exists():
   files.update((root/'COMPLETE',root/'payload_manifest.tsv'));files.update(root/r['path'] for r in w.read_table(root/'payload_manifest.tsv'))
  elif root.exists():
   for directory,subdirs,names in os.walk(str(root),followlinks=False):
    files.update(Path(directory)/name for name in names if name!='.DS_Store' and not (Path(directory)/name).is_symlink())
 generation=(s/'state/CURRENT').read_text().strip();files.update((s/'state/CURRENT',s/'state'/generation/'jobs.tsv',s/'state'/generation/'manifest.tsv',s/'STATUS.md'))
 listing=[dict(path=str(f.relative_to(s)),sha256=w.sha(f)) for f in sorted(files)]
 w.table(s/'tables/delivery_manifest.tsv',listing,['path','sha256']);w.verify_manifest(s,s/'tables/delivery_manifest.tsv')
 print(endpoint,len(jobs),'closed parents;',actual,'new core-seconds;',len(listing),'delivery files verified locally.')
if __name__=='__main__':main()
