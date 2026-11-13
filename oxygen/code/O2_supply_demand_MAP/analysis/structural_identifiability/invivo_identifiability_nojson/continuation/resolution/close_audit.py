"""Close resolution01 with finite scientific outcomes and settled phase bills."""
import argparse,re,os
from pathlib import Path
import workflow as w
import policy as p
import sys
sys.path.append(str(Path(__file__).resolve().parent.parent))

def main():
 parser=argparse.ArgumentParser();parser.add_argument('study',type=Path);a=parser.parse_args();s=a.study.resolve();jobs=w.read_table(s/'tables/jobs.tsv')
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
   if j['mode'] in ('local','global','directions','profile_fresh','resolution_check','resolution_validation','reference_reconstruct') and j['kind']!='engineering':
    audit=w.kv(s/'tables/record_audits'/(j['task']+'_summary.tsv'))
    if audit['pass']!='true' or audit['payload_sha256']!=w.sha(root/'payload_manifest.tsv'):raise ValueError('Missing full native-record audit')
   if j['mode']=='global' and j['kind']!='engineering':
    proof=[r for r in geometry if r['task']==j['task']];populations=list(root.rglob('actual_initial_population.rds'))
    if len(proof)!=len(populations) or not proof:raise ValueError('Missing archived population geometry/RNG audit')
    for row in proof:
     if row['pass'].lower()!='true' or row['population_sha256']!=w.sha(root/row['block']/'actual_initial_population.rds'):raise ValueError('Population audit fingerprint mismatch')
   if j['mode'] in ('fresh','directions','profile_fresh','replay','resolution_check','resolution_validation','reference_reconstruct'):
    proof=w.kv(s/'tables/prediction_audits'/(j['task']+'_summary.tsv'))
    if proof['pass']!='true' or proof['payload_sha256']!=w.sha(root/'payload_manifest.tsv'):raise ValueError('Missing full native predictive-output identity audit')
  else:
   failed+=1
   if j['task'] not in reviewed or int(reviewed[j['task']]['charged_core_seconds'])!=cost:raise ValueError('Unreviewed failed/partial allocation')
 baseline=w.kv(s/'tables/baseline_gate.tsv');strategy=w.read_table(s/'tables/strategy_review.tsv')
 states={r['policy']:r['decision'] for r in strategy}
 stops=w.read_table(s/'tables/resolution_budget_stops.tsv') if (s/'tables/resolution_budget_stops.tsv').exists() else []
 for r in stops:
  if r.get('phase','baseline')=='baseline':
   for policy in states:
    if r['policy'] in (policy,'all') and states[policy]!='accepted':states[policy]='phase_budget_stop'
 if baseline['pass']=='true':
  w.gate(s,'baseline_gate.tsv')
  profile=w.kv(s/'tables/profile_gate.tsv')
  if profile['pass']=='true':w.gate(s,'profile_gate.tsv');endpoint='accepted_density_profiles'
  else:
   terminal=w.kv(s/'tables/resolution_profile_terminal.tsv')
   if terminal['complete']!='true':raise ValueError('Authorized finite profile work remains')
   w.verify_manifest(s,s/terminal['evidence_manifest']);endpoint='density_profiles_unresolved'
 else:
  if set(states)!={'production','data_only'} or not any(v!='accepted' for v in states.values()) or any(v not in ('accepted','resolution_unresolved','phase_budget_stop') for v in states.values()):
   raise ValueError('Authorized finite baseline work remains')
  if any(j['kind']=='profile' for j in jobs):
   pg=w.kv(s/'tables/profile_gate.tsv')
   if pg['free_baseline_invalidated']!='true':raise ValueError('Profile ancestry lacks accepted/invalidated baseline evidence')
  endpoint='baseline_resolution_unresolved'
 ledger=p.assignments(s);phase=p.export_accounting(s,jobs)
 transition=w.kv(s/'policies/resolution01/transition.tsv')
 historical=int(transition['historical_lifetime_core_seconds'])
 if historical!=388524:raise ValueError('Revised historical opening changed')
 old={j['task'] for j in w.read_table(s/'policies/resolution01/legacy_snapshot/jobs.tsv')}
 assigned={r['task'] for r in ledger}
 if len(assigned)!=len(ledger) or {j['task'] for j in jobs} != old|assigned:raise ValueError('Missing or duplicate phase assignment')
 if historical+sum(a['actual'] for a in phase.values())!=301846+actual:raise ValueError('Phase/historical split does not reconcile')
 counts={}
 for a in ledger:
  key=(a['phase'],a['slot_class'],a['slot_key']);counts[key]=counts.get(key,0)+1
  if a['slot_class'] in p.LIMITS and counts[key]>p.LIMITS[a['slot_class']]:raise ValueError('Non-resetting finite allowance exceeded')
 for name,a in phase.items():
  if a['live'] or a['actual']>a['cap']:raise ValueError('Phase is not closed within its cap')
 w.gate(s,'resolution_native_identity_gate.tsv')
 if w.kv(s/'tables/resolution_native_identity_gate.tsv').get('schema')!='2':raise ValueError('Original-parent native identity proof missing')
 w.typed(s/'tables/resolution_terminal.tsv',dict(complete=True,endpoint=endpoint,
    baseline_outcome='accepted' if baseline['pass']=='true' else 'unresolved',
    baseline_actual_core_seconds=phase['baseline']['actual'],baseline_cap_core_seconds=p.CAPS['baseline'],
    profile_actual_core_seconds=phase['profile']['actual'],profile_cap_core_seconds=p.CAPS['profile'],
    revised_historical_core_seconds=historical,model_evaluations=0))
 # Scientific/billing status remain separate, including unsuccessful allocations.
 # The deployed policy overlay and the preserved legacy helper are siblings
 # under different directories on HPC, unlike the local source layout.
 sys.path.append(str(s/'implementation'))
 from final_accounting import finalize
 if finalize(s)!=actual:raise ValueError('Final billing total differs from closure audit')
 w.table(s/'tables/resolution_final_phase_accounting.tsv',[
  dict(phase=name,closed=True,actual_core_seconds=a['actual'],cap_core_seconds=a['cap'],
    unused_cap_core_seconds=a['cap']-a['actual'],released_dependent_escrow_seconds=a['escrow'],
    unspent_closure_reserve_seconds=a['closure']) for name,a in phase.items()],
  ['phase','closed','actual_core_seconds','cap_core_seconds','unused_cap_core_seconds',
   'released_dependent_escrow_seconds','unspent_closure_reserve_seconds'])
 with (s/'reports/final_accounting.md').open('a') as output:
  output.write('\nResolution01 prospective phase accounting: historical closed opening388524core-seconds '
    '(107.92333333core-hours); baseline {:.8f}/220core-hours; profile {:.8f}/360core-hours. '
    'Historical opening plus these two settled phase costs equals lifetime billing exactly. '
    'Unused dependent/closure reservations are released at this finite endpoint and are not actual charges. '
    'The phases cannot borrow from each other or culture.\n'.format(phase['baseline']['actual']/3600,phase['profile']['actual']/3600))
 w.typed(s/'tables/closure_audit.tsv',{'pass':True,'endpoint':endpoint,'owned_parents':len(jobs),'successful_parents':success,'reviewed_unsuccessful_parents':failed,'new_actual_core_seconds':actual,'lifetime_actual_core_seconds':301846+actual,'historical_formal_plus_new_core_seconds':271355+actual,'active_parents':0,'baseline_phase_cap_core_hours':220,'profile_phase_cap_core_hours':360,'model_evaluations':0})
 from resolution_report import write_report
 from recovery_diagnostics import diagnose
 diagnose(s)
 write_report(s,jobs,phase,endpoint,actual,counts)
 files=set()
 for directory in ('manifests','reference','source','prior','releases','metadata_releases','implementation','collections','dispatch','tables','reports','plots','profile_epochs','logs','policies','policy_releases','policy_implementation'):
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
 w.typed(s/'tables/delivery_verified.tsv',{'pass':True,'authority':'HPC full payload and source verification',
   'delivery_manifest_sha256':w.sha(s/'tables/delivery_manifest.tsv'),'files':len(listing),'owned_parents':len(jobs),
   'all_owned_parents_closed':True,'model_evaluations':0})
 print(endpoint,len(jobs),'closed parents;',actual,'new core-seconds;',len(listing),'delivery files verified locally.')
if __name__=='__main__':main()
