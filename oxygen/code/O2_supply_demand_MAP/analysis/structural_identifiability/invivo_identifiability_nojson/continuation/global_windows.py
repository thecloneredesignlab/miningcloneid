"""Bounded global/polish dispatch passes. No aggregate or attempt ceiling."""
import argparse,re
from pathlib import Path
import workflow as w

def successful(j):return j['state']=='COMPLETED' and j['exit_code']=='0:0'
def main():
 p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('action',choices=('start','polish','resume'));p.add_argument('--policy',choices=('production','data_only'));p.add_argument('--batch',type=int);p.add_argument('--reason',default='');a=p.parse_args();s=a.study.resolve()
 rows=w.snapshot(s);by={j['task']:j for j in rows};passing=set()
 reviewed={r['task'] for r in w.read_table(s/'tables/reviewed_failures.tsv')} if (s/'tables/reviewed_failures.tsv').exists() else set()
 if any(not j['job_id'] for j in rows):raise ValueError('Reconcile uncertain receipt before recording another batch')
 if any(w.closed(j) and not successful(j) and j['task'] not in reviewed for j in rows):raise ValueError('Review execution failure before recording another batch')
 if (s/'tables/baseline_policy_gates.tsv').exists():passing={r['policy'] for r in w.read_table(s/'tables/baseline_policy_gates.tsv') if r['pass'].lower()=='true'}
 if a.action=='start':
  if not a.policy or not a.batch or a.batch<1 or not a.reason:raise ValueError('Policy, unbounded batch identity and evidence reason required')
  g=w.kv(s/'tables/baseline_gate.tsv');w.verify_manifest(s,s/g['evidence_manifest'])
  policies=w.read_table(s/'tables/baseline_policy_gates.tsv');policy=next(r for r in policies if r['policy']==a.policy)
  if policy['pass'].lower()=='true':raise ValueError('Policy already accepted; no deficit-driven global batch')
  chains=w.read_table(s/'tables/baseline_chain_summary.tsv')
  if not any(r['policy']==a.policy and r['termination']=='poll_stationary' for r in chains):raise ValueError('No verified local stagnation evidence yet')
  own=min((r for r in chains if r['policy']==a.policy and r['family']=='uniform'),key=lambda r:float(r['verified_best']))
  decisions=s/'tables/search_batches.tsv';old=w.read_table(decisions) if decisions.exists() else []
  if any(r['policy']==a.policy and int(r['batch'])==a.batch for r in old):raise ValueError('Existing batch; reconcile owned tasks, never duplicate')
  old.append(dict(policy=a.policy,batch=str(a.batch),own_root=own['root_id'],before_uniform=own['verified_best'],before_feasible=str(min(float(policy['best_objective']),float(policy['known_feasible_bound']))),reason=a.reason,recorded_utc=w.utc()))
  w.table(decisions,old,['policy','batch','own_root','before_uniform','before_feasible','reason','recorded_utc'])
  for index,mode in enumerate(('hybrid','uniform')):
   task='g_'+a.policy+'_c'+str(a.batch).zfill(4)+'_'+mode+'_w0001'
   options=dict(mode='global',role='baseline_global',policy=a.policy,root_family='uniform',cores=8,wall_seconds=7200,margin_seconds=600,max_calls=200000,init_mode=mode,seed=61000+100*a.batch+(10 if a.policy=='data_only' else 0)+index)
   if mode=='hybrid':options.update(root_id=own['root_id'],selector='selection',selection_name='root_'+own['root_id']+'_selected')
   else:options.update(root_id='full_uniform_'+a.policy+'_c'+str(a.batch).zfill(4),selector='reference')
   w.verify_manifest(s,s/'tables/cache_source_manifest.tsv');w.submit(s,task,**options)
 elif a.action=='polish':
  for j in rows:
   if j['mode']!='global' or j['kind']!='baseline' or not successful(j) or j['policy'] in passing:continue
   task='b_'+j['task']+'_polish_w0001'
   if task in by:continue
   w.verify_manifest(s,s/'tables/cache_source_manifest.tsv')
   w.submit(s,task,policy=j['policy'],root_id=j['root_id'],root_family=j['root_family'],selector='task',source_task=j['task'])
 else:
  latest={}
  for j in rows:
   if j['mode']=='global' and j['kind']=='baseline':latest[re.sub('_w[0-9]+$','',j['task'])]=j
  for component,j in latest.items():
   if not successful(j) or j['policy'] in passing:continue
   z=w.kv(s/'fits'/j['task']/'result.tsv')
   if z['termination'] not in ('call_limit','time_limit'):continue
   gain=z['window_gain']
   # A new uniform branch has no imported objective reference; its first
   # complete native block supplies the comparison for this decision.
   if gain=='NA':
    blocks=w.read_table(s/'fits'/j['task']/'blocks.tsv');improving=float(blocks[0]['objective'])-float(z['objective'])>1e-10
   else:improving=float(gain)>1e-10
   if not improving:
    print('GLOBAL_STALLED',component,'retain complete population and await full local/recovery review.');continue
   if not (s/'fits'/j['task']/'global_checkpoint.rds').exists():raise ValueError('No completed population; record endpoint restart instead of claiming a population resume')
   match=re.fullmatch('(.*)_w([0-9]+)',j['task']);task=match[1]+'_w'+str(int(match[2])+1).zfill(4)
   if task in by:continue
   w.verify_manifest(s,s/'tables/cache_source_manifest.tsv')
   w.submit(s,task,mode='global',role='baseline_global',policy=j['policy'],root_id=j['root_id'],root_family=j['root_family'],selector='global_population',source_task=j['task'],cores=8,wall_seconds=7200,margin_seconds=600,max_calls=200000)
 w.snapshot(s)
if __name__=='__main__':main()
