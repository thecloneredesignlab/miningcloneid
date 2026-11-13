"""One evidence-driven continuation pass, without a total/attempt/sweep cap."""
import argparse,re
from pathlib import Path
import workflow as w

def advance(study):
 s=Path(study).resolve()
 rows=w.snapshot(s);by={j['task']:j for j in rows};latest={}
 for j in rows:
  if j['kind']=='baseline' and j['mode']=='local':latest[re.sub('_w[0-9]+$','',j['task'])]=j
 passing=set()
 if (s/'tables/baseline_policy_gates.tsv').exists():passing={r['policy'] for r in w.read_table(s/'tables/baseline_policy_gates.tsv') if r['pass'].lower()=='true'}
 if passing=={'production','data_only'}:
  print('Both policies accepted; no additional baseline continuation.');return
 for component,j in latest.items():
  root=j['root_id']
  if j['policy'] in passing or not w.closed(j) or j['state']!='COMPLETED' or j['exit_code']!='0:0':continue
  f=by.get(j['task']+'_fresh');d=by.get(j['task']+'_directions')
  if not f or not d or any(x['state']!='COMPLETED' or x['exit_code']!='0:0' for x in (f,d)):continue
  result=w.kv(s/'fits'/j['task']/'result.tsv');diagnostic=w.kv(s/'tests'/d['task']/'result.tsv')
  if float(diagnostic['directional_improvement'])>.002:selector='task';source=d['task'];reason='material independently verified direction improvement; endpoint restart'
  elif result['termination'] in ('call_limit','time_limit'):selector='checkpoint';source=j['task'];reason='resource-window interruption; continue pending poll/checkpoint'
  elif result['termination']=='poll_stationary':
   print('STATIONARY',root,'requires recovery/exploration decision; no blind local restart.');continue
  else:raise ValueError('Unrecognized numerical termination')
  match=re.fullmatch('(.*)_w([0-9]+)',j['task'])
  if not match:raise ValueError('Window identity missing')
  task=match[1]+'_w'+str(int(match[2])+1).zfill(4)
  if task in by:continue
  # The unbounded window number identifies work, not an attempt limit.
  print('DECISION',root,reason)
  w.submit(s,task,policy=j['policy'],root_id=root,root_family=j['root_family'],selector=selector,source_task=source)
 w.snapshot(s)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('study',type=Path)
 advance(p.parse_args().study)
