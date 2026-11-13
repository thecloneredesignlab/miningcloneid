"""One conditional profile dispatch pass, with no window/attempt/core-hour cap.

Run profile_gate.R freeze/review first. This module only dispatches metadata-ready
work; actual R simulations remain inside the pinned SIF allocations.
"""
import argparse,re
from pathlib import Path
import workflow as w

POLICIES=('production','data_only')
BRANCHES=('global_local','ascending','descending')
def successful(j):return j and j['state']=='COMPLETED' and j['exit_code']=='0:0'
def next_id(task):
    match=re.fullmatch('(.*)_w([0-9]+)',task)
    if not match:raise ValueError('Missing unbounded window identity')
    return match[1]+'_w'+str(int(match[2])+1).zfill(4)
def next_local_decision(result,diagnostic,anchor_pass):
    if anchor_pass:return 'accepted'
    if float(diagnostic['directional_improvement'])>.002:return 'endpoint'
    if result['termination'] in ('call_limit','time_limit'):return 'checkpoint'
    if result['termination']=='poll_stationary':return 'stationary'
    raise ValueError('Unrecognized local termination')
def verify_branch_source(study,epoch,task,options,by):
    policy=options['policy'];index=int(options['density_index']);role=options['role'];selector=options['selector']
    if role=='density_check':
        parent=by.get(options['source_task'])
        if not parent or parent['kind']!='profile' or parent['policy']!=policy or parent['mode']!='local' or parent['root_id']!=options['root_id'] or int(parent['density_index'])!=index:raise ValueError('Check must test its own local endpoint')
        return
    branch='global_local' if role in ('density_global','density_global_local') else role[len('density_'):]
    if branch not in BRANCHES:raise ValueError('Unknown profile branch')
    expected='profile_'+epoch+'_'+policy+'_'+(branch+'_r'+str(index).zfill(2) if branch=='global_local' else branch)
    if options['root_id']!=expected:raise ValueError('Profile root/branch identity mismatch')
    if selector=='selection':
        if role=='density_global_local' or (branch=='ascending' and index!=1) or (branch=='descending' and index!=7):raise ValueError('Only branch initialization may use the frozen free reference')
        chosen=study/'tables'/(options['selection_name']+'.rds') if 'selection_name' in options else study/'dispatch'/task/'selection.rds'
        frozen=study/'profile_epochs'/epoch/('baseline_'+policy+'.rds')
        if w.sha(chosen)!=w.sha(frozen):raise ValueError('A mutable/other-branch winner cannot initialize a profile branch')
        return
    parent=by.get(options.get('source_task'))
    if not parent or parent['kind']!='profile' or parent['policy']!=policy or parent['root_id']!=options['root_id']:raise ValueError('Cross-branch profile start rejected')
    source=parent
    if parent['mode']=='profile_fresh':source=by[parent['source_task']]
    before=int(source['density_index'])
    if selector=='checkpoint':
        if source['mode']!='local' or source['role']!=role or before!=index:raise ValueError('Checkpoint must continue the same fixed-density branch')
    elif selector=='global_population':
        if source['mode']!='global' or role!='density_global' or before!=index:raise ValueError('Population restart must retain its density/root')
    elif selector=='task':
        expected_role='density_global' if role=='density_global_local' and source['mode']=='global' else role
        if source['role']!=expected_role:raise ValueError('Other branch winner rejected')
        if before!=index and not (task.endswith('_w0001') and ((branch=='ascending' and before==index-1) or (branch=='descending' and before==index+1))):raise ValueError('Directed continuation must use its adjacent own anchor')
    else:raise ValueError('Unrecognized profile source selector')

def main():
 p=argparse.ArgumentParser();p.add_argument('study',type=Path);a=p.parse_args();s=a.study.resolve()
 rows=w.snapshot(s);by={j['task']:j for j in rows}
 w.gate(s,'baseline_gate.tsv');ep=w.kv(s/'tables/profile_epoch.tsv')
 if ep['input_manifest_sha256']!=w.sha(s/'manifests/input_manifest.tsv') or ep['manifest_sha256']!=w.sha(s/ep['manifest']):raise ValueError('Stale profile epoch')
 w.verify_manifest(s,s/ep['manifest']);epoch=ep['epoch'];prefix='p_'+epoch+'_'
 gate=w.kv(s/'tables/profile_gate.tsv');w.verify_manifest(s,s/gate['evidence_manifest'])
 if gate['free_baseline_invalidated']=='true':raise ValueError('Free baseline invalidated; all new profile dispatch blocked')
 if gate['pass']=='true':print('All profile anchors accepted; no new allocation.');return
 anchors={(r['policy'],int(r['index'])):r for r in w.read_table(s/'tables/profile_anchor_gates.tsv')}
 history=w.read_table(s/'tables/profile_branch_history.tsv') if (s/'tables/profile_branch_history.tsv').exists() else []
 active_baseline=[j for j in rows if j['kind']=='baseline' and not w.closed(j)]
 if active_baseline:raise ValueError('Baseline accounting still active')
 w.verify_manifest(s,s/'tables/cache_source_manifest.tsv')
 # Copy an immutable epoch baseline as the dispatch selection, never a mutable
 # later free/profile winner. Its original native root/provenance remain intact.
 for policy in POLICIES:
  name='profile_'+epoch+'_'+policy+'_free_reference';path=s/'tables'/(name+'.rds')
  payload=(s/'profile_epochs'/epoch/('baseline_'+policy+'.rds')).read_bytes()
  if path.exists() and path.read_bytes()!=payload:raise ValueError('Profile free reference changed')
  if not path.exists():path.write_bytes(payload)
 def submit(task,**opts):
  if task in by:return
  verify_branch_source(s,epoch,task,opts,by)
  w.submit(s,task,kind='profile',root_family='profile',**opts)
  by[task]={'task':task,'state':'SUBMITTED','exit_code':''}
 def root(policy,branch,index):return 'profile_'+epoch+'_'+policy+'_'+(branch+'_r'+str(index).zfill(2) if branch=='global_local' else branch)
 def stem(policy,index,branch):return prefix+policy+'_r'+str(index).zfill(2)+'_'+branch
 for policy in POLICIES:
  for index in range(1,8):
   # Seven finite original density anchors are experimental design, not an
   # optimizer iteration/continuation cap.
   gstem=stem(policy,index,'global');gtask=gstem+'_w0001';rid=root(policy,'global_local',index)
   if gtask not in by:
    submit(gtask,mode='global',role='density_global',policy=policy,root_id=rid,selector='selection',selection_name='profile_'+epoch+'_'+policy+'_free_reference',density_index=index,init_mode='hybrid',seed=71000+1000*(0 if policy=='production' else 1)+index,cores=8,wall_seconds=3600,margin_seconds=600,max_calls=100000)
 # All scientifically ready globals run concurrently; no array/application cap.
 for j in rows:
  if j['kind']!='profile' or j['mode']!='global' or not j['task'].startswith(prefix) or not successful(j):continue
  index=int(j['density_index']);task=j['task']+'_polish_w0001'
  if anchors[(j['policy'],index)]['pass'].lower()!='true':
   submit(task,mode='local',role='density_global_local',policy=j['policy'],root_id=j['root_id'],selector='task',source_task=j['task'],density_index=index,wall_seconds=7200,margin_seconds=600,max_calls=20000)
 # Each directed branch begins at its own end of the original grid, then uses
 # only its own checked predecessor. No profile winner is copied across branches.
 for policy in POLICIES:
  for branch,order in (('ascending',range(1,8)),('descending',range(7,0,-1))):
   previous=None
   for index in order:
    task=stem(policy,index,branch)+'_w0001'
    if task not in by:
     if previous is None:
      if index!=next(iter(order)):break
      opts=dict(selector='selection',selection_name='profile_'+epoch+'_'+policy+'_free_reference')
     else:opts=dict(selector='task',source_task=previous['source_task'])
     submit(task,mode='local',role='density_'+branch,policy=policy,root_id=root(policy,branch,index),density_index=index,wall_seconds=3600,margin_seconds=600,max_calls=10000,**opts)
     break
    checked=[r for r in history if r['policy']==policy and int(r['index'])==index and r['branch']==branch]
    if not checked:break
    best=min(checked,key=lambda r:float(r['objective']))
    if anchors[(policy,index)]['pass'].lower()=='true' or (best['termination']=='poll_stationary' and float(best['directional_gain'])<=.002):previous=best
    else:break
 # One combined fresh/full original-19-coordinate check per local endpoint.
 locals=[j for j in rows if j['kind']=='profile' and j['mode']=='local' and j['task'].startswith(prefix)]
 for j in locals:
  if successful(j):submit(j['task']+'_check',mode='profile_fresh',role='density_check',policy=j['policy'],root_id=j['root_id'],selector='task',source_task=j['task'],density_index=int(j['density_index']),wall_seconds=900,margin_seconds=120,max_calls=500)
 latest={}
 for j in locals:latest[re.sub('_w[0-9]+$','',j['task'])]=j
 for component,j in latest.items():
  check=by.get(j['task']+'_check')
  if not successful(j) or not successful(check):continue
  index=int(j['density_index']);z=w.kv(s/'fits'/j['task']/'result.tsv');d=w.kv(s/'tests'/check['task']/'result.tsv')
  decision=next_local_decision(z,d,anchors[(j['policy'],index)]['pass'].lower()=='true')
  if decision in ('accepted','stationary'):continue
  task=next_id(j['task']);selector='checkpoint' if decision=='checkpoint' else 'task';source=j['task'] if decision=='checkpoint' else check['task']
  submit(task,mode='local',role=j['role'],policy=j['policy'],root_id=j['root_id'],selector=selector,source_task=source,density_index=index,wall_seconds=int(j['wall_seconds']),margin_seconds=600,max_calls=20000 if j['role']=='density_global_local' else 10000)
 # A new population window is justified after its own local polish is fully
 # checked/stationary and branch recovery still fails, with measured global gain.
 latest_global={}
 for j in rows:
  if j['kind']=='profile' and j['mode']=='global' and j['task'].startswith(prefix):latest_global[re.sub('_w[0-9]+$','',j['task'])]=j
 for component,j in latest_global.items():
  if not successful(j) or anchors[(j['policy'],int(j['density_index']))]['pass'].lower()=='true':continue
  component_locals=[r for r in locals if r['source_task']==j['task'] or r['task'].startswith(j['task']+'_polish_')]
  if not component_locals:continue
  local=component_locals[-1];check=by.get(local['task']+'_check')
  if not successful(local) or not successful(check):continue
  result=w.kv(s/'fits'/local['task']/'result.tsv')
  if result['termination']!='poll_stationary':continue
  z=w.kv(s/'fits'/j['task']/'result.tsv')
  if z['window_gain']=='NA':
   blocks=w.read_table(s/'fits'/j['task']/'blocks.tsv');gain=float(blocks[0]['objective'])-float(z['objective'])
  else:gain=float(z['window_gain'])
  if gain<=1e-10:print('PROFILE_GLOBAL_STALLED',j['task']);continue
  if not (s/'fits'/j['task']/'global_checkpoint.rds').exists():raise ValueError('No completed population checkpoint')
  submit(next_id(j['task']),mode='global',role='density_global',policy=j['policy'],root_id=j['root_id'],selector='global_population',source_task=j['task'],density_index=int(j['density_index']),cores=8,wall_seconds=3600,margin_seconds=600,max_calls=100000)
 w.snapshot(s)
if __name__=='__main__':main()
