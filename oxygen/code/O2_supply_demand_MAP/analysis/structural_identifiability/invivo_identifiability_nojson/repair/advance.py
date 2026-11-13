#!/usr/bin/env python3
"""One bounded dependency/control pass; conditional waves need current gates."""
import argparse
import shutil
from pathlib import Path
import workflow as w

def snapshot(study):
    state=w.State(study)
    with state.locked():
        rows=state.load();w.reconcile(study,rows);state.save(rows)
    w.table(study/'tables/jobs.tsv',rows,w.FIELDS)
    return rows

def capture_starts(study,rows):
    for j in rows:
        if j['mode']!='local': continue
        root=study/'fits'/j['task'];start=root/'actual_local_start.rds'
        first=list((root/'evaluators').glob('*/improvement_00000001.rds'))
        if not start.exists() or len(first)!=1: continue
        target=study/'tests/start_identity_snapshots'/j['task']
        if (target/'origins.tsv').exists(): continue
        target.mkdir(parents=True,exist_ok=True)
        shutil.copy2(start,target/'start.rds');shutil.copy2(first[0],target/'first.rds')
        w.table(target/'origins.tsv',[dict(kind='start',path=str(start.relative_to(study)),sha256=w.sha(start)),dict(kind='first',path=str(first[0].relative_to(study)),sha256=w.sha(first[0]))],['kind','path','sha256'])

def gate_pass(study,name):
    p=study/'tables'/name
    return p.exists() and {r['key']:r['value'] for r in w.read_table(p)}.get('pass')=='true'

def publish(study,rows):
    phase={r['role']:r['phase'] for r in w.read_table(study/'manifests/budget.tsv')}
    actual=sum(int(j['actual_core_seconds']) for j in rows if w.closed(j) and j['actual_core_seconds'])
    used=sum(map(w.charge,rows));active=[j for j in rows if not w.closed(j)]
    lines=['# In vivo baseline repair and conditional profile status','',
      'Updated '+w.utc()+'. Current dispatch repair02 with a pure lineage-metadata bridge; existing repair01 allocations unchanged. Exact pinned SIF; culture protected.',
      '', 'Prior lifetime debit222988core-seconds, of which192497belong to the formal phase. New closed actual '+str(actual)+' core-seconds. New actual plus live '+str(used)+'. Lifetime actual plus live '+str(222988+used)+' ('+format((222988+used)/3600,'.8f')+' core-hours); formal actual plus live '+str(192497+used)+' ('+format((192497+used)/3600,'.8f')+' core-hours).',
      '', 'Repair used plus live '+str(sum(w.charge(j) for j in rows if phase[j['role']]=='repair'))+'/129600core-seconds. Profile used plus live '+str(sum(w.charge(j) for j in rows if phase[j['role']]=='profile'))+'/316800core-seconds. Overall formal/lifetime ceilings648000/691200.',
      '', 'Engineering: '+('PASS' if gate_pass(study,'engineering_gate.tsv') else 'PENDING/NOT ACCEPTED')+'. Both baselines: '+('PASS' if gate_pass(study,'baseline_gate.tsv') else 'PENDING/NOT ACCEPTED')+'. Profile: '+('PASS' if gate_pass(study,'profile_gate.tsv') else 'NOT RELEASED')+'.',
      '',str(len(rows))+' owned allocations; '+str(len(rows)-len(active))+' closed; '+str(len(active))+' active/uncertain.',
      '', '| Task | Job | State |','|---|---:|---|']
    lines.extend('| '+j['task']+' | '+j['job_id']+' | '+j['state']+' |' for j in active)
    lines+=['','No background submitter or automatic reserve wave. No density task before both baselines pass a current bound gate. Independent root histories cannot share winning parameters. Existing closed runs remain immutable. No culture operation or commit/push.']
    w.atomic_text(study/'STATUS.md','\n'.join(lines)+'\n')
    print('new_actual_ch',actual/3600,'lifetime_used_plus_live_ch',(222988+used)/3600,'formal_used_plus_live_ch',(192497+used)/3600)
    for j in active: print('ACTIVE',j['task'],j['job_id'],j['state'])

def main():
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('wave',choices=['status','engineering','cheap','extended','reserve','profile','profile_reserve']);a=p.parse_args();s=a.study.resolve()
    rows=snapshot(s);by={j['task']:j for j in rows}
    if any(w.closed(j) and (j['state']!='COMPLETED' or j['exit_code']!='0:0') for j in rows):
        publish(s,rows);print('STOP: unexpected closed allocation failure; no dependent submission.');return 2
    plans=w.read_table(s/'manifests/task_plan.tsv')
    activated={}
    if a.wave in ('reserve','profile_reserve'):
        file=s/'tables'/('reserve_activation.tsv' if a.wave=='reserve' else 'profile_reserve_activation.tsv')
        if file.exists(): activated={r['task']:r for r in w.read_table(file)}
    if a.wave!='status':
        for task in plans:
            if task['stage']!=a.wave or task['task'] in by: continue
            if a.wave in ('reserve','profile_reserve') and task['task'] not in activated: continue
            if a.wave in ('cheap','extended','reserve','profile','profile_reserve') and not gate_pass(s,'engineering_gate.tsv'): continue
            if a.wave in ('profile','profile_reserve') and not gate_pass(s,'baseline_gate.tsv'): continue
            if a.wave=='extended' and not gate_pass(s,'extended_'+task['policy']+'_gate.tsv'): continue
            if a.wave=='extended':
                cheap=[r['task'] for r in plans if r['stage']=='cheap']
                if not all(t in by and by[t]['state']=='COMPLETED' and by[t]['exit_code']=='0:0' for t in cheap): continue
            ready=all(parent in by and by[parent]['state']=='COMPLETED' and by[parent]['exit_code']=='0:0' for parent in (task['source_task'],task['required_check']) if parent)
            if not ready: continue
            if task['selector']=='density_selected':
                selection=s/'tables'/('density_'+task['policy']+'_'+task['density_index']+'_selected.rds')
                if not selection.exists(): continue
            selection=activated.get(task['task'],{}).get('selection_name')
            if task['selector']=='root_selected' and not selection: continue
            w.submit(s,task['task'],selection_name=selection)
            by={j['task']:j for j in w.State(s).load()}
    rows=snapshot(s);capture_starts(s,rows);publish(s,rows)
    return 0

if __name__=='__main__': raise SystemExit(main())
