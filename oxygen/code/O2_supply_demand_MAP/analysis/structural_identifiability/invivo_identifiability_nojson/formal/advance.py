#!/usr/bin/env python3
"""One bounded control pass; no background loop and no automatic reserve wave."""
import argparse
import shutil
from pathlib import Path
import workflow as w

def snapshot(study):
    s=w.State(study)
    with s.locked():
        rows=s.load();w.reconcile(study,rows);s.save(rows)
    w.table(study/'tables/jobs.tsv',rows,w.FIELDS)
    return rows

def publish_status(study,rows):
    actual=sum(int(j['actual_core_seconds']) for j in rows if j['actual_core_seconds'] and w.closed(j))
    committed=sum(map(w.charge,rows));active=[j for j in rows if not w.closed(j)]
    def gate(name,default):
        p=study/'tables'/name
        if not p.exists(): return default
        values={r['key']:r['value'] for r in w.read_table(p)}
        return 'PASS' if values.get('pass')=='true' else 'NOT ACCEPTED'
    lines=['# Formal in vivo execution status','',
        'Updated '+w.utc()+'. Prior pilot debit: 30491 core-seconds. Formal cap: 648000 (180 core-hours). Lifetime cap: 691200 (192 core-hours).',
        '', 'Formal closed actual: '+str(actual)+' core-seconds ('+format(actual/3600,'.8f')+' core-hours). Formal actual plus live reservations: '+str(committed)+' core-seconds. Lifetime actual plus live reservations: '+str(30491+committed)+' core-seconds ('+format((30491+committed)/3600,'.8f')+' core-hours).',
        '', 'Engineering: '+gate('engineering_gate.tsv','PENDING')+'. Both-policy baseline: '+gate('baseline_gate.tsv','PENDING')+'. Complete profile: '+gate('profile_gate.tsv','NOT STARTED')+'.',
        '', str(len(rows))+' owned allocations; '+str(len(rows)-len(active))+' closed; '+str(len(active))+' active/uncertain.',
        '', '| Task | Job | State |','|---|---:|---|']
    lines.extend('| '+j['task']+' | '+j['job_id']+' | '+j['state']+' |' for j in active)
    lines+=['','All work uses the pinned exact SIF, unchanged frozen model, fresh-child scoring, red/xxlarge and no-requeue. No profile may start until both baselines pass their recovery/direction gates. Failed/partial attempts remain charged. Culture and shared production assets are untouched; no commit/push.',
        '', 'This is a bounded interactive control pass, with no background submitter. Slurm jobs survive laptop sleep; orchestration can pause. Read reports/resume_notes.md, reports/execution_log.md and the immutable approved plan for scope and finite reserve caps.']
    w.atomic_text(study/'STATUS.md','\n'.join(lines)+'\n')

def capture_start_snapshots(study,rows):
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

def main():
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('wave',choices=['baseline','profile','status']);a=p.parse_args()
    study=a.study.resolve();rows=snapshot(study);bytask={j['task']:j for j in rows}
    bad=[j for j in rows if w.closed(j) and (j['state']!='COMPLETED' or j['exit_code']!='0:0')]
    if bad:
        for j in bad: print('STOP_UNEXPECTED_ALLOCATION',j['task'],j['state'],j['exit_code'])
        publish_status(study,rows)
        return 2
    if a.wave!='status':
        roles={'baseline_global','baseline_local','baseline_fresh','baseline_directional'} if a.wave=='baseline' else {'density_global','density_global_local','density_continuation_local','density_fresh'}
        for task in w.read_table(study/'manifests/task_plan.tsv'):
            if task['role'] not in roles or task['task'] in bytask: continue
            ready=True
            for key in ('source_task','required_check'):
                parent=task[key]
                if parent and (parent not in bytask or bytask[parent]['state']!='COMPLETED' or bytask[parent]['exit_code']!='0:0'): ready=False
            if task['selector']=='density_selected':
                name='density_'+task['policy']+'_'+task['density_index']+'_selected.rds'
                if not (study/'tables'/name).exists(): ready=False
                prefix='d_'+task['policy']+'_'+task['density_index']
                if not all(prefix+suffix in bytask and bytask[prefix+suffix]['state']=='COMPLETED' for suffix in ('_local','_asc','_desc')): ready=False
            if ready:
                w.submit(study,task['task'])
                bytask={j['task']:j for j in w.State(study).load()}
    rows=snapshot(study)
    capture_start_snapshots(study,rows)
    publish_status(study,rows)
    actual=sum(int(j['actual_core_seconds']) for j in rows if j['actual_core_seconds'] and w.closed(j))
    active=[j for j in rows if not w.closed(j)]
    print('formal_actual_ch',actual/3600,'formal_used_plus_committed_ch',sum(map(w.charge,rows))/3600,'lifetime_used_plus_committed_ch',(30491+sum(map(w.charge,rows)))/3600)
    print('owned_jobs',len(rows),'active',len(active),'closed',len(rows)-len(active))
    for j in active: print('ACTIVE',j['task'],j['job_id'],j['state'])
    return 0

if __name__=='__main__': raise SystemExit(main())
