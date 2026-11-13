#!/usr/bin/env python3
"""One-writer handover; cancel only receipt-identified dispatch controllers."""
from pathlib import Path
import argparse
import shutil
import subprocess
import sys
import time
import uuid
sys.dont_write_bytecode=True
import policy as P


def read_generation(path):
    L=P.L
    for r in L.read_table(path/'manifest.tsv'):
        if L.sha(path/r['path'])!=r['sha256']:
            raise ValueError('Corrupt historical generation')
    return {n:L.read_table(path/(n+'.tsv')) for n in ['jobs','tasks','events','reservations']}


def completed_checkpoint(study,project):
    """Reconstruct a genuine old-policy checkpoint, not the newer A-only record."""
    for path in sorted((study/'state').glob('g_*'),key=lambda p:p.stat().st_mtime,reverse=True):
        jobs=P.L.read_table(path/'jobs.tsv')
        if not any(j['task']=='reserve_1_repair_2' and j['mode']=='check' and P.L.job_ok(j) for j in jobs):
            continue
        data=read_generation(path)
        bg,points,cs=P.legacy_gates(data,study,project)
        if bg['pass'] and all(r['reference_pass'] for r in points):
            return path.name,data,bg,cs
    raise ValueError('No completed reference checkpoint can be reconstructed')


def cancel_dispatchers(state,data):
    L=P.L
    live=[j for j in data['jobs'] if j['role']=='controller' and not L.is_closed(j)]
    if any(not j['job_id'].isdigit() or j['mode']!='controller' for j in live):
        raise ValueError('Unreconciled dispatcher receipt')
    ids=[j['job_id'] for j in live]
    state.event(data,'resolution_handover_controller_stop',' '.join(ids) or 'none');state.save(data)
    if ids:
        subprocess.run([str(L.SLURM/'scancel')]+ids,check=True,timeout=30)
    return ids


def activate(project,study):
    L=P.L;W=P.W;state=L.State(study)
    P.verify_release(study)
    test=L.read_kv(study/'policy/tests/controller_validation.tsv')
    assert test['pass'] and test['policy_sha256']==L.sha(P.RELEASE/'source_manifest.tsv')
    # Immutable old generations may be inspected while the current writer runs.
    proof_name,proof_data,bg,cs=completed_checkpoint(study,project)
    with state.locked():
        data=state.load();L.refresh(state,data)
        if (study/'BLOCKED.txt').exists():
            raise ValueError('Existing hold must be diagnosed before handover')
        active=[r for r in data['events'] if r['event']=='resolution_policy_activated']
        if active:
            raise ValueError('A policy is already active; do not repeat migration')
        successful=[j for j in data['jobs'] if j['mode']=='policy_validate' and L.job_ok(j)]
        assert successful,'HPC mapping/endpoint validation has not passed'
        evidence=P.evidence_results(data,study)
        assert any(r['action']=='endpoint' for r in evidence)
        assert {r['action'] for r in evidence}=={'endpoint','free_map'}
        before_counts=W.counts(data);before_rounds={n:(study/'reserve'/str(n)/'CLOSED').exists() for n in (1,2)}
        assert before_rounds[1],'Consumed first reserve marker is missing'
        snapshot=study/'policy/migration'/('handover_'+uuid.uuid4().hex)
        snapshot.mkdir(parents=True)
        generation=(study/'state/CURRENT').read_text().strip()
        shutil.copytree(study/'state'/generation,snapshot/'original_generation')
        for name in ['point_quality.tsv','baseline_gate.tsv','candidates.tsv','jobs.tsv']:
            shutil.copy2(study/'tables'/name,snapshot/name)
        receipts=[]
        for task in data['tasks']:
            if task['mode'] not in {'global','local'} or task['phase']=='fixture':continue
            result=L.accepted(data,study,task['task_id'])
            if result is None:continue
            d=study/'fits'/task['task_id']
            receipts.append(dict(task_id=task['task_id'],result_sha256=L.sha(d/'result.tsv'),
                fit_sha256=L.sha(d/'accepted_fit.rds'),payload_manifest_sha256=L.sha(d/'payload_manifest.tsv')))
        L.write_table(snapshot/'accepted_receipts.tsv',receipts,['task_id','result_sha256','fit_sha256','payload_manifest_sha256'])
        L.write_kv(snapshot/'before.tsv',dict(utc=L.now(),generation=generation,global_origins=before_counts[0],
            fitting_local_slots=before_counts[1],commitments=L.committed(data),nonfit_commitments=L.committed(data,'nonfit'),
            reserve_1_closed=before_rounds[1],reserve_2_closed=before_rounds[2],
            anchor_proof_generation=proof_name,accepted_payloads=len(receipts)))
        protected={j['job_id']:dict(j) for j in data['jobs'] if j['role'] in {'global','local'} and not L.is_closed(j)}
        stopped=cancel_dispatchers(state,data)
        # The lock prevents every old dispatcher from entering another submission cycle.
        for attempt in range(10):
            L.refresh(state,data)
            if all(L.is_closed(j) for j in data['jobs'] if j['job_id'] in stopped):break
            time.sleep(3)
        else:raise ValueError('Old dispatcher exit is not yet authoritative; migration remains held')
        assert not any(j['role']=='controller' and not L.is_closed(j) for j in data['jobs'])
        assert W.counts(data)==before_counts
        assert {n:(study/'reserve'/str(n)/'CLOSED').exists() for n in (1,2)}==before_rounds
        for jid,j in protected.items():
            current=next(r for r in data['jobs'] if r['job_id']==jid)
            for key in ['task','mode','role','phase','token','cpus','wall_seconds','reserved_ch']:
                assert current[key]==j[key]
        P.save_anchor(state,data,project,bg,cs,proof_name,'migration_reconstructed_completed_checkpoint',proof_data)
        state.event(data,'resolution_policy_activated',str(P.RELEASE.relative_to(study)))
        state.event(data,'resolution_policy_fingerprint',L.sha(P.RELEASE/'source_manifest.tsv'))
        state.event(data,'resolution_migration_snapshot',str(snapshot.relative_to(study)));state.save(data)
        new=P.policy_submit(state,data,project,'controller_resolution_1','controller','controller','control')
        L.write_kv(snapshot/'after.tsv',dict(pass_gate=True,utc=L.now(),policy_id=P.POLICY_ID,
            new_controller_job_id=new['job_id'],stopped_dispatchers=';'.join(stopped),
            unchanged_fitting_attempts=True,unchanged_reserve_markers=True,
            protected_fitting_allocations=len(protected),global_origins=W.counts(data)[0],
            fitting_local_slots=W.counts(data)[1],model_propagations=0,optimizer_calls=0))
        print('Activated',P.POLICY_ID,'controller',new['job_id'],'protected fitting jobs',list(protected))


def main():
    p=argparse.ArgumentParser();p.add_argument('project',type=Path);p.add_argument('study',type=Path);a=p.parse_args()
    P.initialize(a.study);activate(a.project,a.study)


if __name__=='__main__':main()
