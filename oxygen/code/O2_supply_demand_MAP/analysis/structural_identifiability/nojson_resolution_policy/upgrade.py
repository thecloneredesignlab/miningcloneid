#!/usr/bin/env python3
"""Preserve one study and anchor while upgrading only the dispatch implementation."""
from pathlib import Path
import argparse
import shutil
import sys
import time
import uuid
sys.dont_write_bytecode=True
import policy as P
from migrate import cancel_dispatchers


def upgrade(project,study):
    L=P.L;W=P.W;state=L.State(study)
    P.verify_release(study)
    allowed=P.compatible_evidence_manifests(study)
    tests=L.read_kv(study/'policy/tests/controller_quality_validation.tsv')
    assert tests['pass'] and tests['policy_sha256']==L.sha(P.RELEASE/'source_manifest.tsv')
    with state.locked():
        data=state.load();L.refresh(state,data)
        assert not (study/'BLOCKED.txt').exists(),'An existing hold requires diagnosis'
        active=[r['detail'] for r in data['events'] if r['event']=='resolution_policy_activated']
        old=study/active[-1]
        assert old!=P.RELEASE and L.sha(old/'source_manifest.tsv') in allowed
        assert [r['detail'] for r in data['events'] if r['event']=='resolution_policy_fingerprint'][-1]==L.sha(old/'source_manifest.tsv')
        assert any(j['mode']=='policy_validate' and L.job_ok(j) for j in data['jobs'])
        P.evidence_results(data,study)
        anchor_before=P.anchor(data,study);counts_before=W.counts(data)
        rounds={n:(study/'reserve'/str(n)/'CLOSED').exists() for n in (1,2)}
        protected={j['job_id']:dict(j) for j in data['jobs'] if j['role'] in {'global','local'} and not L.is_closed(j)}
        folder=study/'policy/migration'/('quality_branch_'+uuid.uuid4().hex);folder.mkdir(parents=True)
        generation=(study/'state/CURRENT').read_text().strip()
        shutil.copytree(study/'state'/generation,folder/'original_generation')
        L.write_kv(folder/'before.tsv',dict(utc=L.now(),old_release=str(old.relative_to(study)),
            old_manifest_sha256=L.sha(old/'source_manifest.tsv'),new_manifest_sha256=L.sha(P.RELEASE/'source_manifest.tsv'),
            global_origins=counts_before[0],fitting_local_slots=counts_before[1],
            commitments=L.committed(data),anchor_task=anchor_before['task_id'],anchor_q=anchor_before['objective']))
        stopped=cancel_dispatchers(state,data)
        for _ in range(10):
            L.refresh(state,data)
            if all(L.is_closed(j) for j in data['jobs'] if j['job_id'] in stopped):break
            time.sleep(3)
        else:raise ValueError('Predecessor dispatcher exit is not authoritative')
        assert not any(j['role']=='controller' and not L.is_closed(j) for j in data['jobs'])
        assert W.counts(data)==counts_before and P.anchor(data,study)==anchor_before
        assert rounds=={n:(study/'reserve'/str(n)/'CLOSED').exists() for n in (1,2)}
        for jid,j in protected.items():
            current=next(r for r in data['jobs'] if r['job_id']==jid)
            for key in ['task','mode','role','phase','token','cpus','wall_seconds','reserved_ch']:
                assert current[key]==j[key]
        state.event(data,'resolution_policy_activated',str(P.RELEASE.relative_to(study)))
        state.event(data,'resolution_policy_fingerprint',L.sha(P.RELEASE/'source_manifest.tsv'))
        state.event(data,'resolution_quality_branch_adapted','Missing-family near-minimum local evidence is repaired within remaining original shared slots; thresholds and single-fit algorithms unchanged')
        state.save(data)
        j=P.policy_submit(state,data,project,'controller_resolution_quality_1','controller','controller','control')
        L.write_kv(folder/'after.tsv',dict(pass_gate=True,utc=L.now(),new_controller=j['job_id'],
            stopped_controllers=';'.join(stopped),protected_fitting_jobs=';'.join(protected),
            unchanged_anchor=True,unchanged_attempts=True,unchanged_rounds=True,
            reused_identical_R_evidence=True,model_propagations=0,optimizer_calls=0))
        print('Quality-branch dispatcher:',j['job_id'],'protected fitting jobs:',list(protected))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('project',type=Path);p.add_argument('study',type=Path);a=p.parse_args()
    P.initialize(a.study);upgrade(a.project,a.study)
