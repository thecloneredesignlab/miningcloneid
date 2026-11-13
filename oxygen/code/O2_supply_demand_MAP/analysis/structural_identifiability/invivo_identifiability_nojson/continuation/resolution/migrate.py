"""Freeze the approved bridge and activate a sole finite-policy dispatcher."""
from pathlib import Path
import argparse,shutil
import workflow as w
import policy as p

BRIDGE={'b_data_only_anchor_w0003':'21846058','b_production_anchor_w0003':'21849242',
 'b_production_uniform_old_w0003':'21849243','b_production_uniform_new_w0003':'21849244',
 'b_data_only_uniform_old_w0003':'21849245','b_data_only_uniform_new_w0003':'21849435'}

def migrate(s):
    s=Path(s);state=w.State(s);directory=s/'policies/resolution01'
    with state.locked():
        if (s/'tables/resolution_policy.tsv').exists():raise ValueError('Already migrated; reconcile only')
        rows=state.load();w.reconcile(s,rows);state.save(rows)
        by={j['task']:j for j in rows};assert set(BRIDGE)<=set(by)
        if any(not j['job_id'] for j in rows):raise ValueError('Uncertain receipt')
        directory.mkdir(parents=True,exist_ok=True)
        if (directory/'manifest.tsv').exists():raise ValueError('Frozen migration already exists')
        snap=directory/'legacy_snapshot';snap.mkdir()
        for name in ('jobs.tsv','raw_accounting.txt','reviewed_failures.tsv','engineering_gate.tsv',
                     'engineering_gate_evidence.tsv','baseline_gate.tsv','baseline_policy_gates.tsv',
                     'execution_release.tsv','metadata_execution_release.tsv','cache_source_manifest.tsv'):
            src=s/'tables'/name
            if name=='jobs.tsv':w.table(snap/name,rows,w.FIELDS)
            elif src.exists():shutil.copy2(src,snap/name)
        shutil.copytree(str(s/'implementation'),str(snap/'implementation'),ignore=shutil.ignore_patterns('__pycache__'))
        bridge=[];ledger=[]
        for task,job_id in BRIDGE.items():
            j=by[task];assert j['job_id']==job_id and j['mode']=='local' and j['release']=='cont06'
            bridge.append({k:j[k] for k in ('task','job_id','root_id','root_family','source_task','release','spec_sha256')})
            ledger.append(dict(task=task,phase='baseline',slot_class='existing_local',slot_key=j['root_id'],
              source_task=j['source_task'],decision_fit='',reservation_seconds=j['reserved_core_seconds'],
              polish_escrow='0',check_escrow='900',consumed_parent='',policy_version=p.VERSION,recorded_utc=w.utc()))
        w.table(directory/'bridge.tsv',bridge,list(bridge[0]))
        historical=301846+sum(int(j['actual_core_seconds']) for j in rows if w.closed(j) and j['task'] not in BRIDGE)
        assert historical==388524
        w.typed(directory/'transition.tsv',dict(policy_version=p.VERSION,created_utc=w.utc(),
            generation=(s/'state/CURRENT').read_text().strip(),historical_lifetime_core_seconds=historical,
            inherited_parents=6,baseline_cap_core_seconds=p.CAPS['baseline'],profile_cap_core_seconds=p.CAPS['profile'],
            closure_reserve_per_phase_seconds=p.CLOSURE,science_release='cont06',check_release='res01'))
        # Preserve the old host source before guarding its dispatch entry point.
        old=s/'implementation/workflow.py';text=old.read_text();needle='def submit(study,task,**options):\n'
        assert needle in text
        text=text.replace(needle,needle+"    if (Path(study)/'tables/resolution_policy.tsv').exists():raise ValueError('Legacy dispatcher disabled by resolution01; use the finite-policy entry point')\n",1)
        w.table(directory/'manifest.tsv',[dict(path=str(f.relative_to(s)),sha256=w.sha(f))
            for f in sorted(directory.rglob('*')) if f.is_file()],['path','sha256'])
        w.verify_manifest(s,directory/'manifest.tsv')
        w.table(s/'tables/resolution_assignments.tsv',ledger,p.FIELDS)
        w.table(s/'tables/resolution_policy.tsv',[dict(key='policy_version',value=p.VERSION),
          dict(key='manifest',value='policies/resolution01/manifest.tsv'),
          dict(key='manifest_sha256',value=w.sha(directory/'manifest.tsv'))],['key','value'])
        w.atomic_text(old,text)
    w.snapshot(s)
    print('Migrated six inherited jobs once; old dispatcher blocked; compatible cont06 numerical runtime retained.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('study',type=Path)
    migrate(parser.parse_args().study.resolve())
