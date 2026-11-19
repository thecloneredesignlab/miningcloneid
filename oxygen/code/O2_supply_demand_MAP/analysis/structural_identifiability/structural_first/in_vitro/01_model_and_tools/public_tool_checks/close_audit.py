"""Close metadata/evidence for Stage A-Alternative; never executes scientific code."""
from pathlib import Path
import csv
import datetime
import hashlib
import re
import shutil
import sys

def read(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream,delimiter='\t'))

def table(path,rows,fields=None):
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields or list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows)

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main(project,study):
    assert not (study/'COMPLETE').exists()
    jobs=read(study/'accounting/jobs.tsv')
    accounting={}
    for line in (study/'accounting/slurm_final_snapshot.txt').read_text().splitlines():
        bits=line.split('|')
        if len(bits)>=7:
            accounting[bits[0]]=bits
    parent_seconds=0
    for row in jobs:
        jid=row['job_id'];data=accounting[jid]
        assert data[1]=='COMPLETED' and data[2]=='0:0',data
        cpus=int(data[3]);elapsed=int(data[4]);parent_seconds+=cpus*elapsed
        # Be conservative about one-second Slurm parent/extern rounding differences.
        maximum=max(int(v[3])*int(v[4]) for k,v in accounting.items() if k==jid or k.startswith(jid+'.'))
        row.update(state=data[1],exit_code=data[2],allocated_cpus=cpus,elapsed_seconds=elapsed,
                   actual_core_seconds=maximum,reserved_core_seconds=0)
    table(study/'accounting/jobs.tsv',jobs)
    actual=sum(int(r['actual_core_seconds']) for r in jobs)
    assert actual<=24*3600
    ledger=read(study/'accounting/budget_ledger.tsv')
    for row in ledger:
        matching=[r for r in jobs if (r['task'].startswith('method_') if row['phase']=='method_controls' else
                   r['task']=='python_dependency_only_retry' if row['phase']=='one_technical_retry' else
                   r['task'] in ('matlab_runtime_gate','python_fallback_runtime_gate') if row['phase']=='runtime' else False)]
        row['actual_core_seconds']=sum(int(r['actual_core_seconds']) for r in matching)
        row['reserved_core_seconds']=0
        assert row['actual_core_seconds']<=float(row['cap_core_hours'])*3600
    table(study/'accounting/budget_ledger.tsv',ledger)
    table(study/'accounting/charged_summary.tsv',[dict(key=k,value=v) for k,v in {
        'Slurm_parent_allocated_core_seconds':parent_seconds,
        'conservative_budget_core_seconds':actual,'conservative_budget_core_hours':actual/3600,
        'authorized_core_hours':24,'reserved_core_seconds':0,'technical_retries_used':1,
        'technical_retries_remaining':0,'scheduler_submissions':4,'allocated_task_jobs':8}.items()])
    roles=read(study/'inputs/parameter_and_nuisance_audit.tsv')
    assert len(roles)==20 and sum(r['scientific_role']=='mechanism' for r in roles)==14
    table(study/'reports/parameter_roles.tsv',roles)
    mapping=[]
    for row in roles:
        param=row['parameter'];initial=param.startswith('init_')
        mapping.append(dict(parameter=param,role=row['scientific_role'],source=row['source_file'],line=row['source_line'],
            expression=row['source_expression'],
            FISPO='variance_only_outside_mean' if param=='sigma_growth' else 'initial_manifold_pullback_not_validated' if initial else 'analytic_ODE_candidate_runtime_blocked',
            GenSSI='variance_only_outside_mean' if param=='sigma_growth' else 'x0_substitution_candidate_runtime_blocked' if initial else 'symbolic_interface_candidate_runtime_blocked',
            StrikePy='variance_only_outside_mean' if param=='sigma_growth' else 'native_original_initial_family_interface_not_validated' if initial else 'original_joint_nonrational_map_not_tested',
            original_parameter_fixed_or_dropped=False,project_structural_result='not_analyzed'))
    table(study/'reports/model_tool_map.tsv',mapping)
    cases=read(study/'inputs/per_case_backend_status.tsv');base=[r for r in cases if r['tool']=='COMBOS']
    readiness=[]
    for tool in ('FISPO','GenSSI','StrikePy'):
        for row in base:
            readiness.append(dict(tool=tool,case_id=row['case_id'],model_form=row['model_form'],oxygen_pct=row['oxygen_pct'],outputs=row['outputs'],initial_role=row['initial_role'],unknowns=row['unknowns'],
                source_license_gate='pass',runtime='pass_private_dependency_adapted' if tool=='StrikePy' else 'hold_JSON_startup_dependency',
                method_gate='incomplete_controls' if tool=='StrikePy' else 'not_run',
                protocol_map='ODE_result_does_not_establish_D' if row['model_form']=='D' else 'original_joint_map_unvalidated',
                compute_released=False,project_structural_result='not_analyzed'))
    assert len(readiness)==42
    table(study/'reports/project_case_readiness.tsv',readiness)
    manifest=read(study/'controls/control_manifest.tsv');results=[]
    for tool in ('FISPO','GenSSI','StrikePy'):
        for case in manifest:
            cid=case['id'];status='not_run_runtime_hold';application_exit='not_run'
            if tool=='StrikePy':
                status='not_run_upstream_control_gate_and_initial_interface_hold'
                if cid.startswith('C'):
                    folder=study/'controls'/cid
                    application_exit=(folder/'application_status.tsv').read_text().strip().split('\t')[1]
                    if cid=='C01':
                        assert application_exit=='0'
                        assert read(folder/'control_result.tsv')[0]['observed_match']=='True'
                        assert 'NATIVE_METHOD_CONTROL_COMPLETED C01' in (folder/'console.log').read_text()
                        status='passed_published_small_control'
                    else:
                        text=(folder/'console.log').read_text()
                        assert 'same number of dimensions' in text and 'NATIVE_METHOD_CONTROL_COMPLETED' not in text
                        status='invalid_fixture_shape_before_matrix_analysis'
            results.append(dict(tool=tool,control=cid,scope=case['scope'],expected=case['expected'],
                derivative_order_cap=case['derivative_order_cap'],status=status,application_exit_code=application_exit,
                project_result=False))
    table(study/'reports/control_results.tsv',results)
    table(study/'reports/runtime_gate.tsv',[
        dict(configuration='MATLAB_R2022b_nojvm',backend='FISPO_and_GenSSI',gate='failed',application_exit=1,detail='Settings JSON access denied; JSON support libraries loaded; no toolbox entry reached'),
        dict(configuration='Python_original_dependencies',backend='StrikePy',gate='failed',application_exit=1,detail='json import blocked through unused debugger dependency'),
        dict(configuration='Python_private_dependency_only_adaptation',backend='StrikePy',gate='passed_scoped_runtime_probe',application_exit=0,detail='One unused import removed; mathematical AST unchanged; policy trace and symbolic startup passed')])
    accesses=[]
    for p in sorted(study.glob('runtime/**/trace*'))+sorted(study.glob('controls/C*/trace*')):
        if not p.is_file(): continue
        for line in p.read_text().splitlines():
            match=re.search(r'(openat\([^,]+,\s*"|open\("|execve\(")([^"\n]+)"',line)
            ret=re.search(r'=\s*(-?\d+)\s*$',line)
            if not match or not ret or int(ret.group(1))<0: continue
            path=match.group(2);low=path.lower();kind=None
            if '.json' in low: kind='JSON_data_file'
            elif 'json' in Path(low).name and '.so' in low: kind='JSON_support_library'
            elif '/sys/java/' in low: kind='Java_configuration_metadata'
            elif any(t in low for t in ('libjvm','libjava.','libv8','libmozjs')): kind='forbidden_runtime_library'
            elif Path(low).name in ('java','node','nodejs','js') or low.endswith('.js'): kind='forbidden_process_or_script'
            if kind: accesses.append(dict(trace=str(p.relative_to(study)),kind=kind,path=path,successful_syscall=True))
    table(study/'reports/file_access_findings.tsv',accesses,['trace','kind','path','successful_syscall'])
    assert not any('python_retry' in r['trace'] or r['trace'].startswith('controls/') for r in accesses), 'Private configuration policy finding'
    assert not any(r['kind']=='JSON_data_file' for r in accesses), 'Successful JSON data read'
    assert sum(r['kind']=='JSON_support_library' for r in accesses)>=3
    for folder in [study/'runtime/python_retry']+[study/'controls'/c for c in ('C01','C02','C03','C04','C05')]:
        assert not any(p.stat().st_size for p in folder.glob('blocked*'))
    table(study/'reports/acceptance.tsv',[
        dict(gate='A0_source_license_freeze',audit_complete=True,passed=True,result='three licensed public candidates; dependencies pinned'),
        dict(gate='A1_runtime',audit_complete=True,passed='partial',result='MATLAB blocked; private StrikePy configuration passed after one dependency-only retry'),
        dict(gate='A2_method_controls',audit_complete=True,passed=False,result='one published control completed; four invalid fixture attempts; nonrational/negative controls not validated'),
        dict(gate='A3_initial_observer_interface',audit_complete=True,passed=False,result='all14+6 static roles retained; runtime controls held; original Gaussian initial map unvalidated'),
        dict(gate='A4_reuse_decision',audit_complete=True,passed=False,result='delivery complete with HOLD_validation_and_initial_mapping; Stage B and Step 5 not released')])
    protected=[]
    for row in read(study/'inputs/protected_worktree_hashes.tsv'):
        assert digest(project/row['path'])==row['sha256'],row['path']
        protected.append(dict(path=row['path'],sha256=row['sha256'],unchanged=True))
    table(study/'reports/protected_worktree_verification.tsv',protected)
    for row in read(study/'reports/frozen_source_verification.tsv'):
        assert digest(project/row['path'])==row['sha256']
    status=dict(stage_A_alternative_delivery_complete=True,scientific_status='HOLD_validation_and_initial_mapping',
        mechanisms_retained=14,nuisances_retained=6,source_license_candidates=3,upstream_method_calls=5,
        completed_method_controls=1,invalid_fixture_attempts=4,nonrational_controls_validated=0,
        nuisance_interface_controls_executed=0,technical_retries=1,live_jobs=0,
        conservative_allocated_core_hours=actual/3600,Slurm_parent_core_seconds=parent_seconds,
        reserved_core_hours=0,cap_core_hours=24,project_structural_cases=0,production_replays=0,
        synthetic_datasets=0,recovery_fits=0,profiles=0,tumor_fits=0,author_contacts=0,
        Stage_B_released=False,Step_5_released=False,general_solver_created=False,
        closed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    table(study/'reports/final_status.tsv',[dict(key=k,value=v) for k,v in status.items()])
    shutil.copy2(__file__,study/'releases/r03/close_audit.py')
    print('Audit closure tables ready. Conservative use:',actual,'core-seconds; parent allocation:',parent_seconds,'core-seconds; no live reservation.')

if __name__=='__main__':
    main(*(Path(p).resolve() for p in sys.argv[1:]))
