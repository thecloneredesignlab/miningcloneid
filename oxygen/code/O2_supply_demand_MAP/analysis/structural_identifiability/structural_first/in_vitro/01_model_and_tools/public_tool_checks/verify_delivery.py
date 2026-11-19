"""Independent stdlib delivery verifier; imports no upstream scientific software."""
import ast
import csv
import hashlib
from pathlib import Path
import sys

def read(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream,delimiter='\t'))

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path,rows,fields=None):
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields or list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows)

def main(project,study,verify_only=False):
    if verify_only:
        rows=read(study/'collection_manifest.tsv')
        for row in rows:
            p=study/row['path'];assert p.is_file() and p.stat().st_size==int(row['bytes']) and sha(p)==row['sha256'],p
        assert (study/'COMPLETE').read_text().splitlines()[0]==sha(study/'collection_manifest.tsv')
        print('READ_ONLY_MANIFEST_VERIFICATION_PASSED',len(rows),'files')
        return
    assert not (study/'COMPLETE').exists()
    checks=[]
    def check(name,condition):
        assert condition,name
        checks.append(dict(check=name,passed=True))
    for filename in ('source_manifest.tsv','dependency_source_manifest.tsv','runtime_source_manifest.tsv'):
        for row in read(study/'inputs'/filename):
            check(filename+':'+row['path'],sha(study/row['path'])==row['sha256'])
    for row in read(study/'reports/frozen_source_verification.tsv'):
        check('production:'+row['path'],sha(project/row['path'])==row['sha256'])
    for row in read(study/'inputs/protected_worktree_hashes.tsv'):
        check('protected:'+row['path'],sha(project/row['path'])==row['sha256'])
    roles=read(study/'reports/parameter_roles.tsv')
    check('all_twenty_unique_original_roles',len(roles)==20 and len({r['parameter'] for r in roles})==20)
    check('fourteen_mechanisms',sum(r['scientific_role']=='mechanism' for r in roles)==14)
    check('six_nuisances',sum(r['scientific_role']!='mechanism' for r in roles)==6)
    check('no_project_parameter_classification',all(r['structural_result']=='not_analyzed' for r in roles))
    ready=read(study/'reports/project_case_readiness.tsv')
    check('forty_two_backend_case_form_dispositions',len(ready)==42)
    check('no_project_case_released',all(r['compute_released']=='False' for r in ready))
    status={r['key']:r['value'] for r in read(study/'reports/final_status.tsv')}
    check('HOLD_not_structural_completion',status['scientific_status']=='HOLD_validation_and_initial_mapping')
    for key in ('Stage_B_released','Step_5_released','general_solver_created'):
        check(key,status[key]=='False')
    for key in ('project_structural_cases','production_replays','synthetic_datasets','recovery_fits','profiles','tumor_fits','author_contacts','live_jobs'):
        check(key,status[key]=='0')
    jobs=read(study/'accounting/jobs.tsv')
    check('eight_owned_task_jobs',len(jobs)==8 and len({r['job_id'] for r in jobs})==8)
    check('only_xxlarge',all(r['qos']=='xxlarge' for r in jobs))
    check('all_terminal_parent_exit_zero',all(r['state']=='COMPLETED' and r['exit_code']=='0:0' for r in jobs))
    check('no_live_reservations',all(int(r['reserved_core_seconds'])==0 for r in jobs))
    check('independent_24_hour_cap',sum(int(r['actual_core_seconds']) for r in jobs)<=86400)
    check('budget_phase_sum_matches_jobs',sum(int(r['actual_core_seconds']) for r in read(study/'accounting/budget_ledger.tsv'))==sum(int(r['actual_core_seconds']) for r in jobs))
    check('single_technical_retry',sum(r['task']=='python_dependency_only_retry' for r in jobs)==1)
    check('owned_queue_empty',not (study/'accounting/owned_squeue_final.txt').read_text().strip())
    check('native_io_guard_negative_control','NATIVE_IO_GUARD_NEGATIVE_CONTROL_PASSED' in (study/'runtime/guard_control.log').read_text())
    check('MATLAB_JSON_requirement_recorded','factoryFileLocation.json' in (study/'runtime/matlab/blocked_access.tsv').read_text())
    check('original_Python_JSON_import_prevented',(study/'runtime/python/blocked_imports.tsv').read_text().strip()=='BLOCKED\timport\tjson')
    check('private_runtime_completed','PYTHON_RUNTIME_GATE_COMPLETED' in (study/'runtime/python_retry/console.log').read_text())
    for folder in [study/'runtime/python_retry']+[study/'controls'/c for c in ('C01','C02','C03','C04','C05')]:
        check(str(folder.relative_to(study))+':no_blocked_forbidden_operation',not any(p.stat().st_size for p in folder.glob('blocked*')))
        module_file=folder/('final_imported_modules.tsv' if folder.parent.name=='controls' else 'imported_modules.tsv')
        modules=read(module_file)
        check(str(folder.relative_to(study))+':no_loaded_JSON_module',not any(r['module'].split('.')[0] in ('json','_json','simplejson','orjson','ujson') for r in modules))
    orig=study/'sources/dependencies/symbtools/src/symbtools'
    private=study/'runtime/private_prefix/symbtools'
    for p in orig.glob('*.py'):
        if p.name=='time_deriv.py':
            text=p.read_text();check('only_unused_import_removed',text.replace('from ipydex import IPS\n','')==(private/p.name).read_text())
            a=ast.parse(text);b=ast.parse((private/p.name).read_text());a.body=[n for n in a.body if not isinstance(n,ast.ImportFrom) or n.module!='ipydex']
            check('mathematical_AST_unchanged',ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False))
        else:
            check('private_source_unchanged:'+p.name,sha(p)==sha(private/p.name))
    for cid in ('C01','C02','C03','C04','C05'):
        p=study/'controls'/cid/'native_app/strike_goldd.py'
        check(cid+':native_core_byte_identical',sha(p)==sha(study/'sources/StrikePy/StrikePy/strike_goldd.py'))
    results=read(study/'reports/control_results.tsv')
    check('one_completed_benchmark_only',sum(r['status']=='passed_published_small_control' for r in results)==1)
    check('four_invalid_attempts_not_false_failures',sum(r['status']=='invalid_fixture_shape_before_matrix_analysis' for r in results)==4)
    check('corrected_fixture_unexecuted',(study/'releases/r04_prepared_not_executed/method_control.py').is_file() and len(jobs)==8)
    for p in study.rglob('*'):
        if p.is_file():
            check('permitted_artifact:'+str(p.relative_to(study)),p.suffix.lower() not in ('.json','.jsonl','.json5','.js','.java','.jar','.ipynb'))
    for p in [study/'REPORT.md',*list((study/'reports').glob('*.md'))]:
        check('English_authored_report:'+p.name,not any('\u4e00'<=c<='\u9fff' for c in p.read_text()))
    write(study/'reports/delivery_verification.tsv',checks)
    manifest=[]
    for p in sorted(study.rglob('*')):
        if p.is_file() and p.name not in ('collection_manifest.tsv','COMPLETE'):
            manifest.append(dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)))
    write(study/'collection_manifest.tsv',manifest)
    (study/'COMPLETE').write_text(sha(study/'collection_manifest.tsv')+'\n'+str(len(manifest))+' manifest-listed files\n')
    print('DELIVERY_CHECKS_PASSED',len(checks),'checks;',len(manifest),'files; closed with a scientific HOLD.')

if __name__=='__main__':
    main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve(),'--verify-only' in sys.argv[3:])
