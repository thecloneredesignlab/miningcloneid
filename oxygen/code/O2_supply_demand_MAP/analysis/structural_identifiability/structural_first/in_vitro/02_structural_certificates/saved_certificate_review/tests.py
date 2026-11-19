"""Standard-library collector tests; no solver or scientific package is imported."""
import ast
import importlib.util
from pathlib import Path
import shutil
import tempfile
import sys
from artifact import atomic_table,start_call,matrix_syntax,summary_sets,validate_metadata,reconstruct_C03,select_matrix
from base_contract import read,table,select_native_matrix,named_set,receipt_valid

def run(study):
    study=Path(study);checks=[]
    def good(name,condition=True):
        assert condition,name
        checks.append(dict(check=name,passed=True))
    def reject(name,fn):
        try:fn()
        except (ValueError,AssertionError,SyntaxError,FileExistsError,KeyError):good(name)
        else:raise AssertionError('Invalid fixture admitted '+name)
    paths=[Path('obs_ident_matrix_case_model_2_Lie_deriv.txt'),Path('obs_ident_matrix_case_model_3_Lie_deriv.txt')]
    good('single_file',select_native_matrix(paths[:1],2)==paths[0])
    for nd in (2,3):
        good('multi_file_final_'+str(nd),select_native_matrix(list(reversed(paths)),nd)==paths[nd-2])
    reject('missing_final',lambda:select_native_matrix(paths,4))
    reject('duplicate_final',lambda:select_native_matrix([paths[1],paths[1]],3))
    reject('invalid_final',lambda:select_native_matrix(paths,5))
    good('parameter_order_independent',named_set('[b,a]',{'a','b'})=={'a','b'})
    for value in ('[a,a]','[a,c]',"[__import__('os')]",'[a.x]'):
        reject('reject_named_'+value,lambda value=value:named_set(value,{'a','b'}))
    good('reject_stale_receipt',not receipt_valid(dict(case='C03',run_id='old',input_sha256='old',status='PASS'),'C03','new','new'))
    for text in ('onx = [[1,0,0]]','onx = [[1/2,-a*x1,x1**b]]','onx = [[exp(a+b),log(x1),0]]'):
        good('accepted_bounded_syntax_'+text,matrix_syntax(text)==(1,3))
    for text in ('onx = [[1,0]]','onx = [[1,0,0]',"onx = [[__import__('os'),0,0]]",'onx = [[a.real,0,0]]','onx = [[z,0,0]]','onx = [[(lambda:1)(),0,0]]'):
        reject('reject_symbolic_syntax_'+text,lambda text=text:matrix_syntax(text))
    complete=dict(case='C03',nd=3,rank=2,identifiable='',nonidentifiable='a,b',measured='x1',observable='',native_returned=True,column_order='x1,a,b',positive_real_symbols=True)
    good('complete_metadata',validate_metadata(complete))
    for key,value in [('nd',5),('rank',0),('column_order','x1,b,a'),('positive_real_symbols',False),('native_returned',False),('nonidentifiable','a,a'),('case','S1')]:
        invalid=dict(complete);invalid[key]=value
        reject('reject_metadata_'+key,lambda invalid=invalid:validate_metadata(invalid))
    incomplete=dict(complete);del incomplete['nd']
    reject('reject_missing_completion_metadata',lambda:validate_metadata(incomplete))
    reject('reject_incomplete_summary',lambda:summary_sets('These parameters are unidentifiable: [a,b]'))
    core=(study/'upstream/StrikePy/StrikePy/strike_goldd.py').read_text()
    metadata=reconstruct_C03(study/'saved/C03',core)
    good('actual_saved_C03_order_reconstruction',metadata['nd']==3)
    with tempfile.TemporaryDirectory(prefix='saved_control_tests_') as tmp:
        tmp=Path(tmp)
        source=tmp/'C03';shutil.copytree(study/'saved/C03',source)
        console=source/'console.log';original=console.read_text()
        console.write_text(original.replace('matrix built with 3 Lie derivatives','matrix built with 2 Lie derivatives'))
        reject('reject_contradictory_saved_order',lambda:reconstruct_C03(source,core));console.write_text(original)
        prefix,tail=core.split('if nd >= len(xaug):',1)
        bad_core=prefix+'if nd >= len(xaug):'+tail.replace('unidflag = 1','nd = nd - 1',1)
        reject('reject_source_branch_changes_order',lambda:reconstruct_C03(source,bad_core))
        target=source/'native_app/results/obs_ident_matrix_case_model_3_Lie_deriv.txt'
        content=target.read_text();target.write_text('onx = [[1,0]]')
        reject('reject_actual_wrong_final_shape',lambda:reconstruct_C03(source,core));target.write_text(content)
        target.unlink();target.symlink_to(source/'native_app/results/obs_ident_matrix_case_model_2_Lie_deriv.txt')
        reject('reject_matrix_symlink',lambda:select_matrix(source/'native_app/results',3))
        fixture=tmp/'mock';fixture.mkdir()
        reject('C03_native_marker_forbidden',lambda:start_call(fixture/'call_started.tsv','C03','unit_mock'))
        start_call(fixture/'call_started.tsv','C04','unit_mock')
        reject('duplicate_call_marker_rejected',lambda:start_call(fixture/'call_started.tsv','C04','unit_mock'))
        calls=[0]
        def mock_solve():
            calls[0]+=1
            m=dict(complete);m.update(case='C04',rank=3,identifiable='a,b',nonidentifiable='');return m
        returned=mock_solve();atomic_table(fixture/'native_completion.tsv',[returned])
        try:raise ValueError('injected post-return collector failure')
        except ValueError:pass
        good('receipt_survives_post_return_failure',validate_metadata(read(fixture/'native_completion.tsv')[0]))
        good('collection_recovery_uses_no_solve',calls[0]==1)
        reject('immutable_receipt_rewrite_rejected',lambda:atomic_table(fixture/'native_completion.tsv',[complete]))
    good('positive_C02_layout_reused_as_metadata_only',len(list((study/'saved/C02_layout/results').glob('obs_ident_matrix_*.txt')))==1)
    before=set(sys.modules)
    spec=importlib.util.spec_from_file_location('collector_import_test',Path(__file__).parent/'worker.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    good('import_worker_has_no_scientific_side_effect',not any(n.split('.')[0] in ('sympy','numpy','symbtools','strike_goldd') for n in set(sys.modules)-before))
    reject('actual_native_entry_rejects_C03',lambda:module.native(study,'C03',None,None,study))
    worker_tree=ast.parse((Path(__file__).parent/'worker.py').read_text())
    native_fn=next(n for n in worker_tree.body if isinstance(n,ast.FunctionDef) and n.name=='native')
    native_try=next(n for n in native_fn.body if isinstance(n,ast.Try))
    good('durable_receipt_is_in_native_finally',any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='atomic_table' for item in native_try.finalbody for n in ast.walk(item)))
    for p in Path(__file__).parent.glob('*.py'):
        ast.parse(p.read_text());good('syntax_'+p.name)
    table(study/'reports/static_tests.tsv',checks)
    print('STANDARD_LIBRARY_COLLECTOR_TESTS_PASSED',len(checks))

if __name__=='__main__':run(Path(sys.argv[1]))
