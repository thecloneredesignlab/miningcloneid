"""Guarded HPC preflight/native toy control. No culture simulation or fitting."""
import ast
import datetime
import importlib.util
import inspect
import os
from pathlib import Path
import random
import re
import shutil
import sys
from types import SimpleNamespace
from contract import CASES,read,table,sha,static_contract,named_set,safe_matrix,select_native_matrix

study=Path(sys.argv[1]);mode=sys.argv[2];case=sys.argv[3]
assert mode in ('preflight','control') and case in CASES,'Scope rejected before scientific imports'
out=study/('preflight' if mode=='preflight' else 'controls')/case
out.mkdir(parents=True,exist_ok=True)
assert not (out/'result.tsv').exists(),'No stale result reuse'
for row in read(study/'source_manifest.tsv'):
    assert sha(study/row['path'])==row['sha256'],'Source drift '+row['path']

def audit(event,args):
    root=args[0].split('.')[0] if event=='import' else ''
    if root=='ipydex':raise ImportError('Previously approved optional debugger exclusion')
    if root in ('json','_json','simplejson','orjson','ujson','rapidjson','jpype','jnius','execjs','py_mini_racer'):
        with (out/'blocked_imports.tsv').open('a') as stream:stream.write('BLOCKED\timport\t'+args[0]+'\n')
        raise ImportError('Forbidden analysis dependency '+args[0])

sys.addaudithook(audit)
sys.path[:0]=[str(study/'runtime/private_prefix'),'/home/4482173/.local/lib/python3.11/site-packages',*[p for p in os.environ.get('PYTHONPATH','').split(':') if p]]

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module);return module

def flatten(value):
    if isinstance(value,sp.MatrixBase):value=list(value)
    if isinstance(value,np.ndarray):value=value.tolist()
    if isinstance(value,(list,tuple)):
        return [item for element in value for item in flatten(element)]
    return [value]

def packing(model):
    namespace=dict(model=model,np=np,zeros=sp.zeros,w1vector=[],w1vector_dot=[])
    exec(compile((study/'inputs/native_packing_span.py').read_text(),'pinned_native_packing_span','exec'),namespace)
    return namespace['xaug'],namespace['faug']

try:
    import numpy as np
    import sympy as sp
    import symbtools as st
    assert np.__version__=='2.4.6' and sp.__version__=='1.12' and st.__version__=='0.4.1'
    app=out/'native_app';assert not app.exists();(app/'models').mkdir(parents=True)
    (app/'models/__init__.py').write_text('')
    native=study/'upstream/StrikePy/StrikePy'
    shutil.copy2(native/'strike_goldd.py',app/'strike_goldd.py')
    shutil.copytree(native/'functions',app/'functions')
    shutil.copy2(study/'fixtures'/case/'model.py',app/'models/case_model.py')
    shutil.copy2(study/'fixtures'/case/'options.py',app/'options.py')
    static_contract(case,(app/'models/case_model.py').read_text(),(app/'options.py').read_text())
    sys.path.insert(0,str(app));os.chdir(app)
    options=load('options',app/'options.py')
    load('models',app/'models/__init__.py')
    model=load('models.case_model',app/'models/case_model.py')
    core=load('native_strike',app/'strike_goldd.py')
    from functions.elim_and_recalc import elim_and_recalc
    x=model.x1;params=[model.a] if case=='C02' else [model.a,model.b]
    symbols={'x1':x,'a':model.a,'b':model.b}
    expected_symbols=set([x]+params)
    assert all(s.is_positive is True for s in expected_symbols)
    assert set(model.f[0].free_symbols)|set(model.h[0].free_symbols)==expected_symbols
    assert [p[0] for p in model.p]==params and len(set(params))==len(params)
    assert sys.modules['options'] is options and Path(options.__file__).resolve()==app/'options.py'
    xaug,faug=packing(model);q=len(params);n=q+1
    assert xaug.shape==(n,1) and faug.shape==(n,1)
    assert sp.Matrix(xaug)==sp.Matrix([x]+params)
    assert sp.Matrix(faug)==sp.Matrix([model.f[0]]+[0]*q)
    first=sp.Matrix(model.h).jacobian(sp.Matrix(xaug));assert first==sp.Matrix([[1]+[0]*q])
    inspect.signature(elim_and_recalc).bind([],n,sp.eye(n),model.p,model.x,1,[],[],[],[])
    primitive=sp.Matrix([[1,0,0],[0,1,1],[0,2,2]])
    assert primitive.col([0,2])==primitive.extract([0,1,2],[0,2])
    table(out/'interface_checks.tsv',[dict(check=k,passed=True) for k in ('same_symbol_objects','positive_domains','free_symbol_roles','options_resolution','native_column_shapes','native_column_entries','zero_parameter_dynamics','initial_output_Jacobian','helper_signature_binding','helper_column_selection')])
    table(out/'arrays.tsv',[dict(array='xaug',shape=str(xaug.shape),value=str(sp.Matrix(xaug))),dict(array='faug',shape=str(faug.shape),value=str(sp.Matrix(faug))),dict(array='output_Jacobian',shape=str(first.shape),value=str(first))])
    if mode=='preflight':
        bad=SimpleNamespace(x=model.x,p=model.p,f=[[model.f[0]]])
        try:packing(bad)
        except ValueError:pass
        else:raise AssertionError('Old nested dynamics unexpectedly admitted')
        if case=='C02':
            regression=load('published_packing_only',study/'inputs/published_C2M_source.py')
            xx,ff=packing(regression)
            assert xx.shape==(6,1) and ff.shape==(6,1)
            assert sp.Matrix(xx)==sp.Matrix(flatten(regression.x)+flatten(regression.p))
            assert sp.Matrix(ff)==sp.Matrix(flatten(regression.f)+[0]*4)
            table(out/'two_state_regression.tsv',[dict(test='C01_packing_only_no_solver',passed=True,shape='6 by 1')])
        result=dict(case=case,status='PASS',scope='interface_preflight_only',run_id=study.name,input_sha256=sha(study/'fixtures'/case/'model.py'),native_method_calls=0)
    else:
        seed=20261116+int(case[1:]);random.seed(seed);np.random.seed(seed)
        captures={};rank_calls=[];substitutions=[]
        def profile(frame,event,arg):
            if event!='return':return
            filename=Path(frame.f_code.co_filename)
            if frame.f_code.co_name=='strike_goldd' and filename==app/'strike_goldd.py':
                captures.update(frame.f_locals)
            if frame.f_code.co_name=='generic_rank' and filename==study/'runtime/private_prefix/symbtools/core.py':
                values=frame.f_locals;idx=len(rank_calls)
                rank_calls.append(dict(call=idx,seed=values.get('seed'),rank=str(arg),precisions='100,200,300; conditional native 400',native_implementation_sha256=sha(filename)))
                for symbol,value in values.get('rnst',[]):substitutions.append(dict(call=idx,symbol=str(symbol),value=str(value)))
        sys.setprofile(profile)
        try:core.strike_goldd()
        finally:sys.setprofile(None)
        table(out/'rank_calls.tsv',rank_calls,['call','seed','rank','precisions','native_implementation_sha256'])
        table(out/'random_substitutions.tsv',substitutions,['call','symbol','value'])
        required={'p_id','p_un','nd','onx','meas_x','obs_states','rango'}
        assert required.issubset(captures),'Native result frame incomplete'
        actual_id={str(v) for v in flatten(captures['p_id'])};actual_un={str(v) for v in flatten(captures['p_un'])}
        expected={str(p) for p in params};positive=case in ('C02','C04')
        assert actual_id==(expected if positive else set()) and actual_un==(set() if positive else expected),'Native parameter classification mismatch or unresolved'
        assert {str(x)} <= {str(v) for v in flatten(captures['meas_x'])+flatten(captures['obs_states'])}
        nd=int(captures['nd']);assert 1<=nd<=4
        summaries=list((app/'results').glob('id_results_*.txt'));matrices=list((app/'results').glob('obs_ident_matrix_*.txt'))
        assert len(summaries)==1
        final_matrix=select_native_matrix(matrices,nd)
        summary=summaries[0].read_text()
        if positive:assert 'All its parameters are' in summary
        else:
            match=re.search(r'These parameters are unidentifiable:\s*(\[[^\]]*\])',summary)
            assert match and named_set(match.group(1),expected)==expected
        observed=safe_matrix(final_matrix.read_text(),sp,symbols)
        assert observed.shape==(nd+1,n) and observed==sp.Matrix(captures['onx'])
        a=model.a;b=model.b
        if case=='C04':
            jets=[x,-a*x**b,a*a*b*x**(2*b-1),-a**3*b*(2*b-1)*x**(3*b-2),a**4*b*(2*b-1)*(3*b-2)*x**(4*b-3)]
        else:
            k=a if case=='C02' else a+b if case=='C03' else sp.exp(a+b)
            jets=[(-k)**i*x for i in range(5)]
        reference=sp.Matrix(jets[:nd+1]).jacobian(sp.Matrix([x]+params))
        residuals=[]
        for i in range(nd+1):
            for j in range(n):
                residual=sp.simplify(observed[i,j]-reference[i,j]);assert residual==0,'Native row mismatch'
                residuals.append(dict(row=i,column=j,residual=str(residual)))
        table(out/'row_checks.tsv',residuals)
        if positive:
            determinant=sp.simplify(observed[:n,:].det())
            witness=-x if case=='C02' else -a*a*x**(3*b-1)
            assert sp.simplify(determinant-witness)==0
            rank=n
        else:
            assert all(sp.simplify(v)==0 for v in observed[:,1]-observed[:,2])
            witness=-x if case=='C03' else -sp.exp(a+b)*x
            assert sp.simplify(observed[:2,:2].det()-witness)==0
            c=sp.Symbol('c',real=True)
            assert sp.simplify(model.f[0].subs({a:a+c,b:b-c},simultaneous=True)-model.f[0])==0
            determinant=0;rank=2
        assert int(captures['rango'])==rank
        table(out/'native_named_results.tsv',[dict(case=case,identifiable=','.join(sorted(actual_id)),nonidentifiable=','.join(sorted(actual_un)),rank=rank,derivative_order=nd,exact_witness=str(witness),determinant=str(determinant),proof_scope='toy control only; no culture parameter claim')])
        result=dict(case=case,status='PASS',scope='native_toy_control_only',run_id=study.name,input_sha256=sha(study/'fixtures'/case/'model.py'),native_method_calls=1,seed=seed,rank=rank,derivative_order=nd)
    table(out/'result.tsv',[result])
    print('CONTROL_RECEIPT_READY',mode,case,flush=True)
finally:
    with (out/'imported_modules.tsv').open('w') as stream:
        stream.write('module\tpath\n')
        for name,module in sorted(sys.modules.items()):stream.write(name+'\t'+str(getattr(module,'__file__','builtin'))+'\n')
