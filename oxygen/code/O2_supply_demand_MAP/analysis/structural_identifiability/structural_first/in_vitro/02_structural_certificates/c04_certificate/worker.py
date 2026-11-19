"""Exact certificate for the archived C04 toy matrix; no native solve or timers."""
import os
from pathlib import Path
import sys
from protocol import read,table,sha,parse_matrix

study=Path(sys.argv[1]);out=study/'verification';out.mkdir(exist_ok=True)
assert not (out/'certificate.tsv').exists()
for row in read(study/'source_manifest.tsv'):assert sha(study/row['path'])==row['sha256']
def guard(event,args):
    root=args[0].split('.')[0] if event=='import' else ''
    if root=='ipydex':raise ImportError('Interactive debugger excluded')
    bad=root in ('json','_json','simplejson','orjson','ujson','rapidjson','jpype','jnius','execjs','py_mini_racer','symbtools','strike_goldd','native_strike','functions')
    if event=='exec':
        name=str(getattr(args[0],'co_filename',''));bad=bad or name.endswith('strike_goldd.py') or '/symbtools/' in name
    if bad:
        with (out/'blocked_imports.tsv').open('a') as f:f.write('BLOCKED\t'+event+'\t'+str(args[0])+'\n')
        raise ImportError('Forbidden certificate dependency')
sys.addaudithook(guard)
sys.path[:0]=['/home/4482173/.local/lib/python3.11/site-packages',*[p for p in os.environ.get('PYTHONPATH','').split(':') if p]]
try:
    import sympy as sp
    assert sp.__version__=='1.12'
    x,a,b=sp.symbols('x1 a b',positive=True)
    matrix_path=study/'inputs/obs_ident_matrix_case_model_2_Lie_deriv.txt'
    actual=parse_matrix(matrix_path.read_text(),sp,dict(x1=x,a=a,b=b))
    assert actual.shape==(3,3)
    jets=sp.Matrix([x,-a*x**b,a*a*b*x**(2*b-1)])
    expected=jets.jacobian(sp.Matrix([x,a,b]))
    rows=[]
    for i in range(3):
        for j in range(3):
            residual=sp.simplify(actual[i,j]-expected[i,j]);assert residual==0
            rows.append(dict(row=i,column=j,residual=str(residual)))
    table(out/'row_checks.tsv',rows)
    determinant=sp.simplify(actual.det());reference=-a*a*x**(3*b-1)
    residual=sp.simplify(determinant-reference);assert residual==0
    assert sp.simplify(x*jets[2]/jets[1]**2-b)==0
    assert sp.simplify(-jets[1]/x**b-a)==0
    table(out/'rank_certificate.tsv',[dict(augmented_columns='x1,a,b',determinant=str(determinant),reference=str(reference),residual=str(residual),nonzero_domain='x1>0 and a>0; registered b>0 remains unknown real',augmented_rank=3,theorem_scope='analytic observability-identifiability rank condition; local only')])
    table(out/'recovery_identities.tsv',[dict(parameter='b',identity='x1*x_ddot/x_dot^2',residual=0),dict(parameter='a',identity='-x_dot/x1^b',residual=0)])
    table(out/'certificate.tsv',[dict(case='C04',status='PASS_EXACT_SAVED_MATRIX_CERTIFICATE',rank=3,individual_parameters_locally_identifiable='a,b',unknown_initial_state_observable=True,source_matrix_sha256=sha(matrix_path),new_native_calls=0,new_randomized_rank_calls=0,native_backend_status='original_timeout_and_incomplete_return_preserved',scope='toy ODE only; not the original culture model')])
    print('EXACT_SAVED_C04_CERTIFICATE_PASSED',flush=True)
finally:
    with (out/'imported_modules.tsv').open('w') as f:
        f.write('module\tpath\n')
        for name,module in sorted(sys.modules.items()):f.write(name+'\t'+str(getattr(module,'__file__','builtin'))+'\n')
