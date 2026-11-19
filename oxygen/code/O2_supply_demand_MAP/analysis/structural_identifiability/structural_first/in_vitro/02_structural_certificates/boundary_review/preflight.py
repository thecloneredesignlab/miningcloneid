"""Small runtime controls only; no project model/certificate reconstruction."""
from pathlib import Path
import sys
from common import table,read,sha,matrix_records
from runtime import boot

def main():
    study=Path(sys.argv[1]);task=sys.argv[2]
    sp,flint,np,out=boot(study,task)
    from algebra import encode,decode,independent_minor_check
    import ctypes
    import os
    import subprocess
    assert task=='P01'
    checks=[]
    def check(name,condition,detail=''):
        checks.append(dict(check=name,passed=bool(condition),detail=detail));table(out/'preflight_checks.tsv',checks)
        assert condition,name
    flint.ctx.prec=256
    M=flint.arb_mat([[1,2],[3,4]])
    ok,bound=independent_minor_check(M);check('exact_nonsingular_Neumann',ok,str(bound))
    ok,bound=independent_minor_check(flint.arb_mat([[1,2],[2,4]]));check('singular_matrix_rejected',not ok)
    x=flint.arb(flint.fmpq(1,3)).exp();check('dyadic_ball_roundtrip',decode(encode(x)).contains(x))
    check('nonfinite_balls_detectable',not flint.arb('nan').is_finite() and not flint.arb('inf').is_finite())
    check('outward_Neumann_bound_is_finite',independent_minor_check(M)[1].is_finite())
    check('private_pinned_flint',Path(flint.__file__).resolve().is_relative_to(study.resolve()/'runtime/private_prefix'))
    check('single_FLINT_thread',flint.ctx.threads==1)
    marker=study/'inputs/guard_probe.tsv';before=sha(marker)
    try:marker.open('w')
    except PermissionError:check('Python_write_scope_guard',True)
    else:raise AssertionError('Python guard did not reject write')
    check('Python_rejection_preserved_input',sha(marker)==before)
    oldlog=os.environ['CULTURE_GUARD_LOG'];expected=out/'expected_native_guard.tsv'
    os.environ['CULTURE_GUARD_LOG']=str(expected)
    libc=ctypes.CDLL(None,use_errno=True);libc.fopen.argtypes=[ctypes.c_char_p,ctypes.c_char_p];libc.fopen.restype=ctypes.c_void_p
    handle=libc.fopen(os.fsencode(marker),b'w')
    os.environ['CULTURE_GUARD_LOG']=oldlog
    check('native_write_scope_guard',not handle and expected.exists() and 'WRITE_OUTSIDE_PHASE' in expected.read_text())
    check('native_rejection_preserved_input',sha(marker)==before)
    fixture=out/'fixtures/bootstrap_failure';fixture.mkdir(parents=True)
    env=dict(os.environ,CULTURE_ALLOWED_WRITE_ROOT=str(fixture),CULTURE_GUARD_LOG=str(fixture/'blocked_access.tsv'))
    with (fixture/'console.log').open('w') as stream:
        result=subprocess.run([sys.executable,'-B','-S',str(Path(__file__).parent/'failure_fixture.py'),str(study),task+'/fixtures/bootstrap_failure'],env=env,stdout=stream,stderr=subprocess.STDOUT)
    check('known_failure_application_exit_preserved',result.returncode==1)
    check('failure_keeps_module_inventory',(fixture/'imported_modules.tsv').exists())
    check('failure_keeps_bootstrap_receipt',(fixture/'bootstrap_receipt.tsv').exists())
    check('expected_missing_dependency_not_policy_dependency',read(fixture/'expected_failure.tsv')[0]['exception']=='ModuleNotFoundError')
    table(out/'runtime_gate.tsv',[dict(status='PASS',checks=len(checks),new_project_cases=0,new_project_determinants=0)])
    table(out/'application_result.tsv',[dict(status='PASS',task=task,phase='runtime_preflight')])
    print('B6_RUNTIME_PREFLIGHT_PASS',len(checks),flush=True)

if __name__=='__main__':main()
