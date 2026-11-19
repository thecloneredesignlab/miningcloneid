"""Local metadata, import-order and collector tests; no CAS imports."""
import ast
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
from common import *

def import_order(source):
    tree=ast.parse(source)
    allowed={'pathlib','sys','common','runtime'}
    for node in tree.body:
        if isinstance(node,(ast.Import,ast.ImportFrom)):
            names=[a.name.split('.')[0] for a in node.names] if isinstance(node,ast.Import) else [(node.module or '').split('.')[0]]
            if not set(names)<=allowed:return False
    main=next((n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main'),None)
    if main is None:return False
    init=[n.lineno for n in ast.walk(main) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='boot']
    if len(init)!=1:return False
    imports=[n.lineno for n in ast.walk(main) if isinstance(n,(ast.Import,ast.ImportFrom))]
    return all(n>init[0] for n in imports)

def run(project,study):
    project=Path(project).resolve();study=Path(study).resolve();code=Path(__file__).parent;checks=[]
    def check(name,condition=True):
        assert condition,name;checks.append(dict(check=name,passed=True))
    roles_valid(read(study/'inputs/parameter_roles.tsv'));check('original_14_plus_6_roles')
    original=project/'oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/structural_first/in_vitro/02_structural_certificates/finite_state/boundary_audit.py'
    check('old_import_order_regression_detected',not import_order(original.read_text()))
    for name in ('preflight.py','review.py','failure_fixture.py'):
        p=code/name;check('correct_import_order_'+name,import_order(p.read_text()))
        spec=importlib.util.spec_from_file_location('b6_static_'+p.stem,p);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        check('import_did_not_load_flint_'+name,'flint' not in sys.modules)
        class BootstrapSentinel(Exception):pass
        events=[]
        def stub(*args):events.append('boot_before_mathematics');raise BootstrapSentinel()
        module.boot=stub;oldargv=sys.argv;sys.argv=[str(p),str(study),'P01' if name=='preflight.py' else 'P02']
        try:module.main()
        except BootstrapSentinel:check('actual_entry_stub_bootstrap_'+name,events==['boot_before_mathematics'] and 'flint' not in sys.modules)
        finally:sys.argv=oldargv
    for p in code.glob('*.py'):
        tree=ast.parse(p.read_text());check('syntax_'+p.name)
        for n in ast.walk(tree):
            if isinstance(n,(ast.Import,ast.ImportFrom)):
                names=[a.name.split('.')[0] for a in n.names] if isinstance(n,ast.Import) else [(n.module or '').split('.')[0]]
                check('excluded_imports_'+p.name+'_'+str(n.lineno),not set(names)&{'json','jpype','jnius','symbtools','strike_goldd','ipydex'})
            if isinstance(n,ast.Call):
                attr=n.func.attr if isinstance(n.func,ast.Attribute) else n.func.id if isinstance(n.func,ast.Name) else ''
                check('no_timer_or_model_call_'+p.name+'_'+str(n.lineno),attr not in {'alarm','setitimer','generic_rank','output_rows','taylor','operator','cpp_o2simps_simulate_one','ivt_run_lineage'} and not any(k.arg=='timeout' for k in n.keywords))
    rows=[dict(row=i,column=j,mid_man=1,mid_exp=0,rad_man=0,rad_exp=0) for i in range(2) for j in range(2)]
    matrix_records(rows,2);check('complete_matrix_metadata')
    fixtures={'missing':rows[:-1],'duplicate':rows[:3]+[rows[0]],'out_of_range':[dict(r,row=4) if k==0 else r for k,r in enumerate(rows)],
              'negative_radius':[dict(r,rad_man=-1) if k==0 else r for k,r in enumerate(rows)]}
    for name,fixture in fixtures.items():
        try:matrix_records(fixture,2)
        except ValueError:check('reject_matrix_'+name)
        else:raise AssertionError(name)
    for state,code_,app,policy,gate,expected in [('COMPLETED','0:0','PASS','PASS','PASS',True),
        ('COMPLETED','0:0','MISSING','PASS','PASS',False),('COMPLETED','0:0','PASS','MISSING','PASS',False),
        ('OUT_OF_MEMORY','0:9','PASS','PASS','PASS',False),('FAILED','1:0','PASS','PASS','PASS',False)]:
        check('collector_state_'+str(len(checks)),accepted(state,code_,app,policy,gate)==expected)
    check('budget_headroom',reservation(0,0,43200,0)==9000)
    try:reservation(10000,0,43200,0)
    except ValueError:check('budget_exhaustion_rejected')
    else:raise AssertionError('Budget admission')
    for p in code.glob('*.sh'):subprocess.run(['bash','-n',str(p)],check=True);check('shell_'+p.name)
    for p in code.glob('*.R'):
        subprocess.run([str(project/'scripts/agentRrunner.sh'),'-e','invisible(parse(file='+repr(str(p))+'))'],check=True,stdout=subprocess.DEVNULL);check('R_parse_'+p.name)
    table(study/'reports/static_tests.tsv',checks)
    assert not (study/'release/common.py').exists()
    for p in code.iterdir():
        if p.is_file() and p.suffix in ('.py','.R','.sh','.c'):shutil.copy2(p,study/'release'/p.name)
    table(study/'source_manifest.tsv',[dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p))
        for root in ('inputs','runtime','release') for p in sorted((study/root).rglob('*')) if p.is_file()])
    print('B6_LOCAL_STATIC_GATE_PASS',len(checks))

if __name__=='__main__':run(sys.argv[1],sys.argv[2])
