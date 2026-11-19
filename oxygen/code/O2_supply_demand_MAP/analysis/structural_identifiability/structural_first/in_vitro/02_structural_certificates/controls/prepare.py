"""Freeze a new finite rerun and perform standard-library static contract checks."""
import ast
from pathlib import Path
import shutil
import subprocess
import sys
import textwrap
from contract import CASES,OPTIONS,model_text,static_contract,named_set,receipt_valid,read,table,sha

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve()
old=project/'oxygen/results/structural_identifiability/20261116_culture_step04_public_tools_01'
assert not study.exists(),'Never overwrite evidence'
assert (old/'COMPLETE').read_text().splitlines()[0]=='6ea489717a56b25bb1932635b91f8ff8de91a71036bfe0184ea32203088d00ef'
assert sha(old/'collection_manifest.tsv')==(old/'COMPLETE').read_text().splitlines()[0]
for row in read(old/'collection_manifest.tsv'):
    assert sha(old/row['path'])==row['sha256'],row['path']
for name in ('inputs','reports','accounting','upstream','runtime','fixtures','releases/r01','preflight','controls','logs'):
    (study/name).mkdir(parents=True)
shutil.copy2(project/'oxygen/results/analysis_documentation/in_vitro/structural_certificates/culture_structural_step04_four_control_rerun_plan_20261116.md',study/'inputs/approved_plan.md')
shutil.copy2(old/'reports/parameter_roles.tsv',study/'inputs/parameter_roles.tsv')
shutil.copy2(old/'collection_manifest.tsv',study/'inputs/parent_collection_manifest.tsv')
shutil.copy2(old/'controls/published_C2M_source.py',study/'inputs/published_C2M_source.py')
shutil.copytree(old/'sources/StrikePy',study/'upstream/StrikePy')
shutil.copytree(old/'sources/dependencies/symbtools',study/'upstream/symbtools')
shutil.copytree(old/'runtime/private_prefix',study/'runtime/private_prefix')
shutil.copy2(old/'runtime/io_guard.so',study/'runtime/io_guard.so')
shutil.copy2(old/'releases/r01/io_guard.c',study/'runtime/io_guard.c')
roles=read(study/'inputs/parameter_roles.tsv');assert len(roles)==20 and sum(r['scientific_role']=='mechanism' for r in roles)==14
checks=[]
for case in CASES:
    folder=study/'fixtures'/case;folder.mkdir()
    (folder/'model.py').write_text(model_text(case));(folder/'options.py').write_text(OPTIONS)
    checks.append(dict(check=case+':valid_export',passed=static_contract(case,model_text(case))))
    malformed=model_text(case).replace('f=['+__import__('contract').F[case]+']','f=[['+__import__('contract').F[case]+']]')
    bad=[('old_nested_f',malformed,OPTIONS),('extra_symbol',model_text(case).replace('f=[','f=[c+'),OPTIONS),('known_parameter',model_text(case),OPTIONS.replace('prev_ident_pars=[]','prev_ident_pars=[1]')),('wrong_options',model_text(case),OPTIONS.replace('checkObser=1','checkObser=0')),('missing_field',model_text(case).replace('h=[x1]',''),OPTIONS)]
    for label,text,opts in bad:
        try: static_contract(case,text,opts)
        except (ValueError,AssertionError,KeyError,SyntaxError): checks.append(dict(check=case+':reject_'+label,passed=True))
        else: raise AssertionError('Invalid contract admitted '+label)
for cid in ('C01','N01','S1'):
    try:model_text(cid)
    except ValueError:checks.append(dict(check='reject_out_of_scope_'+cid,passed=True))
    else:raise AssertionError(cid)
assert named_set('[b, a]',{'a','b'})==named_set('[a,b]',{'a','b'})
checks.append(dict(check='order_independent_named_results',passed=True))
for text in ('[a,c]','[a,a]',"[__import__('os')]"):
    try:named_set(text,{'a','b'})
    except ValueError:checks.append(dict(check='reject_bad_named_result_'+text,passed=True))
    else:raise AssertionError(text)
assert not receipt_valid(dict(case='C02',input_sha256='old',run_id='old',status='PASS'),'C02','new','new')
assert not receipt_valid({},'C02','new','new')
checks.append(dict(check='reject_stale_and_missing_receipts',passed=True))
core=study/'upstream/StrikePy/StrikePy/strike_goldd.py'
tree=ast.parse(core.read_text());main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='strike_goldd')
packing=next(n for n in main.body if isinstance(n,ast.If) and ast.unparse(n.test)=='len(model.x) == 1')
span=textwrap.dedent(ast.get_source_segment(core.read_text(),packing,padded=True))+'\n'
assert ast.dump(ast.parse(span).body[0],include_attributes=False)==ast.dump(packing,include_attributes=False)
(study/'inputs/native_packing_span.py').write_text(span)
helper=ast.parse((study/'upstream/StrikePy/StrikePy/functions/elim_and_recalc.py').read_text()).body[1]
assert isinstance(helper,ast.FunctionDef)
returns=[n for n in ast.walk(helper) if isinstance(n,ast.Return)]
assert returns and all(isinstance(n.value,ast.Tuple) and len(n.value.elts)==6 for n in returns)
checks.append(dict(check='native_helper_six_return_fields_static',passed=True))
calls=[n for n in ast.walk(main) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='elim_and_recalc']
table(study/'reports/helper_call_audit.tsv',[dict(line=n.lineno,argument_count=len(n.args),status='nominal_signature_compatible' if len(n.args)>=7 else 'unsupported_timeout_fallback_do_not_patch') for n in calls])
table(study/'reports/static_input_checks.tsv',checks)
table(study/'inputs/control_registry.tsv',[dict(case=c,seed=20261116+int(c[1:]),expected_rank=2 if c!='C04' else 3,positive=c in ('C02','C04'),derivative_cap=4,input_sha256=sha(study/'fixtures'/c/'model.py')) for c in CASES])
table(study/'inputs/authorization.tsv',[dict(key=k,value=v) for k,v in dict(scope='C02-C05 finite rerun only',date='2026-11-16 America/New_York',cap_core_hours=3,qos='xxlarge',concurrency_cap='none',technical_retries=1,Stage_B=False,Step_5=False,synthetic_recovery_profiles_tumor=False,author_contact=False,Java_JavaScript_JSON='forbidden').items()])
protected=subprocess.check_output(['git','diff','--name-only'],cwd=project,text=True).splitlines()
table(study/'inputs/protected_worktree_hashes.tsv',[dict(path=p,sha256=sha(project/p)) for p in protected])
table(study/'accounting/jobs.tsv',[],['task','job_id','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds','retry'])
for path in Path(__file__).parent.glob('*'):
    if path.is_file() and path.suffix in ('.py','.sh'):shutil.copy2(path,study/'releases/r01'/path.name)
manifest=[]
for folder in ('upstream','runtime','fixtures','releases/r01','inputs'):
    for p in sorted((study/folder).rglob('*')):
        if p.is_file():manifest.append(dict(path=str(p.relative_to(study)),sha256=sha(p),bytes=p.stat().st_size))
table(study/'source_manifest.tsv',manifest)
print('Preparation passed:',len(checks),'static checks; 20 original roles; closed parent unchanged; no scientific dependency imported locally.')
