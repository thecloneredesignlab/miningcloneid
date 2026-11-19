"""Metadata/AST preparation only; scientific operations are confined to HPC."""
import ast
from pathlib import Path
import shutil
import subprocess
import sys
from protocol import read,table,sha,matrix_tree

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve()
old=project/'oxygen/results/structural_identifiability/20261116_culture_step04_saved_review_01'
assert not study.exists()
assert (old/'COMPLETE').read_text().splitlines()[0]=='222cfface062a37a3a8695a8f51fc2697f34475ce266b8ca05dc4ef4ca03896d'
assert sha(old/'collection_manifest.tsv')==(old/'COMPLETE').read_text().splitlines()[0]
for row in read(old/'collection_manifest.tsv'):assert sha(old/row['path'])==row['sha256']
for name in ('inputs','reports','accounting','runtime','verification','logs','releases/r01'):(study/name).mkdir(parents=True)
for name in ('obs_ident_matrix_case_model_2_Lie_deriv.txt',):shutil.copy2(old/'controls/C04/native_app/results'/name,study/'inputs'/name)
shutil.copy2(old/'fixtures/C04/model.py',study/'inputs/control_model.py')
shutil.copy2(old/'fixtures/C04/options.py',study/'inputs/native_options_historical.py')
shutil.copy2(old/'controls/C04/native_app/strike_goldd.py',study/'inputs/native_producer_source.txt')
shutil.copy2(old/'upstream/StrikePy/LICENSE',study/'inputs/StrikePy_LICENSE')
shutil.copy2(old/'controls/C04/application_status.tsv',study/'inputs/original_application_status.tsv')
shutil.copy2(old/'controls/C04/call_started.tsv',study/'inputs/original_native_call_marker.tsv')
shutil.copy2(old/'inputs/parameter_roles.tsv',study/'inputs/parameter_roles.tsv')
shutil.copy2(old/'collection_manifest.tsv',study/'inputs/parent_collection_manifest.tsv')
shutil.copy2(old/'accounting/budget_summary.tsv',study/'inputs/previous_budget.tsv')
shutil.copy2(project/'oxygen/results/analysis_documentation/in_vitro/design_and_contracts/culture_analysis_no_discretionary_time_limits_20261116.md',study/'inputs/no_discretionary_time_limits.md')
for name in ('io_guard.so','io_guard.c'):shutil.copy2(old/'runtime'/name,study/'runtime'/name)
budget={r['key']:r['value'] for r in read(study/'inputs/previous_budget.tsv')}
assert int(budget['cap_core_hours'])==2 and int(budget['conservative_actual_core_seconds'])==313
remaining=7200-313;minutes=remaining//60
table(study/'inputs/authorization.tsv',[dict(key=k,value=v) for k,v in dict(scope='Reprocess timed-out C04; exact certificate accepted by user; retain native incomplete record',confirmation_date='2026-11-16 America/New_York',new_native_calls=0,new_randomized_rank_calls=0,cap_total_core_seconds=7200,previous_charged_core_seconds=313,remaining_core_seconds=remaining,reservation_core_seconds=minutes*60,budget_derived_allocation_minutes=minutes,discretionary_case_timeout='none',native_derivative_timer='not invoked; historical options are data only',qos='xxlarge',concurrency_cap='none',technical_retry_remaining=1,Stage_B=False,Step_5=False,synthetic_recovery_profiles_tumor=False,Java_JavaScript_JSON='forbidden',author_contact=False).items()])
checks=[]
def good(name,condition=True):
    assert condition,name;checks.append(dict(check=name,passed=True))
tree=ast.parse((study/'inputs/control_model.py').read_text())
assign={n.targets[0].id:n.value for n in tree.body if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name)}
good('unknown_real_exponent_retained',ast.dump(assign['f'])==ast.dump(ast.parse('[-a*x1**b]',mode='eval').body))
good('parameter_columns_a_b',ast.dump(assign['p'])==ast.dump(ast.parse('[[a],[b]]',mode='eval').body))
good('state_output_x1',ast.dump(assign['x'])==ast.dump(assign['h'])==ast.dump(ast.parse('[x1]',mode='eval').body))
declaration=next(n.value for n in tree.body if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Tuple))
good('positive_domain_frozen',any(k.arg=='positive' and ast.literal_eval(k.value) is True for k in declaration.keywords))
good('saved_matrix_valid_AST',len(matrix_tree((study/'inputs/obs_ident_matrix_case_model_2_Lie_deriv.txt').read_text()).elts)==3)
producer=(study/'inputs/native_producer_source.txt').read_text()
good('native_state_parameter_column_order_recorded','xaug = np.append(xaug, model.p, axis=0)' in producer and 'jacobiano = sym.Matrix(model.h).jacobian(xaug)' in producer)
good('native_constant_parameter_dynamics_recorded','faug = np.append(faug, zeros(len(model.p), 1), axis=0)' in producer)
for label,text in [('truncated','onx = [[1,0,0]'),('wrong_columns','onx = [[1,0],[1,0],[1,0]]'),('unknown_symbol','onx = [[c,0,0],[0,1,0],[0,0,1]]'),('code_execution',"onx = [[__import__('os'),0,0],[0,1,0],[0,0,1]]")]:
    try:matrix_tree(text)
    except (ValueError,SyntaxError):good('reject_'+label)
    else:raise AssertionError(label)
roles=read(study/'inputs/parameter_roles.tsv')
good('all_twenty_roles',len(roles)==20 and sum(r['scientific_role']=='mechanism' for r in roles)==14)
good('native_never_reclassified_as_completed',not (old/'controls/C04/native_completion.tsv').exists())
for p in Path(__file__).parent.glob('*.py'):ast.parse(p.read_text());good('syntax_'+p.name)
worker_tree=ast.parse((Path(__file__).parent/'worker.py').read_text())
good('no_solver_or_random_rank_call',not any(isinstance(n,ast.Call) and ((isinstance(n.func,ast.Name) and n.func.id in ('strike_goldd','generic_rank')) or (isinstance(n.func,ast.Attribute) and n.func.attr in ('strike_goldd','generic_rank'))) for n in ast.walk(worker_tree)))
table(study/'reports/static_tests.tsv',checks)
tracked=subprocess.check_output(['git','diff','--name-only'],cwd=project,text=True).splitlines()
table(study/'inputs/protected_worktree_hashes.tsv',[dict(path=p,sha256=sha(project/p)) for p in tracked])
table(study/'accounting/jobs.tsv',[],['job_id','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
for p in Path(__file__).parent.iterdir():
    if p.is_file() and p.suffix in ('.py','.sh'):shutil.copy2(p,study/'releases/r01'/p.name)
manifest=[]
for name in ('inputs','runtime','releases'):
    for p in sorted((study/name).rglob('*')):
        if p.is_file():manifest.append(dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)))
table(study/'source_manifest.tsv',manifest)
print('STATIC_CERTIFICATE_CONTRACT_PASSED',len(checks),'checks; cumulative remaining budget',remaining,'core-seconds; no local scientific imports')
