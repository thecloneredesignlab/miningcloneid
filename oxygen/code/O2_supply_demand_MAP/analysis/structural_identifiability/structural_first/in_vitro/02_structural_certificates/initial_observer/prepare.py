"""Freeze provenance and run only standard-library/static preparation tests."""
import ast
from pathlib import Path
import shutil
import subprocess
import sys
from protocol import (MECHANISMS,NUISANCES,PARENTS,CAP,sha,read,table,verify_archive,
                      validate_roles,validate_schema,validate_collection,allocation_seconds)

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();code=Path(__file__).parent
assert not study.exists(),'Fresh evidence namespace required'
base=project/'oxygen/results/structural_identifiability'
retry_parent=Path(sys.argv[3]).resolve() if len(sys.argv)>3 else None
technical_retries_used=0
checks=[]
def check(name,condition=True):
    assert condition,name;checks.append(dict(check=name,passed=True))
def reject(name,callback):
    try:callback()
    except ValueError:check(name)
    else:raise AssertionError('Invalid metadata accepted: '+name)

parent_rows=[]
for name,digest in PARENTS.items():
    verified=verify_archive(base/name,digest)
    parent_rows.append(dict(namespace=name,manifest_sha256=digest,verified_files=len(verified)))
    check('immutable_parent_'+name)
old=base/'20261116_culture_step04_public_tools_01'
frozen=read(old/'reports/frozen_source_verification.tsv')
for row in frozen:check('production_source:'+row['path'],sha(project/row['path'])==row['sha256'])
if retry_parent:
    verified=verify_archive(retry_parent,'cb3cf649d5293427567aee9449ac10a00ac5e9bfb5be1bb24f2f3249b675e5d0')
    failure={r['key']:r['value'] for r in read(retry_parent/'preparation_status.tsv')}
    assert failure['status']=='FAILED_PREPARATION' and failure['HPC_jobs']=='0' and failure['technical_retry_assigned']=='1'
    technical_retries_used=1
    parent_rows.append(dict(namespace=retry_parent.name,manifest_sha256=sha(retry_parent/'collection_manifest.tsv'),verified_files=len(verified)))
    check('failed_preparation_preserved')
for name in ('inputs/production','reports','accounting','runtime','verification','logs','releases/r01'):(study/name).mkdir(parents=True)
table(study/'inputs/parent_verification.tsv',parent_rows)
table(study/'inputs/frozen_source_verification.tsv',frozen)
if retry_parent:shutil.copy2(retry_parent/'preparation_status.tsv',study/'inputs/failed_preparation_status.tsv')
for p in ('parameter_roles.tsv','model_tool_map.tsv','project_case_readiness.tsv'):
    shutil.copy2(old/'reports'/p,study/'inputs'/p)
roles=read(study/'inputs/parameter_roles.tsv');validate_roles(roles)
check('twenty_original_unknown_roles')
sources=[]
for relative in sorted(set(r['source_file'] for r in roles)):
    dest=study/'inputs/production'/Path(relative).name
    shutil.copy2(project/relative,dest)
    sources.append(dict(original_path=relative,frozen_copy=str(dest.relative_to(study)),sha256=sha(dest)))
table(study/'inputs/production_sources.tsv',sources)
for name in ('culture_structural_step04_initial_observer_validation_plan_20261116.md',
             'culture_analysis_no_discretionary_time_limits_20261116.md',
             'meeting_aligned_synthetic_identifiability_master_plan_20261112.md'):
    shutil.copy2(project/Path({'culture_analysis_no_discretionary_time_limits_20261116.md': 'oxygen/results/analysis_documentation/in_vitro/design_and_contracts/culture_analysis_no_discretionary_time_limits_20261116.md', 'culture_structural_step04_initial_observer_validation_plan_20261116.md': 'oxygen/results/analysis_documentation/in_vitro/structural_certificates/culture_structural_step04_initial_observer_validation_plan_20261116.md', 'meeting_aligned_synthetic_identifiability_master_plan_20261112.md': 'oxygen/results/analysis_documentation/frozen_inputs/in_vitro/design_and_contracts/meeting_aligned_synthetic_identifiability_master_plan_20261112.md'}[name]),study/'inputs'/name)
cert=base/'20261116_culture_step04_C04_certificate_01'
for name in ('io_guard.c','io_guard.so'):shutil.copy2(cert/'runtime'/name,study/'runtime'/name)
shutil.copy2(cert/'reports/control_status.tsv',study/'inputs/inherited_control_status.tsv')
table(study/'inputs/authorization.tsv',[dict(key=k,value=v) for k,v in dict(
 scope='Confirmed Step 4 Stage A remaining initial/observer controls only',
 confirmation_date='2026-11-16 America/New_York',cap_core_seconds=CAP,previous_charged_to_this_package=0,
 technical_retries_max=1,technical_retries_used=technical_retries_used,qos='xxlarge',concurrency_cap='none',
 discretionary_elapsed_limit='none',native_driver_calls=0,random_rank_calls=0,
 language='English',Java_JavaScript_JSON='forbidden',Stage_B=False,Step_5=False,
 fits=False,synthetic_culture_data=False,author_contact=False,old_controllers='unchanged',
 production_SIF='reference_only; no new execution',native_initial_manifold_support='not_verified').items()])
expectations=[
 ('T00','Twenty roles, source integrity and compliant private symbolic runtime'),
 ('T01','Finite three-bin Gaussian normalization, derivatives, inverse and positive minor'),
 ('T02','State-space jets then constrained initial chain rule; toy rate/shape witness'),
 ('T03','Total-only explicit all-time invariance under initial shape changes'),
 ('T04','Column-normalized unknown observer width; full distribution versus mean'),
 ('T05','Unbounded Gaussian variance-sum preserving family; no finite-grid inference'),
 ('T06','Flow post-mixture density normalization versus karyotype column normalization'),
 ('T07','Four initial coordinates shared within cohort across two constant inputs'),
 ('T08','Variance-only growth SD and exact-zero expression activity'),
 ('T09','Shape, missing nuisance, normalization, sharing and stale-collector rejection'),
 ('T10','Static source mapping only; original culture statuses remain not_analyzed')]
table(study/'inputs/case_manifest.tsv',[dict(case=k,expected=v,scope='toy/interface or static source only',max_jet_order=2,project_solver=False) for k,v in expectations])
(study/'inputs/control_expectations.md').write_text('# Frozen interface expectations\n\n'+ '\n'.join(k+': '+v for k,v in expectations)+'\n\nOptional numeric checks are not used; exact symbolic references are required.\n')
(study/'inputs/interface_specification.md').write_text(
 '# Finite initial and observer interface\n\nConstruct H_(k+1)=D_x H_k f in state variables first. Substitute x0=g(eta), then differentiate with respect to original declared unknowns. Preserve D_x H D_eta g and direct parameter dependence.\n\nKaryotype normalizes each Gaussian column. Flow normalizes the mixture by finite-grid spacing times its sum. sigma_growth is a likelihood SD outside the mean map.\n\nOnly registered toys through jet order two are evaluated. Independent references do not import the adapter. No native rank driver, full culture matrix or discrete-time engine is executed.\n')
required=['mu2','sd2','mu4','sd4','kappa']
mapping=[(c,e,('mu2','sd2') if c=='2N' else ('mu4','sd4')) for c in ('2N','4N') for e in ('E1','E2')]
validate_schema(required,required,mapping,['sigma_kary']*2);check('valid_shared_schema')
reject('missing_nuisance_rejected',lambda:validate_schema(required[:-1],required))
reject('duplicate_unknown_rejected',lambda:validate_schema(required+['mu2'],required))
reject('variance_mean_column_rejected',lambda:validate_schema(['sigma_growth'],['sigma_growth']))
bad=list(mapping);bad[1]=('2N','E2',('new_mu','new_sd'))
reject('duplicated_initial_per_oxygen_rejected',lambda:validate_schema(required,required,bad))
reject('duplicated_observer_width_rejected',lambda:validate_schema(['r1','r2'],['r1','r2'],assay_widths=['r1','r2']))
validate_schema(['r1','r2'],['r1','r2'],assay_widths=['r1','r2'],diagnostic=True);check('independent_width_control_explicitly_diagnostic')
records=[dict(case='toy',release_sha256='r',input_sha256='i',order=str(k),path='matrix_'+str(k)) for k in (1,2)]
validate_collection(records,'toy','r','i');check('multiple_legitimate_saved_matrices')
reject('stale_collector_rejected',lambda:validate_collection(records,'toy','changed','i'))
reject('duplicate_collector_record_rejected',lambda:validate_collection(records+records[:1],'toy','r','i'))
reject('missing_collector_record_rejected',lambda:validate_collection([],'toy','r','i'))
check('primary_reservation_preserves_retry_and_settlement',allocation_seconds(0)==21600)
check('retry_reservation_counts_prior_spend',allocation_seconds(100,retry=True)==25080)
check('site_limit_is_separate_constraint',allocation_seconds(0,site_seconds=7200)==7200)
reject('overbudget_admission_rejected',lambda:allocation_seconds(CAP))
reject('live_commitment_overbudget_rejected',lambda:allocation_seconds(0,live=CAP))
for p in code.glob('*.py'):
    tree=ast.parse(p.read_text());check('AST_syntax_'+p.name)
    if p.name not in ('prepare.py',):
        for node in ast.walk(tree):
            if isinstance(node,(ast.Import,ast.ImportFrom)):
                names=[a.name.split('.')[0] for a in node.names] if isinstance(node,ast.Import) else [(node.module or '').split('.')[0]]
                check('imports_'+p.name+'_'+str(node.lineno),not set(names)&{'json','jpype','jnius','symbtools','strike_goldd','numpy','ipydex'})
            if isinstance(node,ast.Call):
                attr=node.func.attr if isinstance(node.func,ast.Attribute) else node.func.id if isinstance(node.func,ast.Name) else ''
                check('bounded_call_'+p.name+'_'+str(node.lineno),attr not in ('generic_rank','strike_goldd','alarm','setitimer','sleep') and not any(k.arg=='timeout' for k in node.keywords))
ref_tree=ast.parse((code/'reference.py').read_text())
check('reference_does_not_import_adapter',not any(isinstance(n,(ast.Import,ast.ImportFrom)) for n in ast.walk(ref_tree)))
job=(code/'job.sh').read_text()
check('no_case_timeout_in_job','/usr/bin/timeout' not in job and '\ntimeout ' not in job and 'maxLietime=' not in job)
subprocess.run(['bash','-n',str(code/'job.sh')],check=True);check('shell_syntax')
check('source_functions_frozen',all((study/r['frozen_copy']).is_file() for r in sources))
tracked=subprocess.check_output(['git','diff','--name-only','-z'],cwd=project).decode().split('\0')
untracked=subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=project).decode().split('\0')
owned=str(code.relative_to(project))+'/'
protected=[p for p in sorted(set(tracked+untracked)) if p and not p.startswith(owned) and (project/p).is_file()]
table(study/'inputs/protected_worktree_hashes.tsv',[dict(path=p,sha256=sha(project/p)) for p in protected])
table(study/'reports/static_tests.tsv',checks)
table(study/'accounting/jobs.tsv',[],['job_id','attempt','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
for p in code.iterdir():
    if p.is_file() and p.suffix in ('.py','.sh'):shutil.copy2(p,study/'releases/r01'/p.name)
manifest=[]
for name in ('inputs','runtime','releases'):
    for p in sorted((study/name).rglob('*')):
        if p.is_file():manifest.append(dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)))
table(study/'source_manifest.tsv',manifest)
print('STATIC_PREPARATION_PASSED',len(checks),'checks;',sum(r['verified_files'] for r in parent_rows),'parent files;',len(frozen),'frozen source/input hashes;',len(protected),'protected worktree files')
