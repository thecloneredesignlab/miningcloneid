"""Metadata-only B0 freeze and allowlisted binary acquisition, without imports."""
import ast
from fractions import Fraction
import hashlib
import io
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from protocol import *

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();code=Path(__file__).parent
assert not study.exists(),'Fresh Stage B namespace required'
previous=project/'oxygen/results/structural_identifiability/20261116_culture_step04_initial_observer_02'
count=verify_archive(previous,'9693b628dc0afa3f7e02a62023588298620f33f4b6ab9a46b01d2802c710f309')
for name in ('inputs/production','runtime/private_prefix','runtime/trace_bundle','reports','accounting','logs','verification','releases/r01'):(study/name).mkdir(parents=True)
checks=[]
def check(name,condition=True):assert condition,name;checks.append(dict(check=name,passed=True))
check('closed_Stage_A_archive',count==69)
frozen=read(previous/'inputs/frozen_source_verification.tsv')
for r in frozen:check('frozen_source:'+r['path'],sha(project/r['path'])==r['sha256'])
table(study/'inputs/frozen_source_verification.tsv',frozen)
roles=read(previous/'reports/parameter_roles.tsv');roles_valid(roles)
for r in roles:
    r['project_structural_result']='not_analyzed';r['structural_result']='not_analyzed';r['Stage_B_released']='False'
    if r['parameter']=='mu_hp':r['source_expression']='mu_hp*h*(N/44)^gamma_mu'
    if r['parameter']=='gamma_mu':
        r['source_expression']='(N/44)^gamma_mu'
        r['handling']='Production source audit indicates ploidy_related; actual SIF dispatch must pass B0'
    if r['parameter']=='sigma_growth':r['source_expression']='Positive growth likelihood SD; variance=sigma_growth^2; outside deterministic mean map'
table(study/'inputs/parameter_roles.tsv',roles)
prod=project/'oxygen/results/structural_identifiability/20261109_seed144_wp0/source/code/oxygen/code/O2_supply_demand_MAP'
for p in [prod/'model/model_O2_supply_demand_MAP.R',prod/'model/model_O2_supply_demand_MAP.cpp',
          prod/'util/o2_supply_demand_map_common_semantics.R',prod/'util/o2_supply_demand_map_shared.R',
          prod/'util/o2_supply_demand_map_invitro_lineage_simulation_utils.R',prod/'util/o2_supply_demand_map_invitro_objective_utils.R']:
    shutil.copy2(p,study/'inputs/production'/p.name)
for name in ('io_guard.so','io_guard.c'):shutil.copy2(previous/'runtime'/name,study/'runtime'/name)
for name in ('culture_structural_step04_stageB_finite_analysis_plan_20261116.md','culture_analysis_no_discretionary_time_limits_20261116.md'):
    shutil.copy2(project/Path({'culture_analysis_no_discretionary_time_limits_20261116.md': 'oxygen/results/analysis_documentation/in_vitro/design_and_contracts/culture_analysis_no_discretionary_time_limits_20261116.md', 'culture_structural_step04_stageB_finite_analysis_plan_20261116.md': 'oxygen/results/analysis_documentation/in_vitro/structural_certificates/culture_structural_step04_stageB_finite_analysis_plan_20261116.md'}[name]),study/'inputs'/name)
fixed=project/'oxygen/results/synthetic_identifiability/20261113_culture_contract_steps01_03/fixed_settings.tsv'
shutil.copy2(fixed,study/'inputs/fixed_settings.tsv')
(study/'inputs/semantic_correction.md').write_text('# Required death-mode correction\n\nThe inherited gamma_mu description used diploid_NULL. The frozen actual cfg and run-parameter inventory both say ploidy_related. New project formulas use mu_hp*h*(N/44)^gamma_mu only after B0 SIF dispatch verification. The old closed Stage A archive remains unchanged; its Gaussian/observer controls did not evaluate the death branch.\n')
table(study/'inputs/authorization.tsv',[dict(key=k,value=v) for k,v in dict(cap_core_seconds=CAP,cap_core_hours=64,technical_retry_max=1,technical_retries_used=0,qos='xxlarge',application_concurrency_cap='none',discretionary_time_limit='none',mechanisms=14,nuisances=6,Step_5=False,fits=False,production_forward_calls=0,kernel_call_max=9,language='English',Java_JavaScript_JSON='forbidden',source_SIF_sha256=SIF_SHA).items()])
points=[]
for wid,q in zip(('W1','W2','W3'),('1/3','1/2','2/3')):
    overrides={'O2_crit':{'W1':'1/3','W2':'1','W3':'2'}[wid],
      'init_sd_2N':{'W1':'2.90','W2':'2.95','W3':'2.975'}[wid],
      'init_sd_4N':{'W1':'4','W2':'4.5','W3':'5'}[wid],
      'sigma_kary':{'W1':'1','W2':'2','W3':'3'}[wid],'sigma_growth':'1/5'}
    for r in roles:
        name=r['parameter'];rule='explicit' if name in overrides else 'geometric' if r['transform']=='log10' else 'log1p' if r['transform']=='log1p' else 'linear'
        points.append(dict(witness=wid,parameter=name,rule=rule,lower=r['lower'],upper=r['upper'],fraction=q,value=overrides.get(name,'')))
table(study/'inputs/witnesses.tsv',points)
table(study/'inputs/cases.tsv',[
 dict(case='B01',oxygens='1/5',output='full',unknown_initial=True,unknown_observer=False,max_order=2,conditional=False),
 dict(case='B02',oxygens='1/5,1,41/2',output='full',unknown_initial=True,unknown_observer=False,max_order=2,conditional=False),
 dict(case='B03',oxygens='1/5,1,41/2',output='kary',unknown_initial=True,unknown_observer=True,max_order=2,conditional=False),
 dict(case='B04',oxygens='1/5,1,41/2',output='live',unknown_initial=True,unknown_observer=False,max_order=4,conditional=False),
 dict(case='B05',oxygens='1/5,1,41/2',output='kary_flow',unknown_initial=True,unknown_observer=True,max_order=2,conditional=True),
 dict(case='B06',oxygens='1/5,1,41/2',output='full',unknown_initial=False,unknown_observer=False,max_order=2,conditional=True)])

# Inspect a pinned compatible wheel in memory; do not execute a package frontend.
url='https://files.pythonhosted.org/packages/db/3f/5f57004f9f43ac6508f874cc97123f11b27b6811b5bb0d320c43f600f0cb/python_flint-0.9.0-cp310-abi3-manylinux2014_x86_64.manylinux_2_17_x86_64.whl'
expected='376b88cacd30612479e839ffdba887599d3f9c8c0e214852bf80bb2b194e4d76'
payload=urllib.request.urlopen(url).read();assert hashlib.sha256(payload).hexdigest()==expected
members=[]
with zipfile.ZipFile(io.BytesIO(payload)) as z:
    for info in z.infolist():
        name=info.filename;p=Path(name)
        safe=not p.is_absolute() and '..' not in p.parts
        forbidden=p.suffix.lower() in ('.json','.json5','.jsonl','.js','.java','.jar','.ipynb')
        selected=safe and not forbidden and not info.is_dir() and (name.startswith(('flint/','python_flint.libs/')) or '.dist-info/' in name and ('license' in name.lower() or p.name in ('METADATA','WHEEL','RECORD')))
        members.append(dict(member=name,selected=selected,forbidden_format=forbidden))
        if selected:
            dest=study/'runtime/private_prefix'/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(z.read(info))
table(study/'inputs/arithmetic_distribution.tsv',[dict(package='python-flint',version='0.9.0',url=url,sha256=expected,archive_saved=False,acquisition='standard-library in-memory allowlist; no package frontend executed')])
table(study/'inputs/arithmetic_archive_members.tsv',members)
check('private_flint_source_present',(study/'runtime/private_prefix/flint/__init__.py').exists())

# Reuse the host's trace binary and its exact dependency set inside the SIF.
remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+study.name
ssh=['ssh','-o','LogLevel=ERROR','-o','BatchMode=yes','4482173@red.moffitt.org']
script='mkdir -p '+shlex.quote(remote+'/runtime/trace_bundle')+'\n'
script+='cp -L /usr/bin/strace '+shlex.quote(remote+'/runtime/trace_bundle/strace')+'\n'
script+='ldd /usr/bin/strace\n'
result=subprocess.run(ssh+['bash -lc '+shlex.quote(script)],capture_output=True,text=True,check=True)
(study/'inputs/host_strace_dependencies.txt').write_text(result.stdout)
import re
libs=sorted(set(re.findall(r'(?:=>\s+)?(/\S+)\s+\(',result.stdout)))
assert libs and '/lib64/ld-linux-x86-64.so.2' in libs
for lib in libs:
    subprocess.run(ssh+['cp -L '+shlex.quote(lib)+' '+shlex.quote(remote+'/runtime/trace_bundle/'+Path(lib).name)],check=True)
subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','-e','ssh -o LogLevel=ERROR -o BatchMode=yes',
 '4482173@red.moffitt.org:'+remote+'/runtime/trace_bundle/',str(study/'runtime/trace_bundle')+'/'],check=True)

roles_valid(roles);check('twenty_roles_preserved')
check('single_and_multi_case_counts',len(read(study/'inputs/cases.tsv'))==6)
check('all_witness_roles',len(points)==60)
check('site_max_from_remaining_budget',reservation(0,0,43200)==43200)
check('existing_live_commitments_counted',reservation(0,180000,43200)==21600)
try:reservation(CAP,0,43200)
except ValueError:check('overbudget_rejected')
else:raise AssertionError('Budget admission failed')
for rows_,task,rel,valid in [([dict(task='x',release_sha256='r',item='m1'),dict(task='x',release_sha256='r',item='m2')],'x','r',True),
 ([dict(task='x',release_sha256='old',item='m')],'x','r',False),([], 'x','r',False)]:
    try:check_payload(rows_,task,rel)
    except ValueError:assert not valid
    else:assert valid
check('collector_valid_missing_stale_contracts')
for p in code.glob('*.py'):
    tree=ast.parse(p.read_text());check('syntax_'+p.name)
    if p.name not in ('prepare.py',):
        for n in ast.walk(tree):
            if isinstance(n,ast.Call):
                attr=n.func.attr if isinstance(n.func,ast.Attribute) else n.func.id if isinstance(n.func,ast.Name) else ''
                check('no_discretionary_timer_'+p.name+'_'+str(n.lineno),attr not in ('alarm','setitimer','generic_rank','strike_goldd') and not any(k.arg=='timeout' for k in n.keywords))
for p in code.glob('*.sh'):subprocess.run(['bash','-n',str(p)],check=True);check('shell_'+p.name)
for p in code.glob('*.R'):
    subprocess.run([str(project/'scripts/agentRrunner.sh'),'-e','invisible(parse(file='+repr(str(p))+'))'],check=True,stdout=subprocess.DEVNULL);check('R_parse_only_'+p.name)
tracked=subprocess.check_output(['git','diff','--name-only','-z'],cwd=project).decode().split('\0')
others=subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=project).decode().split('\0')
owned=str(code.relative_to(project))+'/'
protected=[p for p in sorted(set(tracked+others)) if p and not p.startswith(owned) and (project/p).is_file()]
table(study/'inputs/protected_worktree_hashes.tsv',[dict(path=p,sha256=sha(project/p)) for p in protected])
table(study/'reports/static_B0_tests.tsv',checks)
table(study/'accounting/jobs.tsv',[],['job_id','task','release','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
for p in code.iterdir():
    if p.is_file() and p.suffix in ('.py','.R','.sh'):shutil.copy2(p,study/'releases/r01'/p.name)
manifest=[]
for parent in ('inputs','runtime','releases/r01'):
    for p in sorted((study/parent).rglob('*')):
        if p.is_file():manifest.append(dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)))
table(study/'source_manifest_r01.tsv',manifest)
print('B0_STATIC_GATE_PASSED',len(checks),'checks;',len(frozen),'frozen sources;',len(protected),'protected files; private arithmetic files',sum(r['selected'] for r in members))
