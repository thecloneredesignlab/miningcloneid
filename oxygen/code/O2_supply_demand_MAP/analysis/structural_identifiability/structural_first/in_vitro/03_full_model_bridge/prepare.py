"""Freeze Step5 sources/assumptions and evidence; standard library only."""
from pathlib import Path
import shutil
import subprocess
import sys
from common import *

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();code=Path(__file__).parent
assert not study.exists()
results=project/'oxygen/results/structural_identifiability';old=results/ORIGIN
b6=results/'20261116_culture_step04_stageB_B6_closeout_01'
origins=[(old,'COMPLETE',ORIGIN_SHA),(b6,'CLOSED','4a2f843ad0249211acc74690332254abad5c2788da04d7ee0bb69e2dcf9f19d9')]
for archive,marker,expected in origins:
    assert sha(archive/'collection_manifest.tsv')==expected==(archive/marker).read_text().splitlines()[0]
    for r in read(archive/'collection_manifest.tsv'):
        p=archive/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==int(r['bytes'])
for name in ('inputs/production','inputs/data','inputs/evidence','runtime','release','reports','accounting','verification','logs'):(study/name).mkdir(parents=True)
sources=[]
def copy(p,dest,kind):
    p=Path(p);q=study/dest;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
    q.chmod(0o555 if p.stat().st_mode&0o111 else 0o444)
    sources.append(dict(origin=str(p.relative_to(project)),path=dest,sha256=sha(p),kind=kind));assert sha(q)==sha(p)
for r in read(old/'inputs/frozen_source_verification.tsv'):
    p=project/r['path'];assert sha(p)==r['sha256']
    if '/source/code/' in r['path'] and p.suffix in ('.R','.cpp','.h','.hpp'):
        copy(p,'inputs/production/'+p.name,'frozen_production')
wp0=results/'20261109_seed144_wp0'
for rel,name in [('source/seed144/fit_result.rds','fit_result.rds'),('source/seed144/parameter_table_input.csv','parameter_table_input.csv')]:
    copy(wp0/rel,'inputs/data/'+name,'production_anchor')
for p in (wp0/'source/inputs/oxygen').rglob('*'):
    if p.is_file() and p.name in ('fit_data.Rds','jobs_2N.Rds','jobs_4N.Rds','g0g1_ploidy_density_grid.csv'):copy(p,'inputs/data/'+p.name,'production_observations')
helper=results/'20261112_culture_v2_seven_point_review_nojson/implementation/helpers.R'
external=project/'oxygen/results/synthetic_identifiability/20261113_culture_conditional_d0_seed144_02/manifests/external_manifest.tsv'
row=next(r for r in read(external) if r['path']==str(helper.relative_to(project)))
assert sha(helper)==row['sha256'];copy(helper,'inputs/production/parameter_helpers.R','validated_parameter_mapping_not_dispatcher')
d0=project/'oxygen/results/synthetic_identifiability/20261113_culture_conditional_d0_seed144_02_resume01'
for rel in ('truth/truth_parameters.rds','truth/requested_actual_truth.tsv','design/observation_manifest.tsv'):
    row=next(r for r in read(d0/'truth/manifest.tsv') if r['path']==rel);assert sha(d0/rel)==row['sha256']
    copy(d0/rel,'inputs/evidence/d0_'+Path(rel).name,'existing_engineering_reference_not_new_data')
for rel in ('inputs/parameter_roles.tsv','inputs/witnesses.tsv','inputs/fixed_settings.tsv','reports/project_parameter_status.tsv','reports/structural_results.tsv','reports/nomination_diagnostics.tsv'):
    copy(old/rel,'inputs/'+Path(rel).name,'structural_evidence')
for witness in ('W1','W2','W3'):
    src=old/'verification'/('B04_'+witness)
    for p in src.glob('*.tsv'):
        if p.name in ('candidate_registry.tsv','column_registry.tsv','row_registry.tsv','case_result.tsv','certificate_attempts.tsv','branch_registry.tsv') or p.name.startswith('primary_128_'):
            copy(p,'inputs/evidence/B04_'+witness+'/'+p.name,'saved_live_certificate')
for p in (b6/'runtime').rglob('*'):
    if p.is_file() and (p.is_relative_to(b6/'runtime/private_prefix') or p.is_relative_to(b6/'runtime/trace_bundle')):copy(p,'runtime/'+str(p.relative_to(b6/'runtime')),'validated_arithmetic_or_trace_runtime')
copy(b6/'verification/P01/guard_b6.so','runtime/guard.so','validated_read_write_and_dependency_guard')
contract=project/'oxygen/results/synthetic_identifiability/20261113_culture_contract_steps01_03'
for name in ('unresolved_inputs.tsv','modality_availability.tsv','experimental_schedule.tsv','observation_manifest.tsv','growth_observation_contract.md','flow_observation_contract.md','conditioning_graph.md','parameter_activity_rules.tsv'):
    copy(contract/name,'inputs/contracts/'+name,'existing_experiment_audit')
copy(project/'oxygen/results/analysis_documentation/in_vitro/full_model_bridge/culture_step05_structural_to_full_model_and_synthetic_design_plan_20261116.md','inputs/confirmed_plan.md','confirmed_plan')
table(study/'inputs/evidence_registry.tsv',sources)
roles=read(study/'inputs/parameter_roles.tsv');roles_valid(roles)
mapping=read(study/'inputs/project_parameter_status.tsv');assert len(mapping)==120
table(study/'reports/structural_to_full_model_map.tsv',[dict(**r,production_bridge='mechanism/source correspondence; continuous-to-discrete/observation limits retained',full_D_positive_structural_status='not_proved') for r in mapping])
table(study/'reports/case_status_crosswalk.tsv',[dict(case=c,witness=w,historical_direct_status='UNRESOLVED_REGISTERED_SCOPE',new_disposition='output inclusion/restriction only' if c in ('B02','B06') else 'bounded new observer extension; no relabeling') for c in ('B02','B03','B05','B06') for w in ('W1','W2','W3')])
table(study/'reports/mapping_gate.tsv',[dict(status='PASS_SOURCE_MAP_WITH_EXPLICIT_LIMITS',coordinates=20,mechanisms=14,nuisances=6,case_coordinate_rows=120,unresolved_runs_retained=12,C_to_D_positive_claim=False)])
table(study/'accounting/authorization.tsv',[dict(key=k,value=v) for k,v in dict(cap_core_seconds=CAP,technical_retry_max=1,forward_bundle_max=74,contrast_max=6,observer_calls_max=27,observer_perturbation_passes_max=4,qos='xxlarge',cpus_per_worker=1,memory_GiB=16,concurrency_cap='none',discretionary_timeout='none',generated_observations=0,fits=0,tumor=0,language='English').items()])
table(study/'accounting/retry_ledger.tsv',[dict(key='used',value=0),dict(key='maximum',value=1)])
table(study/'accounting/jobs.tsv',[],['job_id','task','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
tracked=subprocess.check_output(['git','diff','--name-only','-z'],cwd=project).decode().split('\0')
others=subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=project).decode().split('\0')
owned=str(code.relative_to(project))+'/'
table(study/'preservation_baseline.tsv',[dict(project_path=p,sha256=sha(project/p)) for p in sorted(set(tracked+others)) if p and not p.startswith(owned) and (project/p).is_file()])
print('STEP5_MAPPING_FROZEN',len(sources),'source/evidence records;120 map rows;12 unresolved runs retained')
