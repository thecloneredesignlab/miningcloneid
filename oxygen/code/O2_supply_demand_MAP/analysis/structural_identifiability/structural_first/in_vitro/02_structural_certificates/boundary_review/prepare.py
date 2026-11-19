"""Create only the authorized fresh B6 package and verified evidence copies."""
from pathlib import Path
import shutil
import subprocess
import sys
from common import *

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();code=Path(__file__).parent
assert not study.exists(),'Fresh execution namespace required'
old=project/'oxygen/results/structural_identifiability'/ORIGIN
assert sha(old/'collection_manifest.tsv')==ORIGIN_SHA==(old/'COMPLETE').read_text().splitlines()[0]
origin_rows=read(old/'collection_manifest.tsv');origin_map={r['path']:r for r in origin_rows}
for r in origin_rows:
    p=old/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==int(r['bytes'])
assert kv(old/'accounting/technical_retry.tsv')['used']=='1'
assert int(read(old/'accounting/budget_final.tsv')[0]['used_core_seconds'])==530
for name in ('inputs/evidence','inputs/production','runtime','verification','accounting','reports','logs','release'):(study/name).mkdir(parents=True)
reused=[]
def copy_origin(rel,dest=None):
    r=origin_map[rel];p=old/rel;q=study/(dest or 'inputs/evidence/'+rel)
    q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
    assert sha(q)==r['sha256'];q.chmod(0o555 if p.stat().st_mode&0o111 else 0o444)
    reused.append(dict(origin_path=rel,new_path=str(q.relative_to(study)),sha256=r['sha256'],origin_manifest_sha256=ORIGIN_SHA))
for rel in ('inputs/parameter_roles.tsv','inputs/witnesses.tsv','inputs/cases.tsv','inputs/fixed_settings.tsv','inputs/frozen_source_verification.tsv',
            'accounting/jobs.tsv','accounting/budget_final.tsv','accounting/technical_retry.tsv','accounting/stop_record.tsv','package_status.tsv','acceptance.tsv',
            'verification/B0/semantic_gate.tsv','verification/B1/validation_gate.tsv','verification/B1/single_input_invariance.tsv',
            'verification/B2_review/operator_gate.tsv','verification/B2_native/kernel_calls.tsv','accounting/operator_acceptance.tsv'):
    copy_origin(rel)
for case in ('B01','B04'):
    for witness in ('W1','W2','W3'):
        prefix='verification/'+case+'_'+witness+'/'
        a=[r for r in read(old/(prefix+'certificate_attempts.tsv')) if r['status']=='CERTIFIED'];assert len(a)==1
        stem=a[0]['candidate']+'_'+a[0]['precision']
        for name in ('case_result.tsv','certificate_attempts.tsv','candidate_registry.tsv','row_registry.tsv','column_registry.tsv','branch_registry.tsv','application_result.tsv','policy_gate.tsv',
                     stem+'_matrix.tsv',stem+'_reference.tsv',stem+'_determinant.tsv'):
            copy_origin(prefix+name)
for r in origin_rows:
    rel=r['path']
    if rel.startswith('inputs/production/'):copy_origin(rel,rel)
    if rel.startswith('runtime/') and (rel.startswith(('runtime/private_prefix/','runtime/trace_bundle/')) or rel.endswith(('io_guard.c','io_guard.so'))):
        copy_origin(rel,rel)
for rel in ('inputs/parameter_roles.tsv','inputs/fixed_settings.tsv'):copy_origin(rel,rel)
copy_origin('inputs/arithmetic_distribution.tsv','inputs/arithmetic_distribution.tsv')
# Include only the additional pure helper's production source, verified against
# the already frozen production-source registry rather than an unfrozen checkout.
rel='oxygen/results/structural_identifiability/20261109_seed144_wp0/source/code/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_invitro_summary_utils.R'
entry=next(r for r in read(old/'inputs/frozen_source_verification.tsv') if r['path']==rel)
assert sha(project/rel)==entry['sha256']
dest=study/'inputs/production'/Path(rel).name;shutil.copy2(project/rel,dest);dest.chmod(0o444)
table(study/'inputs/additional_source.tsv',[dict(project_path=rel,sha256=entry['sha256'],new_path=str(dest.relative_to(study)))])
shutil.copy2(old/'collection_manifest.tsv',study/'inputs/origin_manifest.tsv')
shutil.copy2(project/'oxygen/results/analysis_documentation/in_vitro/structural_certificates/culture_structural_step04_B6_repair_and_closeout_plan_20261116.md',study/'inputs/confirmed_plan.md')
table(study/'inputs/reused_evidence.tsv',reused)
(study/'inputs/guard_probe.tsv').write_text('purpose\nread_only_guard_fixture\n')
(study/'inputs/guard_probe.tsv').chmod(0o444)
roles_valid(read(study/'inputs/parameter_roles.tsv'))
table(study/'accounting/authorization.tsv',[dict(key=k,value=v) for k,v in dict(cap_core_seconds=CAP,technical_retry_max=1,
    qos='xxlarge',cpu_per_worker=1,memory_GiB=16,application_concurrency_cap='none',discretionary_time_limit='none',
    original_mechanisms=14,original_nuisances=6,Step_5=False,fits=0,new_project_cases=0,new_project_determinants=0,new_G_kernels=0,new_trajectories=0,language='English').items()])
table(study/'accounting/retry_ledger.tsv',[dict(key='used',value=0),dict(key='maximum',value=1)])
table(study/'accounting/original_budget_reference.tsv',[dict(origin=ORIGIN,used_core_seconds=530,used_technical_retries=1,status='closed; no unused authority transferred')])
table(study/'accounting/jobs.tsv',[],['job_id','task','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
tracked=subprocess.check_output(['git','diff','--name-only','-z'],cwd=project).decode().split('\0')
others=subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=project).decode().split('\0')
owned=str(code.relative_to(project))+'/'
paths=[p for p in sorted(set(tracked+others)) if p and not p.startswith(owned) and (project/p).is_file()]
table(study/'preservation_baseline.tsv',[dict(project_path=p,sha256=sha(project/p)) for p in paths])
(study/'inputs/proof_scope.md').write_text('''# Confirmed finite B6 scope

Six saved certificates;12 saved matrices; fixed256-bit Neumann rereview.
Three exact boundary composition arguments and a source-based constant-O
real-arithmetic discrete/passage induction. Nine pure interface fixtures.
No project operator/jet/determinant regeneration; no full structural cases,
trajectories, generated culture observations, fits, profiles or Step5 release.
The closed original failure, stop record, used retry and530-core-second spend
remain unchanged. The new independent cap is14,400 allocated core-seconds.
''')
print('FRESH_B6_EVIDENCE_PREPARED',len(origin_rows),'origin files verified;',len(reused),'copied records;',len(paths),'protected files')
