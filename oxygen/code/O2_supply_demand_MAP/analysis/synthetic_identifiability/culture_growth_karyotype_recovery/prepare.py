"""Freeze the confirmed two-case recovery study. Standard-library metadata only."""
from pathlib import Path
import shutil
import subprocess
import sys
from common import *
project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();assert not study.exists()
old=project/'oxygen/results/synthetic_identifiability/20261116_culture_growth_karyotype_conditional_01'
assert sha(old/'collection_manifest.tsv')=='d8d9e0407a8260f0510ec8e09f0c1d94978c2c2ffaf9a82999ce1e04e5ddbea0'
for r in read(old/'collection_manifest.tsv'):assert sha(old/r['path'])==r['sha256']
for d in ('inputs/production','inputs/private','inputs/public','runtime','release','reports','accounting','verification','logs','controller'):(study/d).mkdir(parents=True)
records=[]
def copy(p,path,kind):
 q=study/path;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q);assert sha(p)==sha(q)
 records.append(dict(origin=str(p.relative_to(project)),path=path,sha256=sha(p),kind=kind))
for p in (old/'inputs/production').iterdir():
 if p.is_file():copy(p,'inputs/production/'+p.name,'unchanged_production')
for p in (old/'runtime').rglob('*'):
 if p.is_file() and ('trace_bundle' in p.parts or p.name=='guard.so'):copy(p,'runtime/'+str(p.relative_to(old/'runtime')),'validated_native_guard_or_trace_runtime')
for src,target in [('verification/I01_R1/truth_parameters.rds','truth_parameters.rds'),('verification/R01/result.rds','reference_result.rds'),('inputs/private/anchor_bundle.rds','source_anchor_bundle.rds')]:copy(old/src,'inputs/private/'+target,'private_truth_or_validation_reference')
for name in ('evaluation_bundle.rds','synthetic_dataset.rds','public_bounds.tsv','README.md'):copy(old/'public'/name,'inputs/public/'+name,'frozen_noisy_public_input')
copy(old/'inputs/parameter_roles.tsv','inputs/parameter_roles.tsv','original14_plus6_identity')
copy(project/'oxygen/results/analysis_documentation/in_vitro/synthetic_recovery/culture_deterministic_recovery_and_blind_noisy_pilot_plan_20261116.md','inputs/confirmed_plan.md','confirmed_protocol')
roles_valid(read(study/'inputs/parameter_roles.tsv'));table(study/'inputs/evidence_registry.tsv',records)
table(study/'accounting/authorization.tsv',[dict(key=k,value=v) for k,v in dict(cap_core_seconds=CAP,technical_retry_max=1,external_forward_max=512,origin_attempt_max=84000,primary_per_case=4,reserve_rounds_per_case=1,reserve_origins_per_case=2,origin_max=12,case_sequence='D19_then_gated_N20',qos='xxlarge',concurrency_cap='none',discretionary_deadline='none',profiles=0,tumor=0).items()])
table(study/'accounting/retry_ledger.tsv',[dict(key='used',value=0),dict(key='maximum',value=1)])
table(study/'accounting/jobs.tsv',[],['job_id','task','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
table(study/'accounting/reserve_ledger.tsv',[dict(case=c,rounds_used=0,origin_E='not_released',origin_F='not_released') for c in ('D','N')])
(study/'reports/authority.md').write_text('''# Confirmed finite recovery authority

User confirmed D19 deterministic validation followed only on its gate by frozen-data N20 blind fitting. Preserve14 mechanisms and6 nuisances. Each case has4 primary origins and at most one2-origin scientific reserve round. Independent400-core-hour cap; at most one technical retry. xxlarge, no application concurrency cap or discretionary research deadline. English artifacts, no Java/JavaScript/JSON, no profiles or tumor. Existing source, closed generator, old tasks/controllers/budgets and unrelated work remain protected. No noisy observation is regenerated.
''')
(study/'reports/block_implementation.md').write_text('''# Registered finite block implementation

The first logical60-generation DE block is split after one generation to honor the400-global-callback cost gate:1 generation then59. All origins use this identical split. The extra200 population reevaluations use the registered1000-call handoff/restoration allowance; global scientific generation count remains360. Native population-linked blocks are not claimed bitwise equivalent to one uninterrupted DE call. Each candidate attempt is charged before propagation; no per-callback directory scan or truth-state cache is reused. Raw normalized coordinates are carried unchanged.
''')
tracked=subprocess.check_output(['git','diff','--name-only','-z'],cwd=project).decode().split('\0');extra=subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=project).decode().split('\0')
owned=str(Path(__file__).resolve().parent.relative_to(project))+'/'
table(study/'preservation_baseline.tsv',[dict(project_path=p,sha256=sha(project/p)) for p in sorted(set(tracked+extra)) if p and not p.startswith(owned) and (project/p).is_file()])
print('RECOVERY_AUTHORITY_AND_INPUTS_FROZEN',len(records))
