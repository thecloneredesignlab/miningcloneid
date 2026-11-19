"""Freeze the confirmed conditional case; metadata only."""
from pathlib import Path
import shutil
import subprocess
import sys
from common import *
project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();assert not study.exists()
old=project/'oxygen/results/structural_identifiability/20261116_culture_step05_bridge_design_01'
assert sha(old/'collection_manifest.tsv')=='48f6549ccd2094b2794b6b579b6caedb12f0a90cf5c89597433fe19b4dcc11cc'
for r in read(old/'collection_manifest.tsv'):assert sha(old/r['path'])==r['sha256']
for d in ('inputs/production','inputs/private/reference','runtime','release','reports','accounting','verification','logs'):(study/d).mkdir(parents=True)
records=[]
def copy(p,q,kind):
 dest=study/q;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
 assert sha(p)==sha(dest);records.append(dict(origin=str(p.relative_to(project)),path=q,sha256=sha(p),kind=kind))
for p in (old/'inputs/production').iterdir():
 if p.is_file():copy(p,'inputs/production/'+p.name,'byte_identical_production')
for p in (old/'runtime').rglob('*'):
 if p.is_file():copy(p,'runtime/'+str(p.relative_to(old/'runtime')),'validated_guard_or_arithmetic_runtime')
for rel,name in [('verification/M01/anchor_bundle.rds','anchor_bundle.rds'),('verification/R01/result.rds','reference/baseline.rds'),('verification/R03/bridge_results.rds','reference/bridge_results.rds'),('verification/O02/observer_perturbations.rds','reference/observer_perturbations.rds'),('inputs/data/parameter_table_input.csv','parameter_table_input.csv')]:copy(old/rel,'inputs/private/'+name,'private_truth_or_saved_diagnostic')
for rel in ['verification/R03/selection_branch_audit.tsv','verification/M01/bridge_points.tsv','verification/M01/observer_points.tsv','verification/O02/clamp_invariance_domain.tsv','verification/R02/anchor_replay_checks.tsv']:
 copy(old/rel,'inputs/private/reference/'+Path(rel).name,'saved_diagnostic_metadata')
copy(old/'inputs/parameter_roles.tsv','inputs/parameter_roles.tsv','original20_identity')
copy(old/'verification/M01/parameter_map20.tsv','inputs/public_bounds.tsv','public_bounds_no_truth_values')
for name in ('structural_review_gate.tsv','case_inference_addendum.tsv','missing_information.tsv','nuisance_configurations.tsv'):
 copy(old/'reports'/name,'inputs/'+name,'retained_scoped_evidence')
copy(project/'oxygen/results/analysis_documentation/in_vitro/synthetic_recovery/culture_growth_karyotype_conditional_synthetic_study_plan_20261116.md','inputs/confirmed_plan.md','user_confirmed_scope')
contract=project/'oxygen/results/synthetic_identifiability/20261113_culture_contract_steps01_03'
for name in ('observation_manifest.tsv','conditioning_graph.md','growth_observation_contract.md'):
 copy(contract/name,'inputs/'+name,'accepted_conditional_contract')
roles_valid(read(study/'inputs/parameter_roles.tsv'))
table(study/'inputs/evidence_registry.tsv',records)
table(study/'accounting/authorization.tsv',[dict(key=k,value=v) for k,v in dict(cap_core_seconds=14400,retry_max=1,forward_attempt_max=4,fixture_pass_max=12,observation_attempt_max=1027,qa_distinct=1024,primary_datasets=1,fits=0,profiles=0,tumor=0,qos='xxlarge',concurrency_cap='none',discretionary_deadline='none').items()])
table(study/'accounting/jobs.tsv',[],['job_id','task','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
table(study/'accounting/retry_ledger.tsv',[dict(key='used',value=0),dict(key='maximum',value=1)])
(study/'reports/conditional_scope_decision.md').write_text('''# Confirmed conditional-generator decision

The user confirmed the English proposal, including its fixed algorithmic count/history inputs, artificial endpoint-Gaussian growth law, provisional audited karyotype assignment, seed144 with boundaries and sigma flatness, and generation-only scope. New conditional design release is limited to validation and one frozen dataset. Broader physical-study HOLD and sparse-discrete identification NOT_PROVED remain. No fit, profile, tumor, new structural solver or real-data optimum search is authorized. Independent cap4 core-hours, at most one technical retry; xxlarge, no application concurrency cap or discretionary deadline.
''')
tracked=subprocess.check_output(['git','diff','--name-only','-z'],cwd=project).decode().split('\0')
untracked=subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=project).decode().split('\0')
owned=str(Path(__file__).resolve().parent.relative_to(project))+'/'
table(study/'preservation_baseline.tsv',[dict(project_path=p,sha256=sha(project/p)) for p in sorted(set(tracked+untracked)) if p and not p.startswith(owned) and (project/p).is_file()])
print('CONDITIONAL_SOURCE_AND_SCOPE_FROZEN',len(records))
