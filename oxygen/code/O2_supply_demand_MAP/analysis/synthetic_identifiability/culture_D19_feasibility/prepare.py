"""Freeze only the new confirmed feasibility study; do not submit or evaluate."""
from pathlib import Path
import shutil,subprocess,sys
from common import *
p=Path(sys.argv[1]).resolve();s=Path(sys.argv[2]).resolve();assert not s.exists()
old=p/'oxygen/results/synthetic_identifiability/20261116_culture_growth_karyotype_recovery_01'
assert sha(old/'collection_manifest.tsv')=='16d38d3b5536d627926843f679620fb601d4ff508b1ca7f413601a012d58f170'
for r in read(old/'collection_manifest.tsv'):assert sha(old/r['path'])==r['sha256']
for d in ('inputs/production','inputs/private','inputs/diagnostic','inputs/public','runtime','release','reports','accounting','verification','logs','controller'):(s/d).mkdir(parents=True)
records=[]
def copy(src,dest,kind):
 q=s/dest;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,q);assert sha(src)==sha(q)
 records.append(dict(source=str(src.relative_to(p)),path=dest,sha256=sha(q),kind=kind))
for f in (old/'inputs/production').iterdir():
 if f.is_file():copy(f,'inputs/production/'+f.name,'unchanged_production')
for f in (old/'runtime').rglob('*'):
 if f.is_file():copy(f,'runtime/'+str(f.relative_to(old/'runtime')),'validated_native_runtime')
for src,dest in [('inputs/public/evaluation_bundle.rds','source_bundle.rds'),('inputs/private/truth_parameters.rds','truth.rds'),('inputs/private/reference_result.rds','reference.rds')]:copy(old/src,'inputs/private/'+dest,'isolated_export_or_positive_control_only')
copy(old/'reports/raw_screen_bank_from_journal.tsv','inputs/diagnostic/raw_bank.tsv','old_exact_unit_vectors_no_response_use')
copy(old/'inputs/parameter_roles.tsv','inputs/parameter_roles.tsv','unchanged20_roles')
copy(p/'oxygen/results/analysis_documentation/in_vitro/synthetic_recovery/culture_D19_feasibility_diagnosis_and_blind_start_plan_20261117.md','inputs/confirmed_plan.md','confirmed32h_scope')
roles_valid(read(s/'inputs/parameter_roles.tsv'));table(s/'inputs/evidence_registry.tsv',records)
table(s/'accounting/authority.tsv',[dict(key=k,value=v) for k,v in dict(max_core_seconds=CAP,retry_max=1,tests_forward_max=128,diagnostic_forward_max=1064,constructors=4,constructor_forward_max=12000,validation_forward_max=2752,forward_total_max=51944,scientific_reserve=0,master_seed=20261118,qos='xxlarge',concurrency_limit='none',application_time_limit='none',D19_fit=0,N20=0,profiles=0,tumor=0).items()])
table(s/'accounting/retry_ledger.tsv',[dict(key='used',value=0),dict(key='maximum',value=1)])
table(s/'accounting/jobs.tsv',[],['job_id','task','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
tracked=subprocess.check_output(['git','diff','--name-only','-z'],cwd=p).decode().split('\0');other=subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=p).decode().split('\0');owned=str(Path(__file__).resolve().parent.relative_to(p))+'/'
table(s/'preservation_baseline.tsv',[dict(path=f,sha256=sha(p/f)) for f in sorted(set(tracked+other)) if f and not f.startswith(owned) and (p/f).is_file()])
(s/'reports/authority.md').write_text('''# Confirmed independent feasibility study

Independent32 allocated core-hours; one technical retry; no scientific reserve. Explain1000 failed old draws, then construct four independent fully feasible NP200 banks, and validate twice. Retain14 mechanisms and6 nuisance roles. No D19 response fitting, N20, profiles, tumor or new observations. English artifacts; no Java, JavaScript or JSON. xxlarge, no application concurrency cap or discretionary elapsed-time cutoff. Old sources, tasks, controllers, budgets and used retry remain unchanged. A passed bank is not a recovery result.
''')
print('PREPARED',s.name,len(records))
