"""Freeze saved evidence and run local standard-library collector tests only."""
from pathlib import Path
import shutil
import subprocess
import sys
from base_contract import read,table,sha
from tests import run

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve()
old=project/'oxygen/results/structural_identifiability/20261116_culture_step04_four_controls_01'
assert not study.exists()
assert (old/'COMPLETE').read_text().splitlines()[0]=='0ff28dac7b66fbec51883a52fc20f8944d7b860e58b0882d4880797da1e61e06'
assert sha(old/'collection_manifest.tsv')==(old/'COMPLETE').read_text().splitlines()[0]
for row in read(old/'collection_manifest.tsv'):assert sha(old/row['path'])==row['sha256'],row['path']
for name in ('inputs','reports','accounting','fixtures','saved','upstream','runtime','releases/r01','review','controls','validation_retry','logs'):(study/name).mkdir(parents=True)
shutil.copy2(project/'oxygen/results/analysis_documentation/in_vitro/structural_certificates/culture_structural_step04_saved_C03_review_C04_C05_plan_20261116.md',study/'inputs/approved_plan.md')
shutil.copy2(old/'collection_manifest.tsv',study/'inputs/parent_collection_manifest.tsv')
shutil.copy2(old/'inputs/parameter_roles.tsv',study/'inputs/parameter_roles.tsv')
shutil.copy2(old/'inputs/parent_collection_manifest.tsv',study/'inputs/ancestor_manifest.tsv')
shutil.copytree(old/'controls/C03',study/'saved/C03')
shutil.copytree(old/'controls/C02/native_app',study/'saved/C02_layout')
shutil.copytree(old/'upstream',study/'upstream',dirs_exist_ok=True)
for case in ('C04','C05'):shutil.copytree(old/'fixtures'/case,study/'fixtures'/case)
shutil.copytree(old/'runtime/private_prefix',study/'runtime/private_prefix')
for name in ('io_guard.so','io_guard.c'):shutil.copy2(old/'runtime'/name,study/'runtime'/name)
shutil.copy2(old/'inputs/native_packing_span_validated.py',study/'inputs/validated_packing_span.py')
roles=read(study/'inputs/parameter_roles.tsv');assert len(roles)==20 and sum(r['scientific_role']=='mechanism' for r in roles)==14
table(study/'inputs/authorization.tsv',[dict(key=k,value=v) for k,v in dict(scope='saved C03 verification; C04/C05 native only after acceptance',date='2026-11-16 America/New_York',cap_core_hours=2,technical_retries=1,new_C03_native_calls=0,max_native_calls_per_C04_C05=1,qos='xxlarge',concurrency_cap='none',Stage_B=False,Step_5=False,synthetic_recovery_profiles_tumor=False,Java_JavaScript_JSON='forbidden',author_contact=False).items()])
table(study/'inputs/control_registry.tsv',[dict(case=c,seed=20261116+int(c[1:]),input_sha256=sha(study/'fixtures'/c/'model.py'),maximum_native_calls=1) for c in ('C04','C05')])
protected=subprocess.check_output(['git','diff','--name-only'],cwd=project,text=True).splitlines()
table(study/'inputs/protected_worktree_hashes.tsv',[dict(path=p,sha256=sha(project/p)) for p in protected])
table(study/'accounting/jobs.tsv',[],['task','job_id','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds','retry'])
run(study)
for p in Path(__file__).parent.iterdir():
    if p.is_file() and p.suffix in ('.py','.sh'):shutil.copy2(p,study/'releases/r01'/p.name)
manifest=[]
for folder in ('inputs','fixtures','saved','upstream','runtime','releases'):
    for p in sorted((study/folder).rglob('*')):
        if p.is_file():manifest.append(dict(path=str(p.relative_to(study)),sha256=sha(p),bytes=p.stat().st_size))
table(study/'source_manifest.tsv',manifest)
print('S0/S1 ready:',len(manifest),'frozen files, 14 mechanisms and 6 nuisances; no scientific package imported locally.')
