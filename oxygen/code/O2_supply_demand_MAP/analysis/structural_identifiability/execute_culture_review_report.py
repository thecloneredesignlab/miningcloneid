#!/usr/bin/env python3
"""Render accepted culture candidates inside the fingerprinted runtime."""
import json
from pathlib import Path
import subprocess
import sys
import hashlib
from datetime import datetime, timezone

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();snapshot=sys.argv[3]
code=Path(__file__).resolve().parent;image=Path('/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif')
settings=json.loads((study/'study.json').read_text())
source_manifest=json.loads((code/'source_manifest.json').read_text())
for name,expected in source_manifest['sources'].items():
    if hashlib.sha256((code/name).read_bytes()).hexdigest()!=expected:
        raise RuntimeError('Frozen reporting source changed: '+name)
digest=hashlib.sha256()
with image.open('rb')as stream:
    for block in iter(lambda:stream.read(8*1024**2),b''):digest.update(block)
if digest.hexdigest()!=settings['image_sha256']:raise RuntimeError('Image fingerprint changed')
subprocess.check_call(['python3',str(code/'prepare_culture_review_report.py'),str(project),str(study),snapshot])
inputs=study/'report_inputs'/snapshot
domains=json.loads((inputs/'domains.json').read_text())
wp0=project/settings['wp0_relative_path'];parent=project/settings['parent_relative_path']
model=wp0/'source/code/oxygen/code/O2_supply_demand_MAP/model'
for entry in domains:
    runtime=study/'runtime'/('report_'+snapshot)/entry['domain']
    for name in['home','rcpp_cache']:(runtime/name).mkdir(parents=True,exist_ok=True)
    subprocess.check_call(['cp','-a',str(study/'runtime/validation_validate/rcpp_cache')+'/.',str(runtime/'rcpp_cache')+'/'])
    command=['apptainer','exec','--cleanenv','--home',str(runtime/'home'),'--bind','/share:/share',
        '--bind',str(wp0/'source')+':'+str(wp0/'source')+':ro',
        '--bind',str(parent)+':'+str(parent)+':ro',
        '--bind',str(study/'implementation')+':'+str(study/'implementation')+':ro',
        '--bind',str(runtime/'rcpp_cache')+':'+str(model/'.rcpp_cache_o2_supply_demand_map'),
        '--env','OMP_NUM_THREADS=1','--env','OPENBLAS_NUM_THREADS=1','--env','MKL_NUM_THREADS=1',
        '--env','RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript',
        '--env','R_ENVIRON_USER=/dev/null','--env','R_PROFILE_USER=/dev/null',
        '--env','O2PROF_IMAGE_SHA256='+settings['image_sha256'],str(image),'bash',str(project/'scripts/agentRrunner.sh')]
    subprocess.check_call(command+[str(code/'test_culture_offspring_identity.R'),str(code/'culture_profile_diagnostics.R')])
    subprocess.check_call(command+[str(code/'report_culture_review.R'),str(project),entry['input_directory'],entry['artifact_directory']])
subprocess.check_call(['python3',str(code/'aggregate_culture_review_diagnostics.py'),str(project),str(study),snapshot])
subprocess.check_call(['python3',str(code/'audit_culture_review_search.py'),str(study),'--snapshot',snapshot])
subprocess.check_call(['python3',str(code/'summarize_culture_review.py'),str(project),str(study),'--snapshot',snapshot,'--accounting'])
record=dict(pass_render=True,generated_utc=datetime.now(timezone.utc).isoformat(),snapshot=snapshot,
    review_state_sha256=hashlib.sha256((inputs/'review_state_snapshot.json').read_bytes()).hexdigest(),
    domains=domains,image_sha256=settings['image_sha256'],optimization_performed=False,
    reporting_sources=source_manifest['sources'],
    description='Diagnostic replays of accepted candidates; original and expanded domains remain separate')
path=study/'reports'/('render_'+snapshot+'.json');path.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
