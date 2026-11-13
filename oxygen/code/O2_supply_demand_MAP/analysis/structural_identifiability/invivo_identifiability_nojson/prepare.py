#!/usr/bin/env python3
"""Prepare the separately authorized tumor pilot without evaluating the model."""
import csv
import hashlib
import shutil
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[6]
RUN_ID = '20261111_invivo_identifiability_202133Z'
STUDY = PROJECT / 'oxygen/results/structural_identifiability' / RUN_ID
SOURCE = PROJECT / 'oxygen/results/structural_identifiability/20261109_seed144_wp0/source/code/oxygen/code/O2_supply_demand_MAP'
FIGURE = Path('/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/data/Figures/Figure4')

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def table(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def main():
    if STUDY.exists():
        raise SystemExit('Run already exists; preparation never overwrites a run')
    for d in ['implementation','source/numerical/model','source/numerical/util','source/seed25','source/inputs','source/paper','manifests','tables','tests','fits','logs','runtime','state','reports','figures']:
        (STUDY / d).mkdir(parents=True, exist_ok=True)
    for p in Path(__file__).parent.iterdir():
        if p.is_file(): shutil.copy2(p, STUDY / 'implementation' / p.name)
    for rel in ['model/model_O2_supply_demand_MAP.R','model/model_O2_supply_demand_MAP.cpp','util/o2_supply_demand_map_shared.R','util/o2_supply_demand_map_common_semantics.R','util/o2_supply_demand_map_fit_invivo_backend.R']:
        shutil.copy2(SOURCE / rel, STUDY / 'source/numerical' / rel)
    for src, dst in [
        ('source_fit_config/fit_config.rds','source_config.rds'),
        ('runtime_fit/seed25/fit_config.rds','runtime_config.rds'),
        ('figure4a_runtime_config_provenance.tsv','figure4a_runtime_config_provenance.tsv'),
        ('invivo_parameter_table_seed25.csv','invivo_parameter_table_seed25.csv'),
        ('invivo_fit_objective_ranking_500seeds.tsv','invivo_fit_objective_ranking_500seeds.tsv')]:
        shutil.copy2(FIGURE / src, STUDY / 'source/paper' / dst)
    rows = [dict(role=r, attempts=n, cores=c, wall_seconds=s, cap_core_hours=h) for r,n,c,s,h in [
        ('replay',3,1,7200,6),('validate_parallel',1,8,7200,16),('validate_serial',2,1,14400,8),
        ('global',2,8,28800,128),('local',2,1,21600,12),('check',10,1,3600,10),('reserve',2,1,21600,12)]]
    table(STUDY/'manifests/budget.tsv', rows, ['role','attempts','cores','wall_seconds','cap_core_hours'])
    settings = dict(run_id=RUN_ID, authorized_core_hours=192, partition='red',qos='xxlarge',requeue=False,
        score_tolerance=1e-8,prediction_absolute_tolerance=1e-8,prediction_relative_tolerance=1e-8,
        fixed_log10_density=5.198970004336019, global_NP=256,global_itermax=1000,
        global_reltol=1e-4,global_steptol=25,local_maxit=15000,
        source_commit_reference='1a5d861c9c4d2af72ea0d07def0023831cc531e4',
        image_sha256='dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372',
        source_identity_limit='Core hashes are archived; full historical dirty backend identity requires audit. Endpoint replay establishes tested compatibility, not every historical byte.',
        authorized_scope='Steps 1-5: source/data audit, production replay, adapter validation, two constrained timing chains. No unrestricted baseline search or bulk profile.',
        authorization='User confirmed independent 192 core-hours on 2026-11-11; culture tasks and budget remain unchanged.')
    table(STUDY/'manifests/settings.tsv',[dict(key=k,type='boolean' if isinstance(v,bool) else 'integer' if isinstance(v,int) else 'double' if isinstance(v,float) else 'string',value=str(v).lower() if isinstance(v,bool) else format(v,'.17g') if isinstance(v,float) else str(v)) for k,v in settings.items()],['key','type','value'])
    print(STUDY)

if __name__ == '__main__': main()
