#!/usr/bin/env python3
"""Bounded, nonnumerical HPC inspection of candidate packages and source files."""
import csv
import hashlib
import sys
from pathlib import Path

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()

study=Path(sys.argv[1])
root=Path('/share/lab_crd/taoli/Project/soft_couping_org')
fit=root/'oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253'
with (study/'manifests/candidate_package_inventory.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['seed','file','original_path','exists','bytes','sha256'])
    for seed in [25,74,366,292,392]:
        for name in ['deoptim_result.rds','fit_config.rds','parameter_table.csv','parameter_table_input.csv','fit_summary.tsv','run_provenance.tsv','run_effective_args.tsv','burden_fit.tsv','terminal_ploidy_fit.tsv','necrosis_fit.tsv']:
            p=fit/('seed'+str(seed))/name
            w.writerow([seed,name,str(p),'true' if p.is_file() else 'false',p.stat().st_size if p.is_file() else 0,sha(p) if p.is_file() else ''])
with (study/'manifests/numerical_source_inventory.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['file','frozen_path','discovery_path','frozen_sha256','current_discovery_sha256','matches_current','historical_identity_limit'])
    for rel in ['model/model_O2_supply_demand_MAP.R','model/model_O2_supply_demand_MAP.cpp','util/o2_supply_demand_map_shared.R','util/o2_supply_demand_map_common_semantics.R','util/o2_supply_demand_map_fit_invivo_backend.R']:
        p=study/'source/numerical'/rel; q=root/'oxygen/code/O2_supply_demand_MAP'/rel
        a=sha(p);b=sha(q)
        w.writerow([rel,str(p),str(q),a,b,str(a==b).lower(),'Full historical dirty-worktree hashes were not archived by the fit; core MD5 records and tested replay are available.'])
print('Inventoried 50 candidate artifacts and five frozen numerical source files; no model evaluation.')
