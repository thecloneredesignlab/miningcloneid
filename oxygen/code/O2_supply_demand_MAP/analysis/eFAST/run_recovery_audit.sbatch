#!/usr/bin/env bash
# Run after successful final aggregation to collect the numerical repair audit.
set -euo pipefail
repo=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
out=${EFAST_OUT_DIR:-$repo/oxygen/results/eFAST/neighborhood10pct}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
echo "host=$(hostname -s) job=${SLURM_JOB_ID:?} start=$(date -Is)"
apptainer exec --cleanenv "$sif" env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 "$code/audit_numerical_recovery.py" --out-dir "$out"
echo "Completed at $(date -Is)"
