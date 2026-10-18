#!/usr/bin/env bash
set -euo pipefail
repo=/share/lab_crd/taoli/Project/soft_couping_eFAST
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/high_precision_review"
root="$repo/oxygen/results/eFAST/neighborhood10pct/pilot/N8193_M4_and_N4097_M8/phase_spectrum_diagnostics/high_precision_review"
sif=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif
scratch=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/efast_precision_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
trap 'exit 143' TERM INT
mkdir -p "$scratch/model" "$scratch/perron"
tar -xzf "$repo/oxygen/results/eFAST/neighborhood10pct/slurm/rcpp_template.tar.gz" -C "$scratch/model"
if [[ "${REVIEW_ROLE:?}" == main ]]; then
  tar -xzf "$root/compiled_review_template.tar.gz" -C "$scratch/perron"
fi
model_cache="$repo/oxygen/code/O2_supply_demand_MAP/model/.rcpp_cache_o2_supply_demand_map"
run_r() {
 apptainer exec --cleanenv --bind "$scratch/model:$model_cache" --bind "$scratch:$scratch" "$sif" env \
  R_HOME= R_ENVIRON_USER=/dev/null R_LIBS_USER=/opt/R/4.4.2/lib64/R/library \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MININGCLONEID_RCPP_REBUILD=FALSE \
  EFAST_REVIEW_CACHE="$scratch/perron" SLURM_JOB_ID="$SLURM_JOB_ID" Rscript --vanilla "$@"
}
cd "$repo"
echo "role=$REVIEW_ROLE host=$(hostname -s) job=$SLURM_JOB_ID start=$(date -Is)"
case "$REVIEW_ROLE" in
 smoke)
  run_r "$code/fixtures.R" "$root"
  ids=$(python3 -c 'import json,sys;print(" ".join(map(str,json.load(open(sys.argv[1])))))' "$root/smoke_points.json")
  for id in $ids; do run_r "$code/point.R" "$root" "$id"; done
  python3 "$code/review.py" smoke-audit --root "$root"
  tar -czf "$root/compiled_review_template.tar.gz" -C "$scratch/perron" .
  ;;
 main) run_r "$code/point.R" "$root" "${SLURM_ARRAY_TASK_ID:?}" ;;
 *) echo 'Unsupported role' >&2; exit 2 ;;
esac
echo "Completed at $(date -Is)"
