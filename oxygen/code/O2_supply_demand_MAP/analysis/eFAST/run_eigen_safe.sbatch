#!/usr/bin/env bash
set -euo pipefail
repo=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
out=${EFAST_OUT_DIR:-$repo/oxygen/results/eFAST/neighborhood10pct}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
scratch=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/efast_eigen_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
trap 'exit 143' TERM INT
mkdir "$scratch/cache"
mkdir "$scratch/perron"
tar -xzf "$out/slurm/rcpp_template.tar.gz" -C "$scratch/cache"
model_cache="$repo/oxygen/code/O2_supply_demand_MAP/model/.rcpp_cache_o2_supply_demand_map"
perron_cache="$code/.rcpp_cache_perron_safe"
if [[ -f "$out/slurm/perron_safe_template.tar.gz" ]]; then
  tar -xzf "$out/slurm/perron_safe_template.tar.gz" -C "$scratch/perron"
fi
echo "role=${EFAST_RECOVERY_ROLE:?} host=$(hostname -s) job=$SLURM_JOB_ID start=$(date -Is)"
case "$EFAST_RECOVERY_ROLE" in
  probe) command=(Rscript --vanilla "$code/probe_eigen_failure.R" "$out" "$out/slurm/eigen_recovery/safe_probe" "${EFAST_PROBE_TASKS:-7639,11881,13526,20442}") ;;
  validate) command=(Rscript --vanilla "$code/validate_perron_safe.R" "$out/slurm/eigen_recovery/safe_probe" "$out/slurm/perron_safe_validation.json") ;;
  repair) command=(python3 "$code/recover_eigen_safe.py" repair --out-dir "$out" --task-id "${SLURM_ARRAY_TASK_ID:?}") ;;
  summary) command=(python3 "$code/neighborhood_slurm.py" seed-summary --out-dir "$out" --seed-number "${SLURM_ARRAY_TASK_ID:?}") ;;
  finalize) command=(python3 "$code/neighborhood_slurm.py" finalize --out-dir "$out") ;;
  tests) command=(python3 -m unittest discover -s "$code" -p "test_neighborhood*.py") ;;
  *) echo "Unknown recovery role" >&2; exit 2 ;;
esac
apptainer exec --cleanenv --bind "$scratch/cache:$model_cache" --bind "$scratch/perron:$perron_cache" "$sif" env \
  R_HOME= R_ENVIRON_USER=/dev/null R_LIBS_USER=/opt/R/4.4.2/lib64/R/library \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MININGCLONEID_RCPP_REBUILD=FALSE \
  EFAST_PERRON_CACHE="$perron_cache" SLURM_JOB_ID="$SLURM_JOB_ID" "${command[@]}"
if [[ "$EFAST_RECOVERY_ROLE" == validate ]]; then
  tar -czf "$out/slurm/perron_safe_template.tar.gz" -C "$scratch/perron" .
  sha256sum "$code/perron_high_precision_safe.cpp" "$out/slurm/perron_safe_template.tar.gz" > "$out/slurm/eigen_recovery/perron_safe_template.sha256"
fi
echo "Completed at $(date -Is)"
