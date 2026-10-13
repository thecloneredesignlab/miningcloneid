#!/usr/bin/env bash
set -euo pipefail
repo=${SOBOL_REPO_ROOT:?}
out=${SOBOL_OUT_DIR:?}
sif=${SOBOL_SIF:?}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
scratch=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/sobol_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
trap 'exit 143' TERM INT
mkdir "$scratch/model_cache" "$scratch/perron_cache"
tar -xzf "$out/slurm/rcpp_template.tar.gz" -C "$scratch/model_cache"
model_cache="$repo/oxygen/code/O2_supply_demand_MAP/model/.rcpp_cache_o2_supply_demand_map"
case ${SOBOL_ROLE:?} in
  evaluate) args=(worker --task-id "${SLURM_ARRAY_TASK_ID:?}") ;;
  summarize) args=(summarize --case-id "${SLURM_ARRAY_TASK_ID:?}") ;;
  finalize) args=(finalize) ;;
  response) args=(response-worker --task-id "${SLURM_ARRAY_TASK_ID:?}") ;;
  finish) args=(finish) ;;
  *) exit 2 ;;
esac
cd "$repo"
echo "role=$SOBOL_ROLE host=$(hostname -s) job=$SLURM_JOB_ID task=${SLURM_ARRAY_TASK_ID:-none} start=$(date -Is)"
apptainer exec --cleanenv --bind "$scratch/model_cache:$model_cache" "$sif" \
  env R_HOME= R_ENVIRON_USER=/dev/null R_LIBS_USER=/opt/R/4.4.2/lib64/R/library \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR="$scratch/matplotlib" \
  MININGCLONEID_RCPP_REBUILD=FALSE SOBOL_PERRON_CACHE="$scratch/perron_cache" \
  SLURM_JOB_ID="$SLURM_JOB_ID" python3 "$code/run.py" "${args[@]}" --out-dir "$out"
echo "Completed at $(date -Is)"
