#!/usr/bin/env bash
set -euo pipefail
out=${SOBOL_OUT_DIR:?}
sif=${SOBOL_SIF:?}
snapshot="$out/slurm/code_snapshot"
code="$snapshot/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol"
scratch=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/sobol_full_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
trap 'exit 143' TERM INT
mkdir "$scratch/model_cache" "$scratch/perron_cache"
expected=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["archive_sha256"])' "$out/slurm/runtime_cache_manifest.json")
actual=$(sha256sum "$out/slurm/rcpp_template.tar.gz"); actual=${actual%% *}
[[ "$actual" = "$expected" ]] || { echo 'Cache archive checksum differs' >&2; exit 2; }
tar -xzf "$out/slurm/rcpp_template.tar.gz" -C "$scratch/model_cache"
model_cache="$snapshot/oxygen/code/O2_supply_demand_MAP/model/.rcpp_cache_o2_supply_demand_map"
case ${SOBOL_ROLE:?} in
  bootstrap|gate|aggregate|finish) args=("$SOBOL_ROLE") ;;
  validate|summarize) args=("$SOBOL_ROLE" --case-id "${SLURM_ARRAY_TASK_ID:?}") ;;
  evaluate) args=(worker --task-id "${SLURM_ARRAY_TASK_ID:?}") ;;
  response) args=(response --task-id "${SLURM_ARRAY_TASK_ID:?}") ;;
  *) exit 2 ;;
esac
cd "$snapshot"
echo "role=$SOBOL_ROLE host=$(hostname -s) job=$SLURM_JOB_ID task=${SLURM_ARRAY_TASK_ID:-none} start=$(date -Is)"
apptainer exec --cleanenv --bind "$scratch:$scratch,$scratch/model_cache:$model_cache" "$sif" \
  env R_HOME= R_ENVIRON_USER=/dev/null R_LIBS_USER=/opt/R/4.4.2/lib64/R/library \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR="$scratch/matplotlib" \
  MININGCLONEID_RCPP_REBUILD=FALSE SOBOL_PERRON_CACHE="$scratch/perron_cache" \
  SOBOL_SCRATCH="$scratch" SOBOL_MODEL_CACHE="$scratch/model_cache" SLURM_JOB_ID="$SLURM_JOB_ID" \
  python3 "$code/full_run.py" "${args[@]}" --out-dir "$out"
echo "Completed at $(date -Is)"
