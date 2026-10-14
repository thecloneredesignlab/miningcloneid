#!/usr/bin/env bash
set -euo pipefail
# All resource requests come from the submitter: xxlarge, 12h, no node or throttle.
repo=${EFAST_REPO_ROOT:?}
out=${EFAST_OUT_DIR:?}
sif=${EFAST_SIF:?}
baseline=${EFAST_BASELINE_DIR:?}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
scratch=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/efast_4097_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
trap 'exit 143' TERM INT
mkdir "$scratch/cache" "$scratch/perron"
tar -xzf "$baseline/slurm/rcpp_template.tar.gz" -C "$scratch/cache"
tar -xzf "$baseline/slurm/perron_safe_template.tar.gz" -C "$scratch/perron"
model_cache="$repo/oxygen/code/O2_supply_demand_MAP/model/.rcpp_cache_o2_supply_demand_map"
perron_cache="$code/.rcpp_cache_perron_safe"
case ${EFAST_ARRAY_ROLE:?} in
  trajectory) args=(worker --task-id "${SLURM_ARRAY_TASK_ID:?}") ;;
  seed_summary) args=(seed-summary --summary-id "${SLURM_ARRAY_TASK_ID:?}") ;;
  finalize) args=(finalize) ;;
  smoke) args=(smoke) ;;
  *) exit 2 ;;
esac
cd "$repo"
echo "role=$EFAST_ARRAY_ROLE host=$(hostname -s) job=$SLURM_JOB_ID array_task=${SLURM_ARRAY_TASK_ID:-none} start=$(date -Is)"
apptainer exec --cleanenv --bind "$scratch/cache:$model_cache" --bind "$scratch/perron:$perron_cache" "$sif" env \
  R_HOME= R_ENVIRON_USER=/dev/null R_LIBS_USER=/opt/R/4.4.2/lib64/R/library \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MININGCLONEID_RCPP_REBUILD=FALSE \
  EFAST_PERRON_CACHE="$perron_cache" SLURM_JOB_ID="$SLURM_JOB_ID" \
  python3 "$code/resolution_4097.py" "${args[@]}" --out-dir "$out"
echo "Completed at $(date -Is)"
