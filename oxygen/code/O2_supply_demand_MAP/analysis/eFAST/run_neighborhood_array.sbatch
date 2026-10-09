#!/usr/bin/env bash
set -euo pipefail
# Resources are passed by submit_neighborhood_slurm.sh. No node or %N limit.
repo=${EFAST_REPO_ROOT:?}
out=${EFAST_OUT_DIR:?}
sif=${EFAST_SIF:?}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
case ${EFAST_ARRAY_ROLE:?} in
  trajectory) args=(worker --task-id "${SLURM_ARRAY_TASK_ID:?}") ;;
  seed_summary) args=(seed-summary --seed-number "${SLURM_ARRAY_TASK_ID:?}") ;;
  finalize) args=(finalize) ;;
  *) exit 2 ;;
esac
cd "$repo"
echo "role=$EFAST_ARRAY_ROLE host=$(hostname -s) job=$SLURM_JOB_ID array_task=${SLURM_ARRAY_TASK_ID:-none} start=$(date -Is)"
apptainer exec --cleanenv "$sif" env R_HOME= R_ENVIRON_USER=/dev/null \
  R_LIBS_USER=/opt/R/4.4.2/lib64/R/library OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 SLURM_JOB_ID="$SLURM_JOB_ID" \
  python3 "$code/neighborhood_slurm.py" "${args[@]}" --out-dir "$out"
echo "Completed at $(date -Is)"
