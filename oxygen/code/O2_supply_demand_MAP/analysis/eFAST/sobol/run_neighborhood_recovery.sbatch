#!/usr/bin/env bash
set -euo pipefail
out=${SOBOL_OUT_DIR:?}
rec=${SOBOL_RECOVERY_DIR:?}
sif=${SOBOL_SIF:?}
snapshot="$rec/slurm/code_snapshot"
code="$snapshot/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol"
scratch=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/sobol_resume_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
trap 'exit 143' TERM INT

# Operational launch helper is versioned separately from the frozen analysis.
source "$(dirname -- "${BASH_SOURCE[0]}")/retry_container_startup.sh"

echo "recovery role=$SOBOL_ROLE job=$SLURM_JOB_ID task=${SLURM_ARRAY_TASK_ID:-none} host=$(hostname -s) start=$(date -Is)"
getent passwd "$(id -u)" || true
case ${SOBOL_ROLE:?} in
  model|response)
    if [[ "$SOBOL_ROLE" = model ]]; then
      python3 - "$rec/plan.json" "${SLURM_ARRAY_TASK_ID:?}" <<'PY'
import json,sys
p=json.load(open(sys.argv[1]))
if int(sys.argv[2]) not in [r['task_id'] for r in p['model_tasks']]:
    raise SystemExit('Model task outside the authorized missing-block manifest')
PY
      export SOBOL_ROLE=evaluate
    fi
    launch bash "$out/slurm/code_snapshot/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/run_full_array.sbatch"
    ;;
  *)
    binding="$out:$out:ro,$rec:$rec"
    case "$SOBOL_ROLE" in case|aggregate|finish) binding="$out:$out";; esac
    args=("$SOBOL_ROLE" --recovery-dir "$rec")
    if [[ "$SOBOL_ROLE" = prepare ]]; then
      args+=(--root "$out" --revision "${RECOVERY_REVISION:?}")
      launch apptainer exec --cleanenv --bind "$binding,$scratch:$scratch" "$sif" \
        env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
        MPLCONFIGDIR="$scratch/matplotlib" python3 -m unittest discover -s "$code" -p 'test*.py'
    elif [[ "$SOBOL_ROLE" = shard || "$SOBOL_ROLE" = case ]]; then
      args+=(--task-id "${SLURM_ARRAY_TASK_ID:?}")
    fi
    launch apptainer exec --cleanenv --bind "$binding,$scratch:$scratch" "$sif" \
      env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
      MPLCONFIGDIR="$scratch/matplotlib" RECOVERY_PLAN_SHA256="${RECOVERY_PLAN_SHA256:-}" \
      python3 "$code/recover_neighborhood.py" "${args[@]}"
    ;;
esac
echo "Completed at $(date -Is)"
