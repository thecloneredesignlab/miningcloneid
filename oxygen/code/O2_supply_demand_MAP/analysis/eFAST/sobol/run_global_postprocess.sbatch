#!/usr/bin/env bash
set -euo pipefail
source_root=${GLOBAL_SOURCE_ROOT:?}
out=${GLOBAL_OUT_DIR:?}
sif=${GLOBAL_SIF:?}
snapshot="$out/slurm/code_snapshot"
scratch=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/sobol_global_${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
trap 'exit 143' TERM INT
fit=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["fit_root"])' "$source_root/slurm/submission_plan.json")
cd "$snapshot"
echo "global post-processing host=$(hostname -s) job=$SLURM_JOB_ID start=$(date -Is)"
code="$snapshot/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol"
apptainer exec --cleanenv --bind "$scratch:$scratch:rw" "$sif" \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$scratch/matplotlib" \
  python3 -m unittest discover -s "$code" -p test_global.py
apptainer exec --cleanenv --bind "$source_root:$source_root:ro,$out:$out:rw,$fit:$fit:ro,$scratch:$scratch:rw" "$sif" \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
  MPLCONFIGDIR="$scratch/matplotlib" SLURM_JOB_ID="$SLURM_JOB_ID" \
  python3 "$snapshot/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/global_postprocess.py" \
  run --out-dir "$out"
echo "Completed at $(date -Is)"
