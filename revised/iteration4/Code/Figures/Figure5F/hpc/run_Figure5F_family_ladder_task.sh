#!/usr/bin/env bash
set -euo pipefail

iteration_root="${FIGURE5F_ITERATION_ROOT:?FIGURE5F_ITERATION_ROOT is required}"
family="${FIGURE5F_FAMILY:?FIGURE5F_FAMILY is required}"
rscript_bin="${FIGURE5F_RSCRIPT_BIN:?FIGURE5F_RSCRIPT_BIN is required}"
task_index="${SLURM_PROCID:?SLURM_PROCID is required}"

case "${family}" in C01|C02|C03) ;; *) echo "Invalid family: ${family}" >&2; exit 64 ;; esac
if (( task_index < 0 || task_index > 1 )); then
  echo "SLURM_PROCID must be 0 or 1; received ${task_index}." >&2
  exit 65
fi

replicate="$((task_index + 1))"
generator="${iteration_root}/Code/Figures/generate_Figure5F_generalized_posterior.R"

echo "Figure 5F family-conditioned ladder: ${family}_R${replicate}; target n_iter=${FIGURE5F_SAMPLE_ITER}"
exec "${rscript_bin}" "${generator}" \
  --mode=sample \
  --n_iter="${FIGURE5F_SAMPLE_ITER:?FIGURE5F_SAMPLE_ITER is required}" \
  --burnin="${FIGURE5F_SAMPLE_BURNIN:-3000}" \
  --thin=1 \
  --ladders_per_family=2 \
  --cores=2 \
  --temperature_cores_per_ladder=5 \
  --checkpoint_every="${FIGURE5F_CHECKPOINT_EVERY:-250}" \
  --ladder_family="${family}" \
  --ladder_replicate="${replicate}" \
  --force=FALSE \
  --resume=TRUE
