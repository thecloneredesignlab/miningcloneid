#!/usr/bin/env bash
set -euo pipefail

mode="${1:?mode is required: pilot or sample}"
case "${mode}" in
  pilot|sample) ;;
  *) echo "Unknown ladder mode: ${mode}" >&2; exit 64 ;;
esac

iteration_root="${FIGURE5F_ITERATION_ROOT:?FIGURE5F_ITERATION_ROOT is required}"
generator="${iteration_root}/Code/Figures/generate_Figure5F_generalized_posterior.R"
rscript_bin="${FIGURE5F_RSCRIPT_BIN:?FIGURE5F_RSCRIPT_BIN is required}"
task_index="${SLURM_PROCID:?SLURM_PROCID is required for independent ladder tasks}"

families=(C01 C01 C02 C02 C03 C03)
replicates=(1 2 1 2 1 2)
if (( task_index < 0 || task_index >= ${#families[@]} )); then
  echo "SLURM_PROCID must be in 0..5; received ${task_index}." >&2
  exit 65
fi

family="${families[task_index]}"
replicate="${replicates[task_index]}"
common_args=(
  "--mode=${mode}"
  "--cores=${FIGURE5F_LADDER_TASKS:-6}"
  "--temperature_cores_per_ladder=${FIGURE5F_TEMPERATURE_CORES_PER_LADDER:-5}"
  "--checkpoint_every=${FIGURE5F_CHECKPOINT_EVERY:-250}"
  "--ladder_family=${family}"
  "--ladder_replicate=${replicate}"
)

echo "Figure 5F independent task ${task_index}: ${family}_R${replicate} (${mode})"
if [[ "${mode}" == "pilot" ]]; then
  exec "${rscript_bin}" "${generator}" \
    "${common_args[@]}" \
    "--n_iter=${FIGURE5F_PILOT_ITER:-600}" \
    "--burnin=${FIGURE5F_PILOT_BURNIN:-300}" \
    --force=TRUE \
    --resume=FALSE
else
  exec "${rscript_bin}" "${generator}" \
    "${common_args[@]}" \
    "--n_iter=${FIGURE5F_SAMPLE_ITER:-6000}" \
    "--burnin=${FIGURE5F_SAMPLE_BURNIN:-3000}" \
    --force=FALSE \
    --resume=TRUE
fi
