#!/usr/bin/env bash
set -euo pipefail

stage="${1:?stage is required: prepare, benchmark, pilot_ladders, pilot_aggregate, calibrate, sample_ladders, or sample_aggregate}"
iteration_root="${FIGURE5F_ITERATION_ROOT:?FIGURE5F_ITERATION_ROOT is required}"
result_root="${FIGURE5F_RESULT_ROOT:?FIGURE5F_RESULT_ROOT is required}"
r_module="${FIGURE5F_R_MODULE:-R/4.4}"
generator="${iteration_root}/Code/Figures/generate_Figure5F_generalized_posterior.R"
calibrator="${iteration_root}/Code/Figures/calibrate_Figure5F_pilot_proposal.R"
preparer="${iteration_root}/Code/Figures/Figure5F/prepare_endpoint_proposal.R"
ladder_task_runner="${iteration_root}/Code/Figures/Figure5F/hpc/run_Figure5F_ladder_task.sh"
pilot_state="${iteration_root}/data/Figures/Figure5/generalized_posterior/figure5f_generalized_posterior_state_draws.rds"

if ! type ml >/dev/null 2>&1; then
  for module_init in /etc/profile.d/z00_lmod.sh /etc/profile.d/lmod.sh; do
    if [[ -r "${module_init}" ]]; then
      # shellcheck disable=SC1090
      source "${module_init}"
      break
    fi
  done
fi
if ! type ml >/dev/null 2>&1; then
  echo "The HPC module command 'ml' is unavailable." >&2
  exit 5
fi
ml "${r_module}"
rscript_bin="$(command -v Rscript || true)"
if [[ -z "${rscript_bin}" ]]; then
  echo "Rscript is unavailable after loading ${r_module}." >&2
  exit 5
fi
echo "Figure 5F R module: ${r_module}"
echo "Figure 5F Rscript: ${rscript_bin}"

for path in "${iteration_root}" "${result_root}" "${generator}" "${calibrator}" "${preparer}" "${ladder_task_runner}"; do
  if [[ ! -e "${path}" ]]; then
    echo "Missing required Figure 5F input: ${path}" >&2
    exit 2
  fi
done

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1
export FIGURE5F_RESULT_ROOT="${result_root}"
export FIGURE5F_PILOT_STATE="${pilot_state}"
export FIGURE5F_RSCRIPT_BIN="${rscript_bin}"

ladder_cores="${FIGURE5F_LADDER_CORES:-6}"
temperature_cores_per_ladder="${FIGURE5F_TEMPERATURE_CORES_PER_LADDER:-5}"
required_sampler_cpus="$((ladder_cores * temperature_cores_per_ladder))"
export FIGURE5F_LADDER_TASKS="${ladder_cores}"
export FIGURE5F_TEMPERATURE_CORES_PER_LADDER="${temperature_cores_per_ladder}"
if [[ "${stage}" == "pilot_ladders" || "${stage}" == "sample_ladders" ]]; then
  allocated_tasks="${SLURM_NTASKS:-0}"
  allocated_cpus_per_task="${SLURM_CPUS_PER_TASK:-0}"
  allocated_cpus="$((allocated_tasks * allocated_cpus_per_task))"
  if (( allocated_tasks != ladder_cores ||
        allocated_cpus_per_task != temperature_cores_per_ladder ||
        allocated_cpus != required_sampler_cpus )); then
    echo "Figure 5F ${stage} requires ${ladder_cores} tasks x ${temperature_cores_per_ladder} CPUs = ${required_sampler_cpus} CPUs; received ${allocated_tasks} x ${allocated_cpus_per_task} = ${allocated_cpus}." >&2
    exit 4
  fi
fi

run_rscript() {
  "${rscript_bin}" "$@"
}

case "${stage}" in
  prepare)
    run_rscript "${preparer}"
    ;;
  benchmark)
    run_rscript "${generator}" \
      --mode=benchmark \
      --benchmark_evals="${FIGURE5F_BENCHMARK_EVALS:-10}" \
      --cores=1 \
      --temperature_cores_per_ladder="${temperature_cores_per_ladder}"
    ;;
  pilot_ladders)
    export FIGURE5F_CHECKPOINT_EVERY="${FIGURE5F_CHECKPOINT_EVERY:-100}"
    srun \
      --ntasks="${ladder_cores}" \
      --cpus-per-task="${temperature_cores_per_ladder}" \
      --cpu-bind=cores \
      --kill-on-bad-exit=1 \
      "${ladder_task_runner}" pilot
    ;;
  pilot_aggregate)
    run_rscript "${generator}" \
      --mode=pilot \
      --n_iter="${FIGURE5F_PILOT_ITER:-600}" \
      --burnin="${FIGURE5F_PILOT_BURNIN:-300}" \
      --checkpoint_every="${FIGURE5F_CHECKPOINT_EVERY:-100}" \
      --cores="${ladder_cores}" \
      --temperature_cores_per_ladder="${temperature_cores_per_ladder}" \
      --aggregate_only=TRUE
    ;;
  calibrate)
    if [[ ! -f "${pilot_state}" ]]; then
      echo "Pilot state is missing: ${pilot_state}" >&2
      exit 3
    fi
    run_rscript "${calibrator}"
    ;;
  sample_ladders)
    export FIGURE5F_CHECKPOINT_EVERY="${FIGURE5F_CHECKPOINT_EVERY:-250}"
    srun \
      --ntasks="${ladder_cores}" \
      --cpus-per-task="${temperature_cores_per_ladder}" \
      --cpu-bind=cores \
      --kill-on-bad-exit=1 \
      "${ladder_task_runner}" sample
    ;;
  sample_aggregate)
    run_rscript "${generator}" \
      --mode=sample \
      --n_iter="${FIGURE5F_SAMPLE_ITER:-6000}" \
      --burnin="${FIGURE5F_SAMPLE_BURNIN:-3000}" \
      --checkpoint_every="${FIGURE5F_CHECKPOINT_EVERY:-250}" \
      --cores="${ladder_cores}" \
      --temperature_cores_per_ladder="${temperature_cores_per_ladder}" \
      --aggregate_only=TRUE
    ;;
  *)
    echo "Unknown Figure 5F stage: ${stage}" >&2
    exit 64
    ;;
esac
