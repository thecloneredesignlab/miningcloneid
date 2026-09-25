#!/usr/bin/env bash
set -euo pipefail

stage="${1:?stage is required: prepare, sample, aggregate, or combine}"
iteration_root="${FIGURE5F_ITERATION_ROOT:?FIGURE5F_ITERATION_ROOT is required}"
result_root="${FIGURE5F_RESULT_ROOT:?FIGURE5F_RESULT_ROOT is required}"
r_module="${FIGURE5F_R_MODULE:-R/4.4}"
rscript_bin="/app/eb/software/R/4.4.2-gfbf-2024a/bin/Rscript"

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
if [[ ! -x "${rscript_bin}" ]]; then
  echo "Required Rscript is unavailable: ${rscript_bin}" >&2
  exit 5
fi
echo "Figure 5F R module: ${r_module}"
echo "Figure 5F Rscript: ${rscript_bin}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1
export FIGURE5F_RESULT_ROOT="${result_root}"
export FIGURE5F_RSCRIPT_BIN="${rscript_bin}"

preparer="${iteration_root}/Code/Figures/prepare_Figure5F_family_conditioning.R"
aggregator="${iteration_root}/Code/Figures/aggregate_Figure5F_family_posterior.R"
combiner="${iteration_root}/Code/Figures/combine_Figure5F_family_posteriors.R"
auditor="${iteration_root}/Code/Figures/audit_Figure5F_prior_posterior_inputs.R"
product_builder="${iteration_root}/Code/Figures/build_Figure5F_generalized_posterior_products.R"
table_builder="${iteration_root}/Code/Figures/build_Figure5F_supplementary_table.R"
figure_builder="${iteration_root}/Code/Figures/draw_Figure5.R"
ladder_runner="${iteration_root}/Code/Figures/Figure5F/hpc/run_Figure5F_family_ladder_task.sh"

case "${stage}" in
  prepare)
    exec "${rscript_bin}" "${preparer}"
    ;;
  sample)
    family="${FIGURE5F_FAMILY:?FIGURE5F_FAMILY is required}"
    if (( ${SLURM_NTASKS:-0} != 2 || ${SLURM_CPUS_PER_TASK:-0} != 5 )); then
      echo "Family sampling requires exactly 2 Slurm tasks x 5 CPUs." >&2
      exit 4
    fi
    exec srun \
      --ntasks=2 \
      --cpus-per-task=5 \
      --cpu-bind=cores \
      --kill-on-bad-exit=1 \
      "${ladder_runner}"
    ;;
  aggregate)
    exec "${rscript_bin}" "${aggregator}" \
      --family="${FIGURE5F_FAMILY:?FIGURE5F_FAMILY is required}" \
      --n_iter="${FIGURE5F_SAMPLE_ITER:?FIGURE5F_SAMPLE_ITER is required}" \
      --burnin="${FIGURE5F_SAMPLE_BURNIN:-3000}" \
      --thin=1
    ;;
  combine)
    "${rscript_bin}" "${combiner}"
    "${rscript_bin}" "${auditor}"
    "${rscript_bin}" "${product_builder}"
    "${rscript_bin}" "${table_builder}"
    repository_root="$(cd "${iteration_root}/../.." && pwd -P)"
    figure5_data_dir="${iteration_root}/data/Figures/Figure5"
    FIGURE5_DRAW_WORKER=1 \
      FIGURE_WORKSPACE_ROOT="${iteration_root}" \
      HYPOXIA_REPO_ROOT="${repository_root}" \
      FIGURE5_DATA_DIR="${figure5_data_dir}" \
      FIGURE5_FIGURE_DIR="${iteration_root}/Figures" \
      FIGURE5_PANEL_DIR="${figure5_data_dir}/panels" \
      JOINT_TSNE_FULL_COORDINATES="${figure5_data_dir}/pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_full_coordinates.csv" \
      JOINT_TSNE_BEST_COORDINATES="${figure5_data_dir}/pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv" \
      "${rscript_bin}" "${figure_builder}"
    echo "Figure 5F release audit, products, supplementary table, and assembled Figure 5 completed."
    ;;
  *)
    echo "Unknown family-conditioned stage: ${stage}" >&2
    exit 64
    ;;
esac
