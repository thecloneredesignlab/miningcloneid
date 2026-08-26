#!/bin/bash
# Unified O2 HPC submitter for in vivo, in vitro, and joint fitting.

set -euo pipefail

O2SD_SHELL_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"

ORIGINAL_SUBMIT_ARGS=("$@")

usage() {
  cat <<'EOF'
Usage:
  bash submit_o2_fit.sh --fitting_mode=invivo [options]
  bash submit_o2_fit.sh --fitting_mode=invitro [options]
  bash submit_o2_fit.sh --fitting_mode=joint --invivo_run_dir=DIR --invitro_run_dir=DIR [options]
  bash submit_o2_fit.sh --fitting_mode=all [options]

Required modes:
  --fitting_mode=invivo|invitro|joint|all

Joint mode behavior:
  Build a pooled parameter-space t-SNE from the specified source runs, cluster
  only the in-vivo best points, select the objective-minimum seed from every
  in-vivo primary cluster, pair each with the single global objective-minimum
  in-vitro seed, then submit one global pair x seed task-table array. This is the
  only supported primary-cluster workflow.

All mode behavior:
  Submit in vivo first, in vitro after it succeeds, then submit the fixed joint
  primary-cluster controller using those generated result directories.

Common options:
  --project_root=/path/to/repo
  Defaults to the repository containing this submit script.
  --config_path=/path/to/O2_supply_demand.yaml
  --out_root=/path/to/oxygen/results
  --r_module=R/4.4
  --dry_run=TRUE|FALSE
  --force_extra_results=TRUE|FALSE
  --log_root=/path/to/results/log

Single-fit options:
  --invivo_run_prefix=name
  --invitro_run_prefix=name
  --invivo_run_dir=/path/to/existing_invivo_run
  --invitro_run_dir=/path/to/existing_invitro_run
  --invivo_total_seeds=500 --invivo_array_tasks=500 --invivo_seeds_per_task=1
  --invitro_total_seeds=500 --invitro_array_tasks=500 --invitro_seeds_per_task=1
  --invivo_qos=xlarge --invivo_time_limit=12:00:00
  --invitro_qos=xxlarge --invitro_time_limit=12:00:00
  --itermax=1000 --itermax_max=1000

Joint options:
  --joint_run_prefix=name
  --joint_job_name=o2_joint_B
  --joint_dependency=JOBID_OR_ARRAY_WILDCARD
  --joint_total_seeds=500 --joint_array_tasks=500 --joint_seeds_per_task=1
  --joint_qos=xxlarge --joint_time_limit=12:00:00
  --postprocess_qos=small --postprocess_time_limit=4:00:00
  --parameter_table=/path/to/invitro_parameter_table.csv
  --fit_objects_dir=/path/to/fit_objects
  --flow_density_path=/path/to/g0g1_ploidy_density_grid.csv
  --joint_soft_coupling_sigma_default=0.65
  --joint_soft_coupling_welsch_c=0.4
  --joint_warmup_sigmaN=0.0304
  --joint_soft_coupling_delta_params=default|all|none|param1,param2
  --prep_qos=small --prep_time_limit=2:00:00 --prep_mem=8G
  --joint_tsne_seed=123
  --joint_cluster_seed=123
  --joint_landscape_max_seeds=N
  --joint_landscape_n_threads=1
  --joint_seed_space_qos=small
  --joint_seed_space_time_limit=1:00:00
  --joint_seed_space_mem=2G
  --joint_seed_space_array_max_concurrent=N
  --seeds_per_pair=200
  --array_max_concurrent=100
  --multi_warmup_task_order=round_robin|pair_major
  --multi_warmup_task_status_filter=all|not_started
  --skip_existing=TRUE|FALSE

After each fitting finishes:
  A dependent postprocess job runs extra_results.R. If extra_results already
  exists, it is skipped unless force_extra_results=TRUE.
EOF
}

check_seed_plan() {
  local label="$1"
  local total="$2"
  local tasks="$3"
  local per_task="$4"
  if (( tasks * per_task != total )); then
    echo "Mismatch for ${label}: ARRAY_TASKS * SEEDS_PER_TASK must equal TOTAL_SEEDS." >&2
    echo "Got total=${total}, array_tasks=${tasks}, seeds_per_task=${per_task}." >&2
    exit 2
  fi
}

record_array_submission() {
  local run_dir="$1"
  local label="$2"
  local job_id="$3"
  local total_seeds="$4"
  local array_tasks="$5"
  local seeds_per_task="$6"
  local qos="$7"
  local walltime="$8"
  local mem="$9"
  local cpus="${10}"
  local sub_script="${11}"
  local runner_script="${12}"
  shift 12
  local sbatch_command
  sbatch_command="$(o2sd_prov_shell_join "$@")"
  o2sd_prov_write_standard "${run_dir}" "${SELF_SCRIPT}" "${SUBMIT_COMMAND_TEXT:-}"
  o2sd_prov_record_sbatch "${run_dir}" "${sbatch_command}" "${job_id}"
  o2sd_prov_write_many "${run_dir}" \
    scripts submit_script "${SELF_SCRIPT}" \
    scripts array_script "${sub_script}" \
    scripts runner_script "${runner_script}" \
    scripts postprocess_script "${POSTPROCESS_SCRIPT}" \
    scripts extra_results_script "${EXTRA_RESULTS_SCRIPT}" \
    input_config config_path "${CONFIG_PATH}" \
    fit fitting_mode "${label}" \
    fit run_prefix "$(basename "${run_dir}")" \
    fit total_seeds "${total_seeds}" \
    fit array_tasks "${array_tasks}" \
    fit seeds_per_task "${seeds_per_task}" \
    optimizer itermax "${ITERMAX:-NA}" \
    optimizer itermax_max "${ITERMAX_MAX:-NA}" \
    optimizer NP "${NP:-NA}" \
    optimizer de_reltol "${DE_RELTOL:-NA}" \
    optimizer de_steptol "${DE_STEPTOL:-NA}" \
    optimizer n_cores "${cpus}" \
    slurm array_job_id "${job_id}" \
    slurm qos "${qos}" \
    slurm walltime "${walltime}" \
    slurm memory "${mem}" \
    slurm cpus "${cpus}" \
    slurm log_root "${LOG_ROOT}"
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h)
        usage
        exit 0
        ;;
      --fitting_mode=*) FITTING_MODE="${arg#*=}" ;;
      --joint_fitting_mode=*)
        echo "--joint_fitting_mode has been removed; --fitting_mode=joint always pairs in-vivo primary clusters with the global-best in-vitro seed." >&2
        exit 2
        ;;
      --internal_stage=*) INTERNAL_STAGE="${arg#*=}" ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --config_path=*) CONFIG_PATH="${arg#*=}" ;;
      --out_root=*) OUT_ROOT="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      --force_extra_results=*) FORCE_EXTRA_RESULTS="${arg#*=}" ;;
      --invivo_run_dir=*) INVIVO_RUN_DIR="${arg#*=}" ;;
      --invitro_run_dir=*) INVITRO_RUN_DIR="${arg#*=}" ;;
      --invivo_run_prefix=*) INVIVO_RUN_PREFIX="${arg#*=}" ;;
      --invitro_run_prefix=*) INVITRO_RUN_PREFIX="${arg#*=}" ;;
      --joint_run_prefix=*) JOINT_RUN_PREFIX="${arg#*=}" ;;
      --multi_warmup_prefix=*) JOINT_RUN_PREFIX="${arg#*=}" ;;
      --joint_job_name=*) JOINT_JOB_NAME="${arg#*=}" ;;
      --joint_dependency=*) JOINT_DEPENDENCY="${arg#*=}" ;;
      --invivo_total_seeds=*) INVIVO_TOTAL_SEEDS="${arg#*=}" ;;
      --invivo_array_tasks=*) INVIVO_ARRAY_TASKS="${arg#*=}" ;;
      --invivo_seeds_per_task=*) INVIVO_SEEDS_PER_TASK="${arg#*=}" ;;
      --invitro_total_seeds=*) INVITRO_TOTAL_SEEDS="${arg#*=}" ;;
      --invitro_array_tasks=*) INVITRO_ARRAY_TASKS="${arg#*=}" ;;
      --invitro_seeds_per_task=*) INVITRO_SEEDS_PER_TASK="${arg#*=}" ;;
      --joint_total_seeds=*) JOINT_TOTAL_SEEDS="${arg#*=}" ;;
      --joint_array_tasks=*) JOINT_ARRAY_TASKS="${arg#*=}" ;;
      --joint_seeds_per_task=*) JOINT_SEEDS_PER_TASK="${arg#*=}" ;;
      --n_cores=*) N_CORES="${arg#*=}" ;;
      --invivo_n_cores=*) INVIVO_N_CORES="${arg#*=}" ;;
      --invitro_n_cores=*) INVITRO_N_CORES="${arg#*=}" ;;
      --joint_n_cores=*) JOINT_N_CORES="${arg#*=}" ;;
      --mem=*) MEM="${arg#*=}" ;;
      --invivo_mem=*) INVIVO_MEM="${arg#*=}" ;;
      --invitro_mem=*) INVITRO_MEM="${arg#*=}" ;;
      --joint_mem=*) JOINT_MEM="${arg#*=}" ;;
      --invivo_qos=*) INVIVO_QOS="${arg#*=}" ;;
      --invitro_qos=*) INVITRO_QOS="${arg#*=}" ;;
      --joint_qos=*) JOINT_QOS="${arg#*=}" ;;
      --postprocess_qos=*) POSTPROCESS_QOS="${arg#*=}" ;;
      --invivo_time_limit=*) INVIVO_TIME_LIMIT="${arg#*=}" ;;
      --invitro_time_limit=*) INVITRO_TIME_LIMIT="${arg#*=}" ;;
      --joint_time_limit=*) JOINT_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_time_limit=*) POSTPROCESS_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_mem=*) POSTPROCESS_MEM="${arg#*=}" ;;
      --parameter_table=*) PARAMETER_TABLE="${arg#*=}" ;;
      --invitro_parameter_table=*) PARAMETER_TABLE="${arg#*=}" ;;
      --parameter_table_invitro=*) PARAMETER_TABLE="${arg#*=}" ;;
      --fit_objects_dir=*) FIT_OBJECTS_DIR="${arg#*=}" ;;
      --flow_density_path=*) FLOW_DENSITY_PATH="${arg#*=}" ;;
      --invivo_best_seed_dir=*|--joint_warmup_invivo_seed_dir=*|--joint_warmup_invivo_best_seed_dir=*|--invitro_best_seed_dir=*|--joint_warmup_invitro_seed_dir=*|--joint_warmup_invitro_best_seed_dir=*|--joint_warmup_vitro_seed_dir=*|--joint_warmup_enable=*|--joint_warmup_seed_label=*|--joint_seed_label=*|--seed_label=*|--joint_soft_coupling_parameters_table=*|--joint_soft_coupling_parameters_table_path=*|--select_required_files=*|--invivo_objective_columns=*|--invitro_objective_columns=*)
        echo "${arg%%=*} has been removed; joint anchors and coupling tables are derived from the specified source runs." >&2
        exit 2
        ;;
      --joint_warmup_sigmaN=*) JOINT_WARMUP_SIGMAN="${arg#*=}" ;;
      --joint_soft_coupling_sigma_default=*) JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${arg#*=}" ;;
      --joint_soft_coupling_welsch_c=*) JOINT_SOFT_COUPLING_WELSCH_C="${arg#*=}" ;;
      --joint_soft_coupling_delta_params=*) JOINT_SOFT_COUPLING_DELTA_PARAMS="${arg#*=}" ;;
      --joint_landscape_max_seeds=*|--multi_warmup_landscape_max_seeds=*|--landscape_max_seeds=*) MULTI_WARMUP_LANDSCAPE_MAX_SEEDS="${arg#*=}" ;;
      --joint_landscape_n_threads=*|--multi_warmup_n_threads=*|--landscape_n_threads=*) MULTI_WARMUP_N_THREADS="${arg#*=}" ;;
      --joint_seed_space_qos=*|--multi_warmup_seed_space_qos=*|--seed_space_qos=*) MULTI_WARMUP_SEED_SPACE_QOS="${arg#*=}" ;;
      --joint_seed_space_time_limit=*|--multi_warmup_seed_space_time_limit=*|--seed_space_time_limit=*) MULTI_WARMUP_SEED_SPACE_TIME_LIMIT="${arg#*=}" ;;
      --joint_seed_space_mem=*|--multi_warmup_seed_space_mem=*|--seed_space_mem=*) MULTI_WARMUP_SEED_SPACE_MEM="${arg#*=}" ;;
      --joint_seed_space_array_max_concurrent=*|--multi_warmup_seed_space_array_max_concurrent=*|--seed_space_array_max_concurrent=*) MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT="${arg#*=}" ;;
      --joint_cluster_seed=*|--multi_warmup_cluster_seed=*|--landscape_cluster_seed=*) MULTI_WARMUP_CLUSTER_SEED="${arg#*=}" ;;
      --joint_tsne_seed=*|--multi_warmup_tsne_seed=*|--landscape_tsne_seed=*) MULTI_WARMUP_TSNE_SEED="${arg#*=}" ;;
      --joint_analysis_root=*) MULTI_WARMUP_ANALYSIS_ROOT="${arg#*=}" ;;
      --multi_warmup_pair_method=*|--pair_method=*|--multi_warmup_reductions=*|--landscape_reductions=*|--multi_warmup_subcluster_seed=*|--landscape_subcluster_seed=*|--multi_warmup_pairing_policy=*|--pairing_policy=*|--multi_warmup_deduplicate_pairs=*|--deduplicate_pairs=*|--multi_warmup_invivo_curve_filter=*|--invivo_curve_filter=*|--multi_warmup_invivo_curve_class=*|--invivo_curve_class=*|--multi_warmup_curve_filter_comparison=*|--curve_filter_comparison=*|--multi_warmup_reference_subcluster_dir=*|--reference_subcluster_dir=*)
        echo "${arg%%=*} has been removed; joint fitting now has one fixed in-vivo-cluster/global-in-vitro-best workflow." >&2
        exit 2
        ;;
      --multi_warmup_monotonicity_qos=*|--multi_warmup_monotonicity_time_limit=*|--multi_warmup_monotonicity_mem=*|--multi_warmup_monotonicity_cpus=*|--multi_warmup_monotonicity_tasks_per_array_task=*|--multi_warmup_validation_qos=*|--multi_warmup_validation_time_limit=*|--multi_warmup_validation_mem=*)
        echo "${arg%%=*} has been removed because the joint workflow no longer runs curve filtering or second-level validation." >&2
        exit 2
        ;;
      --seeds_per_pair=*|--joint_seeds_per_pair=*|--multi_warmup_seeds_per_pair=*) MULTI_WARMUP_SEEDS_PER_PAIR="${arg#*=}" ;;
      --array_max_concurrent=*|--multi_warmup_array_max_concurrent=*) MULTI_WARMUP_ARRAY_MAX_CONCURRENT="${arg#*=}" ;;
      --multi_warmup_task_order=*|--task_order=*) MULTI_WARMUP_TASK_ORDER="${arg#*=}" ;;
      --multi_warmup_task_status_filter=*|--task_status_filter=*) MULTI_WARMUP_TASK_STATUS_FILTER="${arg#*=}" ;;
      --skip_existing=*|--multi_warmup_skip_existing=*) MULTI_WARMUP_SKIP_EXISTING="${arg#*=}" ;;
      --refresh_task_status=*|--multi_warmup_refresh_task_status=*) MULTI_WARMUP_REFRESH_TASK_STATUS="${arg#*=}" ;;
      --prep_qos=*) PREP_QOS="${arg#*=}" ;;
      --prep_time_limit=*) PREP_TIME_LIMIT="${arg#*=}" ;;
      --prep_mem=*) PREP_MEM="${arg#*=}" ;;
      --report_qos=*) REPORT_QOS="${arg#*=}" ;;
      --report_time_limit=*) REPORT_TIME_LIMIT="${arg#*=}" ;;
      --report_mem=*) REPORT_MEM="${arg#*=}" ;;
      --itermax=*) ITERMAX="${arg#*=}" ;;
      --itermax_max=*) ITERMAX_MAX="${arg#*=}" ;;
      --de_reltol=*) DE_RELTOL="${arg#*=}" ;;
      --de_steptol=*) DE_STEPTOL="${arg#*=}" ;;
      --np=*|--NP=*) NP="${arg#*=}" ;;
      --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
      --passage_mode=*)
        echo "--passage_mode has been removed; in vitro always uses the fixed v2 passage implementation." >&2
        exit 2
        ;;
      --log_root=*|--log_dir=*) LOG_ROOT="${arg#*=}" ;;
      *)
        echo "Unknown argument: ${arg}" >&2
        usage >&2
        exit 2
        ;;
    esac
  done
}

ensure_sbatch() {
  if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch not found; run this launcher on the HPC login node or use DRY_RUN=TRUE." >&2
    exit 1
  fi
}

submit_extra_results_job() {
  local label="$1"
  local run_dir="$2"
  local dependency="$3"
  local export_arg="ALL"
  export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
  export_arg+=",RUN_DIR=${run_dir}"
  export_arg+=",EXTRA_RESULTS_SCRIPT=${EXTRA_RESULTS_SCRIPT}"
  export_arg+=",R_MODULE=${R_MODULE}"
  export_arg+=",FORCE_EXTRA_RESULTS=${FORCE_EXTRA_RESULTS}"

  local cmd=(
    sbatch
    --parsable
    "--job-name=${label}_extra"
    "--qos=${POSTPROCESS_QOS}"
    "--time=${POSTPROCESS_TIME_LIMIT}"
    --cpus-per-task=1
    "--mem=${POSTPROCESS_MEM}"
    "--output=${LOG_ROOT}/${label}_extra_%A.out"
    "--error=${LOG_ROOT}/${label}_extra_%A.err"
    "--export=${export_arg}"
  )
  if [[ -n "${dependency}" ]]; then
    cmd+=("--dependency=afterok:${dependency}")
  fi
  cmd+=("${POSTPROCESS_SCRIPT}")

  if truthy "${DRY_RUN}"; then
    print_command "Submit ${label} extra_results" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_${label}_EXTRA_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted ${label} extra_results job: ${LAST_JOB_ID}"
  fi
}

submit_invivo_array() {
  local dependency="${1:-}"
  local run_dir="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
  mkdir -p "${run_dir}"
  local export_arg="ALL"
  export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
  export_arg+=",RUNNER_SCRIPT=${INVIVO_RUNNER_SCRIPT}"
  export_arg+=",CONFIG_PATH=${CONFIG_PATH}"
  export_arg+=",OUT_ROOT=${OUT_ROOT}"
  export_arg+=",RUN_PREFIX=${INVIVO_RUN_PREFIX}"
  export_arg+=",TOTAL_SEEDS=${INVIVO_TOTAL_SEEDS}"
  export_arg+=",ARRAY_TASKS=${INVIVO_ARRAY_TASKS}"
  export_arg+=",SEEDS_PER_TASK=${INVIVO_SEEDS_PER_TASK}"
  export_arg+=",N_CORES=${INVIVO_N_CORES}"
  export_arg+=",ITERMAX=${ITERMAX}"
  export_arg+=",AUTO_VIZ=${AUTO_VIZ}"
  export_arg+=",R_MODULE=${R_MODULE}"
  local cmd=(
    sbatch
    --parsable
    --job-name=o2_ivv_B
    "--qos=${INVIVO_QOS}"
    "--time=${INVIVO_TIME_LIMIT}"
    "--cpus-per-task=${INVIVO_N_CORES}"
    "--mem=${INVIVO_MEM}"
    "--array=1-${INVIVO_ARRAY_TASKS}"
    "--output=${LOG_ROOT}/o2_invivo_%A_%a.out"
    "--error=${LOG_ROOT}/o2_invivo_%A_%a.err"
    "--export=${export_arg}"
  )
  if [[ -n "${dependency}" ]]; then
    cmd+=("--dependency=afterok:${dependency}")
  fi
  cmd+=("${INVIVO_SUB_SCRIPT}")
  if truthy "${DRY_RUN}"; then
    print_command "Submit in vivo" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_INVIVO_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted in vivo array job: ${LAST_JOB_ID}"
  fi
  record_array_submission "${run_dir}" "invivo" "${LAST_JOB_ID}" \
    "${INVIVO_TOTAL_SEEDS}" "${INVIVO_ARRAY_TASKS}" "${INVIVO_SEEDS_PER_TASK}" \
    "${INVIVO_QOS}" "${INVIVO_TIME_LIMIT}" "${INVIVO_MEM}" "${INVIVO_N_CORES}" \
    "${INVIVO_SUB_SCRIPT}" "${INVIVO_RUNNER_SCRIPT}" "${cmd[@]}"
}

submit_invitro_array() {
  local dependency="${1:-}"
  local run_dir="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
  mkdir -p "${run_dir}"
  local export_arg="ALL"
  export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
  export_arg+=",RUNNER_SCRIPT=${INVITRO_RUNNER_SCRIPT}"
  export_arg+=",CONFIG_PATH=${CONFIG_PATH}"
  export_arg+=",OUT_ROOT=${OUT_ROOT}"
  export_arg+=",RUN_PREFIX=${INVITRO_RUN_PREFIX}"
  export_arg+=",TOTAL_SEEDS=${INVITRO_TOTAL_SEEDS}"
  export_arg+=",ARRAY_TASKS=${INVITRO_ARRAY_TASKS}"
  export_arg+=",SEEDS_PER_TASK=${INVITRO_SEEDS_PER_TASK}"
  export_arg+=",N_CORES=${INVITRO_N_CORES}"
  export_arg+=",R_MODULE=${R_MODULE}"
  export_arg+=",PARAMETER_TABLE=${PARAMETER_TABLE}"
  export_arg+=",FIT_OBJECTS_DIR=${FIT_OBJECTS_DIR}"
  export_arg+=",FLOW_DENSITY_PATH=${FLOW_DENSITY_PATH}"
  export_arg+=",ITERMAX=${ITERMAX}"
  export_arg+=",ITERMAX_MAX=${ITERMAX_MAX}"
  export_arg+=",DE_RELTOL=${DE_RELTOL}"
  export_arg+=",DE_STEPTOL=${DE_STEPTOL}"
  export_arg+=",NP=${NP}"
  export_arg+=",AUTO_VIZ=${AUTO_VIZ}"
  local cmd=(
    sbatch
    --parsable
    --job-name=o2_ivt_B
    "--qos=${INVITRO_QOS}"
    "--time=${INVITRO_TIME_LIMIT}"
    "--cpus-per-task=${INVITRO_N_CORES}"
    "--mem=${INVITRO_MEM}"
    "--array=1-${INVITRO_ARRAY_TASKS}"
    "--output=${LOG_ROOT}/o2_invitro_%A_%a.out"
    "--error=${LOG_ROOT}/o2_invitro_%A_%a.err"
    "--export=${export_arg}"
  )
  if [[ -n "${dependency}" ]]; then
    cmd+=("--dependency=afterok:${dependency}")
  fi
  cmd+=("${INVITRO_SUB_SCRIPT}")
  if truthy "${DRY_RUN}"; then
    print_command "Submit in vitro" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_INVITRO_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted in vitro array job: ${LAST_JOB_ID}"
  fi
  record_array_submission "${run_dir}" "invitro" "${LAST_JOB_ID}" \
    "${INVITRO_TOTAL_SEEDS}" "${INVITRO_ARRAY_TASKS}" "${INVITRO_SEEDS_PER_TASK}" \
    "${INVITRO_QOS}" "${INVITRO_TIME_LIMIT}" "${INVITRO_MEM}" "${INVITRO_N_CORES}" \
    "${INVITRO_SUB_SCRIPT}" "${INVITRO_RUNNER_SCRIPT}" "${cmd[@]}"
}

run_multi_warmup_finalize_stage() {
  INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
  INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
  load_r_module

  local multi_root="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
  local progress_log="${multi_root}/multi_warmup_progress.log"
  local manifest="${multi_root}/multi_warmup_manifest.tsv"
  local mode_file="${multi_root}/multi_warmup_seed_plan_mode.tsv"
  local task_table="${multi_root}/multi_warmup_tasks.tsv"
  local analysis_root="${MULTI_WARMUP_ANALYSIS_ROOT:-${multi_root}/joint_primary_clusters}"
  mkdir -p "${multi_root}" "${multi_root}/joint_soft_coupling_tables"
  if [[ "${INTERNAL_STAGE}" == "multi_warmup_finalize_and_submit" && -f "${progress_log}" ]]; then
    touch "${progress_log}"
  else
    : > "${progress_log}"
  fi
  o2sd_prov_write_standard "${multi_root}" "${SELF_SCRIPT}" "${SUBMIT_COMMAND_TEXT:-}"
  o2sd_prov_write_many "${multi_root}" \
    fit fitting_mode "joint" \
    fit run_prefix "${JOINT_RUN_PREFIX}" \
    joint invivo_run_dir "${INVIVO_RUN_DIR}" \
    joint invitro_run_dir "${INVITRO_RUN_DIR}" \
    joint pair_method "invivo_primary_clusters_to_global_invitro_best" \
    joint reduction "tsne" \
    joint manifest_path "${manifest}" \
    joint task_table_path "${task_table}" \
    joint seeds_per_pair "${MULTI_WARMUP_SEEDS_PER_PAIR}" \
    joint joint_soft_coupling_sigma_default "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" \
    joint joint_soft_coupling_welsch_c "${JOINT_SOFT_COUPLING_WELSCH_C}"

  echo "Joint primary-cluster finalizer"
  echo "  joint_root: ${multi_root}"
  echo "  pair_method: invivo_primary_clusters_to_global_invitro_best"
  echo "  invivo_run_dir: ${INVIVO_RUN_DIR}"
  echo "  invitro_run_dir: ${INVITRO_RUN_DIR}"
  echo "  seeds_per_pair: ${MULTI_WARMUP_SEEDS_PER_PAIR}"
  echo "  joint_soft_coupling_sigma_default: ${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}"
  echo "  joint_soft_coupling_welsch_c: ${JOINT_SOFT_COUPLING_WELSCH_C}"

  local landscape_table_dir="${analysis_root}/pooled_invivo_invitro/full_data_in_vivo_clustring/Tables"
  local invivo_landscape_summary="${landscape_table_dir}/pooled_invivo_best_primary_cluster_summary.csv"
  [[ -f "${invivo_landscape_summary}" ]] || { echo "Missing prepared in-vivo primary-cluster summary: ${invivo_landscape_summary}" >&2; exit 1; }
  local seed_table_dir="${analysis_root}/SeedParameterTables"
  local -a seed_plan_cmd=(
    Rscript "${MULTI_WARMUP_SEED_PLAN_SCRIPT}"
    --stage=finalize_pairs
    "--project_root=${PROJECT_ROOT}"
    "--invivo_run_dir=${INVIVO_RUN_DIR}"
    "--invitro_run_dir=${INVITRO_RUN_DIR}"
    "--out_dir=${multi_root}"
    "--analysis_root=${analysis_root}"
    "--tsne_seed=${MULTI_WARMUP_TSNE_SEED}"
    "--cluster_seed=${MULTI_WARMUP_CLUSTER_SEED}"
    "--invitro_best_csv=${seed_table_dir}/invitro_best_params_by_seed.csv"
  )
  print_command "Generate multi-warmup seed plan" "${seed_plan_cmd[@]}" | tee -a "${progress_log}"
  if ! truthy "${DRY_RUN}"; then
    "${seed_plan_cmd[@]}" 2>&1 | tee -a "${progress_log}"
  fi
  if ! truthy "${DRY_RUN}" && [[ ! -f "${manifest}" ]]; then
    echo "Missing generated manifest: ${manifest}" >&2
    exit 1
  fi

  if [[ -f "${manifest}" ]]; then
    local total_pairs
    total_pairs=$(( $(wc -l < "${manifest}") - 1 ))
    if (( total_pairs < 1 )); then
      local plan_mode=""
      if [[ -f "${mode_file}" ]]; then
        plan_mode="$(awk -F $'\t' '$1 == "mode" {print $2; exit}' "${mode_file}")"
      fi
      case "${plan_mode}" in
        invivo_only|invitro_only)
          echo "Single-side cluster-only seed plan completed: mode=${plan_mode}, manifest=${manifest}" | tee -a "${progress_log}"
          return 0
          ;;
        *)
          echo "Generated manifest has no warm-up pairs: ${manifest}" >&2
          exit 1
          ;;
      esac
    fi
    local pair_index=0
    tail -n +2 "${manifest}" | while IFS=$'\t' read -r warmup_label phase invivo_family invivo_rank invivo_seed invivo_seed_dir invitro_family invitro_rank invitro_seed invitro_seed_dir selection_reason joint_run_prefix joint_table; do
      pair_index=$((pair_index + 1))
      echo "Preparing warm-up pair ${pair_index}/${total_pairs}: ${warmup_label}" | tee -a "${progress_log}"
      if [[ -z "${joint_table}" ]]; then
        echo "Manifest row for ${warmup_label} is missing joint_soft_coupling_parameters_table." >&2
        exit 1
      fi
      case "${joint_table}" in
        /*) ;;
        *) joint_table="${multi_root}/${joint_table}" ;;
      esac
      mkdir -p "$(dirname "${joint_table}")"
      local table_cmd=(
        Rscript "${JOINT_WARM_START_SCRIPT}"
        "--invivo-seed-dir=${invivo_seed_dir}"
        "--invitro-seed-dir=${invitro_seed_dir}"
        "--seed-label=${warmup_label}"
        "--out=${joint_table}"
        "--delta-params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}"
      )
      print_command "Generate pair soft-coupling table" "${table_cmd[@]}" | tee -a "${progress_log}"
      if ! truthy "${DRY_RUN}"; then
        "${table_cmd[@]}" 2>&1 | tee -a "${progress_log}"
      fi
      if ! truthy "${DRY_RUN}" && [[ ! -f "${joint_table}" ]]; then
        echo "Pair soft-coupling table was not created for ${warmup_label}: ${joint_table}" >&2
        exit 1
      fi
    done
  fi

  local build_task_cmd=(
    Rscript "${MULTI_WARMUP_TASK_TABLE_SCRIPT}"
    "--multi_warmup_root=${multi_root}"
    "--out=${task_table}"
    "--project_root=${PROJECT_ROOT}"
    "--config_path=${CONFIG_PATH}"
    "--runner_script=${JOINT_RUNNER_SCRIPT}"
    "--seeds_per_pair=${MULTI_WARMUP_SEEDS_PER_PAIR}"
    "--array_tasks=${MULTI_WARMUP_SEEDS_PER_PAIR}"
    "--seeds_per_task=1"
    "--order=${MULTI_WARMUP_TASK_ORDER}"
    "--refresh_status=${MULTI_WARMUP_REFRESH_TASK_STATUS}"
    "--log_root=${LOG_ROOT}"
    "--joint_n_cores=${JOINT_N_CORES}"
    "--r_module=${R_MODULE}"
    "--parameter_table=${PARAMETER_TABLE}"
    "--fit_objects_dir=${FIT_OBJECTS_DIR}"
    "--flow_density_path=${FLOW_DENSITY_PATH}"
    "--itermax=${ITERMAX}"
    "--de_reltol=${DE_RELTOL}"
    "--de_steptol=${DE_STEPTOL}"
    "--NP=${NP}"
    "--auto_viz=${AUTO_VIZ}"
  )
  if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then
    build_task_cmd+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then
    build_task_cmd+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_WELSCH_C}" ]]; then
    build_task_cmd+=("--joint_soft_coupling_welsch_c=${JOINT_SOFT_COUPLING_WELSCH_C}")
  fi
  print_command "Build multi-warmup task table" "${build_task_cmd[@]}" | tee -a "${progress_log}"
  if ! truthy "${DRY_RUN}"; then
    "${build_task_cmd[@]}" 2>&1 | tee -a "${progress_log}"
  fi
  if ! truthy "${DRY_RUN}" && [[ ! -f "${task_table}" ]]; then
    echo "Missing generated task table: ${task_table}" >&2
    exit 1
  fi

  local submit_task_cmd=(
    bash "${MULTI_WARMUP_TASK_SUBMIT_SCRIPT}"
    "--project_root=${PROJECT_ROOT}"
    "--tasks_tsv=${task_table}"
    "--log_root=${LOG_ROOT}"
    "--r_module=${R_MODULE}"
    "--task_status_filter=${MULTI_WARMUP_TASK_STATUS_FILTER}"
    "--skip_existing=${MULTI_WARMUP_SKIP_EXISTING}"
    "--joint_n_cores=${JOINT_N_CORES}"
    "--joint_mem=${JOINT_MEM}"
    "--joint_qos=${JOINT_QOS}"
    "--joint_time_limit=${JOINT_TIME_LIMIT}"
    "--postprocess_qos=${POSTPROCESS_QOS}"
    "--postprocess_time_limit=${POSTPROCESS_TIME_LIMIT}"
    "--postprocess_mem=${POSTPROCESS_MEM}"
    "--report_qos=${REPORT_QOS}"
    "--report_time_limit=${REPORT_TIME_LIMIT}"
    "--report_mem=${REPORT_MEM}"
    "--force_extra_results=${FORCE_EXTRA_RESULTS}"
    "--dry_run=${DRY_RUN}"
  )
  if [[ -n "${MULTI_WARMUP_ARRAY_MAX_CONCURRENT}" ]]; then
    submit_task_cmd+=("--array_max_concurrent=${MULTI_WARMUP_ARRAY_MAX_CONCURRENT}")
  fi
  print_command "Submit multi-warmup task-table array" "${submit_task_cmd[@]}" | tee -a "${progress_log}"
  "${submit_task_cmd[@]}" 2>&1 | tee -a "${progress_log}"
}

submit_landscape_seed_space_array() {
  local dataset="$1"
  local tasks_tsv="$2"
  local task_count
  if [[ -f "${tasks_tsv}" ]]; then
    task_count=$(( $(wc -l < "${tasks_tsv}") - 1 ))
  elif truthy "${DRY_RUN}"; then
    task_count=1
  else
    echo "Missing seed-space task table for ${dataset}: ${tasks_tsv}" >&2
    exit 1
  fi
  if (( task_count < 1 )); then
    echo "No seed-space tasks for ${dataset}: ${tasks_tsv}" >&2
    exit 1
  fi
  local array_spec="1-${task_count}"
  if [[ -n "${MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT}" ]]; then
    array_spec="${array_spec}%${MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT}"
  fi
  local export_arg="ALL"
  export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
  export_arg+=",R_MODULE=${R_MODULE}"
  export_arg+=",TASKS_TSV=${tasks_tsv}"
  export_arg+=",TASK_LOOKUP_COLUMN=task_id"
  export_arg+=",SEED_SPACE_SCRIPT=${MULTI_WARMUP_SEED_PLAN_SCRIPT}"
  export_arg+=",SKIP_EXISTING=${MULTI_WARMUP_SKIP_EXISTING}"
  local cmd=(
    sbatch
    --parsable
    "--job-name=o2mw_ss_${dataset}"
    "--array=${array_spec}"
    --cpus-per-task=1
    "--mem=${MULTI_WARMUP_SEED_SPACE_MEM}"
    "--qos=${MULTI_WARMUP_SEED_SPACE_QOS}"
    "--time=${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}"
    "--output=${LOG_ROOT}/o2mw_seedspace_${dataset}_%A_%a.out"
    "--error=${LOG_ROOT}/o2mw_seedspace_${dataset}_%A_%a.err"
    "--export=${export_arg}"
    "${LANDSCAPE_SEED_SPACE_ARRAY_SCRIPT}"
  )
  if truthy "${DRY_RUN}"; then
    print_command "Submit ${dataset} seed-space array" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_${dataset}_SEED_SPACE_ARRAY"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted ${dataset} seed-space array job: ${LAST_JOB_ID}"
  fi
}

submit_landscape_seed_space_collector() {
  local dataset="$1"
  local tasks_tsv="$2"
  local tables_dir="$3"
  local dependency="$4"
  local inner_cmd="if command -v ml >/dev/null 2>&1; then ml ${R_MODULE}; elif command -v module >/dev/null 2>&1; then module load ${R_MODULE}; fi; $(shell_join cd "${PROJECT_ROOT}") && $(shell_join Rscript "${MULTI_WARMUP_SEED_PLAN_SCRIPT}" "--stage=collect_seed_space_tables" "--dataset=${dataset}" "--tasks_tsv=${tasks_tsv}" "--tables_dir=${tables_dir}")"
  local wrap_cmd
  wrap_cmd="$(shell_join bash -lc "${inner_cmd}")"
  local cmd=(
    sbatch
    --parsable
    "--job-name=o2mw_ss_collect_${dataset}"
    "--dependency=afterok:${dependency}"
    --cpus-per-task=1
    "--mem=${MULTI_WARMUP_SEED_SPACE_MEM}"
    "--qos=${MULTI_WARMUP_SEED_SPACE_QOS}"
    "--time=${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}"
    "--output=${LOG_ROOT}/o2mw_seedspace_collect_${dataset}_%A.out"
    "--error=${LOG_ROOT}/o2mw_seedspace_collect_${dataset}_%A.err"
    "--wrap=${wrap_cmd}"
  )
  if truthy "${DRY_RUN}"; then
    print_command "Submit ${dataset} seed-space collector" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_${dataset}_SEED_SPACE_COLLECT"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted ${dataset} seed-space collector job: ${LAST_JOB_ID}"
  fi
}

submit_multi_warmup_landscape_prepare_job() {
  local prepare_prefix="$1"
  local dependency="$2"
  local multi_root="${OUT_ROOT}/${prepare_prefix}"
  local seed_space_root="${multi_root}/landscape_seed_space"
  local seed_table_dir="${multi_root}/joint_primary_clusters/SeedParameterTables"
  local invivo_best_csv="${seed_table_dir}/invivo_best_params_by_seed.csv"
  local invivo_initial_csv="${seed_table_dir}/invivo_deoptim_initial_population.csv"
  local invitro_best_csv="${seed_table_dir}/invitro_best_params_by_seed.csv"
  local invitro_initial_csv="${seed_table_dir}/invitro_deoptim_initial_population.csv"
  local invivo_tasks="${seed_space_root}/invivo_seed_space_tasks.tsv"
  local invitro_tasks="${seed_space_root}/invitro_seed_space_tasks.tsv"
  local prepare_args=(
    Rscript "${MULTI_WARMUP_SEED_PLAN_SCRIPT}"
    "--stage=prepare_landscape"
    "--project_root=${PROJECT_ROOT}"
    "--invivo_run_dir=${INVIVO_RUN_DIR}"
    "--invitro_run_dir=${INVITRO_RUN_DIR}"
    "--out_dir=${multi_root}"
    "--tsne_seed=${MULTI_WARMUP_TSNE_SEED}"
    "--cluster_seed=${MULTI_WARMUP_CLUSTER_SEED}"
    "--n_threads=${MULTI_WARMUP_N_THREADS}"
    "--render_figures=FALSE"
    "--invivo_best_csv=${invivo_best_csv}"
    "--invivo_initial_csv=${invivo_initial_csv}"
    "--invitro_best_csv=${invitro_best_csv}"
    "--invitro_initial_csv=${invitro_initial_csv}"
  )
  if [[ -n "${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}" ]]; then prepare_args+=("--max_seeds=${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}"); fi
  if [[ -f "${invivo_tasks}" ]]; then prepare_args+=("--invivo_seed_space_tasks=${invivo_tasks}"); fi
  if [[ -f "${invitro_tasks}" ]]; then prepare_args+=("--invitro_seed_space_tasks=${invitro_tasks}"); fi

  local inner_cmd wrap_cmd
  inner_cmd="if command -v ml >/dev/null 2>&1; then ml ${R_MODULE}; elif command -v module >/dev/null 2>&1; then module load ${R_MODULE}; fi; $(shell_join cd "${PROJECT_ROOT}") && $(shell_join "${prepare_args[@]}")"
  wrap_cmd="$(shell_join bash -lc "${inner_cmd}")"
  local cmd=(
    sbatch
    --parsable
    "--job-name=${JOINT_JOB_NAME}_mw_land"
    "--dependency=afterok:${dependency}"
    "--qos=${PREP_QOS}"
    "--time=${PREP_TIME_LIMIT}"
    "--cpus-per-task=${MULTI_WARMUP_N_THREADS}"
    "--mem=${PREP_MEM}"
    "--output=${LOG_ROOT}/o2_multi_warmup_landscape_%A.out"
    "--error=${LOG_ROOT}/o2_multi_warmup_landscape_%A.err"
    "--wrap=${wrap_cmd}"
  )
  if truthy "${DRY_RUN}"; then
    print_command "Submit multi-warm-up landscape prepare" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_MULTI_WARMUP_LANDSCAPE_PREPARE"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted multi-warm-up landscape prepare job: ${LAST_JOB_ID}"
  fi
}

submit_multi_warmup_finalize_job() {
  local dependency="$1"
  local finalize_prefix="${2:-${JOINT_RUN_PREFIX}}"
  local finalize_analysis_root="${3:-${MULTI_WARMUP_ANALYSIS_ROOT}}"
  local multi_root="${OUT_ROOT}/${finalize_prefix}"
  local finalize_log_prefix="${LOG_ROOT}/o2_multi_warmup_finalize"
  local finalize_args=(
    bash "${SELF_SCRIPT}"
    --internal_stage=multi_warmup_finalize_and_submit
    --fitting_mode=joint
    "--project_root=${PROJECT_ROOT}"
    "--config_path=${CONFIG_PATH}"
    "--out_root=${OUT_ROOT}"
    "--log_root=${LOG_ROOT}"
    "--multi_warmup_prefix=${finalize_prefix}"
    "--invivo_run_dir=${INVIVO_RUN_DIR}"
    "--invitro_run_dir=${INVITRO_RUN_DIR}"
    "--r_module=${R_MODULE}"
    "--force_extra_results=${FORCE_EXTRA_RESULTS}"
    "--parameter_table=${PARAMETER_TABLE}"
    "--fit_objects_dir=${FIT_OBJECTS_DIR}"
    "--flow_density_path=${FLOW_DENSITY_PATH}"
    "--itermax=${ITERMAX}"
    "--de_reltol=${DE_RELTOL}"
    "--de_steptol=${DE_STEPTOL}"
    "--NP=${NP}"
    "--auto_viz=${AUTO_VIZ}"
    "--joint_total_seeds=${JOINT_TOTAL_SEEDS}"
    "--joint_array_tasks=${JOINT_ARRAY_TASKS}"
    "--joint_seeds_per_task=${JOINT_SEEDS_PER_TASK}"
    "--joint_n_cores=${JOINT_N_CORES}"
    "--joint_mem=${JOINT_MEM}"
    "--joint_qos=${JOINT_QOS}"
    "--joint_time_limit=${JOINT_TIME_LIMIT}"
    "--postprocess_qos=${POSTPROCESS_QOS}"
    "--postprocess_time_limit=${POSTPROCESS_TIME_LIMIT}"
    "--postprocess_mem=${POSTPROCESS_MEM}"
    "--report_qos=${REPORT_QOS}"
    "--report_time_limit=${REPORT_TIME_LIMIT}"
    "--report_mem=${REPORT_MEM}"
    "--prep_qos=${PREP_QOS}"
    "--prep_time_limit=${PREP_TIME_LIMIT}"
    "--prep_mem=${PREP_MEM}"
    "--joint_soft_coupling_delta_params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}"
    "--joint_cluster_seed=${MULTI_WARMUP_CLUSTER_SEED}"
    "--joint_tsne_seed=${MULTI_WARMUP_TSNE_SEED}"
    "--joint_landscape_n_threads=${MULTI_WARMUP_N_THREADS}"
    "--seeds_per_pair=${MULTI_WARMUP_SEEDS_PER_PAIR}"
    "--multi_warmup_task_order=${MULTI_WARMUP_TASK_ORDER}"
    "--multi_warmup_task_status_filter=${MULTI_WARMUP_TASK_STATUS_FILTER}"
    "--skip_existing=${MULTI_WARMUP_SKIP_EXISTING}"
    "--refresh_task_status=${MULTI_WARMUP_REFRESH_TASK_STATUS}"
    "--multi_warmup_seed_space_qos=${MULTI_WARMUP_SEED_SPACE_QOS}"
    "--multi_warmup_seed_space_time_limit=${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}"
    "--multi_warmup_seed_space_mem=${MULTI_WARMUP_SEED_SPACE_MEM}"
    "--dry_run=FALSE"
  )
  if [[ -n "${finalize_analysis_root}" ]]; then
    finalize_args+=("--joint_analysis_root=${finalize_analysis_root}")
  fi
  if [[ -n "${MULTI_WARMUP_ARRAY_MAX_CONCURRENT}" ]]; then
    finalize_args+=("--array_max_concurrent=${MULTI_WARMUP_ARRAY_MAX_CONCURRENT}")
  fi
  if [[ -n "${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}" ]]; then
    finalize_args+=("--joint_landscape_max_seeds=${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}")
  fi
  if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then
    finalize_args+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then
    finalize_args+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_WELSCH_C}" ]]; then
    finalize_args+=("--joint_soft_coupling_welsch_c=${JOINT_SOFT_COUPLING_WELSCH_C}")
  fi
  local wrap_cmd
  wrap_cmd="$(shell_join bash -lc "$(shell_join "${finalize_args[@]}")")"
  local cmd=(
    sbatch
    --parsable
    "--job-name=${JOINT_JOB_NAME}_mw_fin"
    "--qos=${PREP_QOS}"
    "--time=${PREP_TIME_LIMIT}"
    "--cpus-per-task=${MULTI_WARMUP_N_THREADS}"
    "--mem=${PREP_MEM}"
    "--output=${finalize_log_prefix}_%A.out"
    "--error=${finalize_log_prefix}_%A.err"
  )
  if [[ -n "${dependency}" ]]; then
    cmd+=("--dependency=afterok:${dependency}")
  fi
  cmd+=("--wrap=${wrap_cmd}")
  if truthy "${DRY_RUN}"; then
    print_command "Submit multi-warm-up finalizer" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_MULTI_WARMUP_FINALIZER_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted multi-warm-up finalizer job: ${LAST_JOB_ID}"
  fi
  o2sd_prov_append "${multi_root}" slurm finalizer_job_id "${LAST_JOB_ID}"
}

run_landscape_seed_space_controller_stage() {
  INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
  INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
  load_r_module

  local multi_root="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
  local seed_space_root="${multi_root}/landscape_seed_space"
  local seed_table_dir="${multi_root}/joint_primary_clusters/SeedParameterTables"
  local progress_log="${multi_root}/multi_warmup_progress.log"
  local jobs_tsv="${multi_root}/landscape_seed_space_jobs.tsv"
  mkdir -p "${multi_root}" "${seed_space_root}" "${seed_table_dir}"
  : > "${progress_log}"
  printf "job_type\tdataset\tjob_id\tdependency\ttasks_tsv\tqos\twalltime\n" > "${jobs_tsv}"
  o2sd_prov_write_standard "${multi_root}" "${SELF_SCRIPT}" "${SUBMIT_COMMAND_TEXT:-}"
  o2sd_prov_write_many "${multi_root}" \
    fit fitting_mode "joint" \
    fit run_prefix "${JOINT_RUN_PREFIX}" \
    joint pair_method "invivo_primary_clusters_to_global_invitro_best" \
    multi_warmup invivo_run_dir "${INVIVO_RUN_DIR}" \
    multi_warmup invitro_run_dir "${INVITRO_RUN_DIR}" \
    multi_warmup seed_space_root "${seed_space_root}" \
    multi_warmup seed_table_dir "${seed_table_dir}" \
    multi_warmup n_threads "${MULTI_WARMUP_N_THREADS}"

  local invivo_tasks="${seed_space_root}/invivo_seed_space_tasks.tsv"
  local invitro_tasks="${seed_space_root}/invitro_seed_space_tasks.tsv"
  local invivo_param_table="${PROJECT_ROOT}/oxygen/data/O2_supply_demand/parameter_table_O2.csv"
  local invitro_param_table="${PARAMETER_TABLE}"
  local -a invivo_task_cmd=(
    Rscript "${MULTI_WARMUP_SEED_PLAN_SCRIPT}"
    "--stage=build_seed_space_tasks"
    "--project_root=${PROJECT_ROOT}"
    "--dataset=invivo"
    "--input_dir=${INVIVO_RUN_DIR}"
    "--out_dir=${seed_space_root}"
    "--tables_dir=${seed_table_dir}"
    "--parameter_table_fallback=${invivo_param_table}"
  )
  local -a invitro_task_cmd=(
    Rscript "${MULTI_WARMUP_SEED_PLAN_SCRIPT}"
    "--stage=build_seed_space_tasks"
    "--project_root=${PROJECT_ROOT}"
    "--dataset=invitro"
    "--input_dir=${INVITRO_RUN_DIR}"
    "--out_dir=${seed_space_root}"
    "--tables_dir=${seed_table_dir}"
    "--parameter_table_fallback=${invitro_param_table}"
  )
  if [[ -n "${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}" ]]; then
    invivo_task_cmd+=("--max_seeds=${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}")
    invitro_task_cmd+=("--max_seeds=${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}")
  fi
  print_command "Build in vivo seed-space task list" "${invivo_task_cmd[@]}" | tee -a "${progress_log}"
  if ! truthy "${DRY_RUN}"; then "${invivo_task_cmd[@]}" 2>&1 | tee -a "${progress_log}"; fi
  print_command "Build in vitro seed-space task list" "${invitro_task_cmd[@]}" | tee -a "${progress_log}"
  if ! truthy "${DRY_RUN}"; then "${invitro_task_cmd[@]}" 2>&1 | tee -a "${progress_log}"; fi

  if ! truthy "${DRY_RUN}" && [[ ! -f "${invivo_tasks}" || ! -f "${invitro_tasks}" ]]; then
    echo "Missing generated seed-space task list(s)." >&2
    exit 1
  fi

  submit_landscape_seed_space_array "invivo" "${invivo_tasks}"
  local invivo_array_job="${LAST_JOB_ID}"
  printf "array\tinvivo\t%s\t\t%s\t%s\t%s\n" "${invivo_array_job}" "${invivo_tasks}" "${MULTI_WARMUP_SEED_SPACE_QOS}" "${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}" >> "${jobs_tsv}"
  submit_landscape_seed_space_array "invitro" "${invitro_tasks}"
  local invitro_array_job="${LAST_JOB_ID}"
  printf "array\tinvitro\t%s\t\t%s\t%s\t%s\n" "${invitro_array_job}" "${invitro_tasks}" "${MULTI_WARMUP_SEED_SPACE_QOS}" "${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}" >> "${jobs_tsv}"

  submit_landscape_seed_space_collector "invivo" "${invivo_tasks}" "${seed_table_dir}" "${invivo_array_job}"
  local invivo_collect_job="${LAST_JOB_ID}"
  printf "collect\tinvivo\t%s\tafterok:%s\t%s\t%s\t%s\n" "${invivo_collect_job}" "${invivo_array_job}" "${invivo_tasks}" "${MULTI_WARMUP_SEED_SPACE_QOS}" "${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}" >> "${jobs_tsv}"
  submit_landscape_seed_space_collector "invitro" "${invitro_tasks}" "${seed_table_dir}" "${invitro_array_job}"
  local invitro_collect_job="${LAST_JOB_ID}"
  printf "collect\tinvitro\t%s\tafterok:%s\t%s\t%s\t%s\n" "${invitro_collect_job}" "${invitro_array_job}" "${invitro_tasks}" "${MULTI_WARMUP_SEED_SPACE_QOS}" "${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}" >> "${jobs_tsv}"

  local collect_dependency="${invivo_collect_job%%;*}:${invitro_collect_job%%;*}"
  submit_multi_warmup_landscape_prepare_job "${JOINT_RUN_PREFIX}" "${collect_dependency}"
  local landscape_prepare_job="${LAST_JOB_ID}"
  printf "prepare_landscape\tALL\t%s\tafterok:%s:%s\t\t%s\t%s\n" "${landscape_prepare_job}" "${invivo_collect_job}" "${invitro_collect_job}" "${PREP_QOS}" "${PREP_TIME_LIMIT}" >> "${jobs_tsv}"

  local shared_analysis_root="${multi_root}/joint_primary_clusters"
  local finalize_dependency="${landscape_prepare_job%%;*}"
  submit_multi_warmup_finalize_job "${finalize_dependency}" "${JOINT_RUN_PREFIX}" "${shared_analysis_root}"
  printf "finalize\tALL\t%s\tafterok:%s\t\t%s\t%s\n" "${LAST_JOB_ID}" "${finalize_dependency}" "${PREP_QOS}" "${PREP_TIME_LIMIT}" >> "${jobs_tsv}"
  echo "Landscape seed-space controller submitted dependent jobs."
  echo "  jobs_tsv: ${jobs_tsv}"
}

run_multi_warmup_controller_stage() {
  run_landscape_seed_space_controller_stage
}

submit_multi_warmup_controller_job() {
  local dependency="$1"
  local multi_root="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
  local controller_log_prefix="${LOG_ROOT}/o2_multi_warmup_submit"
  local controller_args=(
    bash "${SELF_SCRIPT}"
    --internal_stage=multi_warmup_prepare_and_submit
    --fitting_mode=joint
    "--project_root=${PROJECT_ROOT}"
    "--config_path=${CONFIG_PATH}"
    "--out_root=${OUT_ROOT}"
    "--log_root=${LOG_ROOT}"
    "--multi_warmup_prefix=${JOINT_RUN_PREFIX}"
    "--invivo_run_dir=${INVIVO_RUN_DIR}"
    "--invitro_run_dir=${INVITRO_RUN_DIR}"
    "--r_module=${R_MODULE}"
    "--force_extra_results=${FORCE_EXTRA_RESULTS}"
    "--parameter_table=${PARAMETER_TABLE}"
    "--fit_objects_dir=${FIT_OBJECTS_DIR}"
    "--flow_density_path=${FLOW_DENSITY_PATH}"
    "--itermax=${ITERMAX}"
    "--de_reltol=${DE_RELTOL}"
    "--de_steptol=${DE_STEPTOL}"
    "--NP=${NP}"
    "--auto_viz=${AUTO_VIZ}"
    "--joint_total_seeds=${JOINT_TOTAL_SEEDS}"
    "--joint_array_tasks=${JOINT_ARRAY_TASKS}"
    "--joint_seeds_per_task=${JOINT_SEEDS_PER_TASK}"
    "--joint_n_cores=${JOINT_N_CORES}"
    "--joint_mem=${JOINT_MEM}"
    "--joint_qos=${JOINT_QOS}"
    "--joint_time_limit=${JOINT_TIME_LIMIT}"
    "--postprocess_qos=${POSTPROCESS_QOS}"
    "--postprocess_time_limit=${POSTPROCESS_TIME_LIMIT}"
    "--postprocess_mem=${POSTPROCESS_MEM}"
    "--report_qos=${REPORT_QOS}"
    "--report_time_limit=${REPORT_TIME_LIMIT}"
    "--report_mem=${REPORT_MEM}"
    "--prep_qos=${PREP_QOS}"
    "--prep_time_limit=${PREP_TIME_LIMIT}"
    "--prep_mem=${PREP_MEM}"
    "--joint_soft_coupling_delta_params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}"
    "--joint_cluster_seed=${MULTI_WARMUP_CLUSTER_SEED}"
    "--joint_tsne_seed=${MULTI_WARMUP_TSNE_SEED}"
    "--joint_landscape_n_threads=${MULTI_WARMUP_N_THREADS}"
    "--seeds_per_pair=${MULTI_WARMUP_SEEDS_PER_PAIR}"
    "--multi_warmup_task_order=${MULTI_WARMUP_TASK_ORDER}"
    "--multi_warmup_task_status_filter=${MULTI_WARMUP_TASK_STATUS_FILTER}"
    "--skip_existing=${MULTI_WARMUP_SKIP_EXISTING}"
    "--refresh_task_status=${MULTI_WARMUP_REFRESH_TASK_STATUS}"
    "--multi_warmup_seed_space_qos=${MULTI_WARMUP_SEED_SPACE_QOS}"
    "--multi_warmup_seed_space_time_limit=${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}"
    "--multi_warmup_seed_space_mem=${MULTI_WARMUP_SEED_SPACE_MEM}"
    "--dry_run=FALSE"
  )
  if [[ -n "${MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT}" ]]; then
    controller_args+=("--multi_warmup_seed_space_array_max_concurrent=${MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT}")
  fi
  if [[ -n "${MULTI_WARMUP_ARRAY_MAX_CONCURRENT}" ]]; then
    controller_args+=("--array_max_concurrent=${MULTI_WARMUP_ARRAY_MAX_CONCURRENT}")
  fi
  if [[ -n "${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}" ]]; then
    controller_args+=("--joint_landscape_max_seeds=${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}")
  fi
  if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then
    controller_args+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then
    controller_args+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}")
  fi
  if [[ -n "${JOINT_SOFT_COUPLING_WELSCH_C}" ]]; then
    controller_args+=("--joint_soft_coupling_welsch_c=${JOINT_SOFT_COUPLING_WELSCH_C}")
  fi

  local wrap_cmd
  wrap_cmd="$(shell_join bash -lc "$(shell_join "${controller_args[@]}")")"
  local cmd=(
    sbatch
    --parsable
    "--job-name=${JOINT_JOB_NAME}_mw"
    "--qos=${PREP_QOS}"
    "--time=${PREP_TIME_LIMIT}"
    --cpus-per-task=1
    "--mem=${PREP_MEM}"
    "--output=${controller_log_prefix}_%A.out"
    "--error=${controller_log_prefix}_%A.err"
  )
  if [[ -n "${dependency}" ]]; then
    cmd+=("--dependency=afterok:${dependency}")
  fi
  cmd+=("--wrap=${wrap_cmd}")

  if truthy "${DRY_RUN}"; then
    print_command "Submit multi-warm-up controller" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_MULTI_WARMUP_CONTROLLER_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted multi-warm-up controller job: ${LAST_JOB_ID}"
  fi
  o2sd_prov_write_standard "${multi_root}" "${SELF_SCRIPT}" "${SUBMIT_COMMAND_TEXT:-}"
  o2sd_prov_record_sbatch "${multi_root}" "$(o2sd_prov_shell_join "${cmd[@]}")" "${LAST_JOB_ID}"
  o2sd_prov_write_many "${multi_root}" \
    fit fitting_mode "joint" \
    fit run_prefix "${JOINT_RUN_PREFIX}" \
    joint invivo_run_dir "${INVIVO_RUN_DIR}" \
    joint invitro_run_dir "${INVITRO_RUN_DIR}" \
    joint pair_method "invivo_primary_clusters_to_global_invitro_best" \
    joint reduction "tsne" \
    joint seeds_per_pair "${MULTI_WARMUP_SEEDS_PER_PAIR}" \
    joint joint_soft_coupling_sigma_default "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" \
    joint joint_soft_coupling_welsch_c "${JOINT_SOFT_COUPLING_WELSCH_C}" \
    slurm controller_job_id "${LAST_JOB_ID}" \
    slurm qos "${PREP_QOS}" \
    slurm walltime "${PREP_TIME_LIMIT}" \
    slurm memory "${PREP_MEM}" \
    slurm cpus "1"
}

submit_multi_warmup_pipeline() {
  local dependency="${1:-}"
  INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
  INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
  submit_multi_warmup_controller_job "${dependency}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HPC_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORKFLOW_ROOT="$(cd "${HPC_ROOT}/.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_INVIVO_RUN_PREFIX="fit_invivo_O2_buffering_500seed"
DEFAULT_INVITRO_RUN_PREFIX="fit_invitro_O2_buffering_500seed"
DEFAULT_JOINT_RUN_PREFIX="fit_joint_O2_buffering_500seed"
DEFAULT_JOINT_JOB_NAME="o2_joint_B"
DEFAULT_TOTAL_SEEDS="500"
DEFAULT_SEEDS_PER_TASK="1"
DEFAULT_N_CORES="22"
DEFAULT_MEM="16G"
DEFAULT_AUTO_VIZ="TRUE"
DEFAULT_INVIVO_QOS="xlarge"
DEFAULT_INVIVO_TIME_LIMIT="12:00:00"
DEFAULT_INVITRO_QOS="xxlarge"
DEFAULT_INVITRO_TIME_LIMIT="12:00:00"
DEFAULT_JOINT_QOS="xxlarge"
DEFAULT_JOINT_TIME_LIMIT="12:00:00"
DEFAULT_POSTPROCESS_QOS="small"
DEFAULT_POSTPROCESS_TIME_LIMIT="4:00:00"
DEFAULT_POSTPROCESS_MEM="8G"
DEFAULT_ITERMAX="1000"
DEFAULT_ITERMAX_MAX="1000"
DEFAULT_DE_RELTOL="1e-4"
DEFAULT_DE_STEPTOL="25"
DEFAULT_NP="80"
DEFAULT_FORCE_EXTRA_RESULTS="FALSE"
DEFAULT_DRY_RUN="FALSE"
DEFAULT_INTERNAL_STAGE="submit"
DEFAULT_PREP_QOS="small"
DEFAULT_PREP_TIME_LIMIT="2:00:00"
DEFAULT_PREP_MEM="8G"
DEFAULT_REPORT_QOS="small"
DEFAULT_REPORT_TIME_LIMIT="4:00:00"
DEFAULT_REPORT_MEM="8G"
DEFAULT_JOINT_WARMUP_SIGMAN=""
DEFAULT_JOINT_SOFT_COUPLING_SIGMA_DEFAULT=""
DEFAULT_JOINT_SOFT_COUPLING_WELSCH_C=""
DEFAULT_JOINT_SOFT_COUPLING_DELTA_PARAMS="default"
DEFAULT_MULTI_WARMUP_LANDSCAPE_MAX_SEEDS=""
DEFAULT_MULTI_WARMUP_N_THREADS="1"
DEFAULT_MULTI_WARMUP_SEED_SPACE_QOS="small"
DEFAULT_MULTI_WARMUP_SEED_SPACE_TIME_LIMIT="1:00:00"
DEFAULT_MULTI_WARMUP_SEED_SPACE_MEM="2G"
DEFAULT_MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT=""
DEFAULT_MULTI_WARMUP_CLUSTER_SEED="123"
DEFAULT_MULTI_WARMUP_TSNE_SEED="123"
DEFAULT_MULTI_WARMUP_ANALYSIS_ROOT=""
DEFAULT_MULTI_WARMUP_SEEDS_PER_PAIR=""
DEFAULT_MULTI_WARMUP_ARRAY_MAX_CONCURRENT=""
DEFAULT_MULTI_WARMUP_TASK_ORDER="round_robin"
DEFAULT_MULTI_WARMUP_TASK_STATUS_FILTER="all"
DEFAULT_MULTI_WARMUP_SKIP_EXISTING="TRUE"
DEFAULT_MULTI_WARMUP_REFRESH_TASK_STATUS="TRUE"

FITTING_MODE="${FITTING_MODE:-}"
INTERNAL_STAGE="${INTERNAL_STAGE:-}"
PROJECT_ROOT="${PROJECT_ROOT:-}"
R_MODULE="${R_MODULE:-}"
CONFIG_PATH="${CONFIG_PATH:-}"
OUT_ROOT="${OUT_ROOT:-}"
INVIVO_RUN_PREFIX="${INVIVO_RUN_PREFIX:-}"
INVITRO_RUN_PREFIX="${INVITRO_RUN_PREFIX:-}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-}"
JOINT_JOB_NAME="${JOINT_JOB_NAME:-}"
JOINT_DEPENDENCY="${JOINT_DEPENDENCY:-}"
INVIVO_TOTAL_SEEDS="${INVIVO_TOTAL_SEEDS:-}"
INVITRO_TOTAL_SEEDS="${INVITRO_TOTAL_SEEDS:-}"
JOINT_TOTAL_SEEDS="${JOINT_TOTAL_SEEDS:-}"
INVIVO_SEEDS_PER_TASK="${INVIVO_SEEDS_PER_TASK:-}"
INVITRO_SEEDS_PER_TASK="${INVITRO_SEEDS_PER_TASK:-}"
JOINT_SEEDS_PER_TASK="${JOINT_SEEDS_PER_TASK:-}"
INVIVO_ARRAY_TASKS="${INVIVO_ARRAY_TASKS:-}"
INVITRO_ARRAY_TASKS="${INVITRO_ARRAY_TASKS:-}"
JOINT_ARRAY_TASKS="${JOINT_ARRAY_TASKS:-}"
N_CORES="${N_CORES:-}"
INVIVO_N_CORES="${INVIVO_N_CORES:-}"
INVITRO_N_CORES="${INVITRO_N_CORES:-}"
JOINT_N_CORES="${JOINT_N_CORES:-}"
MEM="${MEM:-}"
INVIVO_MEM="${INVIVO_MEM:-}"
INVITRO_MEM="${INVITRO_MEM:-}"
JOINT_MEM="${JOINT_MEM:-}"
INVIVO_QOS="${INVIVO_QOS:-}"
INVIVO_TIME_LIMIT="${INVIVO_TIME_LIMIT:-}"
INVITRO_QOS="${INVITRO_QOS:-}"
INVITRO_TIME_LIMIT="${INVITRO_TIME_LIMIT:-}"
JOINT_QOS="${JOINT_QOS:-}"
JOINT_TIME_LIMIT="${JOINT_TIME_LIMIT:-}"
POSTPROCESS_QOS="${POSTPROCESS_QOS:-}"
POSTPROCESS_TIME_LIMIT="${POSTPROCESS_TIME_LIMIT:-}"
POSTPROCESS_MEM="${POSTPROCESS_MEM:-}"
PARAMETER_TABLE="${PARAMETER_TABLE:-}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-}"
ITERMAX="${ITERMAX:-}"
ITERMAX_MAX="${ITERMAX_MAX:-}"
DE_RELTOL="${DE_RELTOL:-}"
DE_STEPTOL="${DE_STEPTOL:-}"
NP="${NP:-}"
AUTO_VIZ="${AUTO_VIZ:-}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-}"
MULTI_WARMUP_LANDSCAPE_MAX_SEEDS="${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS:-}"
MULTI_WARMUP_N_THREADS="${MULTI_WARMUP_N_THREADS:-}"
MULTI_WARMUP_SEED_SPACE_QOS="${MULTI_WARMUP_SEED_SPACE_QOS:-}"
MULTI_WARMUP_SEED_SPACE_TIME_LIMIT="${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT:-}"
MULTI_WARMUP_SEED_SPACE_MEM="${MULTI_WARMUP_SEED_SPACE_MEM:-}"
MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT="${MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT:-}"
MULTI_WARMUP_CLUSTER_SEED="${MULTI_WARMUP_CLUSTER_SEED:-}"
MULTI_WARMUP_TSNE_SEED="${MULTI_WARMUP_TSNE_SEED:-}"
MULTI_WARMUP_ANALYSIS_ROOT="${MULTI_WARMUP_ANALYSIS_ROOT:-}"
MULTI_WARMUP_SEEDS_PER_PAIR="${MULTI_WARMUP_SEEDS_PER_PAIR:-}"
MULTI_WARMUP_ARRAY_MAX_CONCURRENT="${MULTI_WARMUP_ARRAY_MAX_CONCURRENT:-}"
MULTI_WARMUP_TASK_ORDER="${MULTI_WARMUP_TASK_ORDER:-}"
MULTI_WARMUP_TASK_STATUS_FILTER="${MULTI_WARMUP_TASK_STATUS_FILTER:-}"
MULTI_WARMUP_SKIP_EXISTING="${MULTI_WARMUP_SKIP_EXISTING:-}"
MULTI_WARMUP_REFRESH_TASK_STATUS="${MULTI_WARMUP_REFRESH_TASK_STATUS:-}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-}"
DRY_RUN="${DRY_RUN:-}"
PREP_QOS="${PREP_QOS:-}"
PREP_TIME_LIMIT="${PREP_TIME_LIMIT:-}"
PREP_MEM="${PREP_MEM:-}"
REPORT_QOS="${REPORT_QOS:-}"
REPORT_TIME_LIMIT="${REPORT_TIME_LIMIT:-}"
REPORT_MEM="${REPORT_MEM:-}"
LOG_ROOT="${LOG_ROOT:-}"

parse_args "$@"

INTERNAL_STAGE="${INTERNAL_STAGE:-${DEFAULT_INTERNAL_STAGE}}"
PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
INVIVO_RUN_PREFIX="${INVIVO_RUN_PREFIX:-${DEFAULT_INVIVO_RUN_PREFIX}}"
INVITRO_RUN_PREFIX="${INVITRO_RUN_PREFIX:-${DEFAULT_INVITRO_RUN_PREFIX}}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-${DEFAULT_JOINT_RUN_PREFIX}}"
JOINT_JOB_NAME="${JOINT_JOB_NAME:-${DEFAULT_JOINT_JOB_NAME}}"
JOINT_DEPENDENCY="${JOINT_DEPENDENCY:-}"
INVIVO_TOTAL_SEEDS="${INVIVO_TOTAL_SEEDS:-${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}}"
INVITRO_TOTAL_SEEDS="${INVITRO_TOTAL_SEEDS:-${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}}"
JOINT_TOTAL_SEEDS="${JOINT_TOTAL_SEEDS:-${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}}"
INVIVO_SEEDS_PER_TASK="${INVIVO_SEEDS_PER_TASK:-${SEEDS_PER_TASK:-${DEFAULT_SEEDS_PER_TASK}}}"
INVITRO_SEEDS_PER_TASK="${INVITRO_SEEDS_PER_TASK:-${SEEDS_PER_TASK:-${DEFAULT_SEEDS_PER_TASK}}}"
JOINT_SEEDS_PER_TASK="${JOINT_SEEDS_PER_TASK:-${SEEDS_PER_TASK:-${DEFAULT_SEEDS_PER_TASK}}}"
INVIVO_ARRAY_TASKS="${INVIVO_ARRAY_TASKS:-${ARRAY_TASKS:-${INVIVO_TOTAL_SEEDS}}}"
INVITRO_ARRAY_TASKS="${INVITRO_ARRAY_TASKS:-${ARRAY_TASKS:-${INVITRO_TOTAL_SEEDS}}}"
JOINT_ARRAY_TASKS="${JOINT_ARRAY_TASKS:-${ARRAY_TASKS:-${JOINT_TOTAL_SEEDS}}}"
N_CORES="${N_CORES:-${DEFAULT_N_CORES}}"
INVIVO_N_CORES="${INVIVO_N_CORES:-${N_CORES}}"
INVITRO_N_CORES="${INVITRO_N_CORES:-${N_CORES}}"
JOINT_N_CORES="${JOINT_N_CORES:-${N_CORES}}"
MEM="${MEM:-${DEFAULT_MEM}}"
INVIVO_MEM="${INVIVO_MEM:-${MEM}}"
INVITRO_MEM="${INVITRO_MEM:-${MEM}}"
JOINT_MEM="${JOINT_MEM:-${MEM}}"
INVIVO_QOS="${INVIVO_QOS:-${DEFAULT_INVIVO_QOS}}"
INVIVO_TIME_LIMIT="${INVIVO_TIME_LIMIT:-${DEFAULT_INVIVO_TIME_LIMIT}}"
INVITRO_QOS="${INVITRO_QOS:-${DEFAULT_INVITRO_QOS}}"
INVITRO_TIME_LIMIT="${INVITRO_TIME_LIMIT:-${DEFAULT_INVITRO_TIME_LIMIT}}"
JOINT_QOS="${JOINT_QOS:-${DEFAULT_JOINT_QOS}}"
JOINT_TIME_LIMIT="${JOINT_TIME_LIMIT:-${DEFAULT_JOINT_TIME_LIMIT}}"
POSTPROCESS_QOS="${POSTPROCESS_QOS:-${DEFAULT_POSTPROCESS_QOS}}"
POSTPROCESS_TIME_LIMIT="${POSTPROCESS_TIME_LIMIT:-${DEFAULT_POSTPROCESS_TIME_LIMIT}}"
POSTPROCESS_MEM="${POSTPROCESS_MEM:-${DEFAULT_POSTPROCESS_MEM}}"
PARAMETER_TABLE="${PARAMETER_TABLE:-}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-}"
ITERMAX="${ITERMAX:-${DEFAULT_ITERMAX}}"
ITERMAX_MAX="${ITERMAX_MAX:-${DEFAULT_ITERMAX_MAX}}"
DE_RELTOL="${DE_RELTOL:-${DEFAULT_DE_RELTOL}}"
DE_STEPTOL="${DE_STEPTOL:-${DEFAULT_DE_STEPTOL}}"
NP="${NP:-${DEFAULT_NP}}"
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-${DEFAULT_JOINT_WARMUP_SIGMAN}}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-${DEFAULT_JOINT_SOFT_COUPLING_SIGMA_DEFAULT}}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-${DEFAULT_JOINT_SOFT_COUPLING_WELSCH_C}}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-${DEFAULT_JOINT_SOFT_COUPLING_DELTA_PARAMS}}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-${DEFAULT_FORCE_EXTRA_RESULTS}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"
PREP_QOS="${PREP_QOS:-${DEFAULT_PREP_QOS}}"
PREP_TIME_LIMIT="${PREP_TIME_LIMIT:-${DEFAULT_PREP_TIME_LIMIT}}"
PREP_MEM="${PREP_MEM:-${DEFAULT_PREP_MEM}}"
REPORT_QOS="${REPORT_QOS:-${DEFAULT_REPORT_QOS}}"
REPORT_TIME_LIMIT="${REPORT_TIME_LIMIT:-${DEFAULT_REPORT_TIME_LIMIT}}"
REPORT_MEM="${REPORT_MEM:-${DEFAULT_REPORT_MEM}}"
MULTI_WARMUP_LANDSCAPE_MAX_SEEDS="${MULTI_WARMUP_LANDSCAPE_MAX_SEEDS:-${DEFAULT_MULTI_WARMUP_LANDSCAPE_MAX_SEEDS}}"
MULTI_WARMUP_N_THREADS="${MULTI_WARMUP_N_THREADS:-${DEFAULT_MULTI_WARMUP_N_THREADS}}"
MULTI_WARMUP_SEED_SPACE_QOS="${MULTI_WARMUP_SEED_SPACE_QOS:-${DEFAULT_MULTI_WARMUP_SEED_SPACE_QOS}}"
MULTI_WARMUP_SEED_SPACE_TIME_LIMIT="${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT:-${DEFAULT_MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}}"
MULTI_WARMUP_SEED_SPACE_MEM="${MULTI_WARMUP_SEED_SPACE_MEM:-${DEFAULT_MULTI_WARMUP_SEED_SPACE_MEM}}"
MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT="${MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT:-${DEFAULT_MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT}}"
MULTI_WARMUP_CLUSTER_SEED="${MULTI_WARMUP_CLUSTER_SEED:-${DEFAULT_MULTI_WARMUP_CLUSTER_SEED}}"
MULTI_WARMUP_TSNE_SEED="${MULTI_WARMUP_TSNE_SEED:-${DEFAULT_MULTI_WARMUP_TSNE_SEED}}"
MULTI_WARMUP_ANALYSIS_ROOT="${MULTI_WARMUP_ANALYSIS_ROOT:-${DEFAULT_MULTI_WARMUP_ANALYSIS_ROOT}}"
MULTI_WARMUP_SEEDS_PER_PAIR="${MULTI_WARMUP_SEEDS_PER_PAIR:-${DEFAULT_MULTI_WARMUP_SEEDS_PER_PAIR}}"
MULTI_WARMUP_ARRAY_MAX_CONCURRENT="${MULTI_WARMUP_ARRAY_MAX_CONCURRENT:-${DEFAULT_MULTI_WARMUP_ARRAY_MAX_CONCURRENT}}"
MULTI_WARMUP_TASK_ORDER="${MULTI_WARMUP_TASK_ORDER:-${DEFAULT_MULTI_WARMUP_TASK_ORDER}}"
MULTI_WARMUP_TASK_STATUS_FILTER="${MULTI_WARMUP_TASK_STATUS_FILTER:-${DEFAULT_MULTI_WARMUP_TASK_STATUS_FILTER}}"
MULTI_WARMUP_SKIP_EXISTING="${MULTI_WARMUP_SKIP_EXISTING:-${DEFAULT_MULTI_WARMUP_SKIP_EXISTING}}"
MULTI_WARMUP_REFRESH_TASK_STATUS="${MULTI_WARMUP_REFRESH_TASK_STATUS:-${DEFAULT_MULTI_WARMUP_REFRESH_TASK_STATUS}}"

FITTING_MODE="$(normalize_fitting_mode "${FITTING_MODE}")"
if [[ -z "${FITTING_MODE}" ]]; then
  echo "--fitting_mode must be one of invivo, invitro, joint, or all." >&2
  usage >&2
  exit 2
fi

if [[ "${FITTING_MODE}" == "joint" ]]; then
  if is_null_value "${INVIVO_RUN_DIR}" || is_null_value "${INVITRO_RUN_DIR}"; then
    echo "--fitting_mode=joint requires both --invivo_run_dir and --invitro_run_dir." >&2
    exit 2
  fi
fi
if [[ "${FITTING_MODE}" == "joint" || "${FITTING_MODE}" == "all" ]]; then
  if is_null_value "${MULTI_WARMUP_SEEDS_PER_PAIR}"; then
    MULTI_WARMUP_SEEDS_PER_PAIR="${JOINT_TOTAL_SEEDS}"
  fi
  require_positive_int "MULTI_WARMUP_SEEDS_PER_PAIR" "${MULTI_WARMUP_SEEDS_PER_PAIR}"
  JOINT_TOTAL_SEEDS="${MULTI_WARMUP_SEEDS_PER_PAIR}"
  JOINT_ARRAY_TASKS="${MULTI_WARMUP_SEEDS_PER_PAIR}"
  JOINT_SEEDS_PER_TASK="1"
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
if [[ -z "${CONFIG_PATH}" ]]; then
  CONFIG_PATH="${PROJECT_ROOT}/oxygen/config/O2_supply_demand.yaml"
fi
if [[ -z "${OUT_ROOT}" ]]; then
  OUT_ROOT="${PROJECT_ROOT}/oxygen/results"
fi
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
mkdir -p "${OUT_ROOT}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"
if [[ "${FITTING_MODE}" == "all" ]]; then
  INVIVO_RUN_DIR="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
  INVITRO_RUN_DIR="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
fi
if [[ -z "${LOG_ROOT}" ]]; then
  LOG_ROOT="${OUT_ROOT}/log"
fi
mkdir -p "${LOG_ROOT}"
LOG_ROOT="$(cd "${LOG_ROOT}" && pwd)"

HPC_DIR="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/hpc"
HPC_SUBMIT_DIR="${HPC_DIR}/submit"
HPC_ARRAY_WORKER_DIR="${HPC_DIR}/array_workers"
HPC_POSTPROCESS_DIR="${HPC_DIR}/postprocess"
SELF_SCRIPT="${SELF_SCRIPT:-${HPC_SUBMIT_DIR}/submit_o2_fit.sh}"
if [[ ! -f "${SELF_SCRIPT}" ]]; then
  SELF_SCRIPT_CANDIDATE="${HPC_SUBMIT_DIR}/$(basename "${BASH_SOURCE[0]}")"
  if [[ -f "${SELF_SCRIPT_CANDIDATE}" ]]; then
    SELF_SCRIPT="${SELF_SCRIPT_CANDIDATE}"
  fi
fi
PROVENANCE_HELPER="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_shell_utils.sh"
if [[ -f "${PROVENANCE_HELPER}" ]]; then
  # shellcheck source=/dev/null
  source "${PROVENANCE_HELPER}"
else
  echo "Missing provenance helper: ${PROVENANCE_HELPER}" >&2
  exit 1
fi
SUBMIT_COMMAND_TEXT="$(o2sd_prov_shell_join bash "${SELF_SCRIPT}" "${ORIGINAL_SUBMIT_ARGS[@]}")"
INVIVO_SUB_SCRIPT="${INVIVO_SUB_SCRIPT:-${HPC_ARRAY_WORKER_DIR}/submit_fit_seed_array_buffering.sub}"
INVITRO_SUB_SCRIPT="${INVITRO_SUB_SCRIPT:-${HPC_ARRAY_WORKER_DIR}/submit_fit_seed_array_invitro_buffering.sub}"
MULTI_WARMUP_LANDSCAPE_SEED_PLAN_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/build_joint_primary_cluster_pairs.R"
MULTI_WARMUP_SEED_PLAN_SCRIPT="${MULTI_WARMUP_SEED_PLAN_SCRIPT:-${MULTI_WARMUP_LANDSCAPE_SEED_PLAN_SCRIPT}}"
MULTI_WARMUP_TASK_TABLE_SCRIPT="${MULTI_WARMUP_TASK_TABLE_SCRIPT:-${HPC_SUBMIT_DIR}/build_multi_warmup_task_table.R}"
MULTI_WARMUP_TASK_SUBMIT_SCRIPT="${MULTI_WARMUP_TASK_SUBMIT_SCRIPT:-${HPC_SUBMIT_DIR}/submit_multi_warmup_task_table.sh}"
LANDSCAPE_SEED_SPACE_ARRAY_SCRIPT="${LANDSCAPE_SEED_SPACE_ARRAY_SCRIPT:-${HPC_ARRAY_WORKER_DIR}/run_landscape_seed_space_task.sub}"
POSTPROCESS_SCRIPT="${POSTPROCESS_SCRIPT:-${HPC_POSTPROCESS_DIR}/postprocess_extra_results.sh}"
EXTRA_RESULTS_SCRIPT="${EXTRA_RESULTS_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/fit_results/extra_results.R}"
JOINT_WARM_START_SCRIPT="${JOINT_WARM_START_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/warm_start/make_joint_soft_coupling_parameters_table.R}"
INVIVO_RUNNER_SCRIPT="${INVIVO_RUNNER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_fit_model_O2_supply_demand_MAP.sh}"
INVITRO_RUNNER_SCRIPT="${INVITRO_RUNNER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_fit_model_O2_supply_demand_MAP.sh}"
JOINT_RUNNER_SCRIPT="${JOINT_RUNNER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh}"

if [[ -z "${PARAMETER_TABLE}" ]]; then
  PARAMETER_TABLE="${PROJECT_ROOT}/oxygen/data/O2_supply_demand/parameter_table_invitro_buffering.csv"
fi
if [[ -z "${FIT_OBJECTS_DIR}" ]]; then
  FIT_OBJECTS_DIR="${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects"
fi
if [[ -z "${FLOW_DENSITY_PATH}" ]]; then
  FLOW_DENSITY_PATH="${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv"
fi

for path in "${CONFIG_PATH}" "${SELF_SCRIPT}" "${INVIVO_SUB_SCRIPT}" "${INVITRO_SUB_SCRIPT}" \
            "${MULTI_WARMUP_SEED_PLAN_SCRIPT}" "${MULTI_WARMUP_TASK_TABLE_SCRIPT}" "${MULTI_WARMUP_TASK_SUBMIT_SCRIPT}" \
            "${LANDSCAPE_SEED_SPACE_ARRAY_SCRIPT}" \
            "${POSTPROCESS_SCRIPT}" "${EXTRA_RESULTS_SCRIPT}" \
            "${JOINT_WARM_START_SCRIPT}" \
            "${INVIVO_RUNNER_SCRIPT}" "${INVITRO_RUNNER_SCRIPT}" "${JOINT_RUNNER_SCRIPT}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
done
for name in INVIVO_TOTAL_SEEDS INVIVO_ARRAY_TASKS INVIVO_SEEDS_PER_TASK \
            INVITRO_TOTAL_SEEDS INVITRO_ARRAY_TASKS INVITRO_SEEDS_PER_TASK \
            JOINT_TOTAL_SEEDS JOINT_ARRAY_TASKS JOINT_SEEDS_PER_TASK \
            INVIVO_N_CORES INVITRO_N_CORES JOINT_N_CORES ITERMAX ITERMAX_MAX DE_STEPTOL NP MULTI_WARMUP_N_THREADS; do
  require_positive_int "${name}" "${!name}"
done
check_seed_plan INVIVO "${INVIVO_TOTAL_SEEDS}" "${INVIVO_ARRAY_TASKS}" "${INVIVO_SEEDS_PER_TASK}"
check_seed_plan INVITRO "${INVITRO_TOTAL_SEEDS}" "${INVITRO_ARRAY_TASKS}" "${INVITRO_SEEDS_PER_TASK}"
check_seed_plan JOINT "${JOINT_TOTAL_SEEDS}" "${JOINT_ARRAY_TASKS}" "${JOINT_SEEDS_PER_TASK}"


ensure_sbatch

echo "O2 submitter"
echo "  fitting_mode: ${FITTING_MODE}"
echo "  project_root: ${PROJECT_ROOT}"
echo "  out_root: ${OUT_ROOT}"
echo "  log_root: ${LOG_ROOT}"
echo "  r_module: ${R_MODULE}"
echo "  invivo resources: qos=${INVIVO_QOS}, time=${INVIVO_TIME_LIMIT}, mem=${INVIVO_MEM}, cpus=${INVIVO_N_CORES}"
echo "  invitro resources: qos=${INVITRO_QOS}, time=${INVITRO_TIME_LIMIT}, mem=${INVITRO_MEM}, cpus=${INVITRO_N_CORES}"
echo "  optimizer iterations: requested=${ITERMAX}, maximum=${ITERMAX_MAX}"
echo "  joint resources: qos=${JOINT_QOS}, time=${JOINT_TIME_LIMIT}, mem=${JOINT_MEM}, cpus=${JOINT_N_CORES}"
echo "  postprocess resources: qos=${POSTPROCESS_QOS}, time=${POSTPROCESS_TIME_LIMIT}, mem=${POSTPROCESS_MEM}"
echo "  prep resources: qos=${PREP_QOS}, time=${PREP_TIME_LIMIT}, mem=${PREP_MEM}"
if [[ "${FITTING_MODE}" == "joint" || "${FITTING_MODE}" == "all" ]]; then
  echo "  joint workflow: in-vivo primary clusters paired with one global-best in-vitro seed"
  echo "  joint source in vivo: ${INVIVO_RUN_DIR}"
  echo "  joint source in vitro: ${INVITRO_RUN_DIR}"
  echo "  joint seeds per pair: ${MULTI_WARMUP_SEEDS_PER_PAIR}"
  echo "  joint task order: ${MULTI_WARMUP_TASK_ORDER}"
  echo "  joint array_max_concurrent: ${MULTI_WARMUP_ARRAY_MAX_CONCURRENT:-none}"
  echo "  joint landscape threads: ${MULTI_WARMUP_N_THREADS}"
  echo "  joint seed-space resources: qos=${MULTI_WARMUP_SEED_SPACE_QOS}, time=${MULTI_WARMUP_SEED_SPACE_TIME_LIMIT}, mem=${MULTI_WARMUP_SEED_SPACE_MEM}, max_concurrent=${MULTI_WARMUP_SEED_SPACE_ARRAY_MAX_CONCURRENT:-none}"
fi

case "${INTERNAL_STAGE}" in
  submit|multi_warmup_prepare_and_submit|multi_warmup_finalize_and_submit) ;;
  *)
    echo "--internal_stage must be submit, multi_warmup_prepare_and_submit, or multi_warmup_finalize_and_submit, got: ${INTERNAL_STAGE}" >&2
    exit 2
    ;;
esac

if [[ "${INTERNAL_STAGE}" == "multi_warmup_finalize_and_submit" ]]; then
  run_multi_warmup_finalize_stage
  exit 0
fi
if [[ "${INTERNAL_STAGE}" == "multi_warmup_prepare_and_submit" ]]; then
  run_multi_warmup_controller_stage
  exit 0
fi

case "${FITTING_MODE}" in
  invivo)
    INVIVO_RUN_DIR="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
    submit_invivo_array
    INVIVO_JOB_ID="${LAST_JOB_ID}"
    submit_extra_results_job "o2_invivo" "${INVIVO_RUN_DIR}" "${INVIVO_JOB_ID}"
    ;;
  invitro)
    INVITRO_RUN_DIR="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
    submit_invitro_array
    INVITRO_JOB_ID="${LAST_JOB_ID}"
    submit_extra_results_job "o2_invitro" "${INVITRO_RUN_DIR}" "${INVITRO_JOB_ID}"
    ;;
  joint)
    echo "Submitting the fixed in-vivo-primary-cluster/global-in-vitro-best joint workflow."
    submit_multi_warmup_pipeline "${JOINT_DEPENDENCY}"
    ;;
  all)
    echo "Submitting the complete chain: in vivo -> in vitro -> primary clusters -> joint fitting."
    INVIVO_RUN_DIR="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
    submit_invivo_array
    INVIVO_JOB_ID="${LAST_JOB_ID%%;*}"
    submit_extra_results_job "o2_invivo" "${INVIVO_RUN_DIR}" "${INVIVO_JOB_ID}"

    INVITRO_RUN_DIR="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
    submit_invitro_array "${INVIVO_JOB_ID}"
    INVITRO_JOB_ID="${LAST_JOB_ID%%;*}"
    submit_extra_results_job "o2_invitro" "${INVITRO_RUN_DIR}" "${INVITRO_JOB_ID}"

    submit_multi_warmup_pipeline "${INVITRO_JOB_ID}"
    ;;
esac
