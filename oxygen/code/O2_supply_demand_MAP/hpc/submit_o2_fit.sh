#!/bin/bash
# Unified O2 HPC submitter for in vivo, in vitro, and joint fitting.

set -euo pipefail

ORIGINAL_SUBMIT_ARGS=("$@")

usage() {
  cat <<'EOF'
Usage:
  bash submit_o2_fit.sh --fitting_mode=invivo [options]
  bash submit_o2_fit.sh --fitting_mode=invitro [options]
  bash submit_o2_fit.sh --fitting_mode=joint --joint_fitting_mode=JOINT [options]
  bash submit_o2_fit.sh --fitting_mode=joint --joint_fitting_mode=DIRECT [options]
  bash submit_o2_fit.sh --fitting_mode=joint --joint_fitting_mode=MULTI_WARMUP [options]

Required modes:
  --fitting_mode=invivo|invitro|joint
  --joint_fitting_mode=OFF|JOINT|DIRECT|MULTI_WARMUP
  If fitting_mode=joint and joint_fitting_mode is omitted, DIRECT is used.

Joint mode behavior:
  OFF    Do not submit joint fitting. This is forced when fitting_mode is not joint.
  JOINT  Submit or reuse in vivo and in vitro single fits, run extra_results,
         select each best seed, then submit joint fitting from the selected
         single-fit anchors. Provided best_seed_dir skips that side completely;
         provided run_dir skips fitting but still selects the best seed after
         extra_results; missing sides are submitted and selected before joint.
  DIRECT Submit only the current joint fitter directly from the config.
         SINGLE is accepted as a legacy alias for DIRECT.
  MULTI_WARMUP
         Submit or reuse in vivo and in vitro source fits, run extra_results,
         then build a source ratio UMAP/cluster manifest, expand it into a
         global warm-up pair x seed task table, and submit one cross-pair
         joint array. User-specified best_seed_dir is not accepted in this mode.

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
  --select_required_files=best_params.tsv
  --invivo_objective_columns=objective
  --invitro_objective_columns=objective_total,objective

Joint options:
  --joint_run_prefix=name
  --joint_job_name=o2_joint_B
  --joint_total_seeds=500 --joint_array_tasks=500 --joint_seeds_per_task=1
  --joint_qos=xxlarge --joint_time_limit=12:00:00
  --postprocess_qos=small --postprocess_time_limit=4:00:00
  --parameter_table=/path/to/invitro_parameter_table.csv
  --fit_objects_dir=/path/to/fit_objects
  --flow_density_path=/path/to/g0g1_ploidy_density_grid.csv
  --invivo_best_seed_dir=/path/to/invivo/seed50
  --invitro_best_seed_dir=/path/to/invitro/seed350
  --joint_warmup_seed_label=invivo_seed50__invitro_seed350
  --joint_soft_coupling_sigma_default=0.65
  --joint_soft_coupling_welsch_c=0.4
  --joint_soft_coupling_parameters_table=/path/to/joint_soft_coupling_parameters_table.csv
  --joint_warmup_sigmaN=0.0304
  --joint_soft_coupling_delta_params=default|all|none|param1,param2
  --prep_qos=small --prep_time_limit=2:00:00 --prep_mem=8G
  --multi_warmup_top_n=10
  --multi_warmup_umap_seed=1
  --multi_warmup_invivo_k=auto
  --multi_warmup_invitro_anchor_ranks=1
  --multi_warmup_include_phase2=TRUE|FALSE
  --multi_warmup_phase2_invitro_anchor_ranks=auto
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

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

is_null_value() {
  local val
  val="$(echo "${1:-}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')"
  [[ -z "${val}" || "${val}" == "null" || "${val}" == "none" || "${val}" == "na" ]]
}

normalize_fitting_mode() {
  local val
  val="$(echo "${1:-}" | tr '[:upper:]' '[:lower:]')"
  val="${val// /}"
  val="${val//-/}"
  val="${val//_/}"
  case "${val}" in
    invivo) echo "invivo" ;;
    invitro) echo "invitro" ;;
    joint) echo "joint" ;;
    *) echo "" ;;
  esac
}

require_positive_int() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]] || (( value <= 0 )); then
    echo "${name} must be a positive integer, got: ${value}" >&2
    exit 2
  fi
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

print_command() {
  local label="$1"
  shift
  printf "%s:" "${label}"
  printf " %q" "$@"
  printf "\n"
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

load_r_module() {
  if command -v ml >/dev/null 2>&1; then
    ml "${R_MODULE}"
  elif command -v module >/dev/null 2>&1; then
    module load "${R_MODULE}"
  fi
}

sanitize_label() {
  local val
  val="$(printf "%s" "${1:-}" | tr -c '[:alnum:]_.-' '_' | sed 's/^_*//;s/_*$//')"
  if [[ -z "${val}" ]]; then
    val="seed"
  fi
  printf "%s" "${val}"
}

derive_joint_warmup_seed_label() {
  local invivo_label
  local invitro_label
  invivo_label="$(sanitize_label "$(basename "${INVIVO_BEST_SEED_DIR}")")"
  invitro_label="$(sanitize_label "$(basename "${INVITRO_BEST_SEED_DIR}")")"
  printf "invivo_%s__invitro_%s" "${invivo_label}" "${invitro_label}"
}

label_joint_run_prefix() {
  if truthy "${JOINT_WARMUP_ENABLE}" && ! is_null_value "${JOINT_WARMUP_SEED_LABEL}"; then
    if [[ "${JOINT_RUN_PREFIX}" != *"${JOINT_WARMUP_SEED_LABEL}"* ]]; then
      JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX}__${JOINT_WARMUP_SEED_LABEL}"
    fi
  fi
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h)
        usage
        exit 0
        ;;
      --fitting_mode=*) FITTING_MODE="${arg#*=}" ;;
      --joint_fitting_mode=*) JOINT_FITTING_MODE="${arg#*=}" ;;
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
      --invivo_best_seed_dir=*|--joint_warmup_invivo_seed_dir=*|--joint_warmup_invivo_best_seed_dir=*) INVIVO_BEST_SEED_DIR="${arg#*=}"; USER_INVIVO_BEST_SEED_DIR="TRUE" ;;
      --invitro_best_seed_dir=*|--joint_warmup_invitro_seed_dir=*|--joint_warmup_invitro_best_seed_dir=*|--joint_warmup_vitro_seed_dir=*) INVITRO_BEST_SEED_DIR="${arg#*=}"; USER_INVITRO_BEST_SEED_DIR="TRUE" ;;
      --joint_warmup_enable=*) JOINT_WARMUP_ENABLE="${arg#*=}" ;;
      --joint_warmup_seed_label=*|--joint_seed_label=*|--seed_label=*) JOINT_WARMUP_SEED_LABEL="${arg#*=}" ;;
      --joint_warmup_sigmaN=*) JOINT_WARMUP_SIGMAN="${arg#*=}" ;;
      --joint_soft_coupling_sigma_default=*) JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${arg#*=}" ;;
      --joint_soft_coupling_welsch_c=*) JOINT_SOFT_COUPLING_WELSCH_C="${arg#*=}" ;;
      --joint_soft_coupling_parameters_table=*|--joint_soft_coupling_parameters_table_path=*) JOINT_SOFT_COUPLING_PARAMETERS_TABLE="${arg#*=}" ;;
      --joint_soft_coupling_delta_params=*) JOINT_SOFT_COUPLING_DELTA_PARAMS="${arg#*=}" ;;
      --multi_warmup_top_n=*) MULTI_WARMUP_TOP_N="${arg#*=}" ;;
      --multi_warmup_umap_seed=*|--umap_seed=*) MULTI_WARMUP_UMAP_SEED="${arg#*=}" ;;
      --multi_warmup_invivo_k=*) MULTI_WARMUP_INVIVO_K="${arg#*=}" ;;
      --multi_warmup_invitro_k=*) MULTI_WARMUP_INVITRO_K="${arg#*=}" ;;
      --multi_warmup_invitro_anchor_ranks=*) MULTI_WARMUP_INVITRO_ANCHOR_RANKS="${arg#*=}" ;;
      --multi_warmup_include_phase2=*) MULTI_WARMUP_INCLUDE_PHASE2="${arg#*=}" ;;
      --multi_warmup_phase2_invitro_anchor_ranks=*) MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="${arg#*=}" ;;
      --seeds_per_pair=*|--joint_seeds_per_pair=*|--multi_warmup_seeds_per_pair=*) MULTI_WARMUP_SEEDS_PER_PAIR="${arg#*=}" ;;
      --array_max_concurrent=*|--multi_warmup_array_max_concurrent=*) MULTI_WARMUP_ARRAY_MAX_CONCURRENT="${arg#*=}" ;;
      --multi_warmup_task_order=*|--task_order=*) MULTI_WARMUP_TASK_ORDER="${arg#*=}" ;;
      --multi_warmup_task_status_filter=*|--task_status_filter=*) MULTI_WARMUP_TASK_STATUS_FILTER="${arg#*=}" ;;
      --skip_existing=*|--multi_warmup_skip_existing=*) MULTI_WARMUP_SKIP_EXISTING="${arg#*=}" ;;
      --refresh_task_status=*|--multi_warmup_refresh_task_status=*) MULTI_WARMUP_REFRESH_TASK_STATUS="${arg#*=}" ;;
      --select_required_files=*) SELECT_REQUIRED_FILES="${arg#*=}" ;;
      --invivo_objective_columns=*) INVIVO_OBJECTIVE_COLUMNS="${arg#*=}" ;;
      --invitro_objective_columns=*) INVITRO_OBJECTIVE_COLUMNS="${arg#*=}" ;;
      --prep_qos=*) PREP_QOS="${arg#*=}" ;;
      --prep_time_limit=*) PREP_TIME_LIMIT="${arg#*=}" ;;
      --prep_mem=*) PREP_MEM="${arg#*=}" ;;
      --report_qos=*) REPORT_QOS="${arg#*=}" ;;
      --report_time_limit=*) REPORT_TIME_LIMIT="${arg#*=}" ;;
      --report_mem=*) REPORT_MEM="${arg#*=}" ;;
      --itermax=*) ITERMAX="${arg#*=}" ;;
      --de_reltol=*) DE_RELTOL="${arg#*=}" ;;
      --de_steptol=*) DE_STEPTOL="${arg#*=}" ;;
      --np=*|--NP=*) NP="${arg#*=}" ;;
      --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
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
    "${INVIVO_SUB_SCRIPT}"
  )
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
    "${INVITRO_SUB_SCRIPT}"
  )
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

submit_joint_array() {
  local dependency="$1"
  local run_dir="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
  mkdir -p "${run_dir}"
  local export_arg="ALL"
  export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
  export_arg+=",RUNNER_SCRIPT=${JOINT_RUNNER_SCRIPT}"
  export_arg+=",CONFIG_PATH=${CONFIG_PATH}"
  export_arg+=",OUT_ROOT=${OUT_ROOT}"
  export_arg+=",RUN_PREFIX=${JOINT_RUN_PREFIX}"
  export_arg+=",TOTAL_SEEDS=${JOINT_TOTAL_SEEDS}"
  export_arg+=",ARRAY_TASKS=${JOINT_ARRAY_TASKS}"
  export_arg+=",SEEDS_PER_TASK=${JOINT_SEEDS_PER_TASK}"
  export_arg+=",N_CORES=${JOINT_N_CORES}"
  export_arg+=",AUTO_VIZ=${AUTO_VIZ}"
  export_arg+=",R_MODULE=${R_MODULE}"
  export_arg+=",PARAMETER_TABLE=${PARAMETER_TABLE}"
  export_arg+=",FIT_OBJECTS_DIR=${FIT_OBJECTS_DIR}"
  export_arg+=",FLOW_DENSITY_PATH=${FLOW_DENSITY_PATH}"
  export_arg+=",ITERMAX=${ITERMAX}"
  export_arg+=",DE_RELTOL=${DE_RELTOL}"
  export_arg+=",DE_STEPTOL=${DE_STEPTOL}"
  export_arg+=",NP=${NP}"
  export_arg+=",JOINT_FITTING_MODE=${JOINT_FITTING_MODE}"
  export_arg+=",JOINT_WARMUP_ENABLE=${JOINT_WARMUP_ENABLE}"
  export_arg+=",JOINT_WARMUP_SEED_LABEL=${JOINT_WARMUP_SEED_LABEL}"
  export_arg+=",JOINT_WARMUP_INVIVO_SEED_DIR=${INVIVO_BEST_SEED_DIR}"
  export_arg+=",JOINT_WARMUP_INVITRO_SEED_DIR=${INVITRO_BEST_SEED_DIR}"
  export_arg+=",JOINT_WARMUP_SIGMAN=${JOINT_WARMUP_SIGMAN}"
  export_arg+=",JOINT_SOFT_COUPLING_SIGMA_DEFAULT=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}"
  export_arg+=",JOINT_SOFT_COUPLING_WELSCH_C=${JOINT_SOFT_COUPLING_WELSCH_C}"
  export_arg+=",JOINT_SOFT_COUPLING_PARAMETERS_TABLE=${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}"

  local cmd=(
    sbatch
    --parsable
    "--job-name=${JOINT_JOB_NAME}"
    "--qos=${JOINT_QOS}"
    "--time=${JOINT_TIME_LIMIT}"
    "--cpus-per-task=${JOINT_N_CORES}"
    "--mem=${JOINT_MEM}"
    "--array=1-${JOINT_ARRAY_TASKS}"
    "--output=${LOG_ROOT}/o2_joint_fit_%A_%a.out"
    "--error=${LOG_ROOT}/o2_joint_fit_%A_%a.err"
    "--export=${export_arg}"
  )
  if [[ -n "${dependency}" ]]; then
    cmd+=("--dependency=afterok:${dependency}")
  fi
  cmd+=("${JOINT_SUB_SCRIPT}")

  if truthy "${DRY_RUN}"; then
    print_command "Submit joint" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_JOINT_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted joint array job: ${LAST_JOB_ID}"
  fi
  record_array_submission "${run_dir}" "joint" "${LAST_JOB_ID}" \
    "${JOINT_TOTAL_SEEDS}" "${JOINT_ARRAY_TASKS}" "${JOINT_SEEDS_PER_TASK}" \
    "${JOINT_QOS}" "${JOINT_TIME_LIMIT}" "${JOINT_MEM}" "${JOINT_N_CORES}" \
    "${JOINT_SUB_SCRIPT}" "${JOINT_RUNNER_SCRIPT}" "${cmd[@]}"
  o2sd_prov_write_many "${run_dir}" \
    joint joint_fitting_mode "${JOINT_FITTING_MODE}" \
    joint joint_warmup_enable "${JOINT_WARMUP_ENABLE}" \
    joint joint_warmup_seed_label "${JOINT_WARMUP_SEED_LABEL}" \
    joint joint_warmup_invivo_seed_dir "${INVIVO_BEST_SEED_DIR}" \
    joint joint_warmup_invitro_seed_dir "${INVITRO_BEST_SEED_DIR}" \
    joint joint_soft_coupling_sigma_default "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" \
    joint joint_soft_coupling_welsch_c "${JOINT_SOFT_COUPLING_WELSCH_C}" \
    joint joint_soft_coupling_parameters_table "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}"
}

prepare_joint_warm_start_table() {
  if is_null_value "${INVIVO_BEST_SEED_DIR}" && is_null_value "${INVITRO_BEST_SEED_DIR}"; then
    if truthy "${JOINT_WARMUP_ENABLE}"; then
      echo "joint_warmup_enable=TRUE requires --invivo_best_seed_dir and --invitro_best_seed_dir." >&2
      exit 2
    fi
    JOINT_WARMUP_ENABLE="FALSE"
    return
  fi
  if is_null_value "${INVIVO_BEST_SEED_DIR}" || is_null_value "${INVITRO_BEST_SEED_DIR}"; then
    echo "Both --invivo_best_seed_dir and --invitro_best_seed_dir are required for joint warm start." >&2
    exit 2
  fi
  if [[ "${INVIVO_BEST_SEED_DIR}" != /* && -d "${PROJECT_ROOT}/${INVIVO_BEST_SEED_DIR}" ]]; then
    INVIVO_BEST_SEED_DIR="${PROJECT_ROOT}/${INVIVO_BEST_SEED_DIR}"
  fi
  if [[ "${INVITRO_BEST_SEED_DIR}" != /* && -d "${PROJECT_ROOT}/${INVITRO_BEST_SEED_DIR}" ]]; then
    INVITRO_BEST_SEED_DIR="${PROJECT_ROOT}/${INVITRO_BEST_SEED_DIR}"
  fi
  if [[ ! -d "${INVIVO_BEST_SEED_DIR}" ]]; then
    echo "Missing in vivo best seed directory: ${INVIVO_BEST_SEED_DIR}" >&2
    exit 1
  fi
  if [[ ! -d "${INVITRO_BEST_SEED_DIR}" ]]; then
    echo "Missing in vitro best seed directory: ${INVITRO_BEST_SEED_DIR}" >&2
    exit 1
  fi
  INVIVO_BEST_SEED_DIR="$(cd "${INVIVO_BEST_SEED_DIR}" && pwd)"
  INVITRO_BEST_SEED_DIR="$(cd "${INVITRO_BEST_SEED_DIR}" && pwd)"
  JOINT_WARMUP_ENABLE="TRUE"
  if is_null_value "${JOINT_WARMUP_SEED_LABEL}"; then
    JOINT_WARMUP_SEED_LABEL="$(derive_joint_warmup_seed_label)"
  else
    JOINT_WARMUP_SEED_LABEL="$(sanitize_label "${JOINT_WARMUP_SEED_LABEL}")"
  fi
  label_joint_run_prefix

  if [[ -z "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}" ]]; then
    JOINT_SOFT_COUPLING_PARAMETERS_TABLE="${PROJECT_ROOT}/oxygen/data/O2_supply_demand/joint_soft_coupling_parameters_table__${JOINT_WARMUP_SEED_LABEL}.csv"
  fi
  mkdir -p "$(dirname "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}")"
  JOINT_SOFT_COUPLING_PARAMETERS_TABLE="$(cd "$(dirname "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}")" && pwd)/$(basename "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}")"

  local cmd=(
    Rscript "${JOINT_WARM_START_SCRIPT}"
    "--invivo-seed-dir=${INVIVO_BEST_SEED_DIR}"
    "--invitro-seed-dir=${INVITRO_BEST_SEED_DIR}"
    "--seed-label=${JOINT_WARMUP_SEED_LABEL}"
    "--out=${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}"
  )
  if [[ -n "${JOINT_SOFT_COUPLING_DELTA_PARAMS}" ]]; then
    cmd+=("--delta-params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}")
  fi

  if truthy "${DRY_RUN}"; then
    print_command "Generate joint warm-start table" "${cmd[@]}"
  else
    load_r_module
    "${cmd[@]}"
  fi
}

first_line() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    echo "Missing file: ${path}" >&2
    exit 1
  fi
  local line
  line="$(head -n 1 "${path}")"
  line="$(printf "%s" "${line}" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  if [[ -z "${line}" ]]; then
    echo "Empty first line in ${path}" >&2
    exit 1
  fi
  printf "%s" "${line}"
}

resolve_existing_dir() {
  local label="$1"
  local path="$2"
  if [[ "${path}" != /* && -d "${PROJECT_ROOT}/${path}" ]]; then
    path="${PROJECT_ROOT}/${path}"
  fi
  if [[ ! -d "${path}" ]]; then
    echo "Missing ${label}: ${path}" >&2
    exit 1
  fi
  (cd "${path}" && pwd)
}

append_dependency() {
  local job_id="$1"
  if [[ -z "${job_id}" ]]; then
    return
  fi
  if [[ -z "${JOINT_PREP_DEPENDENCY}" ]]; then
    JOINT_PREP_DEPENDENCY="${job_id}"
  else
    JOINT_PREP_DEPENDENCY="${JOINT_PREP_DEPENDENCY}:${job_id}"
  fi
}

select_best_seed() {
  local label="$1"
  local run_dir="$2"
  local objective_columns="$3"
  local prep_dir="${OUT_ROOT}/${JOINT_RUN_PREFIX}/joint_prep"
  local log_path="${prep_dir}/select_best_${label}.log"
  mkdir -p "${prep_dir}"
  local cmd=(
    Rscript "${SELECT_BEST_SCRIPT}"
    "--run_dir=${run_dir}"
    "--objective_columns=${objective_columns}"
    "--required_files=${SELECT_REQUIRED_FILES}"
  )
  echo "Selecting best ${label} seed"
  print_command "Select ${label}" "${cmd[@]}"
  if ! "${cmd[@]}" > "${log_path}" 2>&1; then
    echo "Best-seed selection failed for ${label}. Log: ${log_path}" >&2
    tail -40 "${log_path}" >&2 || true
    exit 1
  fi
  cat "${log_path}"
}

run_joint_prep_stage() {
  if [[ -z "${INVIVO_RUN_DIR}" ]] && is_null_value "${INVIVO_BEST_SEED_DIR}"; then
    echo "Internal joint prep requires invivo_run_dir or invivo_best_seed_dir." >&2
    exit 2
  fi
  if [[ -z "${INVITRO_RUN_DIR}" ]] && is_null_value "${INVITRO_BEST_SEED_DIR}"; then
    echo "Internal joint prep requires invitro_run_dir or invitro_best_seed_dir." >&2
    exit 2
  fi
  load_r_module

  if is_null_value "${INVIVO_BEST_SEED_DIR}"; then
    INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
    select_best_seed "invivo" "${INVIVO_RUN_DIR}" "${INVIVO_OBJECTIVE_COLUMNS}"
    INVIVO_BEST_SEED_DIR="$(first_line "${INVIVO_RUN_DIR}/best_seed_from_summary.dir")"
  else
    INVIVO_BEST_SEED_DIR="$(resolve_existing_dir "in vivo best seed directory" "${INVIVO_BEST_SEED_DIR}")"
    if [[ -z "${INVIVO_RUN_DIR}" ]]; then
      INVIVO_RUN_DIR="$(cd "$(dirname "${INVIVO_BEST_SEED_DIR}")" && pwd)"
    else
      INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
    fi
  fi

  if is_null_value "${INVITRO_BEST_SEED_DIR}"; then
    INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
    select_best_seed "invitro" "${INVITRO_RUN_DIR}" "${INVITRO_OBJECTIVE_COLUMNS}"
    INVITRO_BEST_SEED_DIR="$(first_line "${INVITRO_RUN_DIR}/best_seed_from_summary.dir")"
  else
    INVITRO_BEST_SEED_DIR="$(resolve_existing_dir "in vitro best seed directory" "${INVITRO_BEST_SEED_DIR}")"
    if [[ -z "${INVITRO_RUN_DIR}" ]]; then
      INVITRO_RUN_DIR="$(cd "$(dirname "${INVITRO_BEST_SEED_DIR}")" && pwd)"
    else
      INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
    fi
  fi

  prepare_joint_warm_start_table
  submit_joint_array ""
  JOINT_JOB_ID="${LAST_JOB_ID}"
  submit_extra_results_job "o2_joint" "${OUT_ROOT}/${JOINT_RUN_PREFIX}" "${JOINT_JOB_ID}"
  JOINT_EXTRA_JOB_ID="${LAST_JOB_ID}"

  local manifest_dir="${OUT_ROOT}/${JOINT_RUN_PREFIX}/joint_prep"
  mkdir -p "${manifest_dir}"
  {
    printf "key\tvalue\n"
    printf "project_root\t%s\n" "${PROJECT_ROOT}"
    printf "config_path\t%s\n" "${CONFIG_PATH}"
    printf "out_root\t%s\n" "${OUT_ROOT}"
    printf "log_root\t%s\n" "${LOG_ROOT}"
    printf "invivo_run_dir\t%s\n" "${INVIVO_RUN_DIR}"
    printf "invitro_run_dir\t%s\n" "${INVITRO_RUN_DIR}"
    printf "invivo_best_seed_dir\t%s\n" "${INVIVO_BEST_SEED_DIR}"
    printf "invitro_best_seed_dir\t%s\n" "${INVITRO_BEST_SEED_DIR}"
    printf "joint_run_prefix\t%s\n" "${JOINT_RUN_PREFIX}"
    printf "joint_job_id\t%s\n" "${JOINT_JOB_ID}"
    printf "joint_extra_results_job_id\t%s\n" "${JOINT_EXTRA_JOB_ID}"
    printf "joint_warmup_seed_label\t%s\n" "${JOINT_WARMUP_SEED_LABEL}"
    printf "joint_soft_coupling_parameters_table\t%s\n" "${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}"
  } > "${manifest_dir}/joint_prep_manifest.tsv"
  echo "Joint prep complete."
  echo "  joint_run_prefix: ${JOINT_RUN_PREFIX}"
  echo "  manifest: ${manifest_dir}/joint_prep_manifest.tsv"
}

submit_joint_prep_job() {
  local dependency="$1"
  local prep_dir="${OUT_ROOT}/${JOINT_RUN_PREFIX}/joint_prep"
  mkdir -p "${prep_dir}"
  local cmd=(
    sbatch
    --parsable
    "--job-name=${JOINT_JOB_NAME}_prep"
    "--qos=${PREP_QOS}"
    "--time=${PREP_TIME_LIMIT}"
    --cpus-per-task=1
    "--mem=${PREP_MEM}"
    "--output=${LOG_ROOT}/o2_joint_prep_%A.out"
    "--error=${LOG_ROOT}/o2_joint_prep_%A.err"
  )
  if [[ -n "${dependency}" ]]; then
    cmd+=("--dependency=afterok:${dependency}")
  fi
  cmd+=(
    "${SELF_SCRIPT}"
    --internal_stage=select_and_submit_joint
    --fitting_mode=joint
    --joint_fitting_mode=JOINT
    "--project_root=${PROJECT_ROOT}"
    "--config_path=${CONFIG_PATH}"
    "--out_root=${OUT_ROOT}"
    "--log_root=${LOG_ROOT}"
    "--r_module=${R_MODULE}"
    "--force_extra_results=${FORCE_EXTRA_RESULTS}"
    "--invivo_run_dir=${INVIVO_RUN_DIR}"
    "--invitro_run_dir=${INVITRO_RUN_DIR}"
    "--joint_run_prefix=${JOINT_RUN_PREFIX}"
    "--joint_job_name=${JOINT_JOB_NAME}"
    "--joint_total_seeds=${JOINT_TOTAL_SEEDS}"
    "--joint_array_tasks=${JOINT_ARRAY_TASKS}"
    "--joint_seeds_per_task=${JOINT_SEEDS_PER_TASK}"
    "--joint_qos=${JOINT_QOS}"
    "--joint_time_limit=${JOINT_TIME_LIMIT}"
    "--joint_n_cores=${JOINT_N_CORES}"
    "--joint_mem=${JOINT_MEM}"
    "--postprocess_qos=${POSTPROCESS_QOS}"
    "--postprocess_time_limit=${POSTPROCESS_TIME_LIMIT}"
    "--postprocess_mem=${POSTPROCESS_MEM}"
    "--parameter_table=${PARAMETER_TABLE}"
    "--fit_objects_dir=${FIT_OBJECTS_DIR}"
    "--flow_density_path=${FLOW_DENSITY_PATH}"
    "--itermax=${ITERMAX}"
    "--de_reltol=${DE_RELTOL}"
    "--de_steptol=${DE_STEPTOL}"
    "--NP=${NP}"
    "--auto_viz=${AUTO_VIZ}"
    "--select_required_files=${SELECT_REQUIRED_FILES}"
    "--invivo_objective_columns=${INVIVO_OBJECTIVE_COLUMNS}"
    "--invitro_objective_columns=${INVITRO_OBJECTIVE_COLUMNS}"
    "--joint_warmup_enable=TRUE"
    "--joint_warmup_seed_label=${JOINT_WARMUP_SEED_LABEL}"
    "--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}"
    "--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}"
    "--joint_soft_coupling_welsch_c=${JOINT_SOFT_COUPLING_WELSCH_C}"
    "--joint_soft_coupling_parameters_table=${JOINT_SOFT_COUPLING_PARAMETERS_TABLE}"
    "--joint_soft_coupling_delta_params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}"
    --dry_run=FALSE
  )
  if ! is_null_value "${INVIVO_BEST_SEED_DIR}"; then
    cmd+=("--invivo_best_seed_dir=${INVIVO_BEST_SEED_DIR}")
  fi
  if ! is_null_value "${INVITRO_BEST_SEED_DIR}"; then
    cmd+=("--invitro_best_seed_dir=${INVITRO_BEST_SEED_DIR}")
  fi

  if truthy "${DRY_RUN}"; then
    print_command "Submit joint prep" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_JOINT_PREP_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted joint prep job: ${LAST_JOB_ID}"
  fi
}

submit_best_seed_joint_pipeline() {
  JOINT_PREP_DEPENDENCY=""

  if is_null_value "${INVIVO_BEST_SEED_DIR}"; then
    if is_null_value "${INVIVO_RUN_DIR}"; then
      INVIVO_RUN_DIR="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
      submit_invivo_array
      INVIVO_JOB_ID="${LAST_JOB_ID}"
      submit_extra_results_job "o2_invivo" "${INVIVO_RUN_DIR}" "${INVIVO_JOB_ID}"
      INVIVO_EXTRA_JOB_ID="${LAST_JOB_ID}"
      append_dependency "${INVIVO_EXTRA_JOB_ID}"
    else
      INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
      echo "Skipping in vivo fitting; using existing run directory: ${INVIVO_RUN_DIR}"
      submit_extra_results_job "o2_invivo_existing" "${INVIVO_RUN_DIR}" ""
      INVIVO_EXTRA_JOB_ID="${LAST_JOB_ID}"
      append_dependency "${INVIVO_EXTRA_JOB_ID}"
    fi
  else
    INVIVO_BEST_SEED_DIR="$(resolve_existing_dir "in vivo best seed directory" "${INVIVO_BEST_SEED_DIR}")"
    if is_null_value "${INVIVO_RUN_DIR}"; then
      INVIVO_RUN_DIR="$(cd "$(dirname "${INVIVO_BEST_SEED_DIR}")" && pwd)"
    else
      INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
    fi
    echo "Skipping in vivo fitting and best-seed selection; using: ${INVIVO_BEST_SEED_DIR}"
  fi

  if is_null_value "${INVITRO_BEST_SEED_DIR}"; then
    if is_null_value "${INVITRO_RUN_DIR}"; then
      INVITRO_RUN_DIR="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
      submit_invitro_array
      INVITRO_JOB_ID="${LAST_JOB_ID}"
      submit_extra_results_job "o2_invitro" "${INVITRO_RUN_DIR}" "${INVITRO_JOB_ID}"
      INVITRO_EXTRA_JOB_ID="${LAST_JOB_ID}"
      append_dependency "${INVITRO_EXTRA_JOB_ID}"
    else
      INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
      echo "Skipping in vitro fitting; using existing run directory: ${INVITRO_RUN_DIR}"
      submit_extra_results_job "o2_invitro_existing" "${INVITRO_RUN_DIR}" ""
      INVITRO_EXTRA_JOB_ID="${LAST_JOB_ID}"
      append_dependency "${INVITRO_EXTRA_JOB_ID}"
    fi
  else
    INVITRO_BEST_SEED_DIR="$(resolve_existing_dir "in vitro best seed directory" "${INVITRO_BEST_SEED_DIR}")"
    if is_null_value "${INVITRO_RUN_DIR}"; then
      INVITRO_RUN_DIR="$(cd "$(dirname "${INVITRO_BEST_SEED_DIR}")" && pwd)"
    else
      INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
    fi
    echo "Skipping in vitro fitting and best-seed selection; using: ${INVITRO_BEST_SEED_DIR}"
  fi

  if [[ -n "${JOINT_PREP_DEPENDENCY}" ]]; then
    submit_joint_prep_job "${JOINT_PREP_DEPENDENCY}"
  else
    prepare_joint_warm_start_table
    submit_joint_array ""
    JOINT_JOB_ID="${LAST_JOB_ID}"
    submit_extra_results_job "o2_joint" "${OUT_ROOT}/${JOINT_RUN_PREFIX}" "${JOINT_JOB_ID}"
  fi
}

shell_join() {
  local out=""
  local token
  local quoted
  for token in "$@"; do
    printf -v quoted "%q" "${token}"
    out+="${quoted} "
  done
  printf "%s" "${out% }"
}

run_multi_warmup_controller_stage() {
  INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
  INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
  load_r_module

  local multi_root="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
  local progress_log="${multi_root}/multi_warmup_progress.log"
  local manifest="${multi_root}/multi_warmup_manifest.tsv"
  local task_table="${multi_root}/multi_warmup_tasks.tsv"
  mkdir -p "${multi_root}" "${multi_root}/joint_soft_coupling_tables"
  : > "${progress_log}"
  o2sd_prov_write_standard "${multi_root}" "${SELF_SCRIPT}" "${SUBMIT_COMMAND_TEXT:-}"
  o2sd_prov_write_many "${multi_root}" \
    fit fitting_mode "multi_warmup" \
    fit run_prefix "${JOINT_RUN_PREFIX}" \
    multi_warmup invivo_run_dir "${INVIVO_RUN_DIR}" \
    multi_warmup invitro_run_dir "${INVITRO_RUN_DIR}" \
    multi_warmup manifest_path "${manifest}" \
    multi_warmup task_table_path "${task_table}" \
    multi_warmup invivo_k "${MULTI_WARMUP_INVIVO_K}" \
    multi_warmup invitro_k "${MULTI_WARMUP_INVITRO_K}" \
    multi_warmup umap_seed "${MULTI_WARMUP_UMAP_SEED}" \
    multi_warmup seeds_per_pair "${MULTI_WARMUP_SEEDS_PER_PAIR}" \
    joint joint_soft_coupling_sigma_default "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" \
    joint joint_soft_coupling_welsch_c "${JOINT_SOFT_COUPLING_WELSCH_C}"

  echo "Multi-warm-up controller"
  echo "  multi_warmup_root: ${multi_root}"
  echo "  invivo_run_dir: ${INVIVO_RUN_DIR}"
  echo "  invitro_run_dir: ${INVITRO_RUN_DIR}"
  echo "  seeds_per_pair: ${MULTI_WARMUP_SEEDS_PER_PAIR}"
  echo "  joint_soft_coupling_sigma_default: ${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}"
  echo "  joint_soft_coupling_welsch_c: ${JOINT_SOFT_COUPLING_WELSCH_C}"

  local seed_plan_cmd=(
    Rscript "${MULTI_WARMUP_SEED_PLAN_SCRIPT}"
    "--invivo_run_dir=${INVIVO_RUN_DIR}"
    "--invitro_run_dir=${INVITRO_RUN_DIR}"
    "--out_dir=${multi_root}"
    "--top_n=${MULTI_WARMUP_TOP_N}"
    "--umap_seed=${MULTI_WARMUP_UMAP_SEED}"
    "--invivo_k=${MULTI_WARMUP_INVIVO_K}"
    "--invitro_k=${MULTI_WARMUP_INVITRO_K}"
    "--invitro_anchor_ranks=${MULTI_WARMUP_INVITRO_ANCHOR_RANKS}"
    "--include_phase2=${MULTI_WARMUP_INCLUDE_PHASE2}"
    "--phase2_invitro_anchor_ranks=${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS}"
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
      echo "Generated manifest has no warm-up pairs: ${manifest}" >&2
      exit 1
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

submit_multi_warmup_controller_job() {
  local dependency="$1"
  local multi_root="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
  local controller_log_prefix="${LOG_ROOT}/o2_multi_warmup_submit"
  local controller_args=(
    bash "${SELF_SCRIPT}"
    --internal_stage=multi_warmup_prepare_and_submit
    --fitting_mode=joint
    --joint_fitting_mode=MULTI_WARMUP
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
    "--multi_warmup_top_n=${MULTI_WARMUP_TOP_N}"
    "--multi_warmup_umap_seed=${MULTI_WARMUP_UMAP_SEED}"
    "--multi_warmup_invivo_k=${MULTI_WARMUP_INVIVO_K}"
    "--multi_warmup_invitro_k=${MULTI_WARMUP_INVITRO_K}"
    "--multi_warmup_invitro_anchor_ranks=${MULTI_WARMUP_INVITRO_ANCHOR_RANKS}"
    "--multi_warmup_include_phase2=${MULTI_WARMUP_INCLUDE_PHASE2}"
    "--multi_warmup_phase2_invitro_anchor_ranks=${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS}"
    "--seeds_per_pair=${MULTI_WARMUP_SEEDS_PER_PAIR}"
    "--multi_warmup_task_order=${MULTI_WARMUP_TASK_ORDER}"
    "--multi_warmup_task_status_filter=${MULTI_WARMUP_TASK_STATUS_FILTER}"
    "--skip_existing=${MULTI_WARMUP_SKIP_EXISTING}"
    "--refresh_task_status=${MULTI_WARMUP_REFRESH_TASK_STATUS}"
    "--dry_run=FALSE"
  )
  if [[ -n "${MULTI_WARMUP_ARRAY_MAX_CONCURRENT}" ]]; then
    controller_args+=("--array_max_concurrent=${MULTI_WARMUP_ARRAY_MAX_CONCURRENT}")
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
    fit fitting_mode "multi_warmup" \
    fit run_prefix "${JOINT_RUN_PREFIX}" \
    multi_warmup invivo_run_dir "${INVIVO_RUN_DIR}" \
    multi_warmup invitro_run_dir "${INVITRO_RUN_DIR}" \
    multi_warmup seeds_per_pair "${MULTI_WARMUP_SEEDS_PER_PAIR}" \
    multi_warmup invivo_k "${MULTI_WARMUP_INVIVO_K}" \
    multi_warmup invitro_k "${MULTI_WARMUP_INVITRO_K}" \
    multi_warmup umap_seed "${MULTI_WARMUP_UMAP_SEED}" \
    joint joint_soft_coupling_sigma_default "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" \
    joint joint_soft_coupling_welsch_c "${JOINT_SOFT_COUPLING_WELSCH_C}" \
    slurm controller_job_id "${LAST_JOB_ID}" \
    slurm qos "${PREP_QOS}" \
    slurm walltime "${PREP_TIME_LIMIT}" \
    slurm memory "${PREP_MEM}" \
    slurm cpus "1"
}

submit_multi_warmup_pipeline() {
  if truthy "${USER_INVIVO_BEST_SEED_DIR}" || truthy "${USER_INVITRO_BEST_SEED_DIR}"; then
    echo "MULTI_WARMUP does not accept --invivo_best_seed_dir/--invitro_best_seed_dir. Provide --invivo_run_dir/--invitro_run_dir or omit them to submit source fits." >&2
    exit 2
  fi

  JOINT_PREP_DEPENDENCY=""
  INVIVO_BEST_SEED_DIR=""
  INVITRO_BEST_SEED_DIR=""
  JOINT_WARMUP_ENABLE="FALSE"

  if is_null_value "${INVIVO_RUN_DIR}"; then
    INVIVO_RUN_DIR="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
    submit_invivo_array
    INVIVO_JOB_ID="${LAST_JOB_ID}"
    submit_extra_results_job "o2_invivo" "${INVIVO_RUN_DIR}" "${INVIVO_JOB_ID}"
    INVIVO_EXTRA_JOB_ID="${LAST_JOB_ID}"
    append_dependency "${INVIVO_EXTRA_JOB_ID}"
  else
    INVIVO_RUN_DIR="$(resolve_existing_dir "in vivo run directory" "${INVIVO_RUN_DIR}")"
    echo "Skipping in vivo fitting; using existing run directory: ${INVIVO_RUN_DIR}"
    submit_extra_results_job "o2_invivo_existing" "${INVIVO_RUN_DIR}" ""
    INVIVO_EXTRA_JOB_ID="${LAST_JOB_ID}"
    append_dependency "${INVIVO_EXTRA_JOB_ID}"
  fi

  if is_null_value "${INVITRO_RUN_DIR}"; then
    INVITRO_RUN_DIR="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
    submit_invitro_array
    INVITRO_JOB_ID="${LAST_JOB_ID}"
    submit_extra_results_job "o2_invitro" "${INVITRO_RUN_DIR}" "${INVITRO_JOB_ID}"
    INVITRO_EXTRA_JOB_ID="${LAST_JOB_ID}"
    append_dependency "${INVITRO_EXTRA_JOB_ID}"
  else
    INVITRO_RUN_DIR="$(resolve_existing_dir "in vitro run directory" "${INVITRO_RUN_DIR}")"
    echo "Skipping in vitro fitting; using existing run directory: ${INVITRO_RUN_DIR}"
    submit_extra_results_job "o2_invitro_existing" "${INVITRO_RUN_DIR}" ""
    INVITRO_EXTRA_JOB_ID="${LAST_JOB_ID}"
    append_dependency "${INVITRO_EXTRA_JOB_ID}"
  fi

  submit_multi_warmup_controller_job "${JOINT_PREP_DEPENDENCY}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_INVIVO_RUN_PREFIX="fit_invivo_O2_buffering_500seed"
DEFAULT_INVITRO_RUN_PREFIX="fit_invitro_O2_buffering_500seed"
DEFAULT_JOINT_RUN_PREFIX="fit_joint_O2_buffering_500seed"
DEFAULT_INVIVO_BEST_SEED_REL="oxygen/results/fit_invivo_O2_buffering_500seed/seed50"
DEFAULT_INVITRO_BEST_SEED_REL="oxygen/results/fit_invitro_O2_buffering_500seed/seed350"
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
DEFAULT_ITERMAX="500"
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
DEFAULT_SELECT_REQUIRED_FILES="best_params.tsv"
DEFAULT_INVIVO_OBJECTIVE_COLUMNS="objective"
DEFAULT_INVITRO_OBJECTIVE_COLUMNS="objective_total,objective"
DEFAULT_JOINT_WARMUP_ENABLE="TRUE"
DEFAULT_JOINT_WARMUP_SEED_LABEL=""
DEFAULT_JOINT_WARMUP_SIGMAN=""
DEFAULT_JOINT_SOFT_COUPLING_SIGMA_DEFAULT=""
DEFAULT_JOINT_SOFT_COUPLING_WELSCH_C=""
DEFAULT_JOINT_SOFT_COUPLING_DELTA_PARAMS="default"
DEFAULT_MULTI_WARMUP_TOP_N="10"
DEFAULT_MULTI_WARMUP_UMAP_SEED="1"
DEFAULT_MULTI_WARMUP_INVIVO_K="auto"
DEFAULT_MULTI_WARMUP_INVITRO_K="auto"
DEFAULT_MULTI_WARMUP_INVITRO_ANCHOR_RANKS="1"
DEFAULT_MULTI_WARMUP_INCLUDE_PHASE2="FALSE"
DEFAULT_MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="auto"
DEFAULT_MULTI_WARMUP_SEEDS_PER_PAIR=""
DEFAULT_MULTI_WARMUP_ARRAY_MAX_CONCURRENT=""
DEFAULT_MULTI_WARMUP_TASK_ORDER="round_robin"
DEFAULT_MULTI_WARMUP_TASK_STATUS_FILTER="all"
DEFAULT_MULTI_WARMUP_SKIP_EXISTING="TRUE"
DEFAULT_MULTI_WARMUP_REFRESH_TASK_STATUS="TRUE"

FITTING_MODE="${FITTING_MODE:-}"
JOINT_FITTING_MODE="${JOINT_FITTING_MODE:-}"
INTERNAL_STAGE="${INTERNAL_STAGE:-}"
PROJECT_ROOT="${PROJECT_ROOT:-}"
R_MODULE="${R_MODULE:-}"
CONFIG_PATH="${CONFIG_PATH:-}"
OUT_ROOT="${OUT_ROOT:-}"
INVIVO_RUN_PREFIX="${INVIVO_RUN_PREFIX:-}"
INVITRO_RUN_PREFIX="${INVITRO_RUN_PREFIX:-}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-}"
JOINT_JOB_NAME="${JOINT_JOB_NAME:-}"
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
DE_RELTOL="${DE_RELTOL:-}"
DE_STEPTOL="${DE_STEPTOL:-}"
NP="${NP:-}"
AUTO_VIZ="${AUTO_VIZ:-}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-}"
INVIVO_BEST_SEED_DIR="${INVIVO_BEST_SEED_DIR:-}"
INVITRO_BEST_SEED_DIR="${INVITRO_BEST_SEED_DIR:-}"
USER_INVIVO_BEST_SEED_DIR="${USER_INVIVO_BEST_SEED_DIR:-FALSE}"
USER_INVITRO_BEST_SEED_DIR="${USER_INVITRO_BEST_SEED_DIR:-FALSE}"
JOINT_WARMUP_ENABLE="${JOINT_WARMUP_ENABLE:-}"
JOINT_WARMUP_SEED_LABEL="${JOINT_WARMUP_SEED_LABEL:-}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-}"
JOINT_SOFT_COUPLING_PARAMETERS_TABLE="${JOINT_SOFT_COUPLING_PARAMETERS_TABLE:-}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-}"
MULTI_WARMUP_TOP_N="${MULTI_WARMUP_TOP_N:-}"
MULTI_WARMUP_UMAP_SEED="${MULTI_WARMUP_UMAP_SEED:-}"
MULTI_WARMUP_INVIVO_K="${MULTI_WARMUP_INVIVO_K:-}"
MULTI_WARMUP_INVITRO_K="${MULTI_WARMUP_INVITRO_K:-}"
MULTI_WARMUP_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_INVITRO_ANCHOR_RANKS:-}"
MULTI_WARMUP_INCLUDE_PHASE2="${MULTI_WARMUP_INCLUDE_PHASE2:-}"
MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS:-}"
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
SELECT_REQUIRED_FILES="${SELECT_REQUIRED_FILES:-}"
INVIVO_OBJECTIVE_COLUMNS="${INVIVO_OBJECTIVE_COLUMNS:-}"
INVITRO_OBJECTIVE_COLUMNS="${INVITRO_OBJECTIVE_COLUMNS:-}"
LOG_ROOT="${LOG_ROOT:-}"

parse_args "$@"

INTERNAL_STAGE="${INTERNAL_STAGE:-${DEFAULT_INTERNAL_STAGE}}"
PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
INVIVO_RUN_PREFIX="${INVIVO_RUN_PREFIX:-${DEFAULT_INVIVO_RUN_PREFIX}}"
INVITRO_RUN_PREFIX="${INVITRO_RUN_PREFIX:-${DEFAULT_INVITRO_RUN_PREFIX}}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-${DEFAULT_JOINT_RUN_PREFIX}}"
JOINT_JOB_NAME="${JOINT_JOB_NAME:-${DEFAULT_JOINT_JOB_NAME}}"
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
DE_RELTOL="${DE_RELTOL:-${DEFAULT_DE_RELTOL}}"
DE_STEPTOL="${DE_STEPTOL:-${DEFAULT_DE_STEPTOL}}"
NP="${NP:-${DEFAULT_NP}}"
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-}"
JOINT_WARMUP_ENABLE="${JOINT_WARMUP_ENABLE:-${DEFAULT_JOINT_WARMUP_ENABLE}}"
RAW_JOINT_MODE="${JOINT_FITTING_MODE:-}"
RAW_JOINT_MODE="$(echo "${RAW_JOINT_MODE}" | tr '[:lower:]' '[:upper:]')"
RAW_JOINT_MODE="${RAW_JOINT_MODE//-/_}"
if [[ "${RAW_JOINT_MODE}" == "MULTI_WARMUP" ]]; then
  JOINT_WARMUP_ENABLE="FALSE"
fi
if truthy "${JOINT_WARMUP_ENABLE}"; then
  INVIVO_BEST_SEED_DIR="${INVIVO_BEST_SEED_DIR:-${PROJECT_ROOT}/${DEFAULT_INVIVO_BEST_SEED_REL}}"
  INVITRO_BEST_SEED_DIR="${INVITRO_BEST_SEED_DIR:-${PROJECT_ROOT}/${DEFAULT_INVITRO_BEST_SEED_REL}}"
else
  INVIVO_BEST_SEED_DIR="${INVIVO_BEST_SEED_DIR:-}"
  INVITRO_BEST_SEED_DIR="${INVITRO_BEST_SEED_DIR:-}"
fi
JOINT_WARMUP_SEED_LABEL="${JOINT_WARMUP_SEED_LABEL:-${DEFAULT_JOINT_WARMUP_SEED_LABEL}}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-${DEFAULT_JOINT_WARMUP_SIGMAN}}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-${DEFAULT_JOINT_SOFT_COUPLING_SIGMA_DEFAULT}}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-${DEFAULT_JOINT_SOFT_COUPLING_WELSCH_C}}"
JOINT_SOFT_COUPLING_PARAMETERS_TABLE="${JOINT_SOFT_COUPLING_PARAMETERS_TABLE:-}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-${DEFAULT_JOINT_SOFT_COUPLING_DELTA_PARAMS}}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-${DEFAULT_FORCE_EXTRA_RESULTS}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"
PREP_QOS="${PREP_QOS:-${DEFAULT_PREP_QOS}}"
PREP_TIME_LIMIT="${PREP_TIME_LIMIT:-${DEFAULT_PREP_TIME_LIMIT}}"
PREP_MEM="${PREP_MEM:-${DEFAULT_PREP_MEM}}"
REPORT_QOS="${REPORT_QOS:-${DEFAULT_REPORT_QOS}}"
REPORT_TIME_LIMIT="${REPORT_TIME_LIMIT:-${DEFAULT_REPORT_TIME_LIMIT}}"
REPORT_MEM="${REPORT_MEM:-${DEFAULT_REPORT_MEM}}"
SELECT_REQUIRED_FILES="${SELECT_REQUIRED_FILES:-${DEFAULT_SELECT_REQUIRED_FILES}}"
INVIVO_OBJECTIVE_COLUMNS="${INVIVO_OBJECTIVE_COLUMNS:-${DEFAULT_INVIVO_OBJECTIVE_COLUMNS}}"
INVITRO_OBJECTIVE_COLUMNS="${INVITRO_OBJECTIVE_COLUMNS:-${DEFAULT_INVITRO_OBJECTIVE_COLUMNS}}"
MULTI_WARMUP_TOP_N="${MULTI_WARMUP_TOP_N:-${DEFAULT_MULTI_WARMUP_TOP_N}}"
MULTI_WARMUP_UMAP_SEED="${MULTI_WARMUP_UMAP_SEED:-${DEFAULT_MULTI_WARMUP_UMAP_SEED}}"
MULTI_WARMUP_INVIVO_K="${MULTI_WARMUP_INVIVO_K:-${DEFAULT_MULTI_WARMUP_INVIVO_K}}"
MULTI_WARMUP_INVITRO_K="${MULTI_WARMUP_INVITRO_K:-${DEFAULT_MULTI_WARMUP_INVITRO_K}}"
MULTI_WARMUP_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_INVITRO_ANCHOR_RANKS:-${DEFAULT_MULTI_WARMUP_INVITRO_ANCHOR_RANKS}}"
MULTI_WARMUP_INCLUDE_PHASE2="${MULTI_WARMUP_INCLUDE_PHASE2:-${DEFAULT_MULTI_WARMUP_INCLUDE_PHASE2}}"
MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS:-${DEFAULT_MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS}}"
MULTI_WARMUP_SEEDS_PER_PAIR="${MULTI_WARMUP_SEEDS_PER_PAIR:-${DEFAULT_MULTI_WARMUP_SEEDS_PER_PAIR}}"
MULTI_WARMUP_ARRAY_MAX_CONCURRENT="${MULTI_WARMUP_ARRAY_MAX_CONCURRENT:-${DEFAULT_MULTI_WARMUP_ARRAY_MAX_CONCURRENT}}"
MULTI_WARMUP_TASK_ORDER="${MULTI_WARMUP_TASK_ORDER:-${DEFAULT_MULTI_WARMUP_TASK_ORDER}}"
MULTI_WARMUP_TASK_STATUS_FILTER="${MULTI_WARMUP_TASK_STATUS_FILTER:-${DEFAULT_MULTI_WARMUP_TASK_STATUS_FILTER}}"
MULTI_WARMUP_SKIP_EXISTING="${MULTI_WARMUP_SKIP_EXISTING:-${DEFAULT_MULTI_WARMUP_SKIP_EXISTING}}"
MULTI_WARMUP_REFRESH_TASK_STATUS="${MULTI_WARMUP_REFRESH_TASK_STATUS:-${DEFAULT_MULTI_WARMUP_REFRESH_TASK_STATUS}}"

FITTING_MODE="$(normalize_fitting_mode "${FITTING_MODE}")"
if [[ -z "${FITTING_MODE}" ]]; then
  echo "--fitting_mode must be one of invivo, invitro, or joint." >&2
  usage >&2
  exit 2
fi

if [[ "${FITTING_MODE}" != "joint" ]]; then
  JOINT_FITTING_MODE="OFF"
elif [[ -z "${JOINT_FITTING_MODE}" ]]; then
  JOINT_FITTING_MODE="DIRECT"
fi
JOINT_FITTING_MODE="$(echo "${JOINT_FITTING_MODE}" | tr '[:lower:]' '[:upper:]')"
JOINT_FITTING_MODE="${JOINT_FITTING_MODE//-/_}"
case "${JOINT_FITTING_MODE}" in
  SINGLE) JOINT_FITTING_MODE="DIRECT" ;;
esac
case "${JOINT_FITTING_MODE}" in
  OFF|JOINT|DIRECT|MULTI_WARMUP) ;;
  *)
    echo "--joint_fitting_mode must be OFF, JOINT, DIRECT, or MULTI_WARMUP. SINGLE is accepted as a legacy alias for DIRECT." >&2
    exit 2
    ;;
esac

if [[ "${JOINT_FITTING_MODE}" == "MULTI_WARMUP" ]]; then
  if truthy "${USER_INVIVO_BEST_SEED_DIR}" || truthy "${USER_INVITRO_BEST_SEED_DIR}"; then
    echo "MULTI_WARMUP does not accept --invivo_best_seed_dir/--invitro_best_seed_dir. Provide --invivo_run_dir/--invitro_run_dir or omit them to submit source fits." >&2
    exit 2
  fi
  if is_null_value "${MULTI_WARMUP_SEEDS_PER_PAIR}"; then
    MULTI_WARMUP_SEEDS_PER_PAIR="${JOINT_TOTAL_SEEDS}"
  fi
  require_positive_int "MULTI_WARMUP_SEEDS_PER_PAIR" "${MULTI_WARMUP_SEEDS_PER_PAIR}"
  JOINT_TOTAL_SEEDS="${MULTI_WARMUP_SEEDS_PER_PAIR}"
  JOINT_ARRAY_TASKS="${MULTI_WARMUP_SEEDS_PER_PAIR}"
  JOINT_SEEDS_PER_TASK="1"
  INVIVO_BEST_SEED_DIR=""
  INVITRO_BEST_SEED_DIR=""
  JOINT_WARMUP_ENABLE="FALSE"
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
if [[ -z "${LOG_ROOT}" ]]; then
  LOG_ROOT="${OUT_ROOT}/log"
fi
mkdir -p "${LOG_ROOT}"
LOG_ROOT="$(cd "${LOG_ROOT}" && pwd)"

HPC_DIR="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/hpc"
SELF_SCRIPT="${SELF_SCRIPT:-${HPC_DIR}/submit_o2_fit.sh}"
if [[ ! -f "${SELF_SCRIPT}" ]]; then
  SELF_SCRIPT_CANDIDATE="${HPC_DIR}/$(basename "${BASH_SOURCE[0]}")"
  if [[ -f "${SELF_SCRIPT_CANDIDATE}" ]]; then
    SELF_SCRIPT="${SELF_SCRIPT_CANDIDATE}"
  fi
fi
PROVENANCE_HELPER="${HPC_DIR}/write_run_provenance.sh"
if [[ -f "${PROVENANCE_HELPER}" ]]; then
  # shellcheck source=/dev/null
  source "${PROVENANCE_HELPER}"
else
  echo "Missing provenance helper: ${PROVENANCE_HELPER}" >&2
  exit 1
fi
SUBMIT_COMMAND_TEXT="$(o2sd_prov_shell_join bash "${SELF_SCRIPT}" "${ORIGINAL_SUBMIT_ARGS[@]}")"
INVIVO_SUB_SCRIPT="${INVIVO_SUB_SCRIPT:-${HPC_DIR}/submit_fit_seed_array_buffering.sub}"
INVITRO_SUB_SCRIPT="${INVITRO_SUB_SCRIPT:-${HPC_DIR}/submit_fit_seed_array_invitro_buffering.sub}"
JOINT_SUB_SCRIPT="${JOINT_SUB_SCRIPT:-${HPC_DIR}/submit_fit_seed_array_joint_buffering.sub}"
MULTI_WARMUP_SEED_PLAN_SCRIPT="${MULTI_WARMUP_SEED_PLAN_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/build_multi_warmup_seed_plan.R}"
MULTI_WARMUP_TASK_TABLE_SCRIPT="${MULTI_WARMUP_TASK_TABLE_SCRIPT:-${HPC_DIR}/build_multi_warmup_task_table.R}"
MULTI_WARMUP_TASK_SUBMIT_SCRIPT="${MULTI_WARMUP_TASK_SUBMIT_SCRIPT:-${HPC_DIR}/submit_multi_warmup_task_table.sh}"
POSTPROCESS_SCRIPT="${POSTPROCESS_SCRIPT:-${HPC_DIR}/postprocess_extra_results.sh}"
EXTRA_RESULTS_SCRIPT="${EXTRA_RESULTS_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/extra_results.R}"
SELECT_BEST_SCRIPT="${SELECT_BEST_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/select_best_seed_from_summary.R}"
JOINT_WARM_START_SCRIPT="${JOINT_WARM_START_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/make_joint_soft_coupling_parameters_table.R}"
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

for path in "${CONFIG_PATH}" "${SELF_SCRIPT}" "${INVIVO_SUB_SCRIPT}" "${INVITRO_SUB_SCRIPT}" "${JOINT_SUB_SCRIPT}" \
            "${MULTI_WARMUP_SEED_PLAN_SCRIPT}" "${MULTI_WARMUP_TASK_TABLE_SCRIPT}" "${MULTI_WARMUP_TASK_SUBMIT_SCRIPT}" \
            "${POSTPROCESS_SCRIPT}" "${EXTRA_RESULTS_SCRIPT}" \
            "${SELECT_BEST_SCRIPT}" "${JOINT_WARM_START_SCRIPT}" \
            "${INVIVO_RUNNER_SCRIPT}" "${INVITRO_RUNNER_SCRIPT}" "${JOINT_RUNNER_SCRIPT}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
done

for name in INVIVO_TOTAL_SEEDS INVIVO_ARRAY_TASKS INVIVO_SEEDS_PER_TASK \
            INVITRO_TOTAL_SEEDS INVITRO_ARRAY_TASKS INVITRO_SEEDS_PER_TASK \
            JOINT_TOTAL_SEEDS JOINT_ARRAY_TASKS JOINT_SEEDS_PER_TASK \
            INVIVO_N_CORES INVITRO_N_CORES JOINT_N_CORES ITERMAX DE_STEPTOL NP; do
  require_positive_int "${name}" "${!name}"
done
check_seed_plan INVIVO "${INVIVO_TOTAL_SEEDS}" "${INVIVO_ARRAY_TASKS}" "${INVIVO_SEEDS_PER_TASK}"
check_seed_plan INVITRO "${INVITRO_TOTAL_SEEDS}" "${INVITRO_ARRAY_TASKS}" "${INVITRO_SEEDS_PER_TASK}"
check_seed_plan JOINT "${JOINT_TOTAL_SEEDS}" "${JOINT_ARRAY_TASKS}" "${JOINT_SEEDS_PER_TASK}"

ensure_sbatch

echo "O2 submitter"
echo "  fitting_mode: ${FITTING_MODE}"
echo "  joint_fitting_mode: ${JOINT_FITTING_MODE}"
echo "  project_root: ${PROJECT_ROOT}"
echo "  out_root: ${OUT_ROOT}"
echo "  log_root: ${LOG_ROOT}"
echo "  r_module: ${R_MODULE}"
echo "  invivo resources: qos=${INVIVO_QOS}, time=${INVIVO_TIME_LIMIT}, mem=${INVIVO_MEM}, cpus=${INVIVO_N_CORES}"
echo "  invitro resources: qos=${INVITRO_QOS}, time=${INVITRO_TIME_LIMIT}, mem=${INVITRO_MEM}, cpus=${INVITRO_N_CORES}"
echo "  joint resources: qos=${JOINT_QOS}, time=${JOINT_TIME_LIMIT}, mem=${JOINT_MEM}, cpus=${JOINT_N_CORES}"
echo "  postprocess resources: qos=${POSTPROCESS_QOS}, time=${POSTPROCESS_TIME_LIMIT}, mem=${POSTPROCESS_MEM}"
echo "  prep resources: qos=${PREP_QOS}, time=${PREP_TIME_LIMIT}, mem=${PREP_MEM}"
if [[ "${JOINT_FITTING_MODE}" == "MULTI_WARMUP" ]]; then
  echo "  multi_warmup seeds_per_pair: ${MULTI_WARMUP_SEEDS_PER_PAIR}"
  echo "  multi_warmup task_order: ${MULTI_WARMUP_TASK_ORDER}"
  echo "  multi_warmup array_max_concurrent: ${MULTI_WARMUP_ARRAY_MAX_CONCURRENT:-none}"
  echo "  multi_warmup task_status_filter: ${MULTI_WARMUP_TASK_STATUS_FILTER}"
  echo "  multi_warmup skip_existing: ${MULTI_WARMUP_SKIP_EXISTING}"
fi

case "${INTERNAL_STAGE}" in
  submit|select_and_submit_joint|multi_warmup_prepare_and_submit) ;;
  *)
    echo "--internal_stage must be submit, select_and_submit_joint, or multi_warmup_prepare_and_submit, got: ${INTERNAL_STAGE}" >&2
    exit 2
    ;;
esac

if [[ "${INTERNAL_STAGE}" == "select_and_submit_joint" ]]; then
  run_joint_prep_stage
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
    case "${JOINT_FITTING_MODE}" in
      OFF)
        echo "joint_fitting_mode=OFF; no fitting submitted."
        ;;
      JOINT)
        echo "joint_fitting_mode=JOINT using best-seed selection pipeline."
        submit_best_seed_joint_pipeline
        ;;
      MULTI_WARMUP)
        echo "joint_fitting_mode=MULTI_WARMUP using source-run ratio clustering."
        submit_multi_warmup_pipeline
        ;;
      DIRECT)
        if ! is_null_value "${INVIVO_RUN_DIR}" || ! is_null_value "${INVITRO_RUN_DIR}"; then
          echo "Ignoring invivo_run_dir/invitro_run_dir: current joint fitting reads inputs from the config, not anchor-derived single-fit outputs."
        fi
        prepare_joint_warm_start_table
        submit_joint_array ""
        JOINT_JOB_ID="${LAST_JOB_ID}"
        submit_extra_results_job "o2_joint" "${OUT_ROOT}/${JOINT_RUN_PREFIX}" "${JOINT_JOB_ID}"
        ;;
    esac
    ;;
esac
