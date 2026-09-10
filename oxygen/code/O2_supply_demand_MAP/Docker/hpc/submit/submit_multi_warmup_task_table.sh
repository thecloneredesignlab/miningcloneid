#!/bin/bash
# Submit a global Slurm array from multi_warmup_tasks.tsv.
#
# The task table should contain one row per warm-up pair x joint seed. The
# corresponding array worker reads the row for SLURM_ARRAY_TASK_ID and calls
# the existing joint runner with the row-specific warm-up directories.

set -euo pipefail

O2SD_DOCKER_HPC_ROOT="${O2SD_DOCKER_HPC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck source=../util/o2_supply_demand_map_apptainer_runtime.sh
source "${O2SD_DOCKER_HPC_ROOT}/util/o2_supply_demand_map_apptainer_runtime.sh"

O2SD_SHELL_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../../../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"

ORIGINAL_SUBMIT_ARGS=("$@")

usage() {
  cat <<'EOF'
Usage:
  bash submit_multi_warmup_task_table.sh --tasks_tsv=/path/to/multi_warmup_tasks.tsv [options]
  bash submit_multi_warmup_task_table.sh --multi_warmup_root=/path/to/result_root [options]

Required:
  Either --tasks_tsv=FILE or --multi_warmup_root=DIR

Common options:
  --project_root=DIR
  --multi_warmup_root=DIR
  --seeds_per_pair=200
  --refresh_task_table=TRUE|FALSE
  --joint_soft_coupling_sigma_default=0.65
  --joint_soft_coupling_welsch_c=0.4
  --array_spec=1-1000
  --array_max_concurrent=N
  --task_status_filter=all
  --skip_existing=TRUE|FALSE
  --joint_n_cores=22
  --parameter_table=FILE
  --fit_objects_dir=DIR
  --flow_density_path=FILE
  --itermax=N
  --NP=N
  --auto_viz=TRUE|FALSE
  --joint_mem=32G
  --joint_qos=xxlarge
  --joint_time_limit=12:00:00
  --submit_postprocess=TRUE|FALSE
  --submit_report=TRUE|FALSE
  --force_extra_results=TRUE|FALSE
  --dry_run=TRUE|FALSE

Examples:
  # Submit all rows, but skip seeds that already have fit_summary.tsv and best_params.tsv.
  bash submit_multi_warmup_task_table.sh \
    --tasks_tsv=/share/.../fit_joint_multi_warmup.../multi_warmup_tasks.tsv \
    --array_max_concurrent=100

  # Submit only rows whose task_status is not_started in the table.
  bash submit_multi_warmup_task_table.sh \
    --tasks_tsv=/share/.../multi_warmup_tasks.tsv \
    --task_status_filter=not_started \
    --array_max_concurrent=100

  # Run a specific row or range from the table.
  bash submit_multi_warmup_task_table.sh \
    --tasks_tsv=/share/.../multi_warmup_tasks.tsv \
    --array_spec=37,42,101-125
EOF
}

submit_or_print() {
  local label="$1"
  shift
  {
    printf "%s:" "${label}"
    printf " %q" "$@"
    printf "\n"
  } | tee -a "${PROGRESS_LOG}" >&2
  if truthy "${DRY_RUN}"; then
    echo "DRY_RUN_JOB_ID"
  else
    "$@" | awk '{print $NF}'
  fi
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h) usage; exit 0 ;;
      --tasks_tsv=*|--task_table=*) TASKS_TSV="${arg#*=}" ;;
      --multi_warmup_root=*|--root=*) MULTI_WARMUP_ROOT="${arg#*=}" ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --build_task_table_script=*) BUILD_TASK_TABLE_SCRIPT="${arg#*=}" ;;
      --array_script=*|--task_array_script=*) TASK_ARRAY_SCRIPT="${arg#*=}" ;;
      --postprocess_script=*) POSTPROCESS_SCRIPT="${arg#*=}" ;;
      --collect_script=*) COLLECT_SCRIPT="${arg#*=}" ;;
      --report_script=*) REPORT_SCRIPT="${arg#*=}" ;;
      --seeds_per_pair=*|--joint_seeds_per_pair=*|--multi_warmup_seeds_per_pair=*) SEEDS_PER_PAIR="${arg#*=}" ;;
      --joint_total_seeds=*) SEEDS_PER_PAIR="${arg#*=}" ;;
      --joint_array_tasks=*|--array_tasks=*) TASK_TABLE_ARRAY_TASKS="${arg#*=}" ;;
      --joint_seeds_per_task=*|--seeds_per_task=*) TASK_TABLE_SEEDS_PER_TASK="${arg#*=}" ;;
      --task_order=*|--order=*) TASK_ORDER="${arg#*=}" ;;
      --refresh_task_table=*) REFRESH_TASK_TABLE="${arg#*=}" ;;
      --refresh_task_status=*) REFRESH_TASK_STATUS="${arg#*=}" ;;
      --config_path=*|--config=*) CONFIG_PATH="${arg#*=}" ;;
      --runner_script=*|--joint_runner_script=*) RUNNER_SCRIPT="${arg#*=}" ;;
      --parameter_table=*|--invitro_parameter_table=*|--parameter_table_invitro=*) PARAMETER_TABLE="${arg#*=}" ;;
      --fit_objects_dir=*) FIT_OBJECTS_DIR="${arg#*=}" ;;
      --flow_density_path=*) FLOW_DENSITY_PATH="${arg#*=}" ;;
      --itermax=*) ITERMAX="${arg#*=}" ;;
      --de_reltol=*) DE_RELTOL="${arg#*=}" ;;
      --de_steptol=*) DE_STEPTOL="${arg#*=}" ;;
      --np=*|--NP=*) NP="${arg#*=}" ;;
      --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
      --array_spec=*|--array=*) ARRAY_SPEC="${arg#*=}" ;;
      --array_max_concurrent=*|--max_concurrent=*) ARRAY_MAX_CONCURRENT="${arg#*=}" ;;
      --task_lookup_column=*) TASK_LOOKUP_COLUMN="${arg#*=}" ;;
      --task_status_filter=*) TASK_STATUS_FILTER="${arg#*=}" ;;
      --skip_existing=*) SKIP_EXISTING="${arg#*=}" ;;
      --joint_n_cores=*|--n_cores=*)
        JOINT_N_CORES="${arg#*=}"
        TASK_TABLE_JOINT_N_CORES="${arg#*=}"
        ;;
      --joint_mem=*|--mem=*) JOINT_MEM="${arg#*=}" ;;
      --joint_qos=*|--qos=*) JOINT_QOS="${arg#*=}" ;;
      --joint_time_limit=*|--time=*|--time_limit=*) JOINT_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_qos=*) POSTPROCESS_QOS="${arg#*=}" ;;
      --postprocess_time_limit=*) POSTPROCESS_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_mem=*) POSTPROCESS_MEM="${arg#*=}" ;;
      --report_qos=*) REPORT_QOS="${arg#*=}" ;;
      --report_time_limit=*) REPORT_TIME_LIMIT="${arg#*=}" ;;
      --report_mem=*) REPORT_MEM="${arg#*=}" ;;
      --submit_postprocess=*) SUBMIT_POSTPROCESS="${arg#*=}" ;;
      --submit_report=*) SUBMIT_REPORT="${arg#*=}" ;;
      --force_extra_results=*) FORCE_EXTRA_RESULTS="${arg#*=}" ;;
      --joint_warmup_sigmaN=*) JOINT_WARMUP_SIGMAN="${arg#*=}" ;;
      --joint_soft_coupling_sigma_default=*) JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${arg#*=}" ;;
      --joint_soft_coupling_welsch_c=*) JOINT_SOFT_COUPLING_WELSCH_C="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --log_root=*|--log_dir=*) LOG_ROOT="${arg#*=}" ;;
      --job_name=*) JOB_NAME="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      *) echo "Unknown argument: ${arg}" >&2; usage >&2; exit 2 ;;
    esac
  done
}

tsv_col_index() {
  local col="$1"
  awk -F '\t' -v col="${col}" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == col) { print i; exit } }' "${TASKS_TSV}"
}

first_col_value() {
  local col="$1"
  local idx
  idx="$(tsv_col_index "${col}")"
  if [[ -z "${idx}" ]]; then
    return 1
  fi
  awk -F '\t' -v i="${idx}" 'NR == 2 { print $i; exit }' "${TASKS_TSV}"
}

write_unique_run_dirs() {
  local out_path="$1"
  awk -F '\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) h[$i] = i
      next
    }
    NR > 1 {
      label = $(h["warmup_label"])
      run = $(h["joint_run_dir"])
      if (run != "" && !seen[run]++) print label "\t" run
    }
  ' "${TASKS_TSV}" > "${out_path}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"

TASKS_TSV="${TASKS_TSV:-}"
MULTI_WARMUP_ROOT="${MULTI_WARMUP_ROOT:-}"
PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
BUILD_TASK_TABLE_SCRIPT="${BUILD_TASK_TABLE_SCRIPT:-}"
TASK_ARRAY_SCRIPT="${TASK_ARRAY_SCRIPT:-}"
POSTPROCESS_SCRIPT="${POSTPROCESS_SCRIPT:-}"
COLLECT_SCRIPT="${COLLECT_SCRIPT:-}"
REPORT_SCRIPT="${REPORT_SCRIPT:-}"
SEEDS_PER_PAIR="${SEEDS_PER_PAIR:-}"
TASK_TABLE_ARRAY_TASKS="${TASK_TABLE_ARRAY_TASKS:-}"
TASK_TABLE_SEEDS_PER_TASK="${TASK_TABLE_SEEDS_PER_TASK:-}"
TASK_ORDER="${TASK_ORDER:-round_robin}"
REFRESH_TASK_TABLE="${REFRESH_TASK_TABLE:-FALSE}"
REFRESH_TASK_STATUS="${REFRESH_TASK_STATUS:-TRUE}"
CONFIG_PATH="${CONFIG_PATH:-}"
RUNNER_SCRIPT="${RUNNER_SCRIPT:-}"
PARAMETER_TABLE="${PARAMETER_TABLE:-}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-}"
ITERMAX="${ITERMAX:-}"
DE_RELTOL="${DE_RELTOL:-}"
DE_STEPTOL="${DE_STEPTOL:-}"
NP="${NP:-}"
AUTO_VIZ="${AUTO_VIZ:-}"
ARRAY_SPEC="${ARRAY_SPEC:-}"
ARRAY_MAX_CONCURRENT="${ARRAY_MAX_CONCURRENT:-}"
TASK_LOOKUP_COLUMN="${TASK_LOOKUP_COLUMN:-recommended_sbatch_array_index}"
TASK_STATUS_FILTER="${TASK_STATUS_FILTER:-all}"
SKIP_EXISTING="${SKIP_EXISTING:-TRUE}"
TASK_TABLE_JOINT_N_CORES="${TASK_TABLE_JOINT_N_CORES:-${JOINT_N_CORES:-}}"
JOINT_N_CORES="${JOINT_N_CORES:-22}"
JOINT_MEM="${JOINT_MEM:-32G}"
JOINT_QOS="${JOINT_QOS:-xxlarge}"
JOINT_TIME_LIMIT="${JOINT_TIME_LIMIT:-12:00:00}"
POSTPROCESS_QOS="${POSTPROCESS_QOS:-small}"
POSTPROCESS_TIME_LIMIT="${POSTPROCESS_TIME_LIMIT:-4:00:00}"
POSTPROCESS_MEM="${POSTPROCESS_MEM:-8G}"
REPORT_QOS="${REPORT_QOS:-small}"
REPORT_TIME_LIMIT="${REPORT_TIME_LIMIT:-4:00:00}"
REPORT_MEM="${REPORT_MEM:-8G}"
SUBMIT_POSTPROCESS="${SUBMIT_POSTPROCESS:-TRUE}"
SUBMIT_REPORT="${SUBMIT_REPORT:-TRUE}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-FALSE}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-}"
R_MODULE="${R_MODULE:-R/4.4}"
LOG_ROOT="${LOG_ROOT:-}"
JOB_NAME="${JOB_NAME:-o2mw_tasks}"
DRY_RUN="${DRY_RUN:-FALSE}"

parse_args "$@"

if [[ -z "${TASKS_TSV}" && -z "${MULTI_WARMUP_ROOT}" ]]; then
  echo "Either --tasks_tsv or --multi_warmup_root is required." >&2
  usage >&2
  exit 2
fi
if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found; run on an HPC login node or use --dry_run=TRUE." >&2
  exit 1
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
PROVENANCE_HELPER="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_shell_utils.sh"
if [[ -f "${PROVENANCE_HELPER}" ]]; then
  # shellcheck source=/dev/null
  source "${PROVENANCE_HELPER}"
else
  echo "Missing provenance helper: ${PROVENANCE_HELPER}" >&2
  exit 1
fi

if [[ -z "${BUILD_TASK_TABLE_SCRIPT}" ]]; then
  BUILD_TASK_TABLE_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/submit/build_multi_warmup_task_table.R"
fi
BUILD_TASK_TABLE_SCRIPT="$(cd "$(dirname "${BUILD_TASK_TABLE_SCRIPT}")" && pwd)/$(basename "${BUILD_TASK_TABLE_SCRIPT}")"

if [[ -n "${MULTI_WARMUP_ROOT}" ]]; then
  MULTI_WARMUP_ROOT="$(cd "${MULTI_WARMUP_ROOT}" && pwd)"
  if [[ -z "${TASKS_TSV}" ]]; then
    TASKS_TSV="${MULTI_WARMUP_ROOT}/multi_warmup_tasks.tsv"
  fi
fi
TASKS_TSV="$(cd "$(dirname "${TASKS_TSV}")" && pwd)/$(basename "${TASKS_TSV}")"

if [[ ! -f "${TASKS_TSV}" ]] || truthy "${REFRESH_TASK_TABLE}"; then
  if [[ -z "${MULTI_WARMUP_ROOT}" ]]; then
    MULTI_WARMUP_ROOT="$(cd "$(dirname "${TASKS_TSV}")" && pwd)"
  fi
  if [[ ! -f "${BUILD_TASK_TABLE_SCRIPT}" ]]; then
    echo "Missing task-table builder: ${BUILD_TASK_TABLE_SCRIPT}" >&2
    exit 1
  fi
  build_cmd=(
    Rscript "${BUILD_TASK_TABLE_SCRIPT}"
    "--multi_warmup_root=${MULTI_WARMUP_ROOT}"
    "--out=${TASKS_TSV}"
    "--order=${TASK_ORDER}"
    "--refresh_status=${REFRESH_TASK_STATUS}"
    "--log_root=${LOG_ROOT:-${PROJECT_ROOT}/oxygen/results/log}"
  )
  if [[ -n "${SEEDS_PER_PAIR}" ]]; then build_cmd+=("--seeds_per_pair=${SEEDS_PER_PAIR}"); fi
  if [[ -n "${TASK_TABLE_ARRAY_TASKS}" ]]; then build_cmd+=("--array_tasks=${TASK_TABLE_ARRAY_TASKS}"); fi
  if [[ -n "${TASK_TABLE_SEEDS_PER_TASK}" ]]; then build_cmd+=("--seeds_per_task=${TASK_TABLE_SEEDS_PER_TASK}"); fi
  if [[ -n "${CONFIG_PATH}" ]]; then build_cmd+=("--config_path=${CONFIG_PATH}"); fi
  if [[ -n "${RUNNER_SCRIPT}" ]]; then build_cmd+=("--runner_script=${RUNNER_SCRIPT}"); fi
  if [[ -n "${PARAMETER_TABLE}" ]]; then build_cmd+=("--parameter_table=${PARAMETER_TABLE}"); fi
  if [[ -n "${FIT_OBJECTS_DIR}" ]]; then build_cmd+=("--fit_objects_dir=${FIT_OBJECTS_DIR}"); fi
  if [[ -n "${FLOW_DENSITY_PATH}" ]]; then build_cmd+=("--flow_density_path=${FLOW_DENSITY_PATH}"); fi
  if [[ -n "${TASK_TABLE_JOINT_N_CORES}" ]]; then build_cmd+=("--joint_n_cores=${TASK_TABLE_JOINT_N_CORES}"); fi
  if [[ -n "${R_MODULE}" ]]; then build_cmd+=("--r_module=${R_MODULE}"); fi
  if [[ -n "${ITERMAX}" ]]; then build_cmd+=("--itermax=${ITERMAX}"); fi
  if [[ -n "${DE_RELTOL}" ]]; then build_cmd+=("--de_reltol=${DE_RELTOL}"); fi
  if [[ -n "${DE_STEPTOL}" ]]; then build_cmd+=("--de_steptol=${DE_STEPTOL}"); fi
  if [[ -n "${NP}" ]]; then build_cmd+=("--NP=${NP}"); fi
  if [[ -n "${AUTO_VIZ}" ]]; then build_cmd+=("--auto_viz=${AUTO_VIZ}"); fi
  if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then build_cmd+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}"); fi
  if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then build_cmd+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}"); fi
  if [[ -n "${JOINT_SOFT_COUPLING_WELSCH_C}" ]]; then build_cmd+=("--joint_soft_coupling_welsch_c=${JOINT_SOFT_COUPLING_WELSCH_C}"); fi
  {
    printf "Build task table:"
    printf " %q" "${build_cmd[@]}"
    printf "\n"
  } >&2
  if truthy "${DRY_RUN}"; then
    :
  else
    load_r_module
    "${build_cmd[@]}"
  fi
fi

if [[ ! -f "${TASKS_TSV}" ]]; then
  echo "Missing task table after build: ${TASKS_TSV}" >&2
  exit 1
fi

TASK_COUNT=$(( $(wc -l < "${TASKS_TSV}") - 1 ))
if (( TASK_COUNT < 1 )); then
  echo "Task table has no data rows: ${TASKS_TSV}" >&2
  exit 1
fi

TABLE_ROOT="$(first_col_value "out_root" || true)"
if is_null_value "${MULTI_WARMUP_ROOT}"; then
  if is_null_value "${TABLE_ROOT}"; then
    MULTI_WARMUP_ROOT="$(cd "$(dirname "${TASKS_TSV}")" && pwd)"
  else
    MULTI_WARMUP_ROOT="$(cd "${TABLE_ROOT}" && pwd)"
  fi
else
  MULTI_WARMUP_ROOT="$(cd "${MULTI_WARMUP_ROOT}" && pwd)"
fi

if [[ -z "${TASK_ARRAY_SCRIPT}" ]]; then
  TASK_ARRAY_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/array_workers/run_multi_warmup_task_table_array.sub"
fi
if [[ -z "${POSTPROCESS_SCRIPT}" ]]; then
  POSTPROCESS_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc/postprocess/postprocess_extra_results.sh"
fi
if [[ -z "${COLLECT_SCRIPT}" ]]; then
  COLLECT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/collect_multi_warmup_results.R"
fi
if [[ -z "${REPORT_SCRIPT}" ]]; then
  REPORT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/multi_warmup_results_report.R"
fi
TASK_ARRAY_SCRIPT="$(cd "$(dirname "${TASK_ARRAY_SCRIPT}")" && pwd)/$(basename "${TASK_ARRAY_SCRIPT}")"
POSTPROCESS_SCRIPT="$(cd "$(dirname "${POSTPROCESS_SCRIPT}")" && pwd)/$(basename "${POSTPROCESS_SCRIPT}")"
COLLECT_SCRIPT="$(cd "$(dirname "${COLLECT_SCRIPT}")" && pwd)/$(basename "${COLLECT_SCRIPT}")"
REPORT_SCRIPT="$(cd "$(dirname "${REPORT_SCRIPT}")" && pwd)/$(basename "${REPORT_SCRIPT}")"

for path in "${TASK_ARRAY_SCRIPT}" "${POSTPROCESS_SCRIPT}" "${COLLECT_SCRIPT}" "${REPORT_SCRIPT}" "${BUILD_TASK_TABLE_SCRIPT}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
done

if [[ -z "${ARRAY_SPEC}" ]]; then
  ARRAY_SPEC="1-${TASK_COUNT}"
fi
if [[ -n "${ARRAY_MAX_CONCURRENT}" && "${ARRAY_SPEC}" != *%* ]]; then
  ARRAY_SPEC="${ARRAY_SPEC}%${ARRAY_MAX_CONCURRENT}"
fi

if [[ -z "${LOG_ROOT}" ]]; then
  LOG_ROOT="${PROJECT_ROOT}/oxygen/results/log"
fi
LOG_ROOT="$(mkdir -p "${LOG_ROOT}" && cd "${LOG_ROOT}" && pwd)"

PROGRESS_LOG="${MULTI_WARMUP_ROOT}/multi_warmup_task_table_progress.log"
JOBS_TSV="${MULTI_WARMUP_ROOT}/multi_warmup_task_table_jobs.tsv"
PAIR_LIST="${MULTI_WARMUP_ROOT}/.multi_warmup_task_table_runs.tsv"
: > "${PROGRESS_LOG}"
printf "job_type\twarmup_label\tjob_id\tdependency\trun_dir\tarray_spec\tqos\twalltime\n" > "${JOBS_TSV}"
SUBMIT_COMMAND_TEXT="$(o2sd_prov_shell_join bash "${BASH_SOURCE[0]}" "${ORIGINAL_SUBMIT_ARGS[@]}")"
o2sd_prov_write_standard "${MULTI_WARMUP_ROOT}" "${BASH_SOURCE[0]}" "${SUBMIT_COMMAND_TEXT}"
o2sd_prov_write_many "${MULTI_WARMUP_ROOT}" \
  scripts submit_script "${BASH_SOURCE[0]}" \
  scripts array_script "${TASK_ARRAY_SCRIPT}" \
  scripts postprocess_script "${POSTPROCESS_SCRIPT}" \
  scripts collect_script "${COLLECT_SCRIPT}" \
  scripts report_script "${REPORT_SCRIPT}" \
  multi_warmup task_table_path "${TASKS_TSV}" \
  multi_warmup array_spec "${ARRAY_SPEC}" \
  multi_warmup task_status_filter "${TASK_STATUS_FILTER}" \
  multi_warmup skip_existing "${SKIP_EXISTING}" \
  optimizer n_cores "${JOINT_N_CORES}" \
  optimizer itermax "${ITERMAX}" \
  optimizer NP "${NP}" \
  optimizer de_reltol "${DE_RELTOL}" \
  optimizer de_steptol "${DE_STEPTOL}" \
  slurm qos "${JOINT_QOS}" \
  slurm walltime "${JOINT_TIME_LIMIT}" \
  slurm memory "${JOINT_MEM}"

log_msg "stage=submit_task_table task_table=${TASKS_TSV} task_count=${TASK_COUNT} array_spec=${ARRAY_SPEC}"
log_msg "stage=submit_task_table options skip_existing=${SKIP_EXISTING} task_status_filter=${TASK_STATUS_FILTER}"

export_arg="ALL,PROJECT_ROOT=${PROJECT_ROOT},TASKS_TSV=${TASKS_TSV},TASK_LOOKUP_COLUMN=${TASK_LOOKUP_COLUMN},R_MODULE=${R_MODULE},SKIP_EXISTING=${SKIP_EXISTING},TASK_STATUS_FILTER=${TASK_STATUS_FILTER}"
array_job_id="$(submit_or_print "Submit multi-warmup task-table array" \
  sbatch \
  "--job-name=${JOB_NAME}" \
  "--array=${ARRAY_SPEC}" \
  "--cpus-per-task=${JOINT_N_CORES}" \
  "--mem=${JOINT_MEM}" \
  "--qos=${JOINT_QOS}" \
  "--time=${JOINT_TIME_LIMIT}" \
  "--output=${LOG_ROOT}/o2mw_task_%A_%a.out" \
  "--error=${LOG_ROOT}/o2mw_task_%A_%a.err" \
  "--export=${export_arg}" \
  "${TASK_ARRAY_SCRIPT}")"
o2sd_prov_record_sbatch "${MULTI_WARMUP_ROOT}" "$(o2sd_prov_shell_join \
  sbatch \
  "--job-name=${JOB_NAME}" \
  "--array=${ARRAY_SPEC}" \
  "--cpus-per-task=${JOINT_N_CORES}" \
  "--mem=${JOINT_MEM}" \
  "--qos=${JOINT_QOS}" \
  "--time=${JOINT_TIME_LIMIT}" \
  "--output=${LOG_ROOT}/o2mw_task_%A_%a.out" \
  "--error=${LOG_ROOT}/o2mw_task_%A_%a.err" \
  "--export=${export_arg}" \
  "${TASK_ARRAY_SCRIPT}")" "${array_job_id}"
o2sd_prov_append "${MULTI_WARMUP_ROOT}" slurm array_job_id "${array_job_id}"
printf "array\tALL\t%s\t\t%s\t%s\t%s\t%s\n" "${array_job_id}" "${MULTI_WARMUP_ROOT}" "${ARRAY_SPEC}" "${JOINT_QOS}" "${JOINT_TIME_LIMIT}" >> "${JOBS_TSV}"
log_msg "stage=submitted_array job_id=${array_job_id}"

post_job_ids=()
if truthy "${SUBMIT_POSTPROCESS}"; then
  write_unique_run_dirs "${PAIR_LIST}"
  while IFS=$'\t' read -r warmup_label run_dir; do
    [[ -z "${run_dir}" ]] && continue
    post_inner_cmd="$(shell_join bash "${POSTPROCESS_SCRIPT}" \
      "--project_root=${PROJECT_ROOT}" \
      "--run_dir=${run_dir}" \
      "--r_module=${R_MODULE}" \
      "--force_extra_results=${FORCE_EXTRA_RESULTS}")"
    post_wrap_cmd="$(shell_join bash -lc "${post_inner_cmd}")"
    post_id="$(submit_or_print "Submit extra_results ${warmup_label}" \
      sbatch \
      "--job-name=o2mw_er_${warmup_label}" \
      "--dependency=afterok:${array_job_id}" \
      "--cpus-per-task=1" \
      "--mem=${POSTPROCESS_MEM}" \
      "--qos=${POSTPROCESS_QOS}" \
      "--time=${POSTPROCESS_TIME_LIMIT}" \
      "--output=${LOG_ROOT}/o2mw_er_${warmup_label}_%j.out" \
      "--error=${LOG_ROOT}/o2mw_er_${warmup_label}_%j.err" \
      "--wrap=${post_wrap_cmd}")"
    post_job_ids+=("${post_id}")
    printf "postprocess\t%s\t%s\tafterok:%s\t%s\t\t%s\t%s\n" "${warmup_label}" "${post_id}" "${array_job_id}" "${run_dir}" "${POSTPROCESS_QOS}" "${POSTPROCESS_TIME_LIMIT}" >> "${JOBS_TSV}"
    o2sd_prov_append "${MULTI_WARMUP_ROOT}" slurm "postprocess_${warmup_label}_job_id" "${post_id}"
    log_msg "stage=submitted_postprocess warmup_label=${warmup_label} job_id=${post_id}"
  done < "${PAIR_LIST}"
fi

if truthy "${SUBMIT_REPORT}"; then
  if (( ${#post_job_ids[@]} > 0 )); then
    dependency_ids="$(IFS=:; echo "${post_job_ids[*]}")"
    report_dependency="afterok:${dependency_ids}"
  else
    report_dependency="afterok:${array_job_id}"
  fi
  manifest_path="${MULTI_WARMUP_ROOT}/multi_warmup_manifest.tsv"
  report_inner_cmd="$(shell_join source "${O2SD_DOCKER_HPC_ROOT}/util/o2_supply_demand_map_apptainer_runtime.sh"); $(shell_join cd "${PROJECT_ROOT}") && $(shell_join Rscript "${COLLECT_SCRIPT}" "--multi_warmup_root=${MULTI_WARMUP_ROOT}" "--manifest=${manifest_path}" "--out_dir=${MULTI_WARMUP_ROOT}") && $(shell_join Rscript "${REPORT_SCRIPT}" "--multi_warmup_root=${MULTI_WARMUP_ROOT}" "--out=${MULTI_WARMUP_ROOT}/multi-warm-up_results.html")"
  wrap_cmd="$(shell_join bash -lc "${report_inner_cmd}")"
  report_job_id="$(submit_or_print "Submit final multi-warmup report" \
    sbatch \
    "--job-name=o2mw_report" \
    "--dependency=${report_dependency}" \
    "--cpus-per-task=1" \
    "--mem=${REPORT_MEM}" \
    "--qos=${REPORT_QOS}" \
    "--time=${REPORT_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_report_%j.out" \
    "--error=${LOG_ROOT}/o2mw_report_%j.err" \
    "--wrap=${wrap_cmd}")"
  printf "report\tALL\t%s\t%s\t%s\t\t%s\t%s\n" "${report_job_id}" "${report_dependency}" "${MULTI_WARMUP_ROOT}" "${REPORT_QOS}" "${REPORT_TIME_LIMIT}" >> "${JOBS_TSV}"
  o2sd_prov_append "${MULTI_WARMUP_ROOT}" slurm report_job_id "${report_job_id}"
  log_msg "stage=submitted_report job_id=${report_job_id} dependency=${report_dependency}"
fi

log_msg "stage=done jobs_tsv=${JOBS_TSV}"
echo "Task-table array job: ${array_job_id}"
echo "Jobs table: ${JOBS_TSV}"
