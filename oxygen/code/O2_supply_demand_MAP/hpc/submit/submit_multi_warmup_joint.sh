#!/bin/bash
# Submit one joint array per selected multi-warm-up pair.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash submit_multi_warmup_joint.sh --invivo_run_dir=DIR --invitro_run_dir=DIR [options]

This launcher:
  1. Builds the pre-fitting joint soft-coupling ratio UMAP and manifest.
  2. Generates one joint soft-coupling parameter table per warm-up pair.
  3. Submits one joint seed array per warm-up pair, in manifest order.
  4. Submits per-pair extra_results jobs.
  5. Submits a final collector/report job that writes multi-warm-up_results.html.

Required:
  --invivo_run_dir=DIR
  --invitro_run_dir=DIR

Common options:
  --project_root=DIR
  --config_path=FILE
  --out_root=DIR
  --multi_warmup_prefix=NAME
  --joint_total_seeds=500
  --joint_array_tasks=500
  --joint_seeds_per_task=1
  --joint_qos=xxlarge
  --joint_time_limit=12:00:00
  --joint_soft_coupling_welsch_c=0.4
  --multi_warmup_invivo_top_n=10  (0 disables in vivo source clustering; not both sides)
  --multi_warmup_invitro_top_n=10 (0 disables in vitro source clustering; not both sides)
  --multi_warmup_umap_seed=1
  --dry_run=TRUE|FALSE
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

require_nonnegative_int() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]]; then
    echo "${name} must be a non-negative integer, got: ${value}" >&2
    exit 2
  fi
}

log_msg() {
  local msg="$1"
  local stamp
  stamp="$(date '+%Y-%m-%d %H:%M:%S')"
  printf '[%s] %s\n' "${stamp}" "${msg}" | tee -a "${PROGRESS_LOG}"
}

print_command() {
  local label="$1"
  shift
  printf "%s:" "${label}" | tee -a "${PROGRESS_LOG}"
  printf " %q" "$@" | tee -a "${PROGRESS_LOG}"
  printf "\n" | tee -a "${PROGRESS_LOG}"
}

run_or_print() {
  local label="$1"
  shift
  print_command "${label}" "$@"
  if ! truthy "${DRY_RUN}"; then
    "$@"
  fi
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h) usage; exit 0 ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --config_path=*|--config=*) CONFIG_PATH="${arg#*=}" ;;
      --out_root=*) OUT_ROOT="${arg#*=}" ;;
      --multi_warmup_prefix=*|--joint_run_prefix=*) MULTI_WARMUP_PREFIX="${arg#*=}" ;;
      --invivo_run_dir=*) INVIVO_RUN_DIR="${arg#*=}" ;;
      --invitro_run_dir=*) INVITRO_RUN_DIR="${arg#*=}" ;;
      --joint_total_seeds=*) JOINT_TOTAL_SEEDS="${arg#*=}" ;;
      --joint_array_tasks=*) JOINT_ARRAY_TASKS="${arg#*=}" ;;
      --joint_seeds_per_task=*) JOINT_SEEDS_PER_TASK="${arg#*=}" ;;
      --joint_n_cores=*|--n_cores=*) JOINT_N_CORES="${arg#*=}" ;;
      --joint_mem=*|--mem=*) JOINT_MEM="${arg#*=}" ;;
      --joint_qos=*|--qos=*) JOINT_QOS="${arg#*=}" ;;
      --joint_time_limit=*|--time=*|--time_limit=*) JOINT_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_qos=*) POSTPROCESS_QOS="${arg#*=}" ;;
      --postprocess_time_limit=*) POSTPROCESS_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_mem=*) POSTPROCESS_MEM="${arg#*=}" ;;
      --report_qos=*) REPORT_QOS="${arg#*=}" ;;
      --report_time_limit=*) REPORT_TIME_LIMIT="${arg#*=}" ;;
      --report_mem=*) REPORT_MEM="${arg#*=}" ;;
      --parameter_table=*|--invitro_parameter_table=*) PARAMETER_TABLE="${arg#*=}" ;;
      --fit_objects_dir=*) FIT_OBJECTS_DIR="${arg#*=}" ;;
      --flow_density_path=*) FLOW_DENSITY_PATH="${arg#*=}" ;;
      --itermax=*) ITERMAX="${arg#*=}" ;;
      --de_reltol=*) DE_RELTOL="${arg#*=}" ;;
      --de_steptol=*) DE_STEPTOL="${arg#*=}" ;;
      --np=*|--NP=*) NP="${arg#*=}" ;;
      --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --force_extra_results=*) FORCE_EXTRA_RESULTS="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      --joint_warmup_sigmaN=*) JOINT_WARMUP_SIGMAN="${arg#*=}" ;;
      --joint_soft_coupling_sigma_default=*) JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${arg#*=}" ;;
      --joint_soft_coupling_welsch_c=*) JOINT_SOFT_COUPLING_WELSCH_C="${arg#*=}" ;;
      --joint_soft_coupling_delta_params=*) JOINT_SOFT_COUPLING_DELTA_PARAMS="${arg#*=}" ;;
      --multi_warmup_top_n=*) MULTI_WARMUP_TOP_N="${arg#*=}" ;;
      --multi_warmup_invivo_top_n=*|--invivo_top_n=*) MULTI_WARMUP_INVIVO_TOP_N="${arg#*=}" ;;
      --multi_warmup_invitro_top_n=*|--invitro_top_n=*) MULTI_WARMUP_INVITRO_TOP_N="${arg#*=}" ;;
      --multi_warmup_umap_seed=*|--umap_seed=*) MULTI_WARMUP_UMAP_SEED="${arg#*=}" ;;
      --multi_warmup_invivo_k=*) MULTI_WARMUP_INVIVO_K="${arg#*=}" ;;
      --multi_warmup_invitro_k=*) MULTI_WARMUP_INVITRO_K="${arg#*=}" ;;
      --multi_warmup_invitro_anchor_ranks=*) MULTI_WARMUP_INVITRO_ANCHOR_RANKS="${arg#*=}" ;;
      --multi_warmup_include_phase2=*) MULTI_WARMUP_INCLUDE_PHASE2="${arg#*=}" ;;
      --multi_warmup_phase2_invitro_anchor_ranks=*) MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="${arg#*=}" ;;
      --log_root=*|--log_dir=*) LOG_ROOT="${arg#*=}" ;;
      *) echo "Unknown argument: ${arg}" >&2; usage >&2; exit 2 ;;
    esac
  done
}

load_r_module() {
  if command -v ml >/dev/null 2>&1; then
    ml "${R_MODULE}"
  elif command -v module >/dev/null 2>&1; then
    module load "${R_MODULE}"
  fi
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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
CONFIG_PATH="${CONFIG_PATH:-}"
OUT_ROOT="${OUT_ROOT:-}"
MULTI_WARMUP_PREFIX="${MULTI_WARMUP_PREFIX:-fit_joint_multi_warmup}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-}"
JOINT_TOTAL_SEEDS="${JOINT_TOTAL_SEEDS:-500}"
JOINT_ARRAY_TASKS="${JOINT_ARRAY_TASKS:-500}"
JOINT_SEEDS_PER_TASK="${JOINT_SEEDS_PER_TASK:-1}"
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
PARAMETER_TABLE="${PARAMETER_TABLE:-}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-}"
ITERMAX="${ITERMAX:-500}"
DE_RELTOL="${DE_RELTOL:-1e-4}"
DE_STEPTOL="${DE_STEPTOL:-25}"
NP="${NP:-80}"
AUTO_VIZ="${AUTO_VIZ:-TRUE}"
R_MODULE="${R_MODULE:-R/4.4}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-FALSE}"
DRY_RUN="${DRY_RUN:-FALSE}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-}"
JOINT_SOFT_COUPLING_WELSCH_C="${JOINT_SOFT_COUPLING_WELSCH_C:-}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-default}"
MULTI_WARMUP_TOP_N="${MULTI_WARMUP_TOP_N:-10}"
MULTI_WARMUP_INVIVO_TOP_N="${MULTI_WARMUP_INVIVO_TOP_N:-${MULTI_WARMUP_TOP_N}}"
MULTI_WARMUP_INVITRO_TOP_N="${MULTI_WARMUP_INVITRO_TOP_N:-${MULTI_WARMUP_TOP_N}}"
MULTI_WARMUP_UMAP_SEED="${MULTI_WARMUP_UMAP_SEED:-1}"
MULTI_WARMUP_INVIVO_K="${MULTI_WARMUP_INVIVO_K:-auto}"
MULTI_WARMUP_INVITRO_K="${MULTI_WARMUP_INVITRO_K:-auto}"
MULTI_WARMUP_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_INVITRO_ANCHOR_RANKS:-1}"
MULTI_WARMUP_INCLUDE_PHASE2="${MULTI_WARMUP_INCLUDE_PHASE2:-FALSE}"
MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS="${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS:-auto}"
LOG_ROOT="${LOG_ROOT:-}"

parse_args "$@"

require_nonnegative_int MULTI_WARMUP_INVIVO_TOP_N "${MULTI_WARMUP_INVIVO_TOP_N}"
require_nonnegative_int MULTI_WARMUP_INVITRO_TOP_N "${MULTI_WARMUP_INVITRO_TOP_N}"
if (( MULTI_WARMUP_INVIVO_TOP_N == 0 && MULTI_WARMUP_INVITRO_TOP_N == 0 )); then
  echo "At least one of --invivo_top_n or --invitro_top_n must be greater than 0." >&2
  exit 2
fi
if (( MULTI_WARMUP_INVIVO_TOP_N > 0 )) && is_null_value "${INVIVO_RUN_DIR}"; then
  echo "HPC multi-warm-up submission requires existing --invivo_run_dir when --invivo_top_n > 0." >&2
  echo "Submit the source in vivo fit first, then rerun this launcher." >&2
  exit 2
fi
if (( MULTI_WARMUP_INVITRO_TOP_N > 0 )) && is_null_value "${INVITRO_RUN_DIR}"; then
  echo "HPC multi-warm-up submission requires existing --invitro_run_dir when --invitro_top_n > 0." >&2
  echo "Submit the source in vitro fit first, then rerun this launcher." >&2
  exit 2
fi
if (( JOINT_ARRAY_TASKS * JOINT_SEEDS_PER_TASK != JOINT_TOTAL_SEEDS )); then
  echo "joint_array_tasks * joint_seeds_per_task must equal joint_total_seeds." >&2
  exit 2
fi
if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found; run on an HPC login node or use --dry_run=TRUE." >&2
  exit 1
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
if [[ -z "${CONFIG_PATH}" ]]; then CONFIG_PATH="${PROJECT_ROOT}/oxygen/config/O2_supply_demand.yaml"; fi
if [[ -z "${OUT_ROOT}" ]]; then OUT_ROOT="${PROJECT_ROOT}/oxygen/results"; fi
if [[ -z "${PARAMETER_TABLE}" ]]; then PARAMETER_TABLE="${PROJECT_ROOT}/oxygen/data/O2_supply_demand/parameter_table_invitro_buffering.csv"; fi
if [[ -z "${FIT_OBJECTS_DIR}" ]]; then FIT_OBJECTS_DIR="${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects"; fi
if [[ -z "${FLOW_DENSITY_PATH}" ]]; then FLOW_DENSITY_PATH="${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv"; fi
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
OUT_ROOT="$(mkdir -p "${OUT_ROOT}" && cd "${OUT_ROOT}" && pwd)"
if (( MULTI_WARMUP_INVIVO_TOP_N > 0 )); then INVIVO_RUN_DIR="$(cd "${INVIVO_RUN_DIR}" && pwd)"; else INVIVO_RUN_DIR=""; fi
if (( MULTI_WARMUP_INVITRO_TOP_N > 0 )); then INVITRO_RUN_DIR="$(cd "${INVITRO_RUN_DIR}" && pwd)"; else INVITRO_RUN_DIR=""; fi
PARAMETER_TABLE="$(cd "$(dirname "${PARAMETER_TABLE}")" && pwd)/$(basename "${PARAMETER_TABLE}")"
FIT_OBJECTS_DIR="$(cd "${FIT_OBJECTS_DIR}" && pwd)"
FLOW_DENSITY_PATH="$(cd "$(dirname "${FLOW_DENSITY_PATH}")" && pwd)/$(basename "${FLOW_DENSITY_PATH}")"
if [[ -z "${LOG_ROOT}" ]]; then LOG_ROOT="${OUT_ROOT}/log"; fi
LOG_ROOT="$(mkdir -p "${LOG_ROOT}" && cd "${LOG_ROOT}" && pwd)"

MULTI_WARMUP_ROOT="${OUT_ROOT}/${MULTI_WARMUP_PREFIX}"
mkdir -p "${MULTI_WARMUP_ROOT}"
PROGRESS_LOG="${MULTI_WARMUP_ROOT}/multi_warmup_progress.log"
JOBS_TSV="${MULTI_WARMUP_ROOT}/multi_warmup_jobs.tsv"
: > "${PROGRESS_LOG}"

SEED_PLAN_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/build_multi_warmup_seed_plan.R"
MAKE_TABLE_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/warm_start/make_joint_soft_coupling_parameters_table.R"
COLLECT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/collect_multi_warmup_results.R"
REPORT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/multi_warmup/multi_warmup_results_report.R"
JOINT_ARRAY_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/hpc/array_workers/submit_fit_seed_array_joint_buffering.sub"
JOINT_RUNNER_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh"
POSTPROCESS_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/hpc/postprocess/postprocess_extra_results.sh"

for path in "${SEED_PLAN_SCRIPT}" "${MAKE_TABLE_SCRIPT}" "${COLLECT_SCRIPT}" "${REPORT_SCRIPT}" "${JOINT_ARRAY_SCRIPT}" "${JOINT_RUNNER_SCRIPT}" "${POSTPROCESS_SCRIPT}" "${CONFIG_PATH}" "${PARAMETER_TABLE}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
done

load_r_module
log_msg "stage=seed_plan invivo_run_dir=${INVIVO_RUN_DIR} invitro_run_dir=${INVITRO_RUN_DIR}"
run_or_print "Generate multi-warmup seed plan" \
  Rscript "${SEED_PLAN_SCRIPT}" \
  "--invivo_run_dir=${INVIVO_RUN_DIR}" \
  "--invitro_run_dir=${INVITRO_RUN_DIR}" \
  "--out_dir=${MULTI_WARMUP_ROOT}" \
  "--top_n=${MULTI_WARMUP_TOP_N}" \
  "--invivo_top_n=${MULTI_WARMUP_INVIVO_TOP_N}" \
  "--invitro_top_n=${MULTI_WARMUP_INVITRO_TOP_N}" \
  "--umap_seed=${MULTI_WARMUP_UMAP_SEED}" \
  "--invivo_k=${MULTI_WARMUP_INVIVO_K}" \
  "--invitro_k=${MULTI_WARMUP_INVITRO_K}" \
  "--invitro_anchor_ranks=${MULTI_WARMUP_INVITRO_ANCHOR_RANKS}" \
  "--include_phase2=${MULTI_WARMUP_INCLUDE_PHASE2}" \
  "--phase2_invitro_anchor_ranks=${MULTI_WARMUP_PHASE2_INVITRO_ANCHOR_RANKS}"

MANIFEST="${MULTI_WARMUP_ROOT}/multi_warmup_manifest.tsv"
PLAN_MODE_FILE="${MULTI_WARMUP_ROOT}/multi_warmup_seed_plan_mode.tsv"
if ! truthy "${DRY_RUN}" && [[ ! -f "${MANIFEST}" ]]; then
  echo "Missing generated manifest: ${MANIFEST}" >&2
  exit 1
fi
if truthy "${DRY_RUN}" && [[ ! -f "${MANIFEST}" ]]; then
  log_msg "stage=dry_run complete_after_seed_plan_command"
  exit 0
fi

printf "warmup_label\tpair_index\ttotal_pairs\tjoint_run_dir\tsubmit_status\tjob_id\tpostprocess_job_id\tarray_spec\tqos\twalltime\n" > "${JOBS_TSV}"
total_pairs=$(( $(wc -l < "${MANIFEST}") - 1 ))
if (( total_pairs < 1 )); then
  plan_mode=""
  if [[ -f "${PLAN_MODE_FILE}" ]]; then
    plan_mode="$(awk -F $'\t' '$1 == "mode" {print $2; exit}' "${PLAN_MODE_FILE}")"
  fi
  case "${plan_mode}" in
    invivo_only|invitro_only)
      log_msg "stage=cluster_only mode=${plan_mode} manifest=${MANIFEST}"
      exit 0
      ;;
    *)
      echo "Generated manifest has no warm-up pairs: ${MANIFEST}" >&2
      exit 1
      ;;
  esac
fi
pair_index=0
post_job_ids=()
POSTPROCESS_JOB_IDS_FILE="${MULTI_WARMUP_ROOT}/.postprocess_job_ids.tmp"
: > "${POSTPROCESS_JOB_IDS_FILE}"
log_msg "stage=submit_pairs total_pairs=${total_pairs} qos=${JOINT_QOS} time=${JOINT_TIME_LIMIT}"
mkdir -p "${MULTI_WARMUP_ROOT}/joint_soft_coupling_tables"

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r warmup_label phase invivo_family invivo_rank invivo_seed invivo_seed_dir invitro_family invitro_rank invitro_seed invitro_seed_dir selection_reason joint_run_prefix joint_table; do
  pair_index=$((pair_index + 1))
  log_msg "stage=pair_submit_start pair=${pair_index}/${total_pairs} warmup_label=${warmup_label}"
  if [[ -z "${joint_table}" ]]; then
    echo "Manifest row for ${warmup_label} is missing joint_soft_coupling_parameters_table." >&2
    exit 1
  fi
  case "${joint_table}" in
    /*) ;;
    *) joint_table="${MULTI_WARMUP_ROOT}/${joint_table}" ;;
  esac
  mkdir -p "$(dirname "${joint_table}")"
  run_or_print "Generate pair soft-coupling table" \
    Rscript "${MAKE_TABLE_SCRIPT}" \
    "--invivo-seed-dir=${invivo_seed_dir}" \
    "--invitro-seed-dir=${invitro_seed_dir}" \
    "--seed-label=${warmup_label}" \
    "--out=${joint_table}" \
    "--delta-params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}"
  if ! truthy "${DRY_RUN}" && [[ ! -f "${joint_table}" ]]; then
    echo "Pair soft-coupling table was not created for ${warmup_label}: ${joint_table}" >&2
    exit 1
  fi

  export_arg="ALL,PROJECT_ROOT=${PROJECT_ROOT},RUNNER_SCRIPT=${JOINT_RUNNER_SCRIPT},CONFIG_PATH=${CONFIG_PATH},OUT_ROOT=${MULTI_WARMUP_ROOT},RUN_PREFIX=${joint_run_prefix},TOTAL_SEEDS=${JOINT_TOTAL_SEEDS},ARRAY_TASKS=${JOINT_ARRAY_TASKS},SEEDS_PER_TASK=${JOINT_SEEDS_PER_TASK},N_CORES=${JOINT_N_CORES},R_MODULE=${R_MODULE},PARAMETER_TABLE=${PARAMETER_TABLE},FIT_OBJECTS_DIR=${FIT_OBJECTS_DIR},FLOW_DENSITY_PATH=${FLOW_DENSITY_PATH},ITERMAX=${ITERMAX},DE_RELTOL=${DE_RELTOL},DE_STEPTOL=${DE_STEPTOL},NP=${NP},AUTO_VIZ=${AUTO_VIZ},JOINT_WARMUP_ENABLE=TRUE,JOINT_WARMUP_SEED_LABEL=${warmup_label},JOINT_WARMUP_INVIVO_SEED_DIR=${invivo_seed_dir},JOINT_WARMUP_INVITRO_SEED_DIR=${invitro_seed_dir},JOINT_WARMUP_SIGMAN=${JOINT_WARMUP_SIGMAN},JOINT_SOFT_COUPLING_SIGMA_DEFAULT=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT},JOINT_SOFT_COUPLING_WELSCH_C=${JOINT_SOFT_COUPLING_WELSCH_C},JOINT_SOFT_COUPLING_PARAMETERS_TABLE=${joint_table}"
  job_id="$(submit_or_print "Submit pair joint array ${warmup_label}" \
    sbatch \
    "--job-name=o2mw_${warmup_label}" \
    "--array=1-${JOINT_ARRAY_TASKS}" \
    "--cpus-per-task=${JOINT_N_CORES}" \
    "--mem=${JOINT_MEM}" \
    "--qos=${JOINT_QOS}" \
    "--time=${JOINT_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_${warmup_label}_%A_%a.out" \
    "--error=${LOG_ROOT}/o2mw_${warmup_label}_%A_%a.err" \
    "--export=${export_arg}" \
    "${JOINT_ARRAY_SCRIPT}")"

  post_id="$(submit_or_print "Submit pair extra_results ${warmup_label}" \
    sbatch \
    "--job-name=o2mw_er_${warmup_label}" \
    "--dependency=afterok:${job_id}" \
    "--cpus-per-task=1" \
    "--mem=${POSTPROCESS_MEM}" \
    "--qos=${POSTPROCESS_QOS}" \
    "--time=${POSTPROCESS_TIME_LIMIT}" \
    "--output=${LOG_ROOT}/o2mw_er_${warmup_label}_%j.out" \
    "--error=${LOG_ROOT}/o2mw_er_${warmup_label}_%j.err" \
    "--export=ALL,PROJECT_ROOT=${PROJECT_ROOT},RUN_DIR=${MULTI_WARMUP_ROOT}/${joint_run_prefix},R_MODULE=${R_MODULE},FORCE_EXTRA_RESULTS=${FORCE_EXTRA_RESULTS}" \
    "${POSTPROCESS_SCRIPT}")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${warmup_label}" "${pair_index}" "${total_pairs}" "${MULTI_WARMUP_ROOT}/${joint_run_prefix}" "submitted" "${job_id}" "${post_id}" "1-${JOINT_ARRAY_TASKS}" "${JOINT_QOS}" "${JOINT_TIME_LIMIT}" >> "${JOBS_TSV}"
  echo "${post_id}" >> "${POSTPROCESS_JOB_IDS_FILE}"
  log_msg "stage=pair_submit_done pair=${pair_index}/${total_pairs} warmup_label=${warmup_label} array_job=${job_id} postprocess_job=${post_id}"
done

dependency_ids="$(paste -sd: "${POSTPROCESS_JOB_IDS_FILE}")"
wrap_cmd="if command -v ml >/dev/null 2>&1; then ml ${R_MODULE}; elif command -v module >/dev/null 2>&1; then module load ${R_MODULE}; fi; cd ${PROJECT_ROOT} && Rscript ${COLLECT_SCRIPT} --multi_warmup_root=${MULTI_WARMUP_ROOT} --manifest=${MANIFEST} --out_dir=${MULTI_WARMUP_ROOT} && Rscript ${REPORT_SCRIPT} --multi_warmup_root=${MULTI_WARMUP_ROOT} --out=${MULTI_WARMUP_ROOT}/multi-warm-up_results.html"
report_job_id="$(submit_or_print "Submit final multi-warmup report" \
  sbatch \
  "--job-name=o2mw_report" \
  "--dependency=afterok:${dependency_ids}" \
  "--cpus-per-task=1" \
  "--mem=${REPORT_MEM}" \
  "--qos=${REPORT_QOS}" \
  "--time=${REPORT_TIME_LIMIT}" \
  "--output=${LOG_ROOT}/o2mw_report_%j.out" \
  "--error=${LOG_ROOT}/o2mw_report_%j.err" \
  "--wrap=${wrap_cmd}")"
log_msg "stage=submitted_all final_report_job=${report_job_id} report=${MULTI_WARMUP_ROOT}/multi-warm-up_results.html"
