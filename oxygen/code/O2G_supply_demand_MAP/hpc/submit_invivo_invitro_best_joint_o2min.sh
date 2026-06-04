#!/bin/bash
# Submit in vivo + in vitro fits, select each best seed, then submit warm-started
# soft-coupled joint fitting. The active config is generated from a base YAML
# with o2_min set to 0.5 by default.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT_DEFAULT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

usage() {
  cat <<'EOF'
Usage:
  bash submit_invivo_invitro_best_joint_o2min.sh [options]

Main behavior:
  1. Generate an active config from --base_config with o2_min=0.5.
  2. Submit an interleaved in vivo / in vitro array.
  3. Run dependent extra_results jobs for both single-fit result folders.
  4. Run a dependent prep job that selects each best seed and submits the
     soft-coupled joint fit using those seed directories as warm start.

Common options:
  --project_root=/path/to/repo
  --base_config=/path/to/O2G_supply_demand.yaml
  --config_path=/path/to/generated_config.yaml
  --write_config=TRUE|FALSE
  --o2_min=0.5
  --out_root=/path/to/oxygen/results
  --log_root=/path/to/oxygen/results/log
  --r_module=R/4.4
  --dry_run=TRUE|FALSE

Single-fit options:
  --total_seeds_per_mode=500
  --array_tasks=1000
  --run_prefix_invivo=fit_invivo_O2G_buffering_o2min0p5_500seed
  --run_prefix_invitro=fit_invitro_O2G_buffering_o2min0p5_500seed
  --qos=xxlarge --time_limit=12:00:00 --n_cores=22 --mem=32G

Postprocess / selection options:
  --force_extra_results=TRUE|FALSE
  --select_required_files=best_params.tsv
  --invivo_objective_columns=objective
  --invitro_objective_columns=objective_total,objective

Joint options:
  --joint_run_prefix=fit_joint_O2G_best_o2min0p5_500seed
  --joint_job_name=o2g_joint_best
  --joint_total_seeds=500 --joint_array_tasks=500 --joint_seeds_per_task=1
  --joint_qos=xxlarge --joint_time_limit=12:00:00
  --joint_n_cores=22 --joint_mem=32G
  --joint_soft_coupling_sigma_default=1.5
  --joint_warmup_sigmaN=0.0304
  --joint_soft_coupling_delta_params=default|all|none|param1,param2

The prep stage is normally submitted by this wrapper. Do not use
--internal_stage unless debugging the dependent stage manually.
EOF
}

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
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

require_numeric() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]; then
    echo "${name} must be numeric, got: ${value}" >&2
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

o2_min_label() {
  local label="${1:-0.5}"
  label="${label//./p}"
  label="${label//-/m}"
  label="${label//+/}"
  sanitize_label "${label}"
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

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h)
        usage
        exit 0
        ;;
      --internal_stage=*) INTERNAL_STAGE="${arg#*=}" ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --base_config=*|--config_template=*) BASE_CONFIG="${arg#*=}" ;;
      --config_path=*|--generated_config_path=*) CONFIG_PATH="${arg#*=}" ;;
      --write_config=*) WRITE_CONFIG="${arg#*=}" ;;
      --o2_min=*) O2_MIN="${arg#*=}" ;;
      --out_root=*) OUT_ROOT="${arg#*=}" ;;
      --log_root=*|--log_dir=*) LOG_ROOT="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      --force_extra_results=*) FORCE_EXTRA_RESULTS="${arg#*=}" ;;
      --total_seeds_per_mode=*) TOTAL_SEEDS_PER_MODE="${arg#*=}" ;;
      --array_tasks=*) ARRAY_TASKS="${arg#*=}" ;;
      --run_prefix_invivo=*) RUN_PREFIX_INVIVO="${arg#*=}" ;;
      --run_prefix_invitro=*) RUN_PREFIX_INVITRO="${arg#*=}" ;;
      --qos=*) QOS="${arg#*=}" ;;
      --time_limit=*) TIME_LIMIT="${arg#*=}" ;;
      --n_cores=*) N_CORES="${arg#*=}" ;;
      --mem=*) MEM="${arg#*=}" ;;
      --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
      --glucose=*) GLUCOSE="${arg#*=}" ;;
      --itermax=*) ITERMAX="${arg#*=}" ;;
      --de_reltol=*) DE_RELTOL="${arg#*=}" ;;
      --de_steptol=*) DE_STEPTOL="${arg#*=}" ;;
      --np=*|--NP=*) NP="${arg#*=}" ;;
      --parameter_table=*|--invitro_parameter_table=*) PARAMETER_TABLE="${arg#*=}" ;;
      --fit_objects_dir=*) FIT_OBJECTS_DIR="${arg#*=}" ;;
      --flow_density_path=*) FLOW_DENSITY_PATH="${arg#*=}" ;;
      --postprocess_qos=*) POSTPROCESS_QOS="${arg#*=}" ;;
      --postprocess_time_limit=*) POSTPROCESS_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_mem=*) POSTPROCESS_MEM="${arg#*=}" ;;
      --prep_qos=*) PREP_QOS="${arg#*=}" ;;
      --prep_time_limit=*) PREP_TIME_LIMIT="${arg#*=}" ;;
      --prep_mem=*) PREP_MEM="${arg#*=}" ;;
      --prep_prefix=*) PREP_PREFIX="${arg#*=}" ;;
      --select_required_files=*) SELECT_REQUIRED_FILES="${arg#*=}" ;;
      --invivo_objective_columns=*) INVIVO_OBJECTIVE_COLUMNS="${arg#*=}" ;;
      --invitro_objective_columns=*) INVITRO_OBJECTIVE_COLUMNS="${arg#*=}" ;;
      --joint_run_prefix=*) JOINT_RUN_PREFIX="${arg#*=}" ;;
      --joint_job_name=*) JOINT_JOB_NAME="${arg#*=}" ;;
      --joint_total_seeds=*) JOINT_TOTAL_SEEDS="${arg#*=}" ;;
      --joint_array_tasks=*) JOINT_ARRAY_TASKS="${arg#*=}" ;;
      --joint_seeds_per_task=*) JOINT_SEEDS_PER_TASK="${arg#*=}" ;;
      --joint_qos=*) JOINT_QOS="${arg#*=}" ;;
      --joint_time_limit=*) JOINT_TIME_LIMIT="${arg#*=}" ;;
      --joint_n_cores=*) JOINT_N_CORES="${arg#*=}" ;;
      --joint_mem=*) JOINT_MEM="${arg#*=}" ;;
      --joint_soft_coupling_sigma_default=*) JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${arg#*=}" ;;
      --joint_warmup_sigmaN=*) JOINT_WARMUP_SIGMAN="${arg#*=}" ;;
      --joint_soft_coupling_delta_params=*) JOINT_SOFT_COUPLING_DELTA_PARAMS="${arg#*=}" ;;
      *)
        echo "Unknown argument: ${arg}" >&2
        usage >&2
        exit 2
        ;;
    esac
  done
}

PROJECT_ROOT="${PROJECT_ROOT:-${PROJECT_ROOT_DEFAULT}}"
BASE_CONFIG="${BASE_CONFIG:-}"
CONFIG_PATH="${CONFIG_PATH:-}"
WRITE_CONFIG="${WRITE_CONFIG:-TRUE}"
O2_MIN="${O2_MIN:-0.5}"
OUT_ROOT="${OUT_ROOT:-}"
LOG_ROOT="${LOG_ROOT:-}"
R_MODULE="${R_MODULE:-R/4.4}"
DRY_RUN="${DRY_RUN:-FALSE}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-TRUE}"
INTERNAL_STAGE="${INTERNAL_STAGE:-submit}"

TOTAL_SEEDS_PER_MODE="${TOTAL_SEEDS_PER_MODE:-500}"
ARRAY_TASKS="${ARRAY_TASKS:-}"
RUN_PREFIX_INVIVO="${RUN_PREFIX_INVIVO:-}"
RUN_PREFIX_INVITRO="${RUN_PREFIX_INVITRO:-}"
QOS="${QOS:-xxlarge}"
TIME_LIMIT="${TIME_LIMIT:-12:00:00}"
N_CORES="${N_CORES:-22}"
MEM="${MEM:-32G}"
AUTO_VIZ="${AUTO_VIZ:-TRUE}"
GLUCOSE="${GLUCOSE:-TRUE}"
ITERMAX="${ITERMAX:-500}"
DE_RELTOL="${DE_RELTOL:-1e-4}"
DE_STEPTOL="${DE_STEPTOL:-25}"
NP="${NP:-80}"
PARAMETER_TABLE="${PARAMETER_TABLE:-}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-}"

POSTPROCESS_QOS="${POSTPROCESS_QOS:-small}"
POSTPROCESS_TIME_LIMIT="${POSTPROCESS_TIME_LIMIT:-4:00:00}"
POSTPROCESS_MEM="${POSTPROCESS_MEM:-8G}"
PREP_QOS="${PREP_QOS:-small}"
PREP_TIME_LIMIT="${PREP_TIME_LIMIT:-2:00:00}"
PREP_MEM="${PREP_MEM:-8G}"
PREP_PREFIX="${PREP_PREFIX:-}"
SELECT_REQUIRED_FILES="${SELECT_REQUIRED_FILES:-best_params.tsv}"
INVIVO_OBJECTIVE_COLUMNS="${INVIVO_OBJECTIVE_COLUMNS:-objective}"
INVITRO_OBJECTIVE_COLUMNS="${INVITRO_OBJECTIVE_COLUMNS:-objective_total,objective}"

JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-}"
JOINT_JOB_NAME="${JOINT_JOB_NAME:-o2g_joint_best}"
JOINT_TOTAL_SEEDS="${JOINT_TOTAL_SEEDS:-500}"
JOINT_ARRAY_TASKS="${JOINT_ARRAY_TASKS:-}"
JOINT_SEEDS_PER_TASK="${JOINT_SEEDS_PER_TASK:-1}"
JOINT_QOS="${JOINT_QOS:-xxlarge}"
JOINT_TIME_LIMIT="${JOINT_TIME_LIMIT:-12:00:00}"
JOINT_N_CORES="${JOINT_N_CORES:-22}"
JOINT_MEM="${JOINT_MEM:-32G}"
JOINT_SOFT_COUPLING_SIGMA_DEFAULT="${JOINT_SOFT_COUPLING_SIGMA_DEFAULT:-}"
JOINT_WARMUP_SIGMAN="${JOINT_WARMUP_SIGMAN:-}"
JOINT_SOFT_COUPLING_DELTA_PARAMS="${JOINT_SOFT_COUPLING_DELTA_PARAMS:-default}"

parse_args "$@"

require_numeric O2_MIN "${O2_MIN}"
require_numeric DE_RELTOL "${DE_RELTOL}"
for name in TOTAL_SEEDS_PER_MODE N_CORES ITERMAX DE_STEPTOL NP \
            JOINT_TOTAL_SEEDS JOINT_SEEDS_PER_TASK JOINT_N_CORES; do
  require_positive_int "${name}" "${!name}"
done
if [[ -z "${ARRAY_TASKS}" ]]; then
  ARRAY_TASKS=$((TOTAL_SEEDS_PER_MODE * 2))
fi
if [[ -z "${JOINT_ARRAY_TASKS}" ]]; then
  JOINT_ARRAY_TASKS="${JOINT_TOTAL_SEEDS}"
fi
require_positive_int ARRAY_TASKS "${ARRAY_TASKS}"
require_positive_int JOINT_ARRAY_TASKS "${JOINT_ARRAY_TASKS}"
if (( ARRAY_TASKS != TOTAL_SEEDS_PER_MODE * 2 )); then
  echo "ARRAY_TASKS must equal TOTAL_SEEDS_PER_MODE * 2." >&2
  echo "Got ARRAY_TASKS=${ARRAY_TASKS}, TOTAL_SEEDS_PER_MODE=${TOTAL_SEEDS_PER_MODE}." >&2
  exit 2
fi
if (( JOINT_ARRAY_TASKS * JOINT_SEEDS_PER_TASK != JOINT_TOTAL_SEEDS )); then
  echo "JOINT_ARRAY_TASKS * JOINT_SEEDS_PER_TASK must equal JOINT_TOTAL_SEEDS." >&2
  exit 2
fi
if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then
  require_numeric JOINT_SOFT_COUPLING_SIGMA_DEFAULT "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}"
fi
if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then
  require_numeric JOINT_WARMUP_SIGMAN "${JOINT_WARMUP_SIGMAN}"
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
if [[ -z "${BASE_CONFIG}" ]]; then
  BASE_CONFIG="${PROJECT_ROOT}/oxygen/config/O2G_supply_demand.yaml"
fi
BASE_CONFIG="$(cd "$(dirname "${BASE_CONFIG}")" && pwd)/$(basename "${BASE_CONFIG}")"
O2_LABEL="$(o2_min_label "${O2_MIN}")"
if [[ -z "${CONFIG_PATH}" ]]; then
  CONFIG_PATH="${PROJECT_ROOT}/oxygen/config/generated/O2G_supply_demand_o2min${O2_LABEL}.yaml"
fi
if [[ -z "${OUT_ROOT}" ]]; then
  OUT_ROOT="${PROJECT_ROOT}/oxygen/results"
fi
if [[ -z "${RUN_PREFIX_INVIVO}" ]]; then
  RUN_PREFIX_INVIVO="fit_invivo_O2G_buffering_o2min${O2_LABEL}_${TOTAL_SEEDS_PER_MODE}seed"
fi
if [[ -z "${RUN_PREFIX_INVITRO}" ]]; then
  RUN_PREFIX_INVITRO="fit_invitro_O2G_buffering_o2min${O2_LABEL}_${TOTAL_SEEDS_PER_MODE}seed"
fi
if [[ -z "${JOINT_RUN_PREFIX}" ]]; then
  JOINT_RUN_PREFIX="fit_joint_O2G_best_o2min${O2_LABEL}_${JOINT_TOTAL_SEEDS}seed"
fi
if [[ -z "${PREP_PREFIX}" ]]; then
  PREP_PREFIX="prep_joint_best_o2min${O2_LABEL}_${TOTAL_SEEDS_PER_MODE}seed"
fi

CONFIG_DIR="$(dirname "${CONFIG_PATH}")"
mkdir -p "${CONFIG_DIR}"
CONFIG_PATH="$(cd "${CONFIG_DIR}" && pwd)/$(basename "${CONFIG_PATH}")"
mkdir -p "${OUT_ROOT}" "${LOG_ROOT:-${OUT_ROOT}/log}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"
if [[ -z "${LOG_ROOT}" ]]; then
  LOG_ROOT="${OUT_ROOT}/log"
fi
LOG_ROOT="$(cd "${LOG_ROOT}" && pwd)"

SUB_SCRIPT="${SUB_SCRIPT:-${SCRIPT_DIR}/submit_fit_seed_array_invivo_invitro_interleaved.sub}"
POSTPROCESS_SCRIPT="${POSTPROCESS_SCRIPT:-${SCRIPT_DIR}/postprocess_extra_results.sh}"
SELECT_BEST_SCRIPT="${SELECT_BEST_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/select_best_seed_from_summary.R}"
UNIFIED_SUBMIT="${UNIFIED_SUBMIT:-${SCRIPT_DIR}/submit_o2g_fit.sh}"
EXTRA_RESULTS_SCRIPT="${EXTRA_RESULTS_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R}"
RUNNER_SCRIPT="${RUNNER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh}"
PARAMETER_TABLE="${PARAMETER_TABLE:-${PROJECT_ROOT}/oxygen/data/O2G_supply_demand/parameter_table_invitro_buffering.csv}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv}"

for var_name in SUB_SCRIPT POSTPROCESS_SCRIPT SELECT_BEST_SCRIPT UNIFIED_SUBMIT EXTRA_RESULTS_SCRIPT RUNNER_SCRIPT PARAMETER_TABLE; do
  path_value="${!var_name}"
  if [[ ! -f "${path_value}" ]]; then
    echo "Missing required file (${var_name}): ${path_value}" >&2
    exit 1
  fi
done
if [[ ! -f "${BASE_CONFIG}" ]]; then
  echo "Missing base config: ${BASE_CONFIG}" >&2
  exit 1
fi
if [[ ! -d "${FIT_OBJECTS_DIR}" ]]; then
  echo "Missing fit_objects_dir: ${FIT_OBJECTS_DIR}" >&2
  exit 1
fi

SUB_SCRIPT="$(cd "$(dirname "${SUB_SCRIPT}")" && pwd)/$(basename "${SUB_SCRIPT}")"
POSTPROCESS_SCRIPT="$(cd "$(dirname "${POSTPROCESS_SCRIPT}")" && pwd)/$(basename "${POSTPROCESS_SCRIPT}")"
SELECT_BEST_SCRIPT="$(cd "$(dirname "${SELECT_BEST_SCRIPT}")" && pwd)/$(basename "${SELECT_BEST_SCRIPT}")"
UNIFIED_SUBMIT="$(cd "$(dirname "${UNIFIED_SUBMIT}")" && pwd)/$(basename "${UNIFIED_SUBMIT}")"
EXTRA_RESULTS_SCRIPT="$(cd "$(dirname "${EXTRA_RESULTS_SCRIPT}")" && pwd)/$(basename "${EXTRA_RESULTS_SCRIPT}")"
RUNNER_SCRIPT="$(cd "$(dirname "${RUNNER_SCRIPT}")" && pwd)/$(basename "${RUNNER_SCRIPT}")"
PARAMETER_TABLE="$(cd "$(dirname "${PARAMETER_TABLE}")" && pwd)/$(basename "${PARAMETER_TABLE}")"
FIT_OBJECTS_DIR="$(cd "${FIT_OBJECTS_DIR}" && pwd)"
if [[ -f "${FLOW_DENSITY_PATH}" ]]; then
  FLOW_DENSITY_PATH="$(cd "$(dirname "${FLOW_DENSITY_PATH}")" && pwd)/$(basename "${FLOW_DENSITY_PATH}")"
fi

RUN_DIR_INVIVO="${OUT_ROOT}/${RUN_PREFIX_INVIVO}"
RUN_DIR_INVITRO="${OUT_ROOT}/${RUN_PREFIX_INVITRO}"
PREP_DIR="${OUT_ROOT}/${PREP_PREFIX}"
WRAPPER_SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

write_o2_min_config() {
  if ! truthy "${WRITE_CONFIG}"; then
    if [[ ! -f "${CONFIG_PATH}" ]]; then
      echo "WRITE_CONFIG=FALSE but CONFIG_PATH does not exist: ${CONFIG_PATH}" >&2
      exit 1
    fi
    return
  fi

  local cmd=(
    Rscript -e '
      args <- commandArgs(TRUE)
      if (!requireNamespace("yaml", quietly = TRUE)) stop("Package yaml is required")
      base_config <- normalizePath(args[[1]], mustWork = TRUE)
      out_config <- normalizePath(args[[2]], mustWork = FALSE)
      o2_min <- as.numeric(args[[3]])
      cfg <- yaml::read_yaml(base_config)
      if (is.null(cfg)) cfg <- list()
      base_dir <- dirname(base_config)
      resolve_path <- function(x) {
        if (is.null(x) || !length(x)) return(x)
        txt <- trimws(as.character(x[[1]]))
        if (!nzchar(txt)) return(x)
        if (grepl("^(/|~|[A-Za-z]:[/\\\\])", txt)) {
          return(normalizePath(path.expand(txt), mustWork = FALSE))
        }
        normalizePath(file.path(base_dir, txt), mustWork = FALSE)
      }
      path_keys <- c(
        "out_root", "data_dir", "seeds_file", "parameter_table", "parameters",
        "init_params_tsv", "invitro_parameter_table", "parameter_table_invitro",
        "fit_objects_dir", "flow_density_path",
        "joint_soft_coupling_parameters_table",
        "joint_soft_coupling_parameters_table_path",
        "joint_warmup_invivo_seed_dir", "joint_warmup_invitro_seed_dir"
      )
      for (key in path_keys) {
        if (!is.null(cfg[[key]])) cfg[[key]] <- resolve_path(cfg[[key]])
      }
      cfg$o2_min <- o2_min
      cfg$joint_soft_coupling_enable <- TRUE
      cfg$append_run_prefix_timestamp <- FALSE
      dir.create(dirname(out_config), recursive = TRUE, showWarnings = FALSE)
      yaml::write_yaml(cfg, out_config)
      cat("wrote_config=", out_config, "\n", sep = "")
    ' "${BASE_CONFIG}" "${CONFIG_PATH}" "${O2_MIN}"
  )
  if truthy "${DRY_RUN}"; then
    print_command "Generate o2_min config" "${cmd[@]}"
  else
    load_r_module
    "${cmd[@]}"
    if [[ ! -f "${CONFIG_PATH}" ]]; then
      echo "Config generation failed: ${CONFIG_PATH}" >&2
      exit 1
    fi
  fi
}

submit_postprocess() {
  local label="$1"
  local run_dir="$2"
  local dependency="$3"
  local job_name="o2g_extra_${label}"
  local exports="ALL,PROJECT_ROOT=${PROJECT_ROOT},RUN_DIR=${run_dir},EXTRA_RESULTS_SCRIPT=${EXTRA_RESULTS_SCRIPT},R_MODULE=${R_MODULE},FORCE_EXTRA_RESULTS=${FORCE_EXTRA_RESULTS}"
  local cmd=(
    sbatch
    --parsable
    "--job-name=${job_name}"
    "--qos=${POSTPROCESS_QOS}"
    "--time=${POSTPROCESS_TIME_LIMIT}"
    --cpus-per-task=1
    "--mem=${POSTPROCESS_MEM}"
    "--output=${LOG_ROOT}/${job_name}_%A.out"
    "--error=${LOG_ROOT}/${job_name}_%A.err"
    "--export=${exports}"
    "--dependency=afterok:${dependency}"
    "${POSTPROCESS_SCRIPT}"
  )
  if truthy "${DRY_RUN}"; then
    print_command "Submit ${label} extra_results" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_${label}_EXTRA_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted ${label} extra_results job: ${LAST_JOB_ID}"
  fi
}

submit_main_stage() {
  write_o2_min_config
  if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch not found; run this script on the HPC login node or use DRY_RUN=TRUE." >&2
    exit 1
  fi

  mkdir -p "${RUN_DIR_INVIVO}" "${RUN_DIR_INVITRO}" "${PREP_DIR}"

  local exports="ALL"
  exports+=",PROJECT_ROOT=${PROJECT_ROOT}"
  exports+=",RUNNER_SCRIPT=${RUNNER_SCRIPT}"
  exports+=",CONFIG_PATH=${CONFIG_PATH}"
  exports+=",OUT_ROOT=${OUT_ROOT}"
  exports+=",RUN_PREFIX_INVIVO=${RUN_PREFIX_INVIVO}"
  exports+=",RUN_PREFIX_INVITRO=${RUN_PREFIX_INVITRO}"
  exports+=",TOTAL_SEEDS_PER_MODE=${TOTAL_SEEDS_PER_MODE}"
  exports+=",ARRAY_TASKS=${ARRAY_TASKS}"
  exports+=",N_CORES=${N_CORES}"
  exports+=",AUTO_VIZ=${AUTO_VIZ}"
  exports+=",GLUCOSE=${GLUCOSE}"
  exports+=",R_MODULE=${R_MODULE}"
  exports+=",PARAMETER_TABLE=${PARAMETER_TABLE}"
  exports+=",FIT_OBJECTS_DIR=${FIT_OBJECTS_DIR}"
  exports+=",FLOW_DENSITY_PATH=${FLOW_DENSITY_PATH}"
  exports+=",ITERMAX=${ITERMAX}"
  exports+=",DE_RELTOL=${DE_RELTOL}"
  exports+=",DE_STEPTOL=${DE_STEPTOL}"
  exports+=",NP=${NP}"

  local array_cmd=(
    sbatch
    --parsable
    --job-name=o2g_iv_best
    "--qos=${QOS}"
    "--time=${TIME_LIMIT}"
    "--cpus-per-task=${N_CORES}"
    "--mem=${MEM}"
    "--array=1-${ARRAY_TASKS}"
    "--output=${LOG_ROOT}/o2g_iv_best_%A_%a.out"
    "--error=${LOG_ROOT}/o2g_iv_best_%A_%a.err"
    "--export=${exports}"
    "${SUB_SCRIPT}"
  )

  echo "Submitting in vivo/in vitro fits before best-seed joint"
  echo "  project_root: ${PROJECT_ROOT}"
  echo "  base_config: ${BASE_CONFIG}"
  echo "  active_config: ${CONFIG_PATH}"
  echo "  o2_min: ${O2_MIN}"
  echo "  invivo_run_dir: ${RUN_DIR_INVIVO}"
  echo "  invitro_run_dir: ${RUN_DIR_INVITRO}"
  echo "  joint_run_prefix_base: ${JOINT_RUN_PREFIX}"
  echo "  dry_run: ${DRY_RUN}"

  local mixed_job_id
  if truthy "${DRY_RUN}"; then
    print_command "Submit interleaved in vivo/in vitro array" "${array_cmd[@]}"
    mixed_job_id="DRYRUN_MIXED_JOB"
  else
    mixed_job_id="$("${array_cmd[@]}")"
    echo "Submitted interleaved in vivo/in vitro array job: ${mixed_job_id}"
  fi

  submit_postprocess "invivo" "${RUN_DIR_INVIVO}" "${mixed_job_id}"
  local invivo_extra_job_id="${LAST_JOB_ID}"
  submit_postprocess "invitro" "${RUN_DIR_INVITRO}" "${mixed_job_id}"
  local invitro_extra_job_id="${LAST_JOB_ID}"

  local prep_cmd=(
    sbatch
    --parsable
    --job-name=o2g_best_joint_prep
    "--qos=${PREP_QOS}"
    "--time=${PREP_TIME_LIMIT}"
    --cpus-per-task=1
    "--mem=${PREP_MEM}"
    "--output=${LOG_ROOT}/o2g_best_joint_prep_%A.out"
    "--error=${LOG_ROOT}/o2g_best_joint_prep_%A.err"
    "--dependency=afterok:${invivo_extra_job_id}:${invitro_extra_job_id}"
    --export=ALL
    "${WRAPPER_SCRIPT}"
    --internal_stage=select_and_submit_joint
    "--project_root=${PROJECT_ROOT}"
    "--base_config=${BASE_CONFIG}"
    "--config_path=${CONFIG_PATH}"
    "--write_config=FALSE"
    "--o2_min=${O2_MIN}"
    "--out_root=${OUT_ROOT}"
    "--log_root=${LOG_ROOT}"
    "--r_module=${R_MODULE}"
    "--force_extra_results=${FORCE_EXTRA_RESULTS}"
    "--run_prefix_invivo=${RUN_PREFIX_INVIVO}"
    "--run_prefix_invitro=${RUN_PREFIX_INVITRO}"
    "--prep_prefix=${PREP_PREFIX}"
    "--select_required_files=${SELECT_REQUIRED_FILES}"
    "--invivo_objective_columns=${INVIVO_OBJECTIVE_COLUMNS}"
    "--invitro_objective_columns=${INVITRO_OBJECTIVE_COLUMNS}"
    "--joint_run_prefix=${JOINT_RUN_PREFIX}"
    "--joint_job_name=${JOINT_JOB_NAME}"
    "--joint_total_seeds=${JOINT_TOTAL_SEEDS}"
    "--joint_array_tasks=${JOINT_ARRAY_TASKS}"
    "--joint_seeds_per_task=${JOINT_SEEDS_PER_TASK}"
    "--joint_qos=${JOINT_QOS}"
    "--joint_time_limit=${JOINT_TIME_LIMIT}"
    "--joint_n_cores=${JOINT_N_CORES}"
    "--joint_mem=${JOINT_MEM}"
    "--joint_soft_coupling_delta_params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}"
    "--auto_viz=${AUTO_VIZ}"
    "--glucose=${GLUCOSE}"
    "--itermax=${ITERMAX}"
    "--de_reltol=${DE_RELTOL}"
    "--de_steptol=${DE_STEPTOL}"
    "--NP=${NP}"
    "--parameter_table=${PARAMETER_TABLE}"
    "--fit_objects_dir=${FIT_OBJECTS_DIR}"
    "--flow_density_path=${FLOW_DENSITY_PATH}"
    "--dry_run=FALSE"
  )
  if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then
    prep_cmd+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}")
  fi
  if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then
    prep_cmd+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}")
  fi

  if truthy "${DRY_RUN}"; then
    print_command "Submit best-seed selection and joint submit prep" "${prep_cmd[@]}"
    echo "DRY_RUN complete."
  else
    local prep_job_id
    prep_job_id="$("${prep_cmd[@]}")"
    echo "Submitted best-seed selection and joint submit prep job: ${prep_job_id}"
    {
      printf "key\tvalue\n"
      printf "project_root\t%s\n" "${PROJECT_ROOT}"
      printf "base_config\t%s\n" "${BASE_CONFIG}"
      printf "active_config\t%s\n" "${CONFIG_PATH}"
      printf "o2_min\t%s\n" "${O2_MIN}"
      printf "out_root\t%s\n" "${OUT_ROOT}"
      printf "log_root\t%s\n" "${LOG_ROOT}"
      printf "run_dir_invivo\t%s\n" "${RUN_DIR_INVIVO}"
      printf "run_dir_invitro\t%s\n" "${RUN_DIR_INVITRO}"
      printf "mixed_job_id\t%s\n" "${mixed_job_id}"
      printf "invivo_extra_job_id\t%s\n" "${invivo_extra_job_id}"
      printf "invitro_extra_job_id\t%s\n" "${invitro_extra_job_id}"
      printf "prep_job_id\t%s\n" "${prep_job_id}"
      printf "joint_run_prefix_base\t%s\n" "${JOINT_RUN_PREFIX}"
    } > "${PREP_DIR}/submission_manifest.tsv"
    echo "Submission manifest: ${PREP_DIR}/submission_manifest.tsv"
  fi
}

select_best_seed() {
  local label="$1"
  local run_dir="$2"
  local objective_columns="$3"
  local log_path="${PREP_DIR}/select_best_${label}.log"
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

select_and_submit_joint_stage() {
  if truthy "${DRY_RUN}"; then
    echo "Internal prep stage does not support DRY_RUN=TRUE." >&2
    exit 2
  fi
  if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch not found in prep stage; cannot submit joint job." >&2
    exit 1
  fi
  if [[ ! -f "${CONFIG_PATH}" ]]; then
    echo "Missing active config in prep stage: ${CONFIG_PATH}" >&2
    exit 1
  fi
  if [[ ! -d "${RUN_DIR_INVIVO}" || ! -d "${RUN_DIR_INVITRO}" ]]; then
    echo "Missing single-fit result directories:" >&2
    echo "  invivo: ${RUN_DIR_INVIVO}" >&2
    echo "  invitro: ${RUN_DIR_INVITRO}" >&2
    exit 1
  fi
  mkdir -p "${PREP_DIR}"
  load_r_module

  select_best_seed "invivo" "${RUN_DIR_INVIVO}" "${INVIVO_OBJECTIVE_COLUMNS}"
  select_best_seed "invitro" "${RUN_DIR_INVITRO}" "${INVITRO_OBJECTIVE_COLUMNS}"

  local invivo_best_dir
  local invitro_best_dir
  invivo_best_dir="$(first_line "${RUN_DIR_INVIVO}/best_seed_from_summary.dir")"
  invitro_best_dir="$(first_line "${RUN_DIR_INVITRO}/best_seed_from_summary.dir")"
  if [[ ! -d "${invivo_best_dir}" ]]; then
    echo "Selected in vivo best seed dir does not exist: ${invivo_best_dir}" >&2
    exit 1
  fi
  if [[ ! -d "${invitro_best_dir}" ]]; then
    echo "Selected in vitro best seed dir does not exist: ${invitro_best_dir}" >&2
    exit 1
  fi
  invivo_best_dir="$(cd "${invivo_best_dir}" && pwd)"
  invitro_best_dir="$(cd "${invitro_best_dir}" && pwd)"

  local warmup_seed_label
  warmup_seed_label="invivo_$(sanitize_label "$(basename "${invivo_best_dir}")")__invitro_$(sanitize_label "$(basename "${invitro_best_dir}")")"
  local final_joint_run_prefix="${JOINT_RUN_PREFIX}"
  if [[ "${final_joint_run_prefix}" != *"${warmup_seed_label}"* ]]; then
    final_joint_run_prefix="${final_joint_run_prefix}__${warmup_seed_label}"
  fi

  {
    printf "key\tvalue\n"
    printf "project_root\t%s\n" "${PROJECT_ROOT}"
    printf "active_config\t%s\n" "${CONFIG_PATH}"
    printf "o2_min\t%s\n" "${O2_MIN}"
    printf "run_dir_invivo\t%s\n" "${RUN_DIR_INVIVO}"
    printf "run_dir_invitro\t%s\n" "${RUN_DIR_INVITRO}"
    printf "invivo_best_seed_dir\t%s\n" "${invivo_best_dir}"
    printf "invitro_best_seed_dir\t%s\n" "${invitro_best_dir}"
    printf "joint_warmup_seed_label\t%s\n" "${warmup_seed_label}"
    printf "joint_run_prefix\t%s\n" "${final_joint_run_prefix}"
    printf "select_required_files\t%s\n" "${SELECT_REQUIRED_FILES}"
    printf "invivo_objective_columns\t%s\n" "${INVIVO_OBJECTIVE_COLUMNS}"
    printf "invitro_objective_columns\t%s\n" "${INVITRO_OBJECTIVE_COLUMNS}"
  } > "${PREP_DIR}/best_seed_selection_manifest.tsv"

  local joint_cmd=(
    bash "${UNIFIED_SUBMIT}"
    --fitting_mode=joint
    --joint_fitting_mode=DIRECT
    "--project_root=${PROJECT_ROOT}"
    "--config_path=${CONFIG_PATH}"
    "--out_root=${OUT_ROOT}"
    "--log_root=${LOG_ROOT}"
    "--r_module=${R_MODULE}"
    "--force_extra_results=${FORCE_EXTRA_RESULTS}"
    "--joint_run_prefix=${final_joint_run_prefix}"
    "--joint_job_name=${JOINT_JOB_NAME}"
    "--joint_total_seeds=${JOINT_TOTAL_SEEDS}"
    "--joint_array_tasks=${JOINT_ARRAY_TASKS}"
    "--joint_seeds_per_task=${JOINT_SEEDS_PER_TASK}"
    "--joint_qos=${JOINT_QOS}"
    "--joint_time_limit=${JOINT_TIME_LIMIT}"
    "--joint_n_cores=${JOINT_N_CORES}"
    "--joint_mem=${JOINT_MEM}"
    "--invivo_best_seed_dir=${invivo_best_dir}"
    "--invitro_best_seed_dir=${invitro_best_dir}"
    "--joint_warmup_enable=TRUE"
    "--joint_warmup_seed_label=${warmup_seed_label}"
    "--joint_soft_coupling_delta_params=${JOINT_SOFT_COUPLING_DELTA_PARAMS}"
    "--parameter_table=${PARAMETER_TABLE}"
    "--fit_objects_dir=${FIT_OBJECTS_DIR}"
    "--flow_density_path=${FLOW_DENSITY_PATH}"
    "--auto_viz=${AUTO_VIZ}"
    "--glucose=${GLUCOSE}"
    "--itermax=${ITERMAX}"
    "--de_reltol=${DE_RELTOL}"
    "--de_steptol=${DE_STEPTOL}"
    "--NP=${NP}"
  )
  if [[ -n "${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}" ]]; then
    joint_cmd+=("--joint_soft_coupling_sigma_default=${JOINT_SOFT_COUPLING_SIGMA_DEFAULT}")
  fi
  if [[ -n "${JOINT_WARMUP_SIGMAN}" ]]; then
    joint_cmd+=("--joint_warmup_sigmaN=${JOINT_WARMUP_SIGMAN}")
  fi

  echo "Submitting warm-started joint fit"
  print_command "Submit joint" "${joint_cmd[@]}"
  "${joint_cmd[@]}" | tee "${PREP_DIR}/submit_joint_from_best_seed.log"
  echo "Best-seed joint prep complete."
  echo "  manifest: ${PREP_DIR}/best_seed_selection_manifest.tsv"
  echo "  joint_run_prefix: ${final_joint_run_prefix}"
}

case "${INTERNAL_STAGE}" in
  submit)
    submit_main_stage
    ;;
  select_and_submit_joint)
    select_and_submit_joint_stage
    ;;
  *)
    echo "Unknown --internal_stage: ${INTERNAL_STAGE}" >&2
    exit 2
    ;;
esac
