#!/bin/bash
# Prepare anchors from finished single fits, then submit the joint O2G array.

#SBATCH --job-name=o2g_joint_prep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash postprocess_and_submit_joint_buffering.sh \
    --invivo_run_dir=/path/to/invivo_run \
    --invitro_run_dir=/path/to/invitro_run \
    --joint_run_prefix=fit_joint_name [options]

Options:
  --project_root=/share/lab_crd/lab_crd/taoli/Project/Rescomposite
  --config_path=/path/to/O2G_supply_demand.yaml
  --out_root=/path/to/oxygen/results
  --anchor_mode=auto|manual
  --manual_o2_grid=0,0.1,0.2,0.5,1,2
  --manual_n_grid=44,66,88
  --top_n=3
  --joint_total_seeds=500
  --joint_array_tasks=500
  --joint_seeds_per_task=1
  --joint_qos=xxlarge
  --joint_time_limit=12:00:00
  --joint_n_cores=22
  --joint_mem=16G
  --postprocess_qos=small
  --postprocess_time_limit=4:00:00
  --postprocess_mem=8G
  --r_module=R/4.4
  --force_extra_results=TRUE|FALSE
  --dry_run=TRUE|FALSE
  --help

Behavior:
  1. Ensures in vivo and in vitro extra_results exist, running extra_results.R
     after loading R/4.4 when needed.
  2. Generates joint_anchor.tsv and joint_anchor_config.yaml from the selected
     single-fit top seeds.
  3. Submits the joint seed array using the generated config.
  4. Submits a dependent extra_results postprocess job for the joint run.
EOF
}

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h)
        usage
        exit 0
        ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --config_path=*) CONFIG_PATH="${arg#*=}" ;;
      --out_root=*) OUT_ROOT="${arg#*=}" ;;
      --invivo_run_dir=*) INVIVO_RUN_DIR="${arg#*=}" ;;
      --invitro_run_dir=*) INVITRO_RUN_DIR="${arg#*=}" ;;
      --joint_run_prefix=*) JOINT_RUN_PREFIX="${arg#*=}" ;;
      --anchor_mode=*) ANCHOR_MODE="${arg#*=}" ;;
      --manual_o2_grid=*) MANUAL_O2_GRID="${arg#*=}" ;;
      --manual_n_grid=*) MANUAL_N_GRID="${arg#*=}" ;;
      --top_n=*) TOP_N="${arg#*=}" ;;
      --threshold_n=*) THRESHOLD_N="${arg#*=}" ;;
      --joint_total_seeds=*) JOINT_TOTAL_SEEDS="${arg#*=}" ;;
      --joint_array_tasks=*) JOINT_ARRAY_TASKS="${arg#*=}" ;;
      --joint_seeds_per_task=*) JOINT_SEEDS_PER_TASK="${arg#*=}" ;;
      --joint_qos=*) JOINT_QOS="${arg#*=}" ;;
      --joint_time_limit=*) JOINT_TIME_LIMIT="${arg#*=}" ;;
      --joint_n_cores=*) JOINT_N_CORES="${arg#*=}" ;;
      --joint_mem=*) JOINT_MEM="${arg#*=}" ;;
      --joint_job_name=*) JOINT_JOB_NAME="${arg#*=}" ;;
      --postprocess_qos=*) POSTPROCESS_QOS="${arg#*=}" ;;
      --postprocess_time_limit=*) POSTPROCESS_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_mem=*) POSTPROCESS_MEM="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --force_extra_results=*) FORCE_EXTRA_RESULTS="${arg#*=}" ;;
      --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
      --glucose=*) GLUCOSE="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      *)
        echo "Unknown argument: ${arg}" >&2
        usage >&2
        exit 2
        ;;
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

require_positive_int() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]] || (( value <= 0 )); then
    echo "${name} must be a positive integer, got: ${value}" >&2
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

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/Rescomposite"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_OUT_ROOT=""
DEFAULT_JOINT_RUN_PREFIX="fit_joint_O2G_buffering_500seed"
DEFAULT_ANCHOR_MODE="auto"
DEFAULT_TOP_N="3"
DEFAULT_THRESHOLD_N="44"
DEFAULT_JOINT_TOTAL_SEEDS="500"
DEFAULT_JOINT_ARRAY_TASKS="500"
DEFAULT_JOINT_SEEDS_PER_TASK="1"
DEFAULT_JOINT_N_CORES="22"
DEFAULT_JOINT_MEM="16G"
DEFAULT_JOINT_QOS="xxlarge"
DEFAULT_JOINT_TIME_LIMIT="12:00:00"
DEFAULT_JOINT_JOB_NAME="o2g_joint_B"
DEFAULT_POSTPROCESS_QOS="small"
DEFAULT_POSTPROCESS_TIME_LIMIT="4:00:00"
DEFAULT_POSTPROCESS_MEM="8G"
DEFAULT_AUTO_VIZ="TRUE"
DEFAULT_GLUCOSE="TRUE"
DEFAULT_FORCE_EXTRA_RESULTS="FALSE"
DEFAULT_DRY_RUN="FALSE"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
OUT_ROOT="${OUT_ROOT:-${DEFAULT_OUT_ROOT}}"
CONFIG_PATH="${CONFIG_PATH:-}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-${DEFAULT_JOINT_RUN_PREFIX}}"
ANCHOR_MODE="${ANCHOR_MODE:-${DEFAULT_ANCHOR_MODE}}"
MANUAL_O2_GRID="${MANUAL_O2_GRID:-}"
MANUAL_N_GRID="${MANUAL_N_GRID:-}"
TOP_N="${TOP_N:-${DEFAULT_TOP_N}}"
THRESHOLD_N="${THRESHOLD_N:-${DEFAULT_THRESHOLD_N}}"
JOINT_TOTAL_SEEDS="${JOINT_TOTAL_SEEDS:-${DEFAULT_JOINT_TOTAL_SEEDS}}"
JOINT_ARRAY_TASKS="${JOINT_ARRAY_TASKS:-${DEFAULT_JOINT_ARRAY_TASKS}}"
JOINT_SEEDS_PER_TASK="${JOINT_SEEDS_PER_TASK:-${DEFAULT_JOINT_SEEDS_PER_TASK}}"
JOINT_N_CORES="${JOINT_N_CORES:-${DEFAULT_JOINT_N_CORES}}"
JOINT_MEM="${JOINT_MEM:-${DEFAULT_JOINT_MEM}}"
JOINT_QOS="${JOINT_QOS:-${DEFAULT_JOINT_QOS}}"
JOINT_TIME_LIMIT="${JOINT_TIME_LIMIT:-${DEFAULT_JOINT_TIME_LIMIT}}"
JOINT_JOB_NAME="${JOINT_JOB_NAME:-${DEFAULT_JOINT_JOB_NAME}}"
POSTPROCESS_QOS="${POSTPROCESS_QOS:-${DEFAULT_POSTPROCESS_QOS}}"
POSTPROCESS_TIME_LIMIT="${POSTPROCESS_TIME_LIMIT:-${DEFAULT_POSTPROCESS_TIME_LIMIT}}"
POSTPROCESS_MEM="${POSTPROCESS_MEM:-${DEFAULT_POSTPROCESS_MEM}}"
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
GLUCOSE="${GLUCOSE:-${DEFAULT_GLUCOSE}}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-${DEFAULT_FORCE_EXTRA_RESULTS}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"

parse_args "$@"

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
if [[ -z "${CONFIG_PATH}" ]]; then
  CONFIG_PATH="${PROJECT_ROOT}/oxygen/config/O2G_supply_demand.yaml"
fi
if [[ -z "${OUT_ROOT}" ]]; then
  OUT_ROOT="${PROJECT_ROOT}/oxygen/results"
fi
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
mkdir -p "${OUT_ROOT}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"

if [[ -z "${INVIVO_RUN_DIR}" || -z "${INVITRO_RUN_DIR}" ]]; then
  echo "Both --invivo_run_dir and --invitro_run_dir are required." >&2
  exit 2
fi
if [[ ! -d "${INVIVO_RUN_DIR}" ]]; then
  echo "Missing in vivo run directory: ${INVIVO_RUN_DIR}" >&2
  exit 1
fi
if [[ ! -d "${INVITRO_RUN_DIR}" ]]; then
  echo "Missing in vitro run directory: ${INVITRO_RUN_DIR}" >&2
  exit 1
fi
INVIVO_RUN_DIR="$(cd "${INVIVO_RUN_DIR}" && pwd)"
INVITRO_RUN_DIR="$(cd "${INVITRO_RUN_DIR}" && pwd)"
if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "Missing config: ${CONFIG_PATH}" >&2
  exit 1
fi

ANCHOR_MODE="$(echo "${ANCHOR_MODE}" | tr '[:upper:]' '[:lower:]')"
case "${ANCHOR_MODE}" in
  auto|manual) ;;
  *)
    echo "ANCHOR_MODE must be auto or manual, got: ${ANCHOR_MODE}" >&2
    exit 2
    ;;
esac
require_positive_int TOP_N "${TOP_N}"
require_positive_int JOINT_TOTAL_SEEDS "${JOINT_TOTAL_SEEDS}"
require_positive_int JOINT_ARRAY_TASKS "${JOINT_ARRAY_TASKS}"
require_positive_int JOINT_SEEDS_PER_TASK "${JOINT_SEEDS_PER_TASK}"
require_positive_int JOINT_N_CORES "${JOINT_N_CORES}"
if (( JOINT_ARRAY_TASKS * JOINT_SEEDS_PER_TASK != JOINT_TOTAL_SEEDS )); then
  echo "JOINT_ARRAY_TASKS * JOINT_SEEDS_PER_TASK must equal JOINT_TOTAL_SEEDS." >&2
  exit 2
fi

HPC_DIR="${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc"
POSTPROCESS_SCRIPT="${POSTPROCESS_SCRIPT:-${HPC_DIR}/postprocess_extra_results.sh}"
DERIVE_ANCHOR_SCRIPT="${DERIVE_ANCHOR_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/derive_joint_anchor_from_single_fits.R}"
EXTRA_RESULTS_SCRIPT="${EXTRA_RESULTS_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R}"
JOINT_SUB_SCRIPT="${JOINT_SUB_SCRIPT:-${HPC_DIR}/submit_fit_seed_array_joint_buffering.sub}"
JOINT_RUNNER_SCRIPT="${JOINT_RUNNER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/runner/run_fit_joint_model_O2G_supply_demand_MAP.sh}"

for path in "${POSTPROCESS_SCRIPT}" "${DERIVE_ANCHOR_SCRIPT}" "${EXTRA_RESULTS_SCRIPT}" "${JOINT_SUB_SCRIPT}" "${JOINT_RUNNER_SCRIPT}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required script: ${path}" >&2
    exit 1
  fi
done

JOINT_RUN_DIR="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
mkdir -p "${JOINT_RUN_DIR}"
ANCHOR_TSV="${JOINT_RUN_DIR}/joint_anchor.tsv"
JOINT_CONFIG_PATH="${JOINT_RUN_DIR}/joint_anchor_config.yaml"
MANIFEST_PATH="${JOINT_RUN_DIR}/joint_submission_manifest.tsv"

postprocess_invivo_cmd=(bash "${POSTPROCESS_SCRIPT}" "--project_root=${PROJECT_ROOT}" "--run_dir=${INVIVO_RUN_DIR}" "--extra_results_script=${EXTRA_RESULTS_SCRIPT}" "--r_module=${R_MODULE}" "--force_extra_results=${FORCE_EXTRA_RESULTS}")
postprocess_invitro_cmd=(bash "${POSTPROCESS_SCRIPT}" "--project_root=${PROJECT_ROOT}" "--run_dir=${INVITRO_RUN_DIR}" "--extra_results_script=${EXTRA_RESULTS_SCRIPT}" "--r_module=${R_MODULE}" "--force_extra_results=${FORCE_EXTRA_RESULTS}")

anchor_cmd=(
  Rscript "${DERIVE_ANCHOR_SCRIPT}"
  "--invivo_run_dir=${INVIVO_RUN_DIR}"
  "--invitro_run_dir=${INVITRO_RUN_DIR}"
  "--base_config=${CONFIG_PATH}"
  "--out_dir=${JOINT_RUN_DIR}"
  "--out_config_yaml=${JOINT_CONFIG_PATH}"
  "--out_anchor_tsv=${ANCHOR_TSV}"
  "--run_prefix=${JOINT_RUN_PREFIX}"
  "--out_root=${OUT_ROOT}"
  "--anchor_mode=${ANCHOR_MODE}"
  "--top_n=${TOP_N}"
  "--threshold_N=${THRESHOLD_N}"
  "--extra_results_script=${EXTRA_RESULTS_SCRIPT}"
  "--run_extra_results=TRUE"
  "--force_extra_results=${FORCE_EXTRA_RESULTS}"
)
if [[ "${ANCHOR_MODE}" == "manual" ]]; then
  anchor_cmd+=("--manual_o2_grid=${MANUAL_O2_GRID}" "--manual_n_grid=${MANUAL_N_GRID}")
fi

if truthy "${DRY_RUN}"; then
  echo "DRY_RUN=TRUE; not running postprocess, anchor generation, or sbatch."
  print_command "Postprocess in vivo" "${postprocess_invivo_cmd[@]}"
  print_command "Postprocess in vitro" "${postprocess_invitro_cmd[@]}"
  print_command "Derive anchor" "${anchor_cmd[@]}"
else
  "${postprocess_invivo_cmd[@]}"
  "${postprocess_invitro_cmd[@]}"
  load_r_module
  cd "${PROJECT_ROOT}"
  "${anchor_cmd[@]}"
fi

joint_export="ALL"
joint_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
joint_export+=",RUNNER_SCRIPT=${JOINT_RUNNER_SCRIPT}"
joint_export+=",CONFIG_PATH=${JOINT_CONFIG_PATH}"
joint_export+=",OUT_ROOT=${OUT_ROOT}"
joint_export+=",RUN_PREFIX=${JOINT_RUN_PREFIX}"
joint_export+=",TOTAL_SEEDS=${JOINT_TOTAL_SEEDS}"
joint_export+=",ARRAY_TASKS=${JOINT_ARRAY_TASKS}"
joint_export+=",SEEDS_PER_TASK=${JOINT_SEEDS_PER_TASK}"
joint_export+=",N_CORES=${JOINT_N_CORES}"
joint_export+=",AUTO_VIZ=${AUTO_VIZ}"
joint_export+=",GLUCOSE=${GLUCOSE}"
joint_export+=",R_MODULE=${R_MODULE}"
joint_export+=",ANCHOR_MODE=${ANCHOR_MODE}"
joint_export+=",ANCHOR_TSV=${ANCHOR_TSV}"
joint_export+=",INVIVO_RUN_DIR=${INVIVO_RUN_DIR}"
joint_export+=",INVITRO_RUN_DIR=${INVITRO_RUN_DIR}"
joint_export+=",JOINT_FITTING_MODE=${JOINT_FITTING_MODE:-SINGLE}"

joint_cmd=(
  sbatch
  --parsable
  "--job-name=${JOINT_JOB_NAME}"
  "--qos=${JOINT_QOS}"
  "--time=${JOINT_TIME_LIMIT}"
  "--cpus-per-task=${JOINT_N_CORES}"
  "--mem=${JOINT_MEM}"
  "--array=1-${JOINT_ARRAY_TASKS}"
  "--output=${OUT_ROOT}/o2g_joint_fit_%A_%a.out"
  "--error=${OUT_ROOT}/o2g_joint_fit_%A_%a.err"
  "--export=${joint_export}"
  "${JOINT_SUB_SCRIPT}"
)

if truthy "${DRY_RUN}"; then
  print_command "Submit joint" "${joint_cmd[@]}"
  JOINT_JOB_ID="DRYRUN_JOINT_JOB"
else
  JOINT_JOB_ID="$("${joint_cmd[@]}")"
  echo "Submitted joint array job: ${JOINT_JOB_ID}"
fi

postprocess_export="ALL"
postprocess_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
postprocess_export+=",RUN_DIR=${JOINT_RUN_DIR}"
postprocess_export+=",EXTRA_RESULTS_SCRIPT=${EXTRA_RESULTS_SCRIPT}"
postprocess_export+=",R_MODULE=${R_MODULE}"
postprocess_export+=",FORCE_EXTRA_RESULTS=${FORCE_EXTRA_RESULTS}"

joint_postprocess_cmd=(
  sbatch
  --parsable
  "--dependency=afterok:${JOINT_JOB_ID}"
  "--job-name=${JOINT_JOB_NAME}_extra"
  "--qos=${POSTPROCESS_QOS}"
  "--time=${POSTPROCESS_TIME_LIMIT}"
  --cpus-per-task=1
  "--mem=${POSTPROCESS_MEM}"
  "--output=${OUT_ROOT}/o2g_joint_extra_%A.out"
  "--error=${OUT_ROOT}/o2g_joint_extra_%A.err"
  "--export=${postprocess_export}"
  "${POSTPROCESS_SCRIPT}"
)

if truthy "${DRY_RUN}"; then
  print_command "Submit joint extra_results" "${joint_postprocess_cmd[@]}"
  JOINT_EXTRA_JOB_ID="DRYRUN_JOINT_EXTRA_JOB"
else
  JOINT_EXTRA_JOB_ID="$("${joint_postprocess_cmd[@]}")"
  echo "Submitted joint extra_results postprocess job: ${JOINT_EXTRA_JOB_ID}"
fi

{
  printf "key\tvalue\n"
  printf "project_root\t%s\n" "${PROJECT_ROOT}"
  printf "config_path_base\t%s\n" "${CONFIG_PATH}"
  printf "joint_config_path\t%s\n" "${JOINT_CONFIG_PATH}"
  printf "out_root\t%s\n" "${OUT_ROOT}"
  printf "joint_run_prefix\t%s\n" "${JOINT_RUN_PREFIX}"
  printf "joint_run_dir\t%s\n" "${JOINT_RUN_DIR}"
  printf "invivo_run_dir\t%s\n" "${INVIVO_RUN_DIR}"
  printf "invitro_run_dir\t%s\n" "${INVITRO_RUN_DIR}"
  printf "anchor_mode\t%s\n" "${ANCHOR_MODE}"
  printf "anchor_tsv\t%s\n" "${ANCHOR_TSV}"
  printf "top_n\t%s\n" "${TOP_N}"
  printf "joint_job_id\t%s\n" "${JOINT_JOB_ID}"
  printf "joint_extra_results_job_id\t%s\n" "${JOINT_EXTRA_JOB_ID}"
  printf "r_module\t%s\n" "${R_MODULE}"
} > "${MANIFEST_PATH}"

echo "Joint preparation complete."
echo "  joint_run_dir: ${JOINT_RUN_DIR}"
echo "  anchor_tsv: ${ANCHOR_TSV}"
echo "  joint_config: ${JOINT_CONFIG_PATH}"
echo "  manifest: ${MANIFEST_PATH}"
echo "  joint resources: qos=${JOINT_QOS}, time=${JOINT_TIME_LIMIT}, mem=${JOINT_MEM}, cpus=${JOINT_N_CORES}"
echo "  postprocess resources: qos=${POSTPROCESS_QOS}, time=${POSTPROCESS_TIME_LIMIT}, mem=${POSTPROCESS_MEM}"
