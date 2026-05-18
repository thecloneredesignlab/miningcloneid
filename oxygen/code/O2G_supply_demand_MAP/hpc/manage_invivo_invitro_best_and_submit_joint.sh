#!/bin/bash
# Summarize completed in vivo and in vitro seed arrays, select the best seed
# from each modality, validate both warm-start directories, and submit joint fit.

#SBATCH --job-name=o2g_best_joint
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --qos=xxlarge
#SBATCH --time=12:00:00
#SBATCH --output=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/o2g_best_joint_%j.out
#SBATCH --error=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/o2g_best_joint_%j.err

set -euo pipefail

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid"
DEFAULT_INVIVO_RUN_PREFIX="fit_invivo_O2G_buffering_500seed"
DEFAULT_INVITRO_RUN_PREFIX="fit_invitro_O2G_buffering_500seed"
DEFAULT_JOINT_RUN_PREFIX="fit_joint_O2G_buffering_warmstart_500seed"
DEFAULT_TOTAL_SEEDS="500"
DEFAULT_ARRAY_TASKS="500"
DEFAULT_SEEDS_PER_TASK="1"
DEFAULT_N_CORES="22"
DEFAULT_AUTO_VIZ="TRUE"
DEFAULT_GLUCOSE="TRUE"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_JOINT_QOS="xxlarge"
DEFAULT_JOINT_TIME_LIMIT="12:00:00"
DEFAULT_SELECTION_HORIZON="1000"
DEFAULT_SELECTION_THRESHOLD_N="44"
DEFAULT_SELECTION_COHORT="2N"
DEFAULT_SELECTION_DOSE="ALL"
DEFAULT_NEAR_THRESH="0.05"
DEFAULT_INVITRO_OBJECTIVE_COLUMNS="objective,objective_total"
DEFAULT_DRY_RUN="FALSE"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
CONFIG_PATH="${CONFIG_PATH:-${PROJECT_ROOT}/oxygen/config/O2G_supply_demand.yaml}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/oxygen/results}"
INVIVO_RUN_PREFIX="${INVIVO_RUN_PREFIX:-${DEFAULT_INVIVO_RUN_PREFIX}}"
INVITRO_RUN_PREFIX="${INVITRO_RUN_PREFIX:-${DEFAULT_INVITRO_RUN_PREFIX}}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-${DEFAULT_JOINT_RUN_PREFIX}}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-${OUT_ROOT}/${INVIVO_RUN_PREFIX}}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-${OUT_ROOT}/${INVITRO_RUN_PREFIX}}"
EXTRA_RESULTS_SCRIPT="${EXTRA_RESULTS_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R}"
INVIVO_SELECTOR_SCRIPT="${INVIVO_SELECTOR_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/select_invivo_best_long_ploidy_seed.R}"
INVITRO_SELECTOR_SCRIPT="${INVITRO_SELECTOR_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/select_best_seed_from_summary.R}"
BUILD_SCRIPT="${BUILD_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/build_joint_init_candidates.R}"
JOINT_SUB_SCRIPT="${JOINT_SUB_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/submit_fit_seed_array_joint_buffering_warmstart.sub}"
TOTAL_SEEDS="${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}"
ARRAY_TASKS="${ARRAY_TASKS:-${DEFAULT_ARRAY_TASKS}}"
SEEDS_PER_TASK="${SEEDS_PER_TASK:-${DEFAULT_SEEDS_PER_TASK}}"
N_CORES="${N_CORES:-${DEFAULT_N_CORES}}"
JOINT_N_CORES="${JOINT_N_CORES:-${N_CORES}}"
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
GLUCOSE="${GLUCOSE:-${DEFAULT_GLUCOSE}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
JOINT_QOS="${JOINT_QOS:-${DEFAULT_JOINT_QOS}}"
JOINT_TIME_LIMIT="${JOINT_TIME_LIMIT:-${DEFAULT_JOINT_TIME_LIMIT}}"
SELECTION_HORIZON="${SELECTION_HORIZON:-${DEFAULT_SELECTION_HORIZON}}"
SELECTION_THRESHOLD_N="${SELECTION_THRESHOLD_N:-${DEFAULT_SELECTION_THRESHOLD_N}}"
SELECTION_COHORT="${SELECTION_COHORT:-${DEFAULT_SELECTION_COHORT}}"
SELECTION_DOSE="${SELECTION_DOSE:-${DEFAULT_SELECTION_DOSE}}"
NEAR_THRESH="${NEAR_THRESH:-${DEFAULT_NEAR_THRESH}}"
INVITRO_OBJECTIVE_COLUMNS="${INVITRO_OBJECTIVE_COLUMNS:-${DEFAULT_INVITRO_OBJECTIVE_COLUMNS}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"

if command -v ml >/dev/null 2>&1; then
  ml "${R_MODULE}"
elif command -v module >/dev/null 2>&1; then
  module load "${R_MODULE}"
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found after module loading."
  exit 1
fi
if ! command -v sbatch >/dev/null 2>&1 && [[ "${DRY_RUN}" != "TRUE" && "${DRY_RUN}" != "true" && "${DRY_RUN}" != "1" ]]; then
  echo "sbatch not found; cannot submit joint fitting."
  exit 1
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
mkdir -p "${OUT_ROOT}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"
EXTRA_RESULTS_SCRIPT="$(cd "$(dirname "${EXTRA_RESULTS_SCRIPT}")" && pwd)/$(basename "${EXTRA_RESULTS_SCRIPT}")"
INVIVO_SELECTOR_SCRIPT="$(cd "$(dirname "${INVIVO_SELECTOR_SCRIPT}")" && pwd)/$(basename "${INVIVO_SELECTOR_SCRIPT}")"
INVITRO_SELECTOR_SCRIPT="$(cd "$(dirname "${INVITRO_SELECTOR_SCRIPT}")" && pwd)/$(basename "${INVITRO_SELECTOR_SCRIPT}")"
BUILD_SCRIPT="$(cd "$(dirname "${BUILD_SCRIPT}")" && pwd)/$(basename "${BUILD_SCRIPT}")"
JOINT_SUB_SCRIPT="$(cd "$(dirname "${JOINT_SUB_SCRIPT}")" && pwd)/$(basename "${JOINT_SUB_SCRIPT}")"

if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "Missing config file: ${CONFIG_PATH}"
  exit 1
fi
if [[ ! -d "${INVIVO_RUN_DIR}" ]]; then
  echo "Missing in vivo run directory: ${INVIVO_RUN_DIR}"
  exit 1
fi
if [[ ! -d "${INVITRO_RUN_DIR}" ]]; then
  echo "Missing in vitro run directory: ${INVITRO_RUN_DIR}"
  exit 1
fi
INVIVO_RUN_DIR="$(cd "${INVIVO_RUN_DIR}" && pwd)"
INVITRO_RUN_DIR="$(cd "${INVITRO_RUN_DIR}" && pwd)"
for path in "${EXTRA_RESULTS_SCRIPT}" "${INVIVO_SELECTOR_SCRIPT}" "${INVITRO_SELECTOR_SCRIPT}" "${BUILD_SCRIPT}" "${JOINT_SUB_SCRIPT}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required script: ${path}"
    exit 1
  fi
done

for var_name in TOTAL_SEEDS ARRAY_TASKS SEEDS_PER_TASK JOINT_N_CORES; do
  var_val="${!var_name}"
  if ! [[ "${var_val}" =~ ^[0-9]+$ ]]; then
    echo "${var_name} must be a positive integer, got: ${var_val}"
    exit 1
  fi
  if (( var_val <= 0 )); then
    echo "${var_name} must be > 0, got: ${var_val}"
    exit 1
  fi
done
if (( ARRAY_TASKS * SEEDS_PER_TASK != TOTAL_SEEDS )); then
  echo "Mismatch: ARRAY_TASKS * SEEDS_PER_TASK must equal TOTAL_SEEDS."
  echo "Got ARRAY_TASKS=${ARRAY_TASKS}, SEEDS_PER_TASK=${SEEDS_PER_TASK}, TOTAL_SEEDS=${TOTAL_SEEDS}."
  exit 1
fi

INVIVO_EXTRA_RESULTS_DIR="${INVIVO_EXTRA_RESULTS_DIR:-${INVIVO_RUN_DIR}/extra_results}"
INVITRO_EXTRA_RESULTS_DIR="${INVITRO_EXTRA_RESULTS_DIR:-${INVITRO_RUN_DIR}/extra_results}"
INVIVO_SEED_SUMMARY="${INVIVO_SEED_SUMMARY:-${INVIVO_EXTRA_RESULTS_DIR}/seed_summary.tsv}"
INVITRO_SEED_SUMMARY="${INVITRO_SEED_SUMMARY:-${INVITRO_EXTRA_RESULTS_DIR}/seed_summary.tsv}"
INVIVO_SELECTION_TSV="${INVIVO_SELECTION_TSV:-${INVIVO_RUN_DIR}/best_long_ploidy_gt2_seed.tsv}"
INVIVO_BEST_DIR_FILE="${INVIVO_BEST_DIR_FILE:-${INVIVO_RUN_DIR}/best_long_ploidy_gt2_seed.dir}"
INVITRO_SELECTION_TSV="${INVITRO_SELECTION_TSV:-${INVITRO_RUN_DIR}/best_seed_from_summary.tsv}"
INVITRO_BEST_DIR_FILE="${INVITRO_BEST_DIR_FILE:-${INVITRO_RUN_DIR}/best_seed_from_summary.dir}"
JOINT_RUN_DIR="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
JOINT_INIT_CANDIDATES_TSV="${JOINT_INIT_CANDIDATES_TSV:-${JOINT_RUN_DIR}/joint_init_candidates.tsv}"

read_best_dir_file() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    echo "Missing best-dir file: ${path}" >&2
    return 1
  fi
  tr -d '[:space:]' < "${path}"
}

validate_best_dir() {
  local label="$1"
  local seed_dir="$2"
  if [[ -z "${seed_dir}" || ! -d "${seed_dir}" ]]; then
    echo "${label} best seed directory is missing or invalid: ${seed_dir}"
    exit 1
  fi
  if [[ ! -f "${seed_dir}/best_params_transformed.tsv" && ! -f "${seed_dir}/fit_parameter_stages.tsv" ]]; then
    echo "${label} best seed directory lacks transformed parameters: ${seed_dir}"
    echo "Expected best_params_transformed.tsv or fit_parameter_stages.tsv"
    exit 1
  fi
}

run_extra_results() {
  local label="$1"
  local run_dir="$2"
  local out_dir="$3"
  local summary_path="$4"
  echo "Running ${label} seed summary"
  echo "  run_dir: ${run_dir}"
  echo "  out_dir: ${out_dir}"
  Rscript "${EXTRA_RESULTS_SCRIPT}" \
    --run_dir="${run_dir}" \
    --out_dir="${out_dir}" \
    --near_thresh="${NEAR_THRESH}"
  if [[ ! -f "${summary_path}" ]]; then
    echo "${label} summary did not produce seed_summary.tsv: ${summary_path}"
    exit 1
  fi
}

echo "Starting summary/select manager"
echo "  project_root: ${PROJECT_ROOT}"
echo "  invivo_run_dir: ${INVIVO_RUN_DIR}"
echo "  invitro_run_dir: ${INVITRO_RUN_DIR}"
echo "  joint_run_dir: ${JOINT_RUN_DIR}"

run_extra_results "in vivo" "${INVIVO_RUN_DIR}" "${INVIVO_EXTRA_RESULTS_DIR}" "${INVIVO_SEED_SUMMARY}"

echo "Selecting best in vivo seed"
Rscript "${INVIVO_SELECTOR_SCRIPT}" \
  --invivo_dir="${INVIVO_RUN_DIR}" \
  --out_tsv="${INVIVO_SELECTION_TSV}" \
  --best_dir_file="${INVIVO_BEST_DIR_FILE}" \
  --horizon="${SELECTION_HORIZON}" \
  --threshold_N="${SELECTION_THRESHOLD_N}" \
  --cohort="${SELECTION_COHORT}" \
  --dose="${SELECTION_DOSE}"

INVIVO_BEST_DIR="$(read_best_dir_file "${INVIVO_BEST_DIR_FILE}")"
validate_best_dir "in vivo" "${INVIVO_BEST_DIR}"

run_extra_results "in vitro" "${INVITRO_RUN_DIR}" "${INVITRO_EXTRA_RESULTS_DIR}" "${INVITRO_SEED_SUMMARY}"

echo "Selecting best in vitro seed"
Rscript "${INVITRO_SELECTOR_SCRIPT}" \
  --run_dir="${INVITRO_RUN_DIR}" \
  --summary_tsv="${INVITRO_SEED_SUMMARY}" \
  --out_tsv="${INVITRO_SELECTION_TSV}" \
  --best_dir_file="${INVITRO_BEST_DIR_FILE}" \
  --objective_columns="${INVITRO_OBJECTIVE_COLUMNS}" \
  --required_files="best_params_transformed.tsv"

INVITRO_BEST_DIR="$(read_best_dir_file "${INVITRO_BEST_DIR_FILE}")"
validate_best_dir "in vitro" "${INVITRO_BEST_DIR}"

mkdir -p "${JOINT_RUN_DIR}"

echo "Building joint warm-start candidates"
echo "  invivo_best_dir: ${INVIVO_BEST_DIR}"
echo "  invitro_best_dir: ${INVITRO_BEST_DIR}"
echo "  joint_init_candidates_tsv: ${JOINT_INIT_CANDIDATES_TSV}"

Rscript "${BUILD_SCRIPT}" \
  --config="${CONFIG_PATH}" \
  --out_tsv="${JOINT_INIT_CANDIDATES_TSV}" \
  --glucose="${GLUCOSE}" \
  --invivo_best_dir="${INVIVO_BEST_DIR}" \
  --invitro_best_dir="${INVITRO_BEST_DIR}"

if [[ ! -f "${JOINT_INIT_CANDIDATES_TSV}" ]]; then
  echo "Joint warm-start candidate table was not produced: ${JOINT_INIT_CANDIDATES_TSV}"
  exit 1
fi

joint_export="ALL"
joint_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
joint_export+=",CONFIG_PATH=${CONFIG_PATH}"
joint_export+=",OUT_ROOT=${OUT_ROOT}"
joint_export+=",RUN_PREFIX=${JOINT_RUN_PREFIX}"
joint_export+=",TOTAL_SEEDS=${TOTAL_SEEDS}"
joint_export+=",ARRAY_TASKS=${ARRAY_TASKS}"
joint_export+=",SEEDS_PER_TASK=${SEEDS_PER_TASK}"
joint_export+=",N_CORES=${JOINT_N_CORES}"
joint_export+=",AUTO_VIZ=${AUTO_VIZ}"
joint_export+=",GLUCOSE=${GLUCOSE}"
joint_export+=",JOINT_INIT_CANDIDATES_TSV=${JOINT_INIT_CANDIDATES_TSV}"
joint_export+=",R_MODULE=${R_MODULE}"

joint_cmd=(
  sbatch
  --parsable
  --job-name=o2gj_WS
  --qos="${JOINT_QOS}"
  --time="${JOINT_TIME_LIMIT}"
  --cpus-per-task="${JOINT_N_CORES}"
  --array="1-${ARRAY_TASKS}"
  --output="${OUT_ROOT}/o2g_joint_ws_%A_%a.out"
  --error="${OUT_ROOT}/o2g_joint_ws_%A_%a.err"
  --export="${joint_export}"
  "${JOINT_SUB_SCRIPT}"
)

echo "Submitting warm-start joint array"
echo "  qos: ${JOINT_QOS}"
echo "  time_limit: ${JOINT_TIME_LIMIT}"
echo "  n_cores: ${JOINT_N_CORES}"

if [[ "${DRY_RUN}" == "TRUE" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=TRUE; not submitting joint fitting."
  printf "Joint command:"
  printf " %q" "${joint_cmd[@]}"
  printf "\n"
  exit 0
fi

JOINT_JOB_ID="$("${joint_cmd[@]}")"
echo "Submitted warm-start joint array job: ${JOINT_JOB_ID}"
