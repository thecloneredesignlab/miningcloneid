#!/bin/bash
# Select the best completed in vivo seed with long-term ploidy > 2
# (chromosome count > 44), then submit warm-start joint fitting.

#SBATCH --job-name=o2g_sel_joint
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --qos=xxlarge
#SBATCH --time=12:00:00
#SBATCH --output=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/o2g_select_joint_%j.out
#SBATCH --error=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/o2g_select_joint_%j.err

set -euo pipefail

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid"
DEFAULT_INVIVO_RUN_PREFIX="fit_invivo_O2G_buffering_500seed"
DEFAULT_JOINT_RUN_PREFIX="fit_joint_O2G_buffering_warmstart_500seed"
DEFAULT_TOTAL_SEEDS="500"
DEFAULT_ARRAY_TASKS="500"
DEFAULT_SEEDS_PER_TASK="1"
DEFAULT_N_CORES="22"
DEFAULT_AUTO_VIZ="TRUE"
DEFAULT_GLUCOSE="TRUE"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_QOS="xxlarge"
DEFAULT_TIME_LIMIT="12:00:00"
DEFAULT_SELECTION_HORIZON="1000"
DEFAULT_SELECTION_THRESHOLD_N="44"
DEFAULT_SELECTION_COHORT="2N"
DEFAULT_SELECTION_DOSE="ALL"
DEFAULT_DRY_RUN="FALSE"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
CONFIG_PATH="${CONFIG_PATH:-${PROJECT_ROOT}/oxygen/config/O2G_supply_demand.yaml}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/oxygen/results}"
INVIVO_RUN_PREFIX="${INVIVO_RUN_PREFIX:-${DEFAULT_INVIVO_RUN_PREFIX}}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-${DEFAULT_JOINT_RUN_PREFIX}}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-${OUT_ROOT}/${INVIVO_RUN_PREFIX}}"
SELECTOR_SCRIPT="${SELECTOR_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/select_invivo_best_long_ploidy_seed.R}"
JOINT_LAUNCHER="${JOINT_LAUNCHER:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/submit_fit_seed_array_joint_buffering_warmstart.sh}"
TOTAL_SEEDS="${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}"
ARRAY_TASKS="${ARRAY_TASKS:-${DEFAULT_ARRAY_TASKS}}"
SEEDS_PER_TASK="${SEEDS_PER_TASK:-${DEFAULT_SEEDS_PER_TASK}}"
N_CORES="${N_CORES:-${DEFAULT_N_CORES}}"
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
GLUCOSE="${GLUCOSE:-${DEFAULT_GLUCOSE}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
QOS="${QOS:-${DEFAULT_QOS}}"
TIME_LIMIT="${TIME_LIMIT:-${DEFAULT_TIME_LIMIT}}"
SELECTION_HORIZON="${SELECTION_HORIZON:-${DEFAULT_SELECTION_HORIZON}}"
SELECTION_THRESHOLD_N="${SELECTION_THRESHOLD_N:-${DEFAULT_SELECTION_THRESHOLD_N}}"
SELECTION_COHORT="${SELECTION_COHORT:-${DEFAULT_SELECTION_COHORT}}"
SELECTION_DOSE="${SELECTION_DOSE:-${DEFAULT_SELECTION_DOSE}}"
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
INVIVO_RUN_DIR="$(cd "${INVIVO_RUN_DIR}" && pwd)"
SELECTOR_SCRIPT="$(cd "$(dirname "${SELECTOR_SCRIPT}")" && pwd)/$(basename "${SELECTOR_SCRIPT}")"
JOINT_LAUNCHER="$(cd "$(dirname "${JOINT_LAUNCHER}")" && pwd)/$(basename "${JOINT_LAUNCHER}")"

if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "Missing config file: ${CONFIG_PATH}"
  exit 1
fi
if [[ ! -f "${SELECTOR_SCRIPT}" ]]; then
  echo "Missing selector script: ${SELECTOR_SCRIPT}"
  exit 1
fi
if [[ ! -f "${JOINT_LAUNCHER}" ]]; then
  echo "Missing joint launcher: ${JOINT_LAUNCHER}"
  exit 1
fi

SELECTION_TSV="${SELECTION_TSV:-${INVIVO_RUN_DIR}/best_long_ploidy_gt2_seed.tsv}"
BEST_DIR_FILE="${BEST_DIR_FILE:-${INVIVO_RUN_DIR}/best_long_ploidy_gt2_seed.dir}"

echo "Selecting best in vivo seed for joint warm start"
echo "  invivo_run_dir: ${INVIVO_RUN_DIR}"
echo "  horizon: ${SELECTION_HORIZON}"
echo "  threshold_N: ${SELECTION_THRESHOLD_N}"
echo "  cohort: ${SELECTION_COHORT}"
echo "  dose: ${SELECTION_DOSE}"
echo "  selection_tsv: ${SELECTION_TSV}"
echo "  best_dir_file: ${BEST_DIR_FILE}"

Rscript "${SELECTOR_SCRIPT}" \
  --invivo_dir="${INVIVO_RUN_DIR}" \
  --out_tsv="${SELECTION_TSV}" \
  --best_dir_file="${BEST_DIR_FILE}" \
  --horizon="${SELECTION_HORIZON}" \
  --threshold_N="${SELECTION_THRESHOLD_N}" \
  --cohort="${SELECTION_COHORT}" \
  --dose="${SELECTION_DOSE}"

INVIVO_BEST_DIR="$(tr -d '[:space:]' < "${BEST_DIR_FILE}")"
if [[ -z "${INVIVO_BEST_DIR}" || ! -d "${INVIVO_BEST_DIR}" ]]; then
  echo "Selected INVIVO_BEST_DIR is missing or invalid: ${INVIVO_BEST_DIR}"
  exit 1
fi

echo "Submitting joint warm-start fit"
echo "  invivo_best_dir: ${INVIVO_BEST_DIR}"
echo "  joint_run_prefix: ${JOINT_RUN_PREFIX}"
echo "  qos: ${QOS}"
echo "  time_limit: ${TIME_LIMIT}"

if [[ "${DRY_RUN}" == "TRUE" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=TRUE; not submitting joint fitting."
  echo "INVIVO_BEST_DIR=${INVIVO_BEST_DIR} RUN_PREFIX=${JOINT_RUN_PREFIX} QOS=${QOS} TIME_LIMIT=${TIME_LIMIT} bash ${JOINT_LAUNCHER}"
  exit 0
fi

PROJECT_ROOT="${PROJECT_ROOT}" \
CONFIG_PATH="${CONFIG_PATH}" \
OUT_ROOT="${OUT_ROOT}" \
RUN_PREFIX="${JOINT_RUN_PREFIX}" \
INVIVO_BEST_DIR="${INVIVO_BEST_DIR}" \
TOTAL_SEEDS="${TOTAL_SEEDS}" \
ARRAY_TASKS="${ARRAY_TASKS}" \
SEEDS_PER_TASK="${SEEDS_PER_TASK}" \
N_CORES="${N_CORES}" \
AUTO_VIZ="${AUTO_VIZ}" \
GLUCOSE="${GLUCOSE}" \
R_MODULE="${R_MODULE}" \
QOS="${QOS}" \
TIME_LIMIT="${TIME_LIMIT}" \
bash "${JOINT_LAUNCHER}"
