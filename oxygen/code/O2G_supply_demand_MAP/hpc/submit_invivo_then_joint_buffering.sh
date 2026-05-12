#!/bin/bash
# Submit in vivo-only O2G buffering 500 seeds, then automatically select the
# best long-term ploidy > 2 seed and submit warm-start joint fitting.

set -euo pipefail

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/Constant_WGD"
DEFAULT_INVIVO_RUN_PREFIX="fit_invivo_O2G_buffering_500seed"
DEFAULT_JOINT_RUN_PREFIX="fit_joint_O2G_buffering_warmstart_500seed"
DEFAULT_TOTAL_SEEDS="500"
DEFAULT_ARRAY_TASKS="500"
DEFAULT_SEEDS_PER_TASK="1"
DEFAULT_N_CORES="22"
DEFAULT_AUTO_VIZ="TRUE"
DEFAULT_GLUCOSE="TRUE"
DEFAULT_GLUCOSE_DYNAMIC="FALSE"
DEFAULT_GLUCOSE_STRESS_MODE="coupled_to_O2"
DEFAULT_MISSEG_LOSS_SURVIVAL="buffering"
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
INVIVO_SUB_SCRIPT="${INVIVO_SUB_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/submit_fit_seed_array_buffering.sub}"
CONTINUATION_SCRIPT="${CONTINUATION_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/select_invivo_best_and_submit_joint.sh}"
TOTAL_SEEDS="${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}"
ARRAY_TASKS="${ARRAY_TASKS:-${DEFAULT_ARRAY_TASKS}}"
SEEDS_PER_TASK="${SEEDS_PER_TASK:-${DEFAULT_SEEDS_PER_TASK}}"
N_CORES="${N_CORES:-${DEFAULT_N_CORES}}"
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
GLUCOSE="${GLUCOSE:-${DEFAULT_GLUCOSE}}"
GLUCOSE_DYNAMIC="${GLUCOSE_DYNAMIC:-${DEFAULT_GLUCOSE_DYNAMIC}}"
GLUCOSE_STRESS_MODE="${GLUCOSE_STRESS_MODE:-${DEFAULT_GLUCOSE_STRESS_MODE}}"
MISSEG_LOSS_SURVIVAL="${MISSEG_LOSS_SURVIVAL:-${DEFAULT_MISSEG_LOSS_SURVIVAL}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
QOS="${QOS:-${DEFAULT_QOS}}"
TIME_LIMIT="${TIME_LIMIT:-${DEFAULT_TIME_LIMIT}}"
SELECTION_HORIZON="${SELECTION_HORIZON:-${DEFAULT_SELECTION_HORIZON}}"
SELECTION_THRESHOLD_N="${SELECTION_THRESHOLD_N:-${DEFAULT_SELECTION_THRESHOLD_N}}"
SELECTION_COHORT="${SELECTION_COHORT:-${DEFAULT_SELECTION_COHORT}}"
SELECTION_DOSE="${SELECTION_DOSE:-${DEFAULT_SELECTION_DOSE}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
mkdir -p "${OUT_ROOT}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"
INVIVO_SUB_SCRIPT="$(cd "$(dirname "${INVIVO_SUB_SCRIPT}")" && pwd)/$(basename "${INVIVO_SUB_SCRIPT}")"
CONTINUATION_SCRIPT="$(cd "$(dirname "${CONTINUATION_SCRIPT}")" && pwd)/$(basename "${CONTINUATION_SCRIPT}")"

if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "Missing config file: ${CONFIG_PATH}"
  exit 1
fi
if [[ ! -f "${INVIVO_SUB_SCRIPT}" ]]; then
  echo "Missing in vivo submit script: ${INVIVO_SUB_SCRIPT}"
  exit 1
fi
if [[ ! -f "${CONTINUATION_SCRIPT}" ]]; then
  echo "Missing continuation submit script: ${CONTINUATION_SCRIPT}"
  exit 1
fi

if ! [[ "${TOTAL_SEEDS}" =~ ^[0-9]+$ && "${ARRAY_TASKS}" =~ ^[0-9]+$ && "${SEEDS_PER_TASK}" =~ ^[0-9]+$ ]]; then
  echo "TOTAL_SEEDS, ARRAY_TASKS, and SEEDS_PER_TASK must be positive integers."
  exit 1
fi
if (( TOTAL_SEEDS <= 0 || ARRAY_TASKS <= 0 || SEEDS_PER_TASK <= 0 )); then
  echo "TOTAL_SEEDS, ARRAY_TASKS, and SEEDS_PER_TASK must all be > 0."
  exit 1
fi
if (( ARRAY_TASKS * SEEDS_PER_TASK != TOTAL_SEEDS )); then
  echo "Mismatch: ARRAY_TASKS * SEEDS_PER_TASK must equal TOTAL_SEEDS."
  echo "Got ARRAY_TASKS=${ARRAY_TASKS}, SEEDS_PER_TASK=${SEEDS_PER_TASK}, TOTAL_SEEDS=${TOTAL_SEEDS}."
  exit 1
fi

if ! command -v sbatch >/dev/null 2>&1 && [[ "${DRY_RUN}" != "TRUE" && "${DRY_RUN}" != "true" && "${DRY_RUN}" != "1" ]]; then
  echo "sbatch not found; this launcher must be run on the HPC login node."
  exit 1
fi

INVIVO_RUN_DIR="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
mkdir -p "${INVIVO_RUN_DIR}"

echo "Submitting in vivo-only 500-seed array"
echo "  invivo_run_prefix: ${INVIVO_RUN_PREFIX}"
echo "  invivo_run_dir: ${INVIVO_RUN_DIR}"
echo "  total_seeds: ${TOTAL_SEEDS}"
echo "  array_tasks: ${ARRAY_TASKS}"
echo "  seeds_per_task: ${SEEDS_PER_TASK}"
echo "  qos: ${QOS}"
echo "  time_limit: ${TIME_LIMIT}"

invivo_export="ALL"
invivo_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
invivo_export+=",CONFIG_PATH=${CONFIG_PATH}"
invivo_export+=",OUT_ROOT=${OUT_ROOT}"
invivo_export+=",RUN_PREFIX=${INVIVO_RUN_PREFIX}"
invivo_export+=",TOTAL_SEEDS=${TOTAL_SEEDS}"
invivo_export+=",ARRAY_TASKS=${ARRAY_TASKS}"
invivo_export+=",SEEDS_PER_TASK=${SEEDS_PER_TASK}"
invivo_export+=",N_CORES=${N_CORES}"
invivo_export+=",AUTO_VIZ=${AUTO_VIZ}"
invivo_export+=",GLUCOSE=${GLUCOSE}"
invivo_export+=",GLUCOSE_DYNAMIC=${GLUCOSE_DYNAMIC}"
invivo_export+=",R_MODULE=${R_MODULE}"

invivo_cmd=(
  sbatch
  --parsable
  --qos="${QOS}"
  --time="${TIME_LIMIT}"
  --array="1-${ARRAY_TASKS}"
  --export="${invivo_export}"
  "${INVIVO_SUB_SCRIPT}"
)

continuation_export="ALL"
continuation_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
continuation_export+=",CONFIG_PATH=${CONFIG_PATH}"
continuation_export+=",OUT_ROOT=${OUT_ROOT}"
continuation_export+=",INVIVO_RUN_PREFIX=${INVIVO_RUN_PREFIX}"
continuation_export+=",JOINT_RUN_PREFIX=${JOINT_RUN_PREFIX}"
continuation_export+=",INVIVO_RUN_DIR=${INVIVO_RUN_DIR}"
continuation_export+=",TOTAL_SEEDS=${TOTAL_SEEDS}"
continuation_export+=",ARRAY_TASKS=${ARRAY_TASKS}"
continuation_export+=",SEEDS_PER_TASK=${SEEDS_PER_TASK}"
continuation_export+=",N_CORES=${N_CORES}"
continuation_export+=",AUTO_VIZ=${AUTO_VIZ}"
continuation_export+=",GLUCOSE=${GLUCOSE}"
continuation_export+=",GLUCOSE_DYNAMIC=${GLUCOSE_DYNAMIC}"
continuation_export+=",GLUCOSE_STRESS_MODE=${GLUCOSE_STRESS_MODE}"
continuation_export+=",MISSEG_LOSS_SURVIVAL=${MISSEG_LOSS_SURVIVAL}"
continuation_export+=",R_MODULE=${R_MODULE}"
continuation_export+=",QOS=${QOS}"
continuation_export+=",TIME_LIMIT=${TIME_LIMIT}"
continuation_export+=",SELECTION_HORIZON=${SELECTION_HORIZON}"
continuation_export+=",SELECTION_THRESHOLD_N=${SELECTION_THRESHOLD_N}"
continuation_export+=",SELECTION_COHORT=${SELECTION_COHORT}"
continuation_export+=",SELECTION_DOSE=${SELECTION_DOSE}"

if [[ "${DRY_RUN}" == "TRUE" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=TRUE; not submitting."
  printf "In vivo command:"
  printf " %q" "${invivo_cmd[@]}"
  printf "\n"
  echo "Continuation command will be:"
  echo "  sbatch --parsable --dependency=afterok:<INVIVO_JOB_ID> --qos=${QOS} --time=${TIME_LIMIT} --export=${continuation_export} ${CONTINUATION_SCRIPT}"
  exit 0
fi

INVIVO_JOB_ID="$("${invivo_cmd[@]}")"
INVIVO_JOB_ID_BASE="${INVIVO_JOB_ID%%;*}"
echo "Submitted in vivo array job: ${INVIVO_JOB_ID}"

echo "Submitting dependent selector + joint launcher"
CONTINUATION_JOB_ID="$(
  sbatch \
    --parsable \
    --dependency="afterok:${INVIVO_JOB_ID_BASE}" \
    --qos="${QOS}" \
    --time="${TIME_LIMIT}" \
    --export="${continuation_export}" \
    "${CONTINUATION_SCRIPT}"
)"
echo "Submitted continuation job: ${CONTINUATION_JOB_ID}"
echo "Dependency: afterok:${INVIVO_JOB_ID_BASE}"
