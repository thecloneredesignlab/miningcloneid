#!/bin/bash
# Submit in vivo, in vitro, and joint O2G buffering seed arrays.

set -euo pipefail

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid"
DEFAULT_INVIVO_RUN_PREFIX="fit_invivo_O2G_buffering_200seed"
DEFAULT_INVITRO_RUN_PREFIX="fit_invitro_O2G_buffering_200seed"
DEFAULT_JOINT_RUN_PREFIX="fit_joint_O2G_buffering_500seed"
DEFAULT_INVIVO_TOTAL_SEEDS="200"
DEFAULT_INVITRO_TOTAL_SEEDS="200"
DEFAULT_JOINT_TOTAL_SEEDS="500"
DEFAULT_SEEDS_PER_TASK="1"
DEFAULT_N_CORES="22"
DEFAULT_MEM="16G"
DEFAULT_AUTO_VIZ="TRUE"
DEFAULT_GLUCOSE="TRUE"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_INVIVO_QOS="small"
DEFAULT_INVIVO_TIME_LIMIT="3:00:00"
DEFAULT_INVITRO_QOS="xlarge"
DEFAULT_INVITRO_TIME_LIMIT="4:00:00"
DEFAULT_JOINT_QOS="xxlarge"
DEFAULT_JOINT_TIME_LIMIT="12:00:00"
DEFAULT_ITERMAX="500"
DEFAULT_DE_RELTOL="1e-4"
DEFAULT_DE_STEPTOL="25"
DEFAULT_NP="80"
DEFAULT_DRY_RUN="FALSE"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
CONFIG_PATH="${CONFIG_PATH:-${PROJECT_ROOT}/oxygen/config/O2G_supply_demand.yaml}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/oxygen/results}"
INVIVO_RUN_PREFIX="${INVIVO_RUN_PREFIX:-${DEFAULT_INVIVO_RUN_PREFIX}}"
INVITRO_RUN_PREFIX="${INVITRO_RUN_PREFIX:-${DEFAULT_INVITRO_RUN_PREFIX}}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-${DEFAULT_JOINT_RUN_PREFIX}}"
INVIVO_SUB_SCRIPT="${INVIVO_SUB_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/submit_fit_seed_array_buffering.sub}"
INVITRO_SUB_SCRIPT="${INVITRO_SUB_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/submit_fit_seed_array_invitro_buffering.sub}"
JOINT_SUB_SCRIPT="${JOINT_SUB_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/submit_fit_seed_array_joint_buffering.sub}"
INVIVO_TOTAL_SEEDS="${INVIVO_TOTAL_SEEDS:-${TOTAL_SEEDS:-${DEFAULT_INVIVO_TOTAL_SEEDS}}}"
INVITRO_TOTAL_SEEDS="${INVITRO_TOTAL_SEEDS:-${TOTAL_SEEDS:-${DEFAULT_INVITRO_TOTAL_SEEDS}}}"
JOINT_TOTAL_SEEDS="${JOINT_TOTAL_SEEDS:-${TOTAL_SEEDS:-${DEFAULT_JOINT_TOTAL_SEEDS}}}"
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
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
GLUCOSE="${GLUCOSE:-${DEFAULT_GLUCOSE}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
INVIVO_QOS="${INVIVO_QOS:-${DEFAULT_INVIVO_QOS}}"
INVIVO_TIME_LIMIT="${INVIVO_TIME_LIMIT:-${DEFAULT_INVIVO_TIME_LIMIT}}"
INVITRO_QOS="${INVITRO_QOS:-${DEFAULT_INVITRO_QOS}}"
INVITRO_TIME_LIMIT="${INVITRO_TIME_LIMIT:-${DEFAULT_INVITRO_TIME_LIMIT}}"
JOINT_QOS="${JOINT_QOS:-${DEFAULT_JOINT_QOS}}"
JOINT_TIME_LIMIT="${JOINT_TIME_LIMIT:-${DEFAULT_JOINT_TIME_LIMIT}}"
PARAMETER_TABLE="${PARAMETER_TABLE:-${PROJECT_ROOT}/oxygen/data/O2G_supply_demand/parameter_table_invitro_buffering.csv}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv}"
ITERMAX="${ITERMAX:-${DEFAULT_ITERMAX}}"
DE_RELTOL="${DE_RELTOL:-${DEFAULT_DE_RELTOL}}"
DE_STEPTOL="${DE_STEPTOL:-${DEFAULT_DE_STEPTOL}}"
NP="${NP:-${DEFAULT_NP}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
mkdir -p "${OUT_ROOT}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"
INVIVO_SUB_SCRIPT="$(cd "$(dirname "${INVIVO_SUB_SCRIPT}")" && pwd)/$(basename "${INVIVO_SUB_SCRIPT}")"
INVITRO_SUB_SCRIPT="$(cd "$(dirname "${INVITRO_SUB_SCRIPT}")" && pwd)/$(basename "${INVITRO_SUB_SCRIPT}")"
JOINT_SUB_SCRIPT="$(cd "$(dirname "${JOINT_SUB_SCRIPT}")" && pwd)/$(basename "${JOINT_SUB_SCRIPT}")"

if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "Missing config file: ${CONFIG_PATH}"
  exit 1
fi
if [[ ! -f "${INVIVO_SUB_SCRIPT}" ]]; then
  echo "Missing in vivo submit script: ${INVIVO_SUB_SCRIPT}"
  exit 1
fi
if [[ ! -f "${INVITRO_SUB_SCRIPT}" ]]; then
  echo "Missing in vitro submit script: ${INVITRO_SUB_SCRIPT}"
  exit 1
fi
if [[ ! -f "${JOINT_SUB_SCRIPT}" ]]; then
  echo "Missing joint submit script: ${JOINT_SUB_SCRIPT}"
  exit 1
fi

for var_name in \
  INVIVO_TOTAL_SEEDS INVIVO_ARRAY_TASKS INVIVO_SEEDS_PER_TASK \
  INVITRO_TOTAL_SEEDS INVITRO_ARRAY_TASKS INVITRO_SEEDS_PER_TASK \
  JOINT_TOTAL_SEEDS JOINT_ARRAY_TASKS JOINT_SEEDS_PER_TASK \
  INVIVO_N_CORES INVITRO_N_CORES JOINT_N_CORES ITERMAX DE_STEPTOL NP; do
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

check_seed_plan() {
  local label="$1"
  local total="$2"
  local tasks="$3"
  local per_task="$4"
  if (( tasks * per_task != total )); then
    echo "Mismatch for ${label}: ARRAY_TASKS * SEEDS_PER_TASK must equal TOTAL_SEEDS."
    echo "Got ${label}_ARRAY_TASKS=${tasks}, ${label}_SEEDS_PER_TASK=${per_task}, ${label}_TOTAL_SEEDS=${total}."
    exit 1
  fi
}

check_seed_plan "INVIVO" "${INVIVO_TOTAL_SEEDS}" "${INVIVO_ARRAY_TASKS}" "${INVIVO_SEEDS_PER_TASK}"
check_seed_plan "INVITRO" "${INVITRO_TOTAL_SEEDS}" "${INVITRO_ARRAY_TASKS}" "${INVITRO_SEEDS_PER_TASK}"
check_seed_plan "JOINT" "${JOINT_TOTAL_SEEDS}" "${JOINT_ARRAY_TASKS}" "${JOINT_SEEDS_PER_TASK}"
if ! [[ "${DE_RELTOL}" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]; then
  echo "DE_RELTOL must be a positive numeric value, got: ${DE_RELTOL}"
  exit 1
fi

if ! command -v sbatch >/dev/null 2>&1 && [[ "${DRY_RUN}" != "TRUE" && "${DRY_RUN}" != "true" && "${DRY_RUN}" != "1" ]]; then
  echo "sbatch not found; this launcher must be run on the HPC login node."
  exit 1
fi

INVIVO_RUN_DIR="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
INVITRO_RUN_DIR="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
JOINT_RUN_DIR="${OUT_ROOT}/${JOINT_RUN_PREFIX}"
mkdir -p "${INVIVO_RUN_DIR}" "${INVITRO_RUN_DIR}" "${JOINT_RUN_DIR}"

echo "Submitting in vivo, in vitro, and joint seed arrays"
echo "  project_root: ${PROJECT_ROOT}"
echo "  out_root: ${OUT_ROOT}"
echo "  invivo_run_dir: ${INVIVO_RUN_DIR}"
echo "  invitro_run_dir: ${INVITRO_RUN_DIR}"
echo "  joint_run_dir: ${JOINT_RUN_DIR}"
echo "  invivo seeds: total=${INVIVO_TOTAL_SEEDS}, array_tasks=${INVIVO_ARRAY_TASKS}, seeds_per_task=${INVIVO_SEEDS_PER_TASK}"
echo "  invitro seeds: total=${INVITRO_TOTAL_SEEDS}, array_tasks=${INVITRO_ARRAY_TASKS}, seeds_per_task=${INVITRO_SEEDS_PER_TASK}"
echo "  joint seeds: total=${JOINT_TOTAL_SEEDS}, array_tasks=${JOINT_ARRAY_TASKS}, seeds_per_task=${JOINT_SEEDS_PER_TASK}"
echo "  invivo resources: qos=${INVIVO_QOS}, time=${INVIVO_TIME_LIMIT}, cpus=${INVIVO_N_CORES}, mem=${INVIVO_MEM}"
echo "  invitro resources: qos=${INVITRO_QOS}, time=${INVITRO_TIME_LIMIT}, cpus=${INVITRO_N_CORES}, mem=${INVITRO_MEM}"
echo "  joint resources: qos=${JOINT_QOS}, time=${JOINT_TIME_LIMIT}, cpus=${JOINT_N_CORES}, mem=${JOINT_MEM}"

invivo_export="ALL"
invivo_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
invivo_export+=",CONFIG_PATH=${CONFIG_PATH}"
invivo_export+=",OUT_ROOT=${OUT_ROOT}"
invivo_export+=",RUN_PREFIX=${INVIVO_RUN_PREFIX}"
invivo_export+=",TOTAL_SEEDS=${INVIVO_TOTAL_SEEDS}"
invivo_export+=",ARRAY_TASKS=${INVIVO_ARRAY_TASKS}"
invivo_export+=",SEEDS_PER_TASK=${INVIVO_SEEDS_PER_TASK}"
invivo_export+=",N_CORES=${INVIVO_N_CORES}"
invivo_export+=",AUTO_VIZ=${AUTO_VIZ}"
invivo_export+=",GLUCOSE=${GLUCOSE}"
invivo_export+=",R_MODULE=${R_MODULE}"

invitro_export="ALL"
invitro_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
invitro_export+=",OUT_ROOT=${OUT_ROOT}"
invitro_export+=",RUN_PREFIX=${INVITRO_RUN_PREFIX}"
invitro_export+=",TOTAL_SEEDS=${INVITRO_TOTAL_SEEDS}"
invitro_export+=",ARRAY_TASKS=${INVITRO_ARRAY_TASKS}"
invitro_export+=",SEEDS_PER_TASK=${INVITRO_SEEDS_PER_TASK}"
invitro_export+=",N_CORES=${INVITRO_N_CORES}"
invitro_export+=",R_MODULE=${R_MODULE}"
invitro_export+=",PARAMETER_TABLE=${PARAMETER_TABLE}"
invitro_export+=",FIT_OBJECTS_DIR=${FIT_OBJECTS_DIR}"
invitro_export+=",FLOW_DENSITY_PATH=${FLOW_DENSITY_PATH}"
invitro_export+=",ITERMAX=${ITERMAX}"
invitro_export+=",DE_RELTOL=${DE_RELTOL}"
invitro_export+=",DE_STEPTOL=${DE_STEPTOL}"
invitro_export+=",NP=${NP}"
invitro_export+=",AUTO_VIZ=${AUTO_VIZ}"

joint_export="ALL"
joint_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
joint_export+=",CONFIG_PATH=${CONFIG_PATH}"
joint_export+=",OUT_ROOT=${OUT_ROOT}"
joint_export+=",RUN_PREFIX=${JOINT_RUN_PREFIX}"
joint_export+=",TOTAL_SEEDS=${JOINT_TOTAL_SEEDS}"
joint_export+=",ARRAY_TASKS=${JOINT_ARRAY_TASKS}"
joint_export+=",SEEDS_PER_TASK=${JOINT_SEEDS_PER_TASK}"
joint_export+=",N_CORES=${JOINT_N_CORES}"
joint_export+=",AUTO_VIZ=${AUTO_VIZ}"
joint_export+=",GLUCOSE=${GLUCOSE}"
joint_export+=",R_MODULE=${R_MODULE}"

invivo_cmd=(
  sbatch
  --parsable
  --job-name=o2g_ivv_B
  --qos="${INVIVO_QOS}"
  --time="${INVIVO_TIME_LIMIT}"
  --cpus-per-task="${INVIVO_N_CORES}"
  --mem="${INVIVO_MEM}"
  --array="1-${INVIVO_ARRAY_TASKS}"
  --output="${OUT_ROOT}/o2g_invivo_%A_%a.out"
  --error="${OUT_ROOT}/o2g_invivo_%A_%a.err"
  --export="${invivo_export}"
  "${INVIVO_SUB_SCRIPT}"
)

invitro_cmd=(
  sbatch
  --parsable
  --job-name=o2g_ivt_B
  --qos="${INVITRO_QOS}"
  --time="${INVITRO_TIME_LIMIT}"
  --cpus-per-task="${INVITRO_N_CORES}"
  --mem="${INVITRO_MEM}"
  --array="1-${INVITRO_ARRAY_TASKS}"
  --output="${OUT_ROOT}/o2g_invitro_%A_%a.out"
  --error="${OUT_ROOT}/o2g_invitro_%A_%a.err"
  --export="${invitro_export}"
  "${INVITRO_SUB_SCRIPT}"
)

joint_cmd=(
  sbatch
  --parsable
  --job-name=o2g_joint_B
  --qos="${JOINT_QOS}"
  --time="${JOINT_TIME_LIMIT}"
  --cpus-per-task="${JOINT_N_CORES}"
  --mem="${JOINT_MEM}"
  --array="1-${JOINT_ARRAY_TASKS}"
  --output="${OUT_ROOT}/o2g_joint_fit_%A_%a.out"
  --error="${OUT_ROOT}/o2g_joint_fit_%A_%a.err"
  --export="${joint_export}"
  "${JOINT_SUB_SCRIPT}"
)

if [[ "${DRY_RUN}" == "TRUE" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=TRUE; not submitting."
  printf "In vivo command:"
  printf " %q" "${invivo_cmd[@]}"
  printf "\n"
  printf "In vitro command:"
  printf " %q" "${invitro_cmd[@]}"
  printf "\n"
  printf "Joint command:"
  printf " %q" "${joint_cmd[@]}"
  printf "\n"
  exit 0
fi

INVIVO_JOB_ID="$("${invivo_cmd[@]}")"
echo "Submitted in vivo array job: ${INVIVO_JOB_ID}"

INVITRO_JOB_ID="$("${invitro_cmd[@]}")"
echo "Submitted in vitro array job: ${INVITRO_JOB_ID}"

JOINT_JOB_ID="$("${joint_cmd[@]}")"
echo "Submitted joint array job: ${JOINT_JOB_ID}"
