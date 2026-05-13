#!/bin/bash
# Submit in vivo and in vitro O2G buffering seed arrays in parallel. After both
# arrays complete successfully, submit one manager job that summarizes both runs,
# selects the best seed from each modality, and submits warm-start joint fitting.

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
DEFAULT_MISSEG_LOSS_SURVIVAL="buffering"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_INVIVO_QOS="xxlarge"
DEFAULT_INVIVO_TIME_LIMIT="12:00:00"
DEFAULT_INVITRO_QOS="xxlarge"
DEFAULT_INVITRO_TIME_LIMIT="12:00:00"
DEFAULT_MANAGER_QOS="xxlarge"
DEFAULT_MANAGER_TIME_LIMIT="12:00:00"
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
MANAGER_SCRIPT="${MANAGER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/manage_invivo_invitro_best_and_submit_joint.sh}"
TOTAL_SEEDS="${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}"
ARRAY_TASKS="${ARRAY_TASKS:-${DEFAULT_ARRAY_TASKS}}"
SEEDS_PER_TASK="${SEEDS_PER_TASK:-${DEFAULT_SEEDS_PER_TASK}}"
N_CORES="${N_CORES:-${DEFAULT_N_CORES}}"
INVIVO_N_CORES="${INVIVO_N_CORES:-${N_CORES}}"
INVITRO_N_CORES="${INVITRO_N_CORES:-${N_CORES}}"
JOINT_N_CORES="${JOINT_N_CORES:-${N_CORES}}"
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
GLUCOSE="${GLUCOSE:-${DEFAULT_GLUCOSE}}"
MISSEG_LOSS_SURVIVAL="${MISSEG_LOSS_SURVIVAL:-${DEFAULT_MISSEG_LOSS_SURVIVAL}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
INVIVO_QOS="${INVIVO_QOS:-${DEFAULT_INVIVO_QOS}}"
INVIVO_TIME_LIMIT="${INVIVO_TIME_LIMIT:-${DEFAULT_INVIVO_TIME_LIMIT}}"
INVITRO_QOS="${INVITRO_QOS:-${DEFAULT_INVITRO_QOS}}"
INVITRO_TIME_LIMIT="${INVITRO_TIME_LIMIT:-${DEFAULT_INVITRO_TIME_LIMIT}}"
MANAGER_QOS="${MANAGER_QOS:-${DEFAULT_MANAGER_QOS}}"
MANAGER_TIME_LIMIT="${MANAGER_TIME_LIMIT:-${DEFAULT_MANAGER_TIME_LIMIT}}"
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
MANAGER_SCRIPT="$(cd "$(dirname "${MANAGER_SCRIPT}")" && pwd)/$(basename "${MANAGER_SCRIPT}")"

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
if [[ ! -f "${MANAGER_SCRIPT}" ]]; then
  echo "Missing manager script: ${MANAGER_SCRIPT}"
  exit 1
fi

for var_name in TOTAL_SEEDS ARRAY_TASKS SEEDS_PER_TASK INVIVO_N_CORES INVITRO_N_CORES JOINT_N_CORES ITERMAX DE_STEPTOL NP; do
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

echo "Submitting in vivo and in vitro seed arrays"
echo "  project_root: ${PROJECT_ROOT}"
echo "  out_root: ${OUT_ROOT}"
echo "  invivo_run_dir: ${INVIVO_RUN_DIR}"
echo "  invitro_run_dir: ${INVITRO_RUN_DIR}"
echo "  joint_run_dir: ${JOINT_RUN_DIR}"
echo "  total_seeds: ${TOTAL_SEEDS}"
echo "  array_tasks: ${ARRAY_TASKS}"
echo "  seeds_per_task: ${SEEDS_PER_TASK}"
echo "  invivo qos/time: ${INVIVO_QOS} / ${INVIVO_TIME_LIMIT}"
echo "  invitro qos/time: ${INVITRO_QOS} / ${INVITRO_TIME_LIMIT}"
echo "  manager qos/time: ${MANAGER_QOS} / ${MANAGER_TIME_LIMIT}"
echo "  joint qos/time: ${JOINT_QOS} / ${JOINT_TIME_LIMIT}"

invivo_export="ALL"
invivo_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
invivo_export+=",CONFIG_PATH=${CONFIG_PATH}"
invivo_export+=",OUT_ROOT=${OUT_ROOT}"
invivo_export+=",RUN_PREFIX=${INVIVO_RUN_PREFIX}"
invivo_export+=",TOTAL_SEEDS=${TOTAL_SEEDS}"
invivo_export+=",ARRAY_TASKS=${ARRAY_TASKS}"
invivo_export+=",SEEDS_PER_TASK=${SEEDS_PER_TASK}"
invivo_export+=",N_CORES=${INVIVO_N_CORES}"
invivo_export+=",AUTO_VIZ=${AUTO_VIZ}"
invivo_export+=",GLUCOSE=${GLUCOSE}"
invivo_export+=",R_MODULE=${R_MODULE}"

invitro_export="ALL"
invitro_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
invitro_export+=",OUT_ROOT=${OUT_ROOT}"
invitro_export+=",RUN_PREFIX=${INVITRO_RUN_PREFIX}"
invitro_export+=",TOTAL_SEEDS=${TOTAL_SEEDS}"
invitro_export+=",ARRAY_TASKS=${ARRAY_TASKS}"
invitro_export+=",SEEDS_PER_TASK=${SEEDS_PER_TASK}"
invitro_export+=",N_CORES=${INVITRO_N_CORES}"
invitro_export+=",R_MODULE=${R_MODULE}"
invitro_export+=",MISSEG_LOSS_SURVIVAL=${MISSEG_LOSS_SURVIVAL}"
invitro_export+=",PARAMETER_TABLE=${PARAMETER_TABLE}"
invitro_export+=",FIT_OBJECTS_DIR=${FIT_OBJECTS_DIR}"
invitro_export+=",FLOW_DENSITY_PATH=${FLOW_DENSITY_PATH}"
invitro_export+=",ITERMAX=${ITERMAX}"
invitro_export+=",DE_RELTOL=${DE_RELTOL}"
invitro_export+=",DE_STEPTOL=${DE_STEPTOL}"
invitro_export+=",NP=${NP}"
invitro_export+=",AUTO_VIZ=${AUTO_VIZ}"

manager_export="ALL"
manager_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
manager_export+=",CONFIG_PATH=${CONFIG_PATH}"
manager_export+=",OUT_ROOT=${OUT_ROOT}"
manager_export+=",INVIVO_RUN_PREFIX=${INVIVO_RUN_PREFIX}"
manager_export+=",INVITRO_RUN_PREFIX=${INVITRO_RUN_PREFIX}"
manager_export+=",JOINT_RUN_PREFIX=${JOINT_RUN_PREFIX}"
manager_export+=",INVIVO_RUN_DIR=${INVIVO_RUN_DIR}"
manager_export+=",INVITRO_RUN_DIR=${INVITRO_RUN_DIR}"
manager_export+=",TOTAL_SEEDS=${TOTAL_SEEDS}"
manager_export+=",ARRAY_TASKS=${ARRAY_TASKS}"
manager_export+=",SEEDS_PER_TASK=${SEEDS_PER_TASK}"
manager_export+=",N_CORES=${JOINT_N_CORES}"
manager_export+=",AUTO_VIZ=${AUTO_VIZ}"
manager_export+=",GLUCOSE=${GLUCOSE}"
manager_export+=",MISSEG_LOSS_SURVIVAL=${MISSEG_LOSS_SURVIVAL}"
manager_export+=",R_MODULE=${R_MODULE}"
manager_export+=",JOINT_QOS=${JOINT_QOS}"
manager_export+=",JOINT_TIME_LIMIT=${JOINT_TIME_LIMIT}"
manager_export+=",JOINT_N_CORES=${JOINT_N_CORES}"

invivo_cmd=(
  sbatch
  --parsable
  --job-name=o2g_ivv_B
  --qos="${INVIVO_QOS}"
  --time="${INVIVO_TIME_LIMIT}"
  --cpus-per-task="${INVIVO_N_CORES}"
  --array="1-${ARRAY_TASKS}"
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
  --array="1-${ARRAY_TASKS}"
  --output="${OUT_ROOT}/o2g_invitro_%A_%a.out"
  --error="${OUT_ROOT}/o2g_invitro_%A_%a.err"
  --export="${invitro_export}"
  "${INVITRO_SUB_SCRIPT}"
)

if [[ "${DRY_RUN}" == "TRUE" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=TRUE; not submitting."
  printf "In vivo command:"
  printf " %q" "${invivo_cmd[@]}"
  printf "\n"
  printf "In vitro command:"
  printf " %q" "${invitro_cmd[@]}"
  printf "\n"
  echo "Manager command will be:"
  echo "  sbatch --parsable --dependency=afterok:<INVIVO_JOB_ID>:<INVITRO_JOB_ID> --job-name=o2g_best_joint --qos=${MANAGER_QOS} --time=${MANAGER_TIME_LIMIT} --output=${OUT_ROOT}/o2g_best_joint_%j.out --error=${OUT_ROOT}/o2g_best_joint_%j.err --export=${manager_export} ${MANAGER_SCRIPT}"
  exit 0
fi

INVIVO_JOB_ID="$("${invivo_cmd[@]}")"
INVIVO_JOB_ID_BASE="${INVIVO_JOB_ID%%;*}"
echo "Submitted in vivo array job: ${INVIVO_JOB_ID}"

INVITRO_JOB_ID="$("${invitro_cmd[@]}")"
INVITRO_JOB_ID_BASE="${INVITRO_JOB_ID%%;*}"
echo "Submitted in vitro array job: ${INVITRO_JOB_ID}"

echo "Submitting dependent summary/select manager + joint launcher"
MANAGER_JOB_ID="$(
  sbatch \
    --parsable \
    --dependency="afterok:${INVIVO_JOB_ID_BASE}:${INVITRO_JOB_ID_BASE}" \
    --job-name=o2g_best_joint \
    --qos="${MANAGER_QOS}" \
    --time="${MANAGER_TIME_LIMIT}" \
    --output="${OUT_ROOT}/o2g_best_joint_%j.out" \
    --error="${OUT_ROOT}/o2g_best_joint_%j.err" \
    --export="${manager_export}" \
    "${MANAGER_SCRIPT}"
)"
echo "Submitted manager job: ${MANAGER_JOB_ID}"
echo "Dependency: afterok:${INVIVO_JOB_ID_BASE}:${INVITRO_JOB_ID_BASE}"
