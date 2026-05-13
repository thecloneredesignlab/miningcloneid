#!/bin/bash
# Submit one SLURM task per seed to rebuild viz/ and report/ for an existing
# warm-start joint O2G run.

set -euo pipefail

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid"
DEFAULT_RUN_DIR="${DEFAULT_PROJECT_ROOT}/oxygen/results/fit_joint_O2G_buffering_warmstart_500seed"
DEFAULT_SUB_SCRIPT="${DEFAULT_PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/rerender_joint_warmstart_500seed_array.sub"
DEFAULT_VIZ_SCRIPT="${DEFAULT_PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/vis/viz_invivo_model_O2G_supply_demand_MAP_results.R"
DEFAULT_INVITRO_VIZ_SCRIPT="${DEFAULT_PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/vis/viz_invitro_model_O2G_supply_demand_MAP_results.R"
DEFAULT_REPORT_SCRIPT="${DEFAULT_PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/report/render_fit_report.R"
DEFAULT_DATA_DIR="${DEFAULT_PROJECT_ROOT}/data/InVivoData_Gemcitabine"
DEFAULT_TOTAL_SEEDS="500"
DEFAULT_SEED_OFFSET="1"
DEFAULT_REPORT_DT="1"
DEFAULT_TOP_N="6"
DEFAULT_N_CORES="1"
DEFAULT_REPORT_BASENAME="fit_report"
DEFAULT_RENDER_PDF="FALSE"
DEFAULT_DELETE_LEGACY_REPROT="TRUE"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_QOS="xxlarge"
DEFAULT_TIME_LIMIT="04:00:00"
DEFAULT_MEM="16G"
DEFAULT_CPUS_PER_TASK="1"
DEFAULT_DRY_RUN="FALSE"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
RUN_DIR="${RUN_DIR:-${DEFAULT_RUN_DIR}}"
SUB_SCRIPT="${SUB_SCRIPT:-${DEFAULT_SUB_SCRIPT}}"
VIZ_SCRIPT="${VIZ_SCRIPT:-${DEFAULT_VIZ_SCRIPT}}"
INVITRO_VIZ_SCRIPT="${INVITRO_VIZ_SCRIPT:-${DEFAULT_INVITRO_VIZ_SCRIPT}}"
REPORT_SCRIPT="${REPORT_SCRIPT:-${DEFAULT_REPORT_SCRIPT}}"
DATA_DIR="${DATA_DIR:-${DEFAULT_DATA_DIR}}"
TOTAL_SEEDS="${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}"
SEED_OFFSET="${SEED_OFFSET:-${DEFAULT_SEED_OFFSET}}"
REPORT_DT="${REPORT_DT:-${DEFAULT_REPORT_DT}}"
TOP_N="${TOP_N:-${DEFAULT_TOP_N}}"
N_CORES="${N_CORES:-${DEFAULT_N_CORES}}"
REPORT_BASENAME="${REPORT_BASENAME:-${DEFAULT_REPORT_BASENAME}}"
RENDER_PDF="${RENDER_PDF:-${DEFAULT_RENDER_PDF}}"
DELETE_LEGACY_REPROT="${DELETE_LEGACY_REPROT:-${DEFAULT_DELETE_LEGACY_REPROT}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
QOS="${QOS:-${DEFAULT_QOS}}"
TIME_LIMIT="${TIME_LIMIT:-${DEFAULT_TIME_LIMIT}}"
MEM="${MEM:-${DEFAULT_MEM}}"
CPUS_PER_TASK="${CPUS_PER_TASK:-${DEFAULT_CPUS_PER_TASK}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"

for numeric_name in TOTAL_SEEDS SEED_OFFSET N_CORES CPUS_PER_TASK; do
  numeric_value="${!numeric_name}"
  if ! [[ "${numeric_value}" =~ ^[0-9]+$ ]]; then
    echo "${numeric_name} must be a positive integer, got: ${numeric_value}"
    exit 1
  fi
  if (( numeric_value <= 0 )); then
    echo "${numeric_name} must be > 0, got: ${numeric_value}"
    exit 1
  fi
done

if [[ ! -d "${RUN_DIR}" ]]; then
  echo "Missing run directory: ${RUN_DIR}"
  exit 1
fi
if [[ ! -f "${SUB_SCRIPT}" ]]; then
  echo "Missing SLURM submit script: ${SUB_SCRIPT}"
  exit 1
fi
if [[ ! -f "${VIZ_SCRIPT}" ]]; then
  echo "Missing in vivo viz script: ${VIZ_SCRIPT}"
  exit 1
fi
if [[ ! -f "${INVITRO_VIZ_SCRIPT}" ]]; then
  echo "Missing in vitro viz script: ${INVITRO_VIZ_SCRIPT}"
  exit 1
fi
if [[ ! -f "${REPORT_SCRIPT}" ]]; then
  echo "Missing report script: ${REPORT_SCRIPT}"
  exit 1
fi
if [[ ! -d "${DATA_DIR}" ]]; then
  echo "Missing data directory: ${DATA_DIR}"
  exit 1
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
RUN_DIR="$(cd "${RUN_DIR}" && pwd)"
SUB_SCRIPT="$(cd "$(dirname "${SUB_SCRIPT}")" && pwd)/$(basename "${SUB_SCRIPT}")"
VIZ_SCRIPT="$(cd "$(dirname "${VIZ_SCRIPT}")" && pwd)/$(basename "${VIZ_SCRIPT}")"
INVITRO_VIZ_SCRIPT="$(cd "$(dirname "${INVITRO_VIZ_SCRIPT}")" && pwd)/$(basename "${INVITRO_VIZ_SCRIPT}")"
REPORT_SCRIPT="$(cd "$(dirname "${REPORT_SCRIPT}")" && pwd)/$(basename "${REPORT_SCRIPT}")"
DATA_DIR="$(cd "${DATA_DIR}" && pwd)"

LOG_ROOT="${LOG_ROOT:-$(dirname "${RUN_DIR}")}"
mkdir -p "${LOG_ROOT}"
LOG_ROOT="$(cd "${LOG_ROOT}" && pwd)"

export_arg="ALL"
export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
export_arg+=",RUN_DIR=${RUN_DIR}"
export_arg+=",VIZ_SCRIPT=${VIZ_SCRIPT}"
export_arg+=",INVITRO_VIZ_SCRIPT=${INVITRO_VIZ_SCRIPT}"
export_arg+=",REPORT_SCRIPT=${REPORT_SCRIPT}"
export_arg+=",DATA_DIR=${DATA_DIR}"
export_arg+=",TOTAL_SEEDS=${TOTAL_SEEDS}"
export_arg+=",SEED_OFFSET=${SEED_OFFSET}"
export_arg+=",REPORT_DT=${REPORT_DT}"
export_arg+=",TOP_N=${TOP_N}"
export_arg+=",N_CORES=${N_CORES}"
export_arg+=",REPORT_BASENAME=${REPORT_BASENAME}"
export_arg+=",RENDER_PDF=${RENDER_PDF}"
export_arg+=",DELETE_LEGACY_REPROT=${DELETE_LEGACY_REPROT}"
export_arg+=",R_MODULE=${R_MODULE}"

echo "Submitting O2G joint warm-start viz/report rebuild"
echo "  run_dir: ${RUN_DIR}"
echo "  seed tasks: ${SEED_OFFSET}-$((SEED_OFFSET + TOTAL_SEEDS - 1))"
echo "  array: 1-${TOTAL_SEEDS}"
echo "  in vivo viz: ${VIZ_SCRIPT}"
echo "  in vitro viz: ${INVITRO_VIZ_SCRIPT}"
echo "  report: ${REPORT_SCRIPT}"
echo "  report_subdir: report"
echo "  delete_legacy_reprot: ${DELETE_LEGACY_REPROT}"
echo "  qos: ${QOS}"
echo "  time_limit: ${TIME_LIMIT}"
echo "  mem: ${MEM}"
echo "  cpus_per_task: ${CPUS_PER_TASK}"

cmd=(
  sbatch
  "--qos=${QOS}"
  "--time=${TIME_LIMIT}"
  "--mem=${MEM}"
  "--cpus-per-task=${CPUS_PER_TASK}"
  "--array=1-${TOTAL_SEEDS}"
  "--output=${LOG_ROOT}/o2g_joint_rerender_%A_%a.out"
  "--error=${LOG_ROOT}/o2g_joint_rerender_%A_%a.err"
  "--export=${export_arg}"
  "${SUB_SCRIPT}"
)

if [[ "${DRY_RUN}" == "TRUE" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=TRUE; not submitting."
  printf "%q " "${cmd[@]}"
  printf "\n"
  exit 0
fi

if ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found; this launcher must be run on the HPC login node."
  exit 1
fi

"${cmd[@]}"
