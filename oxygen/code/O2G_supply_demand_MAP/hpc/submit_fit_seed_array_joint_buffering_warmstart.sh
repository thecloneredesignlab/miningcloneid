#!/bin/bash
# Build joint warm-start candidates from completed single-modality fits, then
# submit the 500-seed warm-start joint fitting SLURM array.

set -euo pipefail

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid"
DEFAULT_CONFIG_PATH="${DEFAULT_PROJECT_ROOT}/oxygen/config/O2G_supply_demand.yaml"
DEFAULT_OUT_ROOT="${DEFAULT_PROJECT_ROOT}/oxygen/results"
DEFAULT_RUN_PREFIX="fit_joint_O2G_buffering_warmstart_500seed"
DEFAULT_BUILD_SCRIPT="${DEFAULT_PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/build_joint_init_candidates.R"
DEFAULT_SUB_SCRIPT="${DEFAULT_PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc/submit_fit_seed_array_joint_buffering_warmstart.sub"
DEFAULT_TOTAL_SEEDS="500"
DEFAULT_ARRAY_TASKS="50"
DEFAULT_SEEDS_PER_TASK="10"
DEFAULT_N_CORES="20"
DEFAULT_AUTO_VIZ="TRUE"
DEFAULT_GLUCOSE="TRUE"
DEFAULT_GLUCOSE_DYNAMIC="FALSE"
DEFAULT_GLUCOSE_STRESS_MODE="coupled_to_O2"
DEFAULT_MISSEG_LOSS_SURVIVAL="buffering"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_DRY_RUN="FALSE"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
CONFIG_PATH="${CONFIG_PATH:-${DEFAULT_CONFIG_PATH}}"
OUT_ROOT="${OUT_ROOT:-${DEFAULT_OUT_ROOT}}"
RUN_PREFIX="${RUN_PREFIX:-${DEFAULT_RUN_PREFIX}}"
INVIVO_BEST_DIR="${INVIVO_BEST_DIR:-}"
INVITRO_BEST_DIR="${INVITRO_BEST_DIR:-}"
BUILD_SCRIPT="${BUILD_SCRIPT:-${DEFAULT_BUILD_SCRIPT}}"
SUB_SCRIPT="${SUB_SCRIPT:-${DEFAULT_SUB_SCRIPT}}"
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
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"

RUN_DIR="${OUT_ROOT}/${RUN_PREFIX}"
JOINT_INIT_CANDIDATES_TSV="${JOINT_INIT_CANDIDATES_TSV:-${RUN_DIR}/joint_init_candidates.tsv}"

if command -v ml >/dev/null 2>&1; then
  ml "${R_MODULE}"
elif command -v module >/dev/null 2>&1; then
  module load "${R_MODULE}"
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found after module loading."
  exit 1
fi
if ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found; this launcher must be run on the HPC login node."
  exit 1
fi

if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "Missing config file: ${CONFIG_PATH}"
  exit 1
fi
if [[ ! -f "${BUILD_SCRIPT}" ]]; then
  echo "Missing warm-start build script: ${BUILD_SCRIPT}"
  exit 1
fi
if [[ ! -f "${SUB_SCRIPT}" ]]; then
  echo "Missing SLURM submit script: ${SUB_SCRIPT}"
  exit 1
fi

mkdir -p "${RUN_DIR}"

echo "Building joint warm-start candidates"
echo "  invivo_best_dir: ${INVIVO_BEST_DIR:-<read from config>}"
echo "  invitro_best_dir: ${INVITRO_BEST_DIR:-<read from config>}"
echo "  config: ${CONFIG_PATH}"
echo "  out_tsv: ${JOINT_INIT_CANDIDATES_TSV}"
echo "  glucose: ${GLUCOSE}"
echo "  glucose_dynamic: ${GLUCOSE_DYNAMIC}"
echo "  glucose_stress_mode: ${GLUCOSE_STRESS_MODE}"
echo "  misseg_loss_survival: ${MISSEG_LOSS_SURVIVAL}"

build_args=(
  "${BUILD_SCRIPT}"
  --config="${CONFIG_PATH}"
  --out_tsv="${JOINT_INIT_CANDIDATES_TSV}"
  --glucose="${GLUCOSE}"
  --glucose_dynamic="${GLUCOSE_DYNAMIC}"
  --glucose_stress_mode="${GLUCOSE_STRESS_MODE}"
  --misseg_loss_survival="${MISSEG_LOSS_SURVIVAL}"
)
if [[ -n "${INVIVO_BEST_DIR}" ]]; then
  build_args+=("--invivo_best_dir=${INVIVO_BEST_DIR}")
fi
if [[ -n "${INVITRO_BEST_DIR}" ]]; then
  build_args+=("--invitro_best_dir=${INVITRO_BEST_DIR}")
fi

Rscript "${build_args[@]}"

echo "Warm-start candidates ready: ${JOINT_INIT_CANDIDATES_TSV}"
echo "Submitting warm-start joint array"
echo "  run_prefix: ${RUN_PREFIX}"
echo "  total_seeds: ${TOTAL_SEEDS}"
echo "  array_tasks: ${ARRAY_TASKS}"
echo "  seeds_per_task: ${SEEDS_PER_TASK}"
echo "  n_cores: ${N_CORES}"
echo "  auto_viz: ${AUTO_VIZ}"

export_arg="ALL"
export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
export_arg+=",CONFIG_PATH=${CONFIG_PATH}"
export_arg+=",OUT_ROOT=${OUT_ROOT}"
export_arg+=",RUN_PREFIX=${RUN_PREFIX}"
export_arg+=",TOTAL_SEEDS=${TOTAL_SEEDS}"
export_arg+=",ARRAY_TASKS=${ARRAY_TASKS}"
export_arg+=",SEEDS_PER_TASK=${SEEDS_PER_TASK}"
export_arg+=",N_CORES=${N_CORES}"
export_arg+=",AUTO_VIZ=${AUTO_VIZ}"
export_arg+=",GLUCOSE=${GLUCOSE}"
export_arg+=",GLUCOSE_DYNAMIC=${GLUCOSE_DYNAMIC}"
export_arg+=",GLUCOSE_STRESS_MODE=${GLUCOSE_STRESS_MODE}"
export_arg+=",MISSEG_LOSS_SURVIVAL=${MISSEG_LOSS_SURVIVAL}"
export_arg+=",JOINT_INIT_CANDIDATES_TSV=${JOINT_INIT_CANDIDATES_TSV}"
export_arg+=",R_MODULE=${R_MODULE}"

if [[ "${DRY_RUN}" == "TRUE" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=TRUE; not submitting."
  echo "sbatch --export=${export_arg} ${SUB_SCRIPT}"
  exit 0
fi

sbatch --export="${export_arg}" "${SUB_SCRIPT}"
