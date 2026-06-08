#!/bin/bash
# Submit one interleaved 400-task array for O2G fit_invivo and fit_invitro.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../../../.." && pwd)}"
SUB_SCRIPT="${SUB_SCRIPT:-${SCRIPT_DIR}/submit_fit_seed_array_invivo_invitro_interleaved.sub}"
POSTPROCESS_SCRIPT="${POSTPROCESS_SCRIPT:-${SCRIPT_DIR}/postprocess_extra_results.sh}"
EXTRA_RESULTS_SCRIPT="${EXTRA_RESULTS_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R}"
RUNNER_SCRIPT="${RUNNER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh}"

CONFIG_PATH="${CONFIG_PATH:-${PROJECT_ROOT}/oxygen/config/O2G_supply_demand.yaml}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/oxygen/results}"
RUN_PREFIX_INVIVO="${RUN_PREFIX_INVIVO:-fit_invivo_O2G_buffering_200seed}"
RUN_PREFIX_INVITRO="${RUN_PREFIX_INVITRO:-fit_invitro_O2G_buffering_200seed}"

TOTAL_SEEDS_PER_MODE="${TOTAL_SEEDS_PER_MODE:-200}"
ARRAY_TASKS="${ARRAY_TASKS:-400}"
N_CORES="${N_CORES:-22}"
MEM="${MEM:-32G}"
AUTO_VIZ="${AUTO_VIZ:-TRUE}"
GLUCOSE="${GLUCOSE:-FALSE}"
R_MODULE="${R_MODULE:-R/4.4}"
ITERMAX="${ITERMAX:-500}"
DE_RELTOL="${DE_RELTOL:-1e-4}"
DE_STEPTOL="${DE_STEPTOL:-25}"
NP="${NP:-80}"
PARAMETER_TABLE="${PARAMETER_TABLE:-${PROJECT_ROOT}/oxygen/data/O2G_supply_demand/parameter_table_invitro_buffering.csv}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-TRUE}"
DRY_RUN="${DRY_RUN:-FALSE}"

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

print_command() {
  local label="$1"
  shift
  printf "%s:" "${label}"
  printf " %q" "$@"
  printf "\n"
}

require_positive_int() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]] || (( value <= 0 )); then
    echo "${name} must be a positive integer, got: ${value}" >&2
    exit 2
  fi
}

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
SUB_SCRIPT="$(cd "$(dirname "${SUB_SCRIPT}")" && pwd)/$(basename "${SUB_SCRIPT}")"
POSTPROCESS_SCRIPT="$(cd "$(dirname "${POSTPROCESS_SCRIPT}")" && pwd)/$(basename "${POSTPROCESS_SCRIPT}")"
EXTRA_RESULTS_SCRIPT="$(cd "$(dirname "${EXTRA_RESULTS_SCRIPT}")" && pwd)/$(basename "${EXTRA_RESULTS_SCRIPT}")"
RUNNER_SCRIPT="$(cd "$(dirname "${RUNNER_SCRIPT}")" && pwd)/$(basename "${RUNNER_SCRIPT}")"
CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
mkdir -p "${OUT_ROOT}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"
PARAMETER_TABLE="$(cd "$(dirname "${PARAMETER_TABLE}")" && pwd)/$(basename "${PARAMETER_TABLE}")"
FIT_OBJECTS_DIR="$(cd "${FIT_OBJECTS_DIR}" && pwd)"
FLOW_DENSITY_PATH="$(cd "$(dirname "${FLOW_DENSITY_PATH}")" && pwd)/$(basename "${FLOW_DENSITY_PATH}")"

for path in "${SUB_SCRIPT}" "${POSTPROCESS_SCRIPT}" "${EXTRA_RESULTS_SCRIPT}" "${RUNNER_SCRIPT}" \
            "${CONFIG_PATH}" "${PARAMETER_TABLE}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
done
if [[ ! -d "${FIT_OBJECTS_DIR}" ]]; then
  echo "Missing fit_objects_dir: ${FIT_OBJECTS_DIR}" >&2
  exit 1
fi

for name in TOTAL_SEEDS_PER_MODE ARRAY_TASKS N_CORES ITERMAX DE_STEPTOL NP; do
  require_positive_int "${name}" "${!name}"
done
if (( ARRAY_TASKS != TOTAL_SEEDS_PER_MODE * 2 )); then
  echo "ARRAY_TASKS must equal TOTAL_SEEDS_PER_MODE * 2." >&2
  echo "Got ARRAY_TASKS=${ARRAY_TASKS}, TOTAL_SEEDS_PER_MODE=${TOTAL_SEEDS_PER_MODE}." >&2
  exit 2
fi
if ! [[ "${DE_RELTOL}" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]; then
  echo "DE_RELTOL must be a positive numeric value, got: ${DE_RELTOL}" >&2
  exit 2
fi
if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found; run this script on the HPC login node or use DRY_RUN=TRUE." >&2
  exit 1
fi

RUN_DIR_INVIVO="${OUT_ROOT}/${RUN_PREFIX_INVIVO}"
RUN_DIR_INVITRO="${OUT_ROOT}/${RUN_PREFIX_INVITRO}"
mkdir -p "${RUN_DIR_INVIVO}" "${RUN_DIR_INVITRO}"

EXPORTS="ALL"
EXPORTS+=",PROJECT_ROOT=${PROJECT_ROOT}"
EXPORTS+=",RUNNER_SCRIPT=${RUNNER_SCRIPT}"
EXPORTS+=",CONFIG_PATH=${CONFIG_PATH}"
EXPORTS+=",OUT_ROOT=${OUT_ROOT}"
EXPORTS+=",RUN_PREFIX_INVIVO=${RUN_PREFIX_INVIVO}"
EXPORTS+=",RUN_PREFIX_INVITRO=${RUN_PREFIX_INVITRO}"
EXPORTS+=",TOTAL_SEEDS_PER_MODE=${TOTAL_SEEDS_PER_MODE}"
EXPORTS+=",ARRAY_TASKS=${ARRAY_TASKS}"
EXPORTS+=",N_CORES=${N_CORES}"
EXPORTS+=",AUTO_VIZ=${AUTO_VIZ}"
EXPORTS+=",GLUCOSE=${GLUCOSE}"
EXPORTS+=",R_MODULE=${R_MODULE}"
EXPORTS+=",PARAMETER_TABLE=${PARAMETER_TABLE}"
EXPORTS+=",FIT_OBJECTS_DIR=${FIT_OBJECTS_DIR}"
EXPORTS+=",FLOW_DENSITY_PATH=${FLOW_DENSITY_PATH}"
EXPORTS+=",ITERMAX=${ITERMAX}"
EXPORTS+=",DE_RELTOL=${DE_RELTOL}"
EXPORTS+=",DE_STEPTOL=${DE_STEPTOL}"
EXPORTS+=",NP=${NP}"

array_cmd=(
  sbatch
  --parsable
  --job-name=o2g_iv_mix
  --qos=xxlarge
  --time=12:00:00
  "--cpus-per-task=${N_CORES}"
  "--mem=${MEM}"
  "--array=1-${ARRAY_TASKS}"
  "--output=${OUT_ROOT}/o2g_invivo_invitro_mix_%A_%a.out"
  "--error=${OUT_ROOT}/o2g_invivo_invitro_mix_%A_%a.err"
  "--export=${EXPORTS}"
  "${SUB_SCRIPT}"
)

echo "Submitting interleaved invivo/invitro array"
echo "  project_root: ${PROJECT_ROOT}"
echo "  out_root: ${OUT_ROOT}"
echo "  config_path: ${CONFIG_PATH}"
echo "  invivo -> ${RUN_DIR_INVIVO}"
echo "  invitro -> ${RUN_DIR_INVITRO}"
echo "  array: 1-${ARRAY_TASKS} (odd=invivo, even=invitro)"
echo "  resources: qos=xxlarge, time=12:00:00, cpus=${N_CORES}, mem=${MEM}"
echo "  force_extra_results: ${FORCE_EXTRA_RESULTS}"
echo "  dry_run: ${DRY_RUN}"

if truthy "${DRY_RUN}"; then
  print_command "Submit interleaved array" "${array_cmd[@]}"
  MIXED_JOB_ID="DRYRUN_MIXED_JOB"
else
  MIXED_JOB_ID="$("${array_cmd[@]}")"
  echo "Submitted interleaved invivo/invitro array job: ${MIXED_JOB_ID}"
fi

submit_postprocess() {
  local label="$1"
  local run_dir="$2"
  local job_name="o2g_extra_${label}"
  local exports="ALL,PROJECT_ROOT=${PROJECT_ROOT},RUN_DIR=${run_dir},EXTRA_RESULTS_SCRIPT=${EXTRA_RESULTS_SCRIPT},R_MODULE=${R_MODULE},FORCE_EXTRA_RESULTS=${FORCE_EXTRA_RESULTS}"
  local cmd=(
    sbatch
    --parsable
    "--job-name=${job_name}"
    --qos=small
    --time=4:00:00
    --cpus-per-task=1
    --mem=8G
    "--output=${OUT_ROOT}/${job_name}_%A.out"
    "--error=${OUT_ROOT}/${job_name}_%A.err"
    "--export=${exports}"
    "--dependency=afterok:${MIXED_JOB_ID}"
    "${POSTPROCESS_SCRIPT}"
  )
  if truthy "${DRY_RUN}"; then
    print_command "Submit ${label} extra_results" "${cmd[@]}"
  else
    local post_job_id
    post_job_id="$("${cmd[@]}")"
    echo "Submitted ${label} extra_results job: ${post_job_id}"
  fi
}

submit_postprocess "invivo" "${RUN_DIR_INVIVO}"
submit_postprocess "invitro" "${RUN_DIR_INVITRO}"
