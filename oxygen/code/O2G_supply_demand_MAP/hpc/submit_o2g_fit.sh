#!/bin/bash
# Unified O2G HPC submitter for in vivo, in vitro, and joint fitting.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash submit_o2g_fit.sh --fitting_mode=invivo [options]
  bash submit_o2g_fit.sh --fitting_mode=invitro [options]
  bash submit_o2g_fit.sh --fitting_mode=joint --joint_fitting_mode=JOINT [options]
  bash submit_o2g_fit.sh --fitting_mode=joint --joint_fitting_mode=SINGLE \
    --invivo_run_dir=/path/to/invivo_run --invitro_run_dir=/path/to/invitro_run [options]

Required modes:
  --fitting_mode=invivo|invitro|joint
  --joint_fitting_mode=OFF|JOINT|SINGLE

Joint mode behavior:
  OFF    Do not submit joint fitting. This is forced when fitting_mode is not joint.
  JOINT  Submit in vivo and in vitro single fits first, run extra_results for each,
         derive anchors, then submit joint fitting. invivo_run_dir and
         invitro_run_dir must be NULL/empty in this mode.
  SINGLE Only submit joint fitting from existing in vivo and in vitro result dirs.
         invivo_run_dir and invitro_run_dir must be specified.

Common options:
  --project_root=/share/lab_crd/lab_crd/taoli/Project/Rescomposite
  --config_path=/path/to/O2G_supply_demand.yaml
  --out_root=/path/to/oxygen/results
  --r_module=R/4.4
  --dry_run=TRUE|FALSE
  --force_extra_results=TRUE|FALSE

Single-fit options:
  --invivo_run_prefix=name
  --invitro_run_prefix=name
  --invivo_total_seeds=500 --invivo_array_tasks=500 --invivo_seeds_per_task=1
  --invitro_total_seeds=500 --invitro_array_tasks=500 --invitro_seeds_per_task=1
  --invivo_qos=xlarge --invivo_time_limit=12:00:00
  --invitro_qos=xxlarge --invitro_time_limit=12:00:00

Joint options:
  --joint_run_prefix=name
  --anchor_mode=auto|manual
  --manual_o2_grid=0,0.1,0.2,0.5,1,2
  --manual_n_grid=44,66,88
  --joint_total_seeds=500 --joint_array_tasks=500 --joint_seeds_per_task=1
  --joint_qos=xxlarge --joint_time_limit=12:00:00
  --prep_qos=small --prep_time_limit=2:00:00
  --postprocess_qos=small --postprocess_time_limit=4:00:00
  Manual grids may use commas or semicolons; commas are converted before Slurm export.

After each fitting finishes:
  A dependent postprocess job runs extra_results.R. If extra_results already
  exists, it is skipped unless force_extra_results=TRUE.
EOF
}

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

is_null_value() {
  local val
  val="$(echo "${1:-}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')"
  [[ -z "${val}" || "${val}" == "null" || "${val}" == "none" || "${val}" == "na" ]]
}

normalize_fitting_mode() {
  local val
  val="$(echo "${1:-}" | tr '[:upper:]' '[:lower:]')"
  val="${val// /}"
  val="${val//-/}"
  val="${val//_/}"
  case "${val}" in
    invivo) echo "invivo" ;;
    invitro) echo "invitro" ;;
    joint) echo "joint" ;;
    *) echo "" ;;
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

check_seed_plan() {
  local label="$1"
  local total="$2"
  local tasks="$3"
  local per_task="$4"
  if (( tasks * per_task != total )); then
    echo "Mismatch for ${label}: ARRAY_TASKS * SEEDS_PER_TASK must equal TOTAL_SEEDS." >&2
    echo "Got total=${total}, array_tasks=${tasks}, seeds_per_task=${per_task}." >&2
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

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h)
        usage
        exit 0
        ;;
      --fitting_mode=*) FITTING_MODE="${arg#*=}" ;;
      --joint_fitting_mode=*) JOINT_FITTING_MODE="${arg#*=}" ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --config_path=*) CONFIG_PATH="${arg#*=}" ;;
      --out_root=*) OUT_ROOT="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      --force_extra_results=*) FORCE_EXTRA_RESULTS="${arg#*=}" ;;
      --invivo_run_dir=*) INVIVO_RUN_DIR="${arg#*=}" ;;
      --invitro_run_dir=*) INVITRO_RUN_DIR="${arg#*=}" ;;
      --invivo_run_prefix=*) INVIVO_RUN_PREFIX="${arg#*=}" ;;
      --invitro_run_prefix=*) INVITRO_RUN_PREFIX="${arg#*=}" ;;
      --joint_run_prefix=*) JOINT_RUN_PREFIX="${arg#*=}" ;;
      --invivo_total_seeds=*) INVIVO_TOTAL_SEEDS="${arg#*=}" ;;
      --invivo_array_tasks=*) INVIVO_ARRAY_TASKS="${arg#*=}" ;;
      --invivo_seeds_per_task=*) INVIVO_SEEDS_PER_TASK="${arg#*=}" ;;
      --invitro_total_seeds=*) INVITRO_TOTAL_SEEDS="${arg#*=}" ;;
      --invitro_array_tasks=*) INVITRO_ARRAY_TASKS="${arg#*=}" ;;
      --invitro_seeds_per_task=*) INVITRO_SEEDS_PER_TASK="${arg#*=}" ;;
      --joint_total_seeds=*) JOINT_TOTAL_SEEDS="${arg#*=}" ;;
      --joint_array_tasks=*) JOINT_ARRAY_TASKS="${arg#*=}" ;;
      --joint_seeds_per_task=*) JOINT_SEEDS_PER_TASK="${arg#*=}" ;;
      --n_cores=*) N_CORES="${arg#*=}" ;;
      --invivo_n_cores=*) INVIVO_N_CORES="${arg#*=}" ;;
      --invitro_n_cores=*) INVITRO_N_CORES="${arg#*=}" ;;
      --joint_n_cores=*) JOINT_N_CORES="${arg#*=}" ;;
      --mem=*) MEM="${arg#*=}" ;;
      --invivo_mem=*) INVIVO_MEM="${arg#*=}" ;;
      --invitro_mem=*) INVITRO_MEM="${arg#*=}" ;;
      --joint_mem=*) JOINT_MEM="${arg#*=}" ;;
      --invivo_qos=*) INVIVO_QOS="${arg#*=}" ;;
      --invitro_qos=*) INVITRO_QOS="${arg#*=}" ;;
      --joint_qos=*) JOINT_QOS="${arg#*=}" ;;
      --prep_qos=*) PREP_QOS="${arg#*=}" ;;
      --postprocess_qos=*) POSTPROCESS_QOS="${arg#*=}" ;;
      --invivo_time_limit=*) INVIVO_TIME_LIMIT="${arg#*=}" ;;
      --invitro_time_limit=*) INVITRO_TIME_LIMIT="${arg#*=}" ;;
      --joint_time_limit=*) JOINT_TIME_LIMIT="${arg#*=}" ;;
      --prep_time_limit=*) PREP_TIME_LIMIT="${arg#*=}" ;;
      --postprocess_time_limit=*) POSTPROCESS_TIME_LIMIT="${arg#*=}" ;;
      --prep_mem=*) PREP_MEM="${arg#*=}" ;;
      --postprocess_mem=*) POSTPROCESS_MEM="${arg#*=}" ;;
      --anchor_mode=*) ANCHOR_MODE="${arg#*=}" ;;
      --manual_o2_grid=*) MANUAL_O2_GRID="${arg#*=}" ;;
      --manual_n_grid=*) MANUAL_N_GRID="${arg#*=}" ;;
      --top_n=*) TOP_N="${arg#*=}" ;;
      --parameter_table=*) PARAMETER_TABLE="${arg#*=}" ;;
      --fit_objects_dir=*) FIT_OBJECTS_DIR="${arg#*=}" ;;
      --flow_density_path=*) FLOW_DENSITY_PATH="${arg#*=}" ;;
      --itermax=*) ITERMAX="${arg#*=}" ;;
      --de_reltol=*) DE_RELTOL="${arg#*=}" ;;
      --de_steptol=*) DE_STEPTOL="${arg#*=}" ;;
      --np=*|--NP=*) NP="${arg#*=}" ;;
      --auto_viz=*) AUTO_VIZ="${arg#*=}" ;;
      --glucose=*) GLUCOSE="${arg#*=}" ;;
      *)
        echo "Unknown argument: ${arg}" >&2
        usage >&2
        exit 2
        ;;
    esac
  done
}

ensure_sbatch() {
  if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch not found; run this launcher on the HPC login node or use DRY_RUN=TRUE." >&2
    exit 1
  fi
}

submit_extra_results_job() {
  local label="$1"
  local run_dir="$2"
  local dependency="$3"
  local export_arg="ALL"
  export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
  export_arg+=",RUN_DIR=${run_dir}"
  export_arg+=",EXTRA_RESULTS_SCRIPT=${EXTRA_RESULTS_SCRIPT}"
  export_arg+=",R_MODULE=${R_MODULE}"
  export_arg+=",FORCE_EXTRA_RESULTS=${FORCE_EXTRA_RESULTS}"

  local cmd=(
    sbatch
    --parsable
    "--job-name=${label}_extra"
    "--qos=${POSTPROCESS_QOS}"
    "--time=${POSTPROCESS_TIME_LIMIT}"
    --cpus-per-task=1
    "--mem=${POSTPROCESS_MEM}"
    "--output=${OUT_ROOT}/${label}_extra_%A.out"
    "--error=${OUT_ROOT}/${label}_extra_%A.err"
    "--export=${export_arg}"
  )
  if [[ -n "${dependency}" ]]; then
    cmd+=("--dependency=afterok:${dependency}")
  fi
  cmd+=("${POSTPROCESS_SCRIPT}")

  if truthy "${DRY_RUN}"; then
    print_command "Submit ${label} extra_results" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_${label}_EXTRA_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted ${label} extra_results job: ${LAST_JOB_ID}"
  fi
}

submit_invivo_array() {
  local run_dir="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
  mkdir -p "${run_dir}"
  local export_arg="ALL"
  export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
  export_arg+=",RUNNER_SCRIPT=${INVIVO_RUNNER_SCRIPT}"
  export_arg+=",CONFIG_PATH=${CONFIG_PATH}"
  export_arg+=",OUT_ROOT=${OUT_ROOT}"
  export_arg+=",RUN_PREFIX=${INVIVO_RUN_PREFIX}"
  export_arg+=",TOTAL_SEEDS=${INVIVO_TOTAL_SEEDS}"
  export_arg+=",ARRAY_TASKS=${INVIVO_ARRAY_TASKS}"
  export_arg+=",SEEDS_PER_TASK=${INVIVO_SEEDS_PER_TASK}"
  export_arg+=",N_CORES=${INVIVO_N_CORES}"
  export_arg+=",AUTO_VIZ=${AUTO_VIZ}"
  export_arg+=",GLUCOSE=${GLUCOSE}"
  export_arg+=",R_MODULE=${R_MODULE}"
  local cmd=(
    sbatch
    --parsable
    --job-name=o2g_ivv_B
    "--qos=${INVIVO_QOS}"
    "--time=${INVIVO_TIME_LIMIT}"
    "--cpus-per-task=${INVIVO_N_CORES}"
    "--mem=${INVIVO_MEM}"
    "--array=1-${INVIVO_ARRAY_TASKS}"
    "--output=${OUT_ROOT}/o2g_invivo_%A_%a.out"
    "--error=${OUT_ROOT}/o2g_invivo_%A_%a.err"
    "--export=${export_arg}"
    "${INVIVO_SUB_SCRIPT}"
  )
  if truthy "${DRY_RUN}"; then
    print_command "Submit in vivo" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_INVIVO_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted in vivo array job: ${LAST_JOB_ID}"
  fi
}

submit_invitro_array() {
  local run_dir="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
  mkdir -p "${run_dir}"
  local export_arg="ALL"
  export_arg+=",PROJECT_ROOT=${PROJECT_ROOT}"
  export_arg+=",RUNNER_SCRIPT=${INVITRO_RUNNER_SCRIPT}"
  export_arg+=",CONFIG_PATH=${CONFIG_PATH}"
  export_arg+=",OUT_ROOT=${OUT_ROOT}"
  export_arg+=",RUN_PREFIX=${INVITRO_RUN_PREFIX}"
  export_arg+=",TOTAL_SEEDS=${INVITRO_TOTAL_SEEDS}"
  export_arg+=",ARRAY_TASKS=${INVITRO_ARRAY_TASKS}"
  export_arg+=",SEEDS_PER_TASK=${INVITRO_SEEDS_PER_TASK}"
  export_arg+=",N_CORES=${INVITRO_N_CORES}"
  export_arg+=",R_MODULE=${R_MODULE}"
  export_arg+=",PARAMETER_TABLE=${PARAMETER_TABLE}"
  export_arg+=",FIT_OBJECTS_DIR=${FIT_OBJECTS_DIR}"
  export_arg+=",FLOW_DENSITY_PATH=${FLOW_DENSITY_PATH}"
  export_arg+=",ITERMAX=${ITERMAX}"
  export_arg+=",DE_RELTOL=${DE_RELTOL}"
  export_arg+=",DE_STEPTOL=${DE_STEPTOL}"
  export_arg+=",NP=${NP}"
  export_arg+=",AUTO_VIZ=${AUTO_VIZ}"
  local cmd=(
    sbatch
    --parsable
    --job-name=o2g_ivt_B
    "--qos=${INVITRO_QOS}"
    "--time=${INVITRO_TIME_LIMIT}"
    "--cpus-per-task=${INVITRO_N_CORES}"
    "--mem=${INVITRO_MEM}"
    "--array=1-${INVITRO_ARRAY_TASKS}"
    "--output=${OUT_ROOT}/o2g_invitro_%A_%a.out"
    "--error=${OUT_ROOT}/o2g_invitro_%A_%a.err"
    "--export=${export_arg}"
    "${INVITRO_SUB_SCRIPT}"
  )
  if truthy "${DRY_RUN}"; then
    print_command "Submit in vitro" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_INVITRO_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted in vitro array job: ${LAST_JOB_ID}"
  fi
}

submit_joint_prep() {
  local dependency="$1"
  local manual_o2_export="${MANUAL_O2_GRID//,/;}"
  local manual_n_export="${MANUAL_N_GRID//,/;}"
  local prep_export="ALL"
  prep_export+=",PROJECT_ROOT=${PROJECT_ROOT}"
  prep_export+=",CONFIG_PATH=${CONFIG_PATH}"
  prep_export+=",OUT_ROOT=${OUT_ROOT}"
  prep_export+=",INVIVO_RUN_DIR=${INVIVO_RUN_DIR}"
  prep_export+=",INVITRO_RUN_DIR=${INVITRO_RUN_DIR}"
  prep_export+=",JOINT_FITTING_MODE=${JOINT_FITTING_MODE}"
  prep_export+=",JOINT_RUN_PREFIX=${JOINT_RUN_PREFIX}"
  prep_export+=",ANCHOR_MODE=${ANCHOR_MODE}"
  prep_export+=",MANUAL_O2_GRID=${manual_o2_export}"
  prep_export+=",MANUAL_N_GRID=${manual_n_export}"
  prep_export+=",TOP_N=${TOP_N}"
  prep_export+=",JOINT_TOTAL_SEEDS=${JOINT_TOTAL_SEEDS}"
  prep_export+=",JOINT_ARRAY_TASKS=${JOINT_ARRAY_TASKS}"
  prep_export+=",JOINT_SEEDS_PER_TASK=${JOINT_SEEDS_PER_TASK}"
  prep_export+=",JOINT_N_CORES=${JOINT_N_CORES}"
  prep_export+=",JOINT_MEM=${JOINT_MEM}"
  prep_export+=",JOINT_QOS=${JOINT_QOS}"
  prep_export+=",JOINT_TIME_LIMIT=${JOINT_TIME_LIMIT}"
  prep_export+=",R_MODULE=${R_MODULE}"
  prep_export+=",AUTO_VIZ=${AUTO_VIZ}"
  prep_export+=",GLUCOSE=${GLUCOSE}"
  prep_export+=",FORCE_EXTRA_RESULTS=${FORCE_EXTRA_RESULTS}"
  prep_export+=",POSTPROCESS_SCRIPT=${POSTPROCESS_SCRIPT}"
  prep_export+=",DERIVE_ANCHOR_SCRIPT=${DERIVE_ANCHOR_SCRIPT}"
  prep_export+=",EXTRA_RESULTS_SCRIPT=${EXTRA_RESULTS_SCRIPT}"
  prep_export+=",JOINT_SUB_SCRIPT=${JOINT_SUB_SCRIPT}"
  prep_export+=",JOINT_RUNNER_SCRIPT=${JOINT_RUNNER_SCRIPT}"
  prep_export+=",POSTPROCESS_QOS=${POSTPROCESS_QOS}"
  prep_export+=",POSTPROCESS_TIME_LIMIT=${POSTPROCESS_TIME_LIMIT}"
  prep_export+=",POSTPROCESS_MEM=${POSTPROCESS_MEM}"

  local cmd=(
    sbatch
    --parsable
    --job-name=o2g_joint_prep
    "--qos=${PREP_QOS}"
    "--time=${PREP_TIME_LIMIT}"
    --cpus-per-task=1
    "--mem=${PREP_MEM}"
    "--output=${OUT_ROOT}/o2g_joint_prep_%A.out"
    "--error=${OUT_ROOT}/o2g_joint_prep_%A.err"
    "--export=${prep_export}"
  )
  if [[ -n "${dependency}" ]]; then
    cmd+=("--dependency=afterok:${dependency}")
  fi
  cmd+=("${JOINT_PREP_SCRIPT}")

  if truthy "${DRY_RUN}"; then
    print_command "Submit joint prep" "${cmd[@]}"
    LAST_JOB_ID="DRYRUN_JOINT_PREP_JOB"
  else
    LAST_JOB_ID="$("${cmd[@]}")"
    echo "Submitted joint preparation job: ${LAST_JOB_ID}"
  fi
}

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/Rescomposite"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_INVIVO_RUN_PREFIX="fit_invivo_O2G_buffering_500seed"
DEFAULT_INVITRO_RUN_PREFIX="fit_invitro_O2G_buffering_500seed"
DEFAULT_JOINT_RUN_PREFIX="fit_joint_O2G_buffering_500seed"
DEFAULT_TOTAL_SEEDS="500"
DEFAULT_SEEDS_PER_TASK="1"
DEFAULT_N_CORES="22"
DEFAULT_MEM="16G"
DEFAULT_AUTO_VIZ="TRUE"
DEFAULT_GLUCOSE="TRUE"
DEFAULT_INVIVO_QOS="xlarge"
DEFAULT_INVIVO_TIME_LIMIT="12:00:00"
DEFAULT_INVITRO_QOS="xxlarge"
DEFAULT_INVITRO_TIME_LIMIT="12:00:00"
DEFAULT_JOINT_QOS="xxlarge"
DEFAULT_JOINT_TIME_LIMIT="12:00:00"
DEFAULT_PREP_QOS="small"
DEFAULT_PREP_TIME_LIMIT="2:00:00"
DEFAULT_PREP_MEM="8G"
DEFAULT_POSTPROCESS_QOS="small"
DEFAULT_POSTPROCESS_TIME_LIMIT="4:00:00"
DEFAULT_POSTPROCESS_MEM="8G"
DEFAULT_ITERMAX="500"
DEFAULT_DE_RELTOL="1e-4"
DEFAULT_DE_STEPTOL="25"
DEFAULT_NP="80"
DEFAULT_ANCHOR_MODE="auto"
DEFAULT_TOP_N="3"
DEFAULT_FORCE_EXTRA_RESULTS="FALSE"
DEFAULT_DRY_RUN="FALSE"

FITTING_MODE="${FITTING_MODE:-}"
JOINT_FITTING_MODE="${JOINT_FITTING_MODE:-OFF}"
PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
CONFIG_PATH="${CONFIG_PATH:-}"
OUT_ROOT="${OUT_ROOT:-}"
INVIVO_RUN_PREFIX="${INVIVO_RUN_PREFIX:-${DEFAULT_INVIVO_RUN_PREFIX}}"
INVITRO_RUN_PREFIX="${INVITRO_RUN_PREFIX:-${DEFAULT_INVITRO_RUN_PREFIX}}"
JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX:-${DEFAULT_JOINT_RUN_PREFIX}}"
INVIVO_TOTAL_SEEDS="${INVIVO_TOTAL_SEEDS:-${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}}"
INVITRO_TOTAL_SEEDS="${INVITRO_TOTAL_SEEDS:-${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}}"
JOINT_TOTAL_SEEDS="${JOINT_TOTAL_SEEDS:-${TOTAL_SEEDS:-${DEFAULT_TOTAL_SEEDS}}}"
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
INVIVO_QOS="${INVIVO_QOS:-${DEFAULT_INVIVO_QOS}}"
INVIVO_TIME_LIMIT="${INVIVO_TIME_LIMIT:-${DEFAULT_INVIVO_TIME_LIMIT}}"
INVITRO_QOS="${INVITRO_QOS:-${DEFAULT_INVITRO_QOS}}"
INVITRO_TIME_LIMIT="${INVITRO_TIME_LIMIT:-${DEFAULT_INVITRO_TIME_LIMIT}}"
JOINT_QOS="${JOINT_QOS:-${DEFAULT_JOINT_QOS}}"
JOINT_TIME_LIMIT="${JOINT_TIME_LIMIT:-${DEFAULT_JOINT_TIME_LIMIT}}"
PREP_QOS="${PREP_QOS:-${DEFAULT_PREP_QOS}}"
PREP_TIME_LIMIT="${PREP_TIME_LIMIT:-${DEFAULT_PREP_TIME_LIMIT}}"
PREP_MEM="${PREP_MEM:-${DEFAULT_PREP_MEM}}"
POSTPROCESS_QOS="${POSTPROCESS_QOS:-${DEFAULT_POSTPROCESS_QOS}}"
POSTPROCESS_TIME_LIMIT="${POSTPROCESS_TIME_LIMIT:-${DEFAULT_POSTPROCESS_TIME_LIMIT}}"
POSTPROCESS_MEM="${POSTPROCESS_MEM:-${DEFAULT_POSTPROCESS_MEM}}"
PARAMETER_TABLE="${PARAMETER_TABLE:-}"
FIT_OBJECTS_DIR="${FIT_OBJECTS_DIR:-}"
FLOW_DENSITY_PATH="${FLOW_DENSITY_PATH:-}"
ITERMAX="${ITERMAX:-${DEFAULT_ITERMAX}}"
DE_RELTOL="${DE_RELTOL:-${DEFAULT_DE_RELTOL}}"
DE_STEPTOL="${DE_STEPTOL:-${DEFAULT_DE_STEPTOL}}"
NP="${NP:-${DEFAULT_NP}}"
AUTO_VIZ="${AUTO_VIZ:-${DEFAULT_AUTO_VIZ}}"
GLUCOSE="${GLUCOSE:-${DEFAULT_GLUCOSE}}"
ANCHOR_MODE="${ANCHOR_MODE:-${DEFAULT_ANCHOR_MODE}}"
MANUAL_O2_GRID="${MANUAL_O2_GRID:-}"
MANUAL_N_GRID="${MANUAL_N_GRID:-}"
TOP_N="${TOP_N:-${DEFAULT_TOP_N}}"
INVIVO_RUN_DIR="${INVIVO_RUN_DIR:-}"
INVITRO_RUN_DIR="${INVITRO_RUN_DIR:-}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-${DEFAULT_FORCE_EXTRA_RESULTS}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"

parse_args "$@"

FITTING_MODE="$(normalize_fitting_mode "${FITTING_MODE}")"
if [[ -z "${FITTING_MODE}" ]]; then
  echo "--fitting_mode must be one of invivo, invitro, or joint." >&2
  usage >&2
  exit 2
fi

JOINT_FITTING_MODE="$(echo "${JOINT_FITTING_MODE}" | tr '[:lower:]' '[:upper:]')"
if [[ "${FITTING_MODE}" != "joint" ]]; then
  JOINT_FITTING_MODE="OFF"
fi
case "${JOINT_FITTING_MODE}" in
  OFF|JOINT|SINGLE) ;;
  *)
    echo "--joint_fitting_mode must be OFF, JOINT, or SINGLE." >&2
    exit 2
    ;;
esac

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

HPC_DIR="${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/hpc"
INVIVO_SUB_SCRIPT="${INVIVO_SUB_SCRIPT:-${HPC_DIR}/submit_fit_seed_array_buffering.sub}"
INVITRO_SUB_SCRIPT="${INVITRO_SUB_SCRIPT:-${HPC_DIR}/submit_fit_seed_array_invitro_buffering.sub}"
JOINT_SUB_SCRIPT="${JOINT_SUB_SCRIPT:-${HPC_DIR}/submit_fit_seed_array_joint_buffering.sub}"
JOINT_PREP_SCRIPT="${JOINT_PREP_SCRIPT:-${HPC_DIR}/postprocess_and_submit_joint_buffering.sh}"
POSTPROCESS_SCRIPT="${POSTPROCESS_SCRIPT:-${HPC_DIR}/postprocess_extra_results.sh}"
DERIVE_ANCHOR_SCRIPT="${DERIVE_ANCHOR_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/derive_joint_anchor_from_single_fits.R}"
EXTRA_RESULTS_SCRIPT="${EXTRA_RESULTS_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R}"
INVIVO_RUNNER_SCRIPT="${INVIVO_RUNNER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh}"
INVITRO_RUNNER_SCRIPT="${INVITRO_RUNNER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh}"
JOINT_RUNNER_SCRIPT="${JOINT_RUNNER_SCRIPT:-${PROJECT_ROOT}/oxygen/code/O2G_supply_demand_MAP/runner/run_fit_joint_model_O2G_supply_demand_MAP.sh}"

if [[ -z "${PARAMETER_TABLE}" ]]; then
  PARAMETER_TABLE="${PROJECT_ROOT}/oxygen/data/O2G_supply_demand/parameter_table_invitro_buffering.csv"
fi
if [[ -z "${FIT_OBJECTS_DIR}" ]]; then
  FIT_OBJECTS_DIR="${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects"
fi
if [[ -z "${FLOW_DENSITY_PATH}" ]]; then
  FLOW_DENSITY_PATH="${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv"
fi

for path in "${CONFIG_PATH}" "${INVIVO_SUB_SCRIPT}" "${INVITRO_SUB_SCRIPT}" "${JOINT_SUB_SCRIPT}" \
            "${JOINT_PREP_SCRIPT}" "${POSTPROCESS_SCRIPT}" "${DERIVE_ANCHOR_SCRIPT}" "${EXTRA_RESULTS_SCRIPT}" \
            "${INVIVO_RUNNER_SCRIPT}" "${INVITRO_RUNNER_SCRIPT}" "${JOINT_RUNNER_SCRIPT}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
done

for name in INVIVO_TOTAL_SEEDS INVIVO_ARRAY_TASKS INVIVO_SEEDS_PER_TASK \
            INVITRO_TOTAL_SEEDS INVITRO_ARRAY_TASKS INVITRO_SEEDS_PER_TASK \
            JOINT_TOTAL_SEEDS JOINT_ARRAY_TASKS JOINT_SEEDS_PER_TASK \
            INVIVO_N_CORES INVITRO_N_CORES JOINT_N_CORES ITERMAX DE_STEPTOL NP TOP_N; do
  require_positive_int "${name}" "${!name}"
done
check_seed_plan INVIVO "${INVIVO_TOTAL_SEEDS}" "${INVIVO_ARRAY_TASKS}" "${INVIVO_SEEDS_PER_TASK}"
check_seed_plan INVITRO "${INVITRO_TOTAL_SEEDS}" "${INVITRO_ARRAY_TASKS}" "${INVITRO_SEEDS_PER_TASK}"
check_seed_plan JOINT "${JOINT_TOTAL_SEEDS}" "${JOINT_ARRAY_TASKS}" "${JOINT_SEEDS_PER_TASK}"

ensure_sbatch

echo "O2G submitter"
echo "  fitting_mode: ${FITTING_MODE}"
echo "  joint_fitting_mode: ${JOINT_FITTING_MODE}"
echo "  project_root: ${PROJECT_ROOT}"
echo "  out_root: ${OUT_ROOT}"
echo "  r_module: ${R_MODULE}"
echo "  invivo resources: qos=${INVIVO_QOS}, time=${INVIVO_TIME_LIMIT}, mem=${INVIVO_MEM}, cpus=${INVIVO_N_CORES}"
echo "  invitro resources: qos=${INVITRO_QOS}, time=${INVITRO_TIME_LIMIT}, mem=${INVITRO_MEM}, cpus=${INVITRO_N_CORES}"
echo "  joint resources: qos=${JOINT_QOS}, time=${JOINT_TIME_LIMIT}, mem=${JOINT_MEM}, cpus=${JOINT_N_CORES}"
echo "  prep resources: qos=${PREP_QOS}, time=${PREP_TIME_LIMIT}, mem=${PREP_MEM}"
echo "  postprocess resources: qos=${POSTPROCESS_QOS}, time=${POSTPROCESS_TIME_LIMIT}, mem=${POSTPROCESS_MEM}"

case "${FITTING_MODE}" in
  invivo)
    INVIVO_RUN_DIR="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
    submit_invivo_array
    INVIVO_JOB_ID="${LAST_JOB_ID}"
    submit_extra_results_job "o2g_invivo" "${INVIVO_RUN_DIR}" "${INVIVO_JOB_ID}"
    ;;
  invitro)
    INVITRO_RUN_DIR="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
    submit_invitro_array
    INVITRO_JOB_ID="${LAST_JOB_ID}"
    submit_extra_results_job "o2g_invitro" "${INVITRO_RUN_DIR}" "${INVITRO_JOB_ID}"
    ;;
  joint)
    case "${JOINT_FITTING_MODE}" in
      OFF)
        echo "joint_fitting_mode=OFF; no fitting submitted."
        ;;
      JOINT)
        if ! is_null_value "${INVIVO_RUN_DIR}" || ! is_null_value "${INVITRO_RUN_DIR}"; then
          echo "joint_fitting_mode=JOINT requires invivo_run_dir and invitro_run_dir to be NULL/empty." >&2
          exit 2
        fi
        INVIVO_RUN_DIR="${OUT_ROOT}/${INVIVO_RUN_PREFIX}"
        INVITRO_RUN_DIR="${OUT_ROOT}/${INVITRO_RUN_PREFIX}"
        submit_invivo_array
        INVIVO_JOB_ID="${LAST_JOB_ID}"
        submit_invitro_array
        INVITRO_JOB_ID="${LAST_JOB_ID}"
        submit_extra_results_job "o2g_invivo" "${INVIVO_RUN_DIR}" "${INVIVO_JOB_ID}"
        INVIVO_EXTRA_JOB_ID="${LAST_JOB_ID}"
        submit_extra_results_job "o2g_invitro" "${INVITRO_RUN_DIR}" "${INVITRO_JOB_ID}"
        INVITRO_EXTRA_JOB_ID="${LAST_JOB_ID}"
        submit_joint_prep "${INVIVO_EXTRA_JOB_ID}:${INVITRO_EXTRA_JOB_ID}"
        ;;
      SINGLE)
        if is_null_value "${INVIVO_RUN_DIR}" || is_null_value "${INVITRO_RUN_DIR}"; then
          echo "joint_fitting_mode=SINGLE requires --invivo_run_dir and --invitro_run_dir." >&2
          exit 2
        fi
        if [[ ! -d "${INVIVO_RUN_DIR}" || ! -d "${INVITRO_RUN_DIR}" ]]; then
          echo "joint_fitting_mode=SINGLE requires existing result directories." >&2
          echo "  invivo_run_dir: ${INVIVO_RUN_DIR}" >&2
          echo "  invitro_run_dir: ${INVITRO_RUN_DIR}" >&2
          exit 1
        fi
        INVIVO_RUN_DIR="$(cd "${INVIVO_RUN_DIR}" && pwd)"
        INVITRO_RUN_DIR="$(cd "${INVITRO_RUN_DIR}" && pwd)"
        submit_joint_prep ""
        ;;
    esac
    ;;
esac
