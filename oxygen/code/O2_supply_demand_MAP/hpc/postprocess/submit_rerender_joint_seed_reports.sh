#!/usr/bin/env bash
# Re-render figures and HTML reports for existing joint-fit seed directories.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash submit_rerender_joint_seed_reports.sh [options]

Primary options:
  --root=/share/.../oxygen/results/fit_joint_...
  --project_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling
  --child_prefix=fit_joint_tsne_vi_
  --data_dir=/share/.../data/InVivoData_Gemcitabine
  --r_module=R/4.4
  --qos=xxlarge
  --partition=NAME
  --account=NAME
  --time=04:00:00
  --mem=12G
  --cpus=1
  --report_dt=1
  --top_n=6
  --max_concurrent=N       Optional Slurm array concurrency cap. Default: none.
  --dry_run=TRUE|FALSE
  --help

Behavior:
  Scans ROOT/fit_joint_tsne_vi_*/seed*/fit_summary.tsv, keeps fit_joint seed
  directories, then submits one Slurm array task per seed. Each task runs:
    1. viz_invivo_model_O2_supply_demand_MAP_results.R -> seed/viz/invivo
    2. viz_invitro_model_O2_supply_demand_MAP_results.R -> seed/viz/invitro
    3. render_fit_report.R -> seed/report
EOF
}

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling"
ROOT=""
CHILD_PREFIX="fit_joint_tsne_vi_"
DATA_DIR=""
R_MODULE="R/4.4"
QOS="xxlarge"
PARTITION=""
ACCOUNT=""
TIME_LIMIT="04:00:00"
MEM="12G"
CPUS="1"
REPORT_DT="1"
TOP_N="6"
MAX_CONCURRENT=""
DRY_RUN="FALSE"
JOB_NAME="joint_seed_report"

for arg in "$@"; do
  case "${arg}" in
    --help|-h)
      usage
      exit 0
      ;;
    --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
    --root=*) ROOT="${arg#*=}" ;;
    --child_prefix=*) CHILD_PREFIX="${arg#*=}" ;;
    --data_dir=*) DATA_DIR="${arg#*=}" ;;
    --r_module=*) R_MODULE="${arg#*=}" ;;
    --qos=*) QOS="${arg#*=}" ;;
    --partition=*) PARTITION="${arg#*=}" ;;
    --account=*) ACCOUNT="${arg#*=}" ;;
    --time=*) TIME_LIMIT="${arg#*=}" ;;
    --mem=*) MEM="${arg#*=}" ;;
    --cpus=*) CPUS="${arg#*=}" ;;
    --report_dt=*) REPORT_DT="${arg#*=}" ;;
    --top_n=*) TOP_N="${arg#*=}" ;;
    --max_concurrent=*) MAX_CONCURRENT="${arg#*=}" ;;
    --dry_run=*) DRY_RUN="${arg#*=}" ;;
    --job_name=*) JOB_NAME="${arg#*=}" ;;
    *)
      echo "Unknown argument: ${arg}" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "${ROOT}" ]]; then
  ROOT="${PROJECT_ROOT}/oxygen/results/fit_joint_O2_buffering_100seed_tsne_invivo_subclusters_to_invitro_best_full_20260701_2f1772d"
fi
if [[ -z "${DATA_DIR}" ]]; then
  DATA_DIR="${PROJECT_ROOT}/data/InVivoData_Gemcitabine"
fi

if [[ ! -d "${PROJECT_ROOT}" ]]; then
  echo "PROJECT_ROOT does not exist: ${PROJECT_ROOT}" >&2
  exit 1
fi
if [[ ! -d "${ROOT}" ]]; then
  echo "ROOT does not exist: ${ROOT}" >&2
  exit 1
fi
if [[ ! -d "${DATA_DIR}" ]]; then
  echo "DATA_DIR does not exist: ${DATA_DIR}" >&2
  exit 1
fi

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="${ROOT}/rerender_joint_seed_reports_logs_${TIMESTAMP}"
mkdir -p "${LOG_DIR}"

TASK_FILE="${LOG_DIR}/joint_seed_fit_dirs.txt"
find "${ROOT}" \
  -mindepth 3 \
  -maxdepth 3 \
  -type f \
  -name fit_summary.tsv \
  -path "${ROOT}/${CHILD_PREFIX}*/seed*/fit_summary.tsv" \
  -print |
  while IFS= read -r summary_path; do
    if grep -q $'^fit_mode\tfit_joint$' "${summary_path}"; then
      dirname "${summary_path}"
    fi
  done |
  sort -V > "${TASK_FILE}"

N_TASKS="$(wc -l < "${TASK_FILE}" | tr -d '[:space:]')"
if [[ "${N_TASKS}" -lt 1 ]]; then
  echo "No joint seed fit directories found under ${ROOT}/${CHILD_PREFIX}*/seed*" >&2
  exit 1
fi

WORKER_SCRIPT="${LOG_DIR}/rerender_joint_seed_report_worker.sh"
cat > "${WORKER_SCRIPT}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

load_r_module() {
  if command -v ml >/dev/null 2>&1; then
    ml "${R_MODULE}"
  elif command -v module >/dev/null 2>&1; then
    module load "${R_MODULE}"
  fi
}

run_logged() {
  local label="$1"
  local log_path="$2"
  shift 2

  mkdir -p "$(dirname "${log_path}")"
  {
    echo "[$(date '+%F %T')] ${label} start"
    printf 'Command:'
    printf ' %q' "$@"
    printf '\n'
  } > "${log_path}"

  if "$@" >> "${log_path}" 2>&1; then
    echo "[$(date '+%F %T')] ${label} done" >> "${log_path}"
  else
    local status="$?"
    echo "[$(date '+%F %T')] ${label} failed with status ${status}" >> "${log_path}"
    tail -80 "${log_path}" >&2 || true
    return "${status}"
  fi
}

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID is not set." >&2
  exit 1
fi

FIT_DIR="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${TASK_FILE}")"
if [[ -z "${FIT_DIR}" || ! -d "${FIT_DIR}" ]]; then
  echo "Invalid FIT_DIR for task ${SLURM_ARRAY_TASK_ID}: ${FIT_DIR}" >&2
  exit 1
fi

load_r_module
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found after loading ${R_MODULE}." >&2
  exit 1
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

INVIVO_VIZ_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/vis/viz_invivo_model_O2_supply_demand_MAP_results.R"
INVITRO_VIZ_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/vis/viz_invitro_model_O2_supply_demand_MAP_results.R"
REPORT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/report/render_fit_report.R"

for script in "${INVIVO_VIZ_SCRIPT}" "${INVITRO_VIZ_SCRIPT}" "${REPORT_SCRIPT}"; do
  if [[ ! -f "${script}" ]]; then
    echo "Missing script: ${script}" >&2
    exit 1
  fi
done

echo "Task ${SLURM_ARRAY_TASK_ID}: ${FIT_DIR}"

run_logged "in vivo viz" "${FIT_DIR}/viz_status.log" \
  Rscript "${INVIVO_VIZ_SCRIPT}" \
    --fit_dir="${FIT_DIR}" \
    --out_dir="${FIT_DIR}/viz/invivo" \
    --data_dir="${DATA_DIR}" \
    --report_dt="${REPORT_DT}" \
    --top_n="${TOP_N}" \
    --n_cores=1

run_logged "in vitro viz" "${FIT_DIR}/invitro_viz_status.log" \
  Rscript "${INVITRO_VIZ_SCRIPT}" \
    --fit_dir="${FIT_DIR}" \
    --out_dir="${FIT_DIR}/viz/invitro"

run_logged "HTML report" "${FIT_DIR}/report_status.log" \
  Rscript "${REPORT_SCRIPT}" \
    --fit_dir="${FIT_DIR}" \
    --render_pdf=FALSE
EOF
chmod +x "${WORKER_SCRIPT}"

ARRAY_SPEC="1-${N_TASKS}"
if [[ -n "${MAX_CONCURRENT}" ]]; then
  ARRAY_SPEC="${ARRAY_SPEC}%${MAX_CONCURRENT}"
fi

SBATCH_ARGS=(
  --job-name="${JOB_NAME}"
  --ntasks=1
  --cpus-per-task="${CPUS}"
  --mem="${MEM}"
  --time="${TIME_LIMIT}"
  --array="${ARRAY_SPEC}"
  --output="${LOG_DIR}/%x_%A_%a.out"
  --error="${LOG_DIR}/%x_%A_%a.err"
  --export="ALL,PROJECT_ROOT=${PROJECT_ROOT},DATA_DIR=${DATA_DIR},TASK_FILE=${TASK_FILE},R_MODULE=${R_MODULE},REPORT_DT=${REPORT_DT},TOP_N=${TOP_N}"
)

if [[ -n "${QOS}" ]]; then
  SBATCH_ARGS+=(--qos="${QOS}")
fi
if [[ -n "${PARTITION}" ]]; then
  SBATCH_ARGS+=(--partition="${PARTITION}")
fi
if [[ -n "${ACCOUNT}" ]]; then
  SBATCH_ARGS+=(--account="${ACCOUNT}")
fi

echo "Root: ${ROOT}"
echo "Child prefix: ${CHILD_PREFIX}"
echo "Task file: ${TASK_FILE}"
echo "Tasks: ${N_TASKS}"
echo "Log dir: ${LOG_DIR}"
echo "Worker: ${WORKER_SCRIPT}"
echo "Array: ${ARRAY_SPEC}"

if truthy "${DRY_RUN}"; then
  echo "Dry run; not submitting."
  printf 'sbatch'
  printf ' %q' "${SBATCH_ARGS[@]}" "${WORKER_SCRIPT}"
  printf '\n'
  exit 0
fi

sbatch "${SBATCH_ARGS[@]}" "${WORKER_SCRIPT}"
