#!/bin/bash
# Submit fixed-O2 simulation tasks as a Slurm array.
#
# The task list is built by simulation/fix_o2_simulation.R, so the HPC array
# uses the same seed x O2 x ploidy x simulation_id logic as local batch mode.
# Each Slurm array task requests one core and 2G memory by default.

set -euo pipefail

ORIGINAL_SUBMIT_ARGS=("$@")

usage() {
  cat <<'EOF'
Usage:
  bash submit_fix_o2_simulation_array.sh \
    --fit_dir=oxygen/results/fit_invivo_O2_buffering_500seed \
    --simulation=invivo \
    --seeds=1:500 \
    --o2_values=0,0.5,1,2,5 \
    --initial_ploidy_values=2,4 \
    --time_days=1000 \
    --n_sim=10

Required simulation options:
  --fit_dir=DIR or --run_dir=DIR   Parent result dir containing seedXX/best_params.tsv.
  --simulation=MODE                invivo, invitro, or joint.
  --o2_values=CSV or --o2=CSV      Fixed O2 values.
  --initial_ploidy_values=CSV      Initial ploidy values. --initial_ploidy=CSV also works.
  --time_days=N                    Simulation horizon in days.
  --n_sim=N                        Number of replicate simulations per seed/O2/ploidy.

Common optional simulation options:
  --seeds=CSV/RANGE                Seed filter, e.g. 1:500, 25,322, or seed25,seed322.
  --out_dir=DIR                    Batch output root; defaults to PROJECT_ROOT/oxygen/results/O2_fixed_simulation.
  --dt=N
  --save_every_days=N
  --report_dt=N
  --initial_cells=N
  --joint_scope=shared_invivo|invitro_effective
  --Crowding=TRUE|FALSE
  --O2_growth=TRUE|FALSE
  --start_with=MODE
  --ploidy_O2_death=MODE

HPC options:
  --project_root=DIR
  --simulation_script=FILE
  --array_script=FILE
  --tasks_tsv=FILE                 Existing or generated task table path.
  --refresh_task_list=TRUE|FALSE   Rebuild task table; defaults TRUE.
  --array_spec=SPEC                Slurm array spec; defaults to all task IDs.
  --array_max_concurrent=N         Adds %N to array spec.
  --cpus_per_task=N                Defaults to 1.
  --mem=SIZE                       Defaults to 2G.
  --qos=NAME                       Defaults to small.
  --time_limit=HH:MM:SS            Defaults to 4:00:00.
  --job_name=NAME                  Defaults to o2fix_sim.
  --log_root=DIR                   Defaults to PROJECT_ROOT/oxygen/results/log.
  --r_module=MODULE                Defaults to R/4.4.
  --skip_existing=TRUE|FALSE       Worker skips tasks with all required outputs; defaults FALSE.
  --dry_run=TRUE|FALSE             Build task list and print sbatch command only.
EOF
}

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

load_r_module() {
  if command -v ml >/dev/null 2>&1; then
    ml "${R_MODULE}"
  elif command -v module >/dev/null 2>&1; then
    module load "${R_MODULE}"
  fi
}

shell_join() {
  local out=""
  local token
  local quoted
  for token in "$@"; do
    printf -v quoted "%q" "${token}"
    out+="${quoted} "
  done
  printf "%s" "${out% }"
}

write_env_line() {
  local name="$1"
  local value="${2-}"
  local quoted
  printf -v quoted "%q" "${value}"
  printf "%s=%s\n" "${name}" "${quoted}" >> "${TASK_ENV_FILE}"
}

submit_or_print() {
  {
    printf "Submit fixed-O2 simulation array:"
    printf " %q" "$@"
    printf "\n"
  } | tee -a "${PROGRESS_LOG}" >&2
  if truthy "${DRY_RUN}"; then
    echo "DRY_RUN_JOB_ID"
  else
    "$@" | awk '{print $NF}'
  fi
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h) usage; exit 0 ;;
      --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
      --simulation_script=*) SIM_SCRIPT="${arg#*=}" ;;
      --array_script=*) ARRAY_SCRIPT="${arg#*=}" ;;
      --tasks_tsv=*|--task_list=*) TASKS_TSV="${arg#*=}" ;;
      --refresh_task_list=*) REFRESH_TASK_LIST="${arg#*=}" ;;
      --fit_dir=*) FIT_DIR="${arg#*=}" ;;
      --run_dir=*) RUN_DIR="${arg#*=}" ;;
      --simulation=*) SIMULATION="${arg#*=}" ;;
      --seeds=*) SEEDS="${arg#*=}" ;;
      --o2_values=*) O2_VALUES="${arg#*=}" ;;
      --o2=*) O2_VALUES="${arg#*=}" ;;
      --initial_ploidy_values=*) INITIAL_PLOIDY_VALUES="${arg#*=}" ;;
      --initial_ploidy=*) INITIAL_PLOIDY_VALUES="${arg#*=}" ;;
      --time_days=*) TIME_DAYS="${arg#*=}" ;;
      --n_sim=*) N_SIM="${arg#*=}" ;;
      --out_dir=*) OUT_ROOT="${arg#*=}" ;;
      --dt=*) DT="${arg#*=}" ;;
      --save_every_days=*) SAVE_EVERY_DAYS="${arg#*=}" ;;
      --report_dt=*) REPORT_DT="${arg#*=}" ;;
      --initial_cells=*) INITIAL_CELLS="${arg#*=}" ;;
      --joint_scope=*) JOINT_SCOPE="${arg#*=}" ;;
      --Crowding=*) CROWDING="${arg#*=}" ;;
      --O2_growth=*) O2_GROWTH="${arg#*=}" ;;
      --start_with=*) START_WITH="${arg#*=}" ;;
      --ploidy_O2_death=*) PLOIDY_O2_DEATH="${arg#*=}" ;;
      --array_spec=*|--array=*) ARRAY_SPEC="${arg#*=}" ;;
      --array_max_concurrent=*|--max_concurrent=*) ARRAY_MAX_CONCURRENT="${arg#*=}" ;;
      --cpus_per_task=*|--cpus-per-task=*) CPUS_PER_TASK="${arg#*=}" ;;
      --mem=*) MEM="${arg#*=}" ;;
      --qos=*) QOS="${arg#*=}" ;;
      --time_limit=*|--time=*) TIME_LIMIT="${arg#*=}" ;;
      --job_name=*) JOB_NAME="${arg#*=}" ;;
      --log_root=*|--log_dir=*) LOG_ROOT="${arg#*=}" ;;
      --r_module=*) R_MODULE="${arg#*=}" ;;
      --skip_existing=*) SKIP_EXISTING="${arg#*=}" ;;
      --dry_run=*) DRY_RUN="${arg#*=}" ;;
      *) echo "Unknown argument: ${arg}" >&2; usage >&2; exit 2 ;;
    esac
  done
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
SIM_SCRIPT="${SIM_SCRIPT:-}"
ARRAY_SCRIPT="${ARRAY_SCRIPT:-}"
TASKS_TSV="${TASKS_TSV:-}"
REFRESH_TASK_LIST="${REFRESH_TASK_LIST:-TRUE}"
FIT_DIR="${FIT_DIR:-}"
RUN_DIR="${RUN_DIR:-}"
SIMULATION="${SIMULATION:-}"
SEEDS="${SEEDS:-}"
O2_VALUES="${O2_VALUES:-}"
INITIAL_PLOIDY_VALUES="${INITIAL_PLOIDY_VALUES:-}"
TIME_DAYS="${TIME_DAYS:-}"
N_SIM="${N_SIM:-}"
OUT_ROOT="${OUT_ROOT:-}"
DT="${DT:-}"
SAVE_EVERY_DAYS="${SAVE_EVERY_DAYS:-}"
REPORT_DT="${REPORT_DT:-}"
INITIAL_CELLS="${INITIAL_CELLS:-}"
JOINT_SCOPE="${JOINT_SCOPE:-}"
CROWDING="${CROWDING:-}"
O2_GROWTH="${O2_GROWTH:-}"
START_WITH="${START_WITH:-}"
PLOIDY_O2_DEATH="${PLOIDY_O2_DEATH:-}"
ARRAY_SPEC="${ARRAY_SPEC:-}"
ARRAY_MAX_CONCURRENT="${ARRAY_MAX_CONCURRENT:-}"
CPUS_PER_TASK="${CPUS_PER_TASK:-1}"
MEM="${MEM:-2G}"
QOS="${QOS:-small}"
TIME_LIMIT="${TIME_LIMIT:-4:00:00}"
JOB_NAME="${JOB_NAME:-o2fix_sim}"
LOG_ROOT="${LOG_ROOT:-}"
R_MODULE="${R_MODULE:-R/4.4}"
SKIP_EXISTING="${SKIP_EXISTING:-FALSE}"
DRY_RUN="${DRY_RUN:-FALSE}"

parse_args "$@"

if [[ -z "${RUN_DIR}" && -z "${FIT_DIR}" ]]; then
  echo "Provide --fit_dir or --run_dir pointing to a parent directory with seedXX subdirectories." >&2
  exit 2
fi
for required_name in SIMULATION O2_VALUES INITIAL_PLOIDY_VALUES TIME_DAYS N_SIM; do
  if [[ -z "${!required_name}" ]]; then
    echo "Missing required option: ${required_name}" >&2
    usage >&2
    exit 2
  fi
done
if ! [[ "${CPUS_PER_TASK}" =~ ^[0-9]+$ ]] || (( CPUS_PER_TASK < 1 )); then
  echo "cpus_per_task must be a positive integer, got: ${CPUS_PER_TASK}" >&2
  exit 2
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
if [[ -z "${SIM_SCRIPT}" ]]; then
  SIM_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/simulation/fix_o2_simulation.R"
fi
if [[ -z "${ARRAY_SCRIPT}" ]]; then
  ARRAY_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/hpc/run_fix_o2_simulation_array.sub"
fi
if [[ -z "${OUT_ROOT}" ]]; then
  OUT_ROOT="${PROJECT_ROOT}/oxygen/results/O2_fixed_simulation"
elif [[ "${OUT_ROOT}" != /* ]]; then
  OUT_ROOT="${PROJECT_ROOT}/${OUT_ROOT}"
fi
if [[ -z "${LOG_ROOT}" ]]; then
  LOG_ROOT="${PROJECT_ROOT}/oxygen/results/log"
elif [[ "${LOG_ROOT}" != /* ]]; then
  LOG_ROOT="${PROJECT_ROOT}/${LOG_ROOT}"
fi

SIM_SCRIPT="$(cd "$(dirname "${SIM_SCRIPT}")" && pwd)/$(basename "${SIM_SCRIPT}")"
ARRAY_SCRIPT="$(cd "$(dirname "${ARRAY_SCRIPT}")" && pwd)/$(basename "${ARRAY_SCRIPT}")"
mkdir -p "${OUT_ROOT}" "${LOG_ROOT}"
OUT_ROOT="$(cd "${OUT_ROOT}" && pwd)"
LOG_ROOT="$(cd "${LOG_ROOT}" && pwd)"

if [[ -z "${TASKS_TSV}" ]]; then
  TASKS_TSV="${OUT_ROOT}/${SIMULATION}/task_list.tsv"
elif [[ "${TASKS_TSV}" != /* ]]; then
  TASKS_TSV="${PROJECT_ROOT}/${TASKS_TSV}"
fi
mkdir -p "$(dirname "${TASKS_TSV}")"
TASKS_TSV="$(cd "$(dirname "${TASKS_TSV}")" && pwd)/$(basename "${TASKS_TSV}")"

TASK_ROOT="$(cd "$(dirname "${TASKS_TSV}")" && pwd)"
STATUS_DIR="${TASK_ROOT}/hpc_task_status"
TASK_ENV_FILE="${TASK_ROOT}/hpc_submit_config.sh"
PROGRESS_LOG="${TASK_ROOT}/hpc_submit_progress.log"
JOBS_TSV="${TASK_ROOT}/hpc_submit_jobs.tsv"
: > "${PROGRESS_LOG}"
printf "job_type\tjob_id\tarray_spec\tqos\twalltime\tmem\tcpus_per_task\ttask_table\n" > "${JOBS_TSV}"

if [[ ! -f "${SIM_SCRIPT}" ]]; then
  echo "Missing simulation script: ${SIM_SCRIPT}" >&2
  exit 1
fi
if [[ ! -f "${ARRAY_SCRIPT}" ]]; then
  echo "Missing array script: ${ARRAY_SCRIPT}" >&2
  exit 1
fi

build_cmd=(
  Rscript "${SIM_SCRIPT}"
  "--build_task_list_only=TRUE"
  "--n_core=1"
  "--simulation=${SIMULATION}"
  "--o2_values=${O2_VALUES}"
  "--initial_ploidy_values=${INITIAL_PLOIDY_VALUES}"
  "--time_days=${TIME_DAYS}"
  "--n_sim=${N_SIM}"
  "--out_dir=${OUT_ROOT}"
)
if [[ -n "${FIT_DIR}" ]]; then build_cmd+=("--fit_dir=${FIT_DIR}"); fi
if [[ -n "${RUN_DIR}" ]]; then build_cmd+=("--run_dir=${RUN_DIR}"); fi
if [[ -n "${SEEDS}" ]]; then build_cmd+=("--seeds=${SEEDS}"); fi
if [[ -n "${DT}" ]]; then build_cmd+=("--dt=${DT}"); fi
if [[ -n "${SAVE_EVERY_DAYS}" ]]; then build_cmd+=("--save_every_days=${SAVE_EVERY_DAYS}"); fi
if [[ -n "${REPORT_DT}" ]]; then build_cmd+=("--report_dt=${REPORT_DT}"); fi
if [[ -n "${INITIAL_CELLS}" ]]; then build_cmd+=("--initial_cells=${INITIAL_CELLS}"); fi
if [[ -n "${JOINT_SCOPE}" ]]; then build_cmd+=("--joint_scope=${JOINT_SCOPE}"); fi
if [[ -n "${CROWDING}" ]]; then build_cmd+=("--Crowding=${CROWDING}"); fi
if [[ -n "${O2_GROWTH}" ]]; then build_cmd+=("--O2_growth=${O2_GROWTH}"); fi
if [[ -n "${START_WITH}" ]]; then build_cmd+=("--start_with=${START_WITH}"); fi
if [[ -n "${PLOIDY_O2_DEATH}" ]]; then build_cmd+=("--ploidy_O2_death=${PLOIDY_O2_DEATH}"); fi

if [[ ! -f "${TASKS_TSV}" ]] || truthy "${REFRESH_TASK_LIST}"; then
  {
    printf "Build fixed-O2 task list:"
    printf " %q" "${build_cmd[@]}"
    printf "\n"
  } | tee -a "${PROGRESS_LOG}" >&2
  load_r_module
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript not found after loading ${R_MODULE}." >&2
    exit 1
  fi
  (cd "${PROJECT_ROOT}" && "${build_cmd[@]}") | tee -a "${PROGRESS_LOG}"
fi

if [[ ! -f "${TASKS_TSV}" ]]; then
  echo "Missing task table after build: ${TASKS_TSV}" >&2
  exit 1
fi
TASK_COUNT=$(( $(wc -l < "${TASKS_TSV}") - 1 ))
if (( TASK_COUNT < 1 )); then
  echo "Task table has no data rows: ${TASKS_TSV}" >&2
  exit 1
fi

if [[ -z "${ARRAY_SPEC}" ]]; then
  ARRAY_SPEC="1-${TASK_COUNT}"
fi
if [[ -n "${ARRAY_MAX_CONCURRENT}" && "${ARRAY_SPEC}" != *%* ]]; then
  ARRAY_SPEC="${ARRAY_SPEC}%${ARRAY_MAX_CONCURRENT}"
fi

: > "${TASK_ENV_FILE}"
write_env_line "PROJECT_ROOT" "${PROJECT_ROOT}"
write_env_line "SIM_SCRIPT" "${SIM_SCRIPT}"
write_env_line "OUT_ROOT" "${OUT_ROOT}"
write_env_line "SIMULATION" "${SIMULATION}"
write_env_line "TIME_DAYS" "${TIME_DAYS}"
write_env_line "DT" "${DT}"
write_env_line "SAVE_EVERY_DAYS" "${SAVE_EVERY_DAYS}"
write_env_line "REPORT_DT" "${REPORT_DT}"
write_env_line "INITIAL_CELLS" "${INITIAL_CELLS}"
write_env_line "JOINT_SCOPE" "${JOINT_SCOPE}"
write_env_line "CROWDING" "${CROWDING}"
write_env_line "O2_GROWTH" "${O2_GROWTH}"
write_env_line "START_WITH" "${START_WITH}"
write_env_line "PLOIDY_O2_DEATH" "${PLOIDY_O2_DEATH}"

SUBMIT_COMMAND_TEXT="$(shell_join bash "${BASH_SOURCE[0]}" "${ORIGINAL_SUBMIT_ARGS[@]}")"
{
  printf "key\tvalue\n"
  printf "project_root\t%s\n" "${PROJECT_ROOT}"
  printf "simulation_script\t%s\n" "${SIM_SCRIPT}"
  printf "array_script\t%s\n" "${ARRAY_SCRIPT}"
  printf "task_table\t%s\n" "${TASKS_TSV}"
  printf "task_count\t%s\n" "${TASK_COUNT}"
  printf "out_root\t%s\n" "${OUT_ROOT}"
  printf "simulation\t%s\n" "${SIMULATION}"
  printf "o2_values\t%s\n" "${O2_VALUES}"
  printf "initial_ploidy_values\t%s\n" "${INITIAL_PLOIDY_VALUES}"
  printf "time_days\t%s\n" "${TIME_DAYS}"
  printf "n_sim\t%s\n" "${N_SIM}"
  printf "array_spec\t%s\n" "${ARRAY_SPEC}"
  printf "cpus_per_task\t%s\n" "${CPUS_PER_TASK}"
  printf "mem\t%s\n" "${MEM}"
  printf "qos\t%s\n" "${QOS}"
  printf "time_limit\t%s\n" "${TIME_LIMIT}"
  printf "r_module\t%s\n" "${R_MODULE}"
  printf "skip_existing\t%s\n" "${SKIP_EXISTING}"
  printf "submit_command\t%s\n" "${SUBMIT_COMMAND_TEXT}"
} > "${TASK_ROOT}/hpc_submit_manifest.tsv"

export_arg="ALL,PROJECT_ROOT=${PROJECT_ROOT},TASKS_TSV=${TASKS_TSV},TASK_ENV_FILE=${TASK_ENV_FILE},TASK_LOOKUP_COLUMN=task_id,R_MODULE=${R_MODULE},SKIP_EXISTING=${SKIP_EXISTING},STATUS_DIR=${STATUS_DIR},DRY_RUN=FALSE"

if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found; run on an HPC login node or use --dry_run=TRUE." >&2
  exit 1
fi

array_job_id="$(submit_or_print \
  sbatch \
  "--job-name=${JOB_NAME}" \
  "--array=${ARRAY_SPEC}" \
  "--cpus-per-task=${CPUS_PER_TASK}" \
  "--mem=${MEM}" \
  "--qos=${QOS}" \
  "--time=${TIME_LIMIT}" \
  "--output=${LOG_ROOT}/o2fix_sim_%A_%a.out" \
  "--error=${LOG_ROOT}/o2fix_sim_%A_%a.err" \
  "--export=${export_arg}" \
  "${ARRAY_SCRIPT}")"

printf "array\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "${array_job_id}" "${ARRAY_SPEC}" "${QOS}" "${TIME_LIMIT}" "${MEM}" "${CPUS_PER_TASK}" "${TASKS_TSV}" >> "${JOBS_TSV}"

echo "Submitted fixed-O2 simulation array"
echo "  job_id: ${array_job_id}"
echo "  task_count: ${TASK_COUNT}"
echo "  array_spec: ${ARRAY_SPEC}"
echo "  task_table: ${TASKS_TSV}"
echo "  status_dir: ${STATUS_DIR}"
echo "  logs: ${LOG_ROOT}/o2fix_sim_%A_%a.out"
