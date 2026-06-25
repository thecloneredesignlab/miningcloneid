#!/bin/bash
# Submit fixed-O2 analytical-vs-simulation scatter plotting to Slurm.
#
# This launcher is meant to run on the HPC login node. It submits a Slurm job
# and lets Slurm choose the compute node unless --node is explicitly provided.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash oxygen/code/O2_supply_demand_MAP/hpc/submit_fixo2_noemi_scatters.sh \
    --simulation_ids=1,2,3

HPC options:
  --project_root=DIR              Defaults to this repository root.
  --node=NAME                     Optional compute node for Slurm --nodelist.
  --qos=NAME                      Defaults xxlarge.
  --time_limit=HH:MM:SS           Defaults 04:00:00.
  --mem=SIZE                      Defaults 32G.
  --cpus_per_task=N               Defaults 8.
  --n_workers=N                   Defaults to cpus_per_task.
  --job_name=NAME                 Defaults fixo2_scat.
  --log_root=DIR                  Defaults PROJECT_ROOT/oxygen/results/log.
  --r_module=MODULE               Defaults R/4.4.
  --dependency=SPEC               Optional Slurm dependency, e.g. afterok:12345.
  --dry_run=TRUE|FALSE            Print the sbatch command without submitting.

Plotting options:
  --simulation_dir=DIR            Defaults ~/oxygen/results/O2_fixed_simulation.
  --analysis_dir=DIR              Defaults ~/oxygen/results/analysis/FixO2_invivo_500seed.
  --fit_dir=DIR                   Defaults ~/oxygen/results/fit_invitro_O2_buffering_500seed.
  --run_dir=DIR                   Fitting seed directory root used to generate analytical solutions. Defaults fit_dir.
  --out_dir=DIR                   Defaults ~/oxygen/results/analysis/FixO2_invivo_500seed.
  --simulation_mode=MODE          Defaults invivo.
  --simulation_ids=CSV            Defaults 1,2,3.
  --seed_ids=CSV                  Optional seed filter for analytical generation.
  --time_points=CSV               Defaults 25,50,100,200,300,500,700,1000.
  --o2_values=CSV                 Defaults 0,0.1,0.5,1,2,5.
  --initial_ploidy_values=CSV     Defaults 2,4.
  --progress_every=N              Defaults 200.
  --analytical_methods=CSV        Defaults eigen,expm.
  --objective_transform=MODE      identity or log10. Defaults identity.
  --recompute=TRUE|FALSE          Defaults FALSE.
  --recompute_analytical=TRUE|FALSE
                                 Rebuild generated analytical trajectory cache. Defaults FALSE.
  --cache_all_times=TRUE|FALSE    Build/use all-time simulation metric cache. Defaults TRUE.
  --analytical_cache_table=FILE   Optional generated analytical trajectory cache path.
  --all_time_simulation_metric_table=FILE
                                 Optional all-time metric cache path.
EOF
}

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
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
  printf "%s=%s\n" "${name}" "${quoted}" >> "${ARG_ENV_FILE}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
RUN_SCRIPT="${SCRIPT_DIR}/run_fixo2_noemi_scatters.sub"
PLOT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/plot_fixo2_noemi_scatters.R"

NODE=""
QOS="xxlarge"
TIME_LIMIT="04:00:00"
MEM="32G"
CPUS_PER_TASK="8"
N_WORKERS=""
JOB_NAME="fixo2_scat"
LOG_ROOT=""
R_MODULE="R/4.4"
DEPENDENCY=""
DRY_RUN="FALSE"

SIMULATION_DIR="~/oxygen/results/O2_fixed_simulation"
ANALYSIS_DIR="~/oxygen/results/analysis/FixO2_invivo_500seed"
FIT_DIR="~/oxygen/results/fit_invitro_O2_buffering_500seed"
RUN_DIR=""
OUT_DIR="~/oxygen/results/analysis/FixO2_invivo_500seed"
SIMULATION_MODE="invivo"
SIMULATION_IDS="1,2,3"
SEED_IDS=""
TIME_POINTS="25,50,100,200,300,500,700,1000"
O2_VALUES="0,0.1,0.5,1,2,5"
INITIAL_PLOIDY_VALUES="2,4"
PROGRESS_EVERY="200"
ANALYTICAL_METHODS="eigen,expm"
OBJECTIVE_TRANSFORM="identity"
RECOMPUTE="FALSE"
RECOMPUTE_ANALYTICAL="FALSE"
CACHE_ALL_TIMES="TRUE"
ANALYTICAL_CACHE_TABLE=""
ALL_TIME_SIMULATION_METRIC_TABLE=""
SIMULATION_METRIC_TABLE=""
SIMULATION_SUMMARY_TABLE=""
SCATTER_DATA_TABLE=""

for arg in "$@"; do
  case "${arg}" in
    --help|-h) usage; exit 0 ;;
    --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
    --node=*) NODE="${arg#*=}" ;;
    --qos=*) QOS="${arg#*=}" ;;
    --time_limit=*|--time=*) TIME_LIMIT="${arg#*=}" ;;
    --mem=*) MEM="${arg#*=}" ;;
    --cpus_per_task=*|--cpus-per-task=*) CPUS_PER_TASK="${arg#*=}" ;;
    --n_workers=*|--n-workers=*) N_WORKERS="${arg#*=}" ;;
    --job_name=*) JOB_NAME="${arg#*=}" ;;
    --log_root=*) LOG_ROOT="${arg#*=}" ;;
    --r_module=*) R_MODULE="${arg#*=}" ;;
    --dependency=*) DEPENDENCY="${arg#*=}" ;;
    --dry_run=*) DRY_RUN="${arg#*=}" ;;
    --simulation_dir=*) SIMULATION_DIR="${arg#*=}" ;;
    --analysis_dir=*) ANALYSIS_DIR="${arg#*=}" ;;
    --fit_dir=*) FIT_DIR="${arg#*=}" ;;
    --run_dir=*) RUN_DIR="${arg#*=}" ;;
    --out_dir=*) OUT_DIR="${arg#*=}" ;;
    --simulation_mode=*) SIMULATION_MODE="${arg#*=}" ;;
    --simulation_ids=*) SIMULATION_IDS="${arg#*=}" ;;
    --seed_ids=*) SEED_IDS="${arg#*=}" ;;
    --time_points=*) TIME_POINTS="${arg#*=}" ;;
    --o2_values=*|--o2=*) O2_VALUES="${arg#*=}" ;;
    --initial_ploidy_values=*|--initial_ploidy=*) INITIAL_PLOIDY_VALUES="${arg#*=}" ;;
    --progress_every=*) PROGRESS_EVERY="${arg#*=}" ;;
    --analytical_methods=*) ANALYTICAL_METHODS="${arg#*=}" ;;
    --objective_transform=*) OBJECTIVE_TRANSFORM="${arg#*=}" ;;
    --recompute=*) RECOMPUTE="${arg#*=}" ;;
    --recompute_analytical=*) RECOMPUTE_ANALYTICAL="${arg#*=}" ;;
    --cache_all_times=*) CACHE_ALL_TIMES="${arg#*=}" ;;
    --analytical_cache_table=*) ANALYTICAL_CACHE_TABLE="${arg#*=}" ;;
    --all_time_simulation_metric_table=*) ALL_TIME_SIMULATION_METRIC_TABLE="${arg#*=}" ;;
    --simulation_metric_table=*) SIMULATION_METRIC_TABLE="${arg#*=}" ;;
    --simulation_summary_table=*) SIMULATION_SUMMARY_TABLE="${arg#*=}" ;;
    --scatter_data_table=*) SCATTER_DATA_TABLE="${arg#*=}" ;;
    *) echo "Unknown argument: ${arg}" >&2; usage >&2; exit 1 ;;
  esac
done

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
RUN_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/hpc/run_fixo2_noemi_scatters.sub"
PLOT_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/plot_fixo2_noemi_scatters.R"
LOG_ROOT="${LOG_ROOT:-${PROJECT_ROOT}/oxygen/results/log}"
N_WORKERS="${N_WORKERS:-${CPUS_PER_TASK}}"
mkdir -p "${LOG_ROOT}"

if [[ ! -f "${RUN_SCRIPT}" ]]; then
  echo "Missing Slurm worker script: ${RUN_SCRIPT}" >&2
  exit 1
fi
if [[ ! -f "${PLOT_SCRIPT}" ]]; then
  echo "Missing plotting script: ${PLOT_SCRIPT}" >&2
  exit 1
fi
if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found; run this launcher on an HPC login node or use --dry_run=TRUE." >&2
  exit 1
fi

timestamp="$(date +%Y%m%d_%H%M%S)"
ARG_ENV_FILE="${LOG_ROOT}/fixo2_noemi_scatters_${timestamp}.env"
: > "${ARG_ENV_FILE}"
write_env_line PROJECT_ROOT "${PROJECT_ROOT}"
write_env_line PLOT_SCRIPT "${PLOT_SCRIPT}"
write_env_line R_MODULE "${R_MODULE}"
write_env_line SIMULATION_DIR "${SIMULATION_DIR}"
write_env_line ANALYSIS_DIR "${ANALYSIS_DIR}"
write_env_line FIT_DIR "${FIT_DIR}"
write_env_line RUN_DIR "${RUN_DIR:-${FIT_DIR}}"
write_env_line OUT_DIR "${OUT_DIR}"
write_env_line SIMULATION_MODE "${SIMULATION_MODE}"
write_env_line SIMULATION_IDS "${SIMULATION_IDS}"
write_env_line SEED_IDS "${SEED_IDS}"
write_env_line TIME_POINTS "${TIME_POINTS}"
write_env_line O2_VALUES "${O2_VALUES}"
write_env_line INITIAL_PLOIDY_VALUES "${INITIAL_PLOIDY_VALUES}"
write_env_line PROGRESS_EVERY "${PROGRESS_EVERY}"
write_env_line ANALYTICAL_METHODS "${ANALYTICAL_METHODS}"
write_env_line OBJECTIVE_TRANSFORM "${OBJECTIVE_TRANSFORM}"
write_env_line RECOMPUTE "${RECOMPUTE}"
write_env_line RECOMPUTE_ANALYTICAL "${RECOMPUTE_ANALYTICAL}"
write_env_line CACHE_ALL_TIMES "${CACHE_ALL_TIMES}"
write_env_line N_WORKERS "${N_WORKERS}"
write_env_line ANALYTICAL_CACHE_TABLE "${ANALYTICAL_CACHE_TABLE}"
write_env_line ALL_TIME_SIMULATION_METRIC_TABLE "${ALL_TIME_SIMULATION_METRIC_TABLE}"
write_env_line SIMULATION_METRIC_TABLE "${SIMULATION_METRIC_TABLE}"
write_env_line SIMULATION_SUMMARY_TABLE "${SIMULATION_SUMMARY_TABLE}"
write_env_line SCATTER_DATA_TABLE "${SCATTER_DATA_TABLE}"

cmd=(
  sbatch
  "--job-name=${JOB_NAME}"
  "--cpus-per-task=${CPUS_PER_TASK}"
  "--mem=${MEM}"
  "--qos=${QOS}"
  "--time=${TIME_LIMIT}"
  "--output=${LOG_ROOT}/fixo2_noemi_scatters_%j.out"
  "--error=${LOG_ROOT}/fixo2_noemi_scatters_%j.err"
  "--export=ALL,ARG_ENV_FILE=${ARG_ENV_FILE}"
)

if [[ -n "${NODE}" ]]; then
  cmd+=("--nodelist=${NODE}")
fi
if [[ -n "${DEPENDENCY}" ]]; then
  cmd+=("--dependency=${DEPENDENCY}")
fi

cmd+=("${RUN_SCRIPT}")

echo "Submit fixed-O2 Noemi scatter plotting:"
echo "  command: $(shell_join "${cmd[@]}")"
echo "  arg_env_file: ${ARG_ENV_FILE}"
echo "  logs: ${LOG_ROOT}/fixo2_noemi_scatters_%j.out"

if truthy "${DRY_RUN}"; then
  echo "DRY_RUN_JOB_ID"
else
  "${cmd[@]}" | awk '{print $NF}'
fi
