#!/usr/bin/env bash
# Submit the full O2 parameter-landscape analysis to Slurm.
#
# The batch job computes fixed-O2 mode labels, estimates parameter contributions
# to Mode1/Mode2 separation, then regenerates the UMAP outputs.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash submit_parameter_landscape_full.sh [options]

Options:
  --project_root=/share/.../miningcloneid_soft_coupling
  --result_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/parameter_landscape
  --invivo_input=oxygen/results/fit_invivo_O2_buffering_500seed
  --invitro_input=oxygen/results/fit_invitro_O2_buffering_500seed
  --r_module=R/4.4.2-gfbf-2024a
  --job_name=o2_param_landscape
  --cpus=8
  --mem=96G
  --time=2-00:00:00
  --qos=xxlarge
  --partition=NAME
  --account=NAME
  --max_seeds=N
  --mode_reference_o2=2
  --mode_reference_o2_values=0,0.1,0.5,1,2,5
  --attractor_o2_grid=0,0.05,...,5
  --summary_o2=0,0.1,0.5,1,2,5
  --attractor_feature_o2_values=0,0.1,0.5,1,2,5
  --overwrite_modes=TRUE|FALSE
  --overwrite_parameter_contribution=TRUE|FALSE
  --mode_contribution_bootstrap=100
  --run_parameter_contribution=TRUE|FALSE
  --run_umap=TRUE|FALSE
  --dry_run=TRUE|FALSE
  --help

Default behavior:
  1. Generate UMAP tables, including fixed-O2 mode tables for in vivo.
  2. Estimate parameter contributions to fixed-O2 Mode1/Mode2 separation.
  3. Generate in vivo, in vitro, and pooled in vivo/in vitro UMAPs.
  4. Render the parameter_landscape_clustering_umap_cluster_report.html report.

By default, mode tables are written for reference O2 values 0,0.1,0.5,1,2,5,
and UMAP shape grouping uses --mode_reference_o2=2.
Existing fixed-O2 mode outputs are reused unless --overwrite_modes=TRUE.
Existing mode-contribution outputs are reused unless --overwrite_parameter_contribution=TRUE.
EOF
}

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

resolve_hpc_path() {
  local base="$1"
  local path="$2"
  case "${path}" in
    "~") path="${HOME}" ;;
    "~/"*) path="${HOME}/${path#"~/"}" ;;
  esac
  case "${path}" in
    /*) printf "%s" "${path}" ;;
    *) printf "%s/%s" "${base}" "${path}" ;;
  esac
}

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling"
DEFAULT_RESULT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/parameter_landscape"
DEFAULT_INVIVO_INPUT="oxygen/results/fit_invivo_O2_buffering_500seed"
DEFAULT_INVITRO_INPUT="oxygen/results/fit_invitro_O2_buffering_500seed"
DEFAULT_R_MODULE="R/4.4.2-gfbf-2024a"
DEFAULT_JOB_NAME="o2_param_landscape"
DEFAULT_CPUS="8"
DEFAULT_MEM="96G"
DEFAULT_TIME="2-00:00:00"
DEFAULT_QOS="xxlarge"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
RESULT_ROOT="${RESULT_ROOT:-${DEFAULT_RESULT_ROOT}}"
INVIVO_INPUT="${INVIVO_INPUT:-${DEFAULT_INVIVO_INPUT}}"
INVITRO_INPUT="${INVITRO_INPUT:-${DEFAULT_INVITRO_INPUT}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
JOB_NAME="${JOB_NAME:-${DEFAULT_JOB_NAME}}"
CPUS="${CPUS:-${DEFAULT_CPUS}}"
MEM="${MEM:-${DEFAULT_MEM}}"
TIME_LIMIT="${TIME_LIMIT:-${DEFAULT_TIME}}"
QOS="${QOS:-${DEFAULT_QOS}}"
PARTITION="${PARTITION:-}"
ACCOUNT="${ACCOUNT:-}"
MAX_SEEDS="${MAX_SEEDS:-}"
MODE_REFERENCE_O2="${MODE_REFERENCE_O2:-2}"
MODE_REFERENCE_O2_VALUES="${MODE_REFERENCE_O2_VALUES:-0,0.1,0.5,1,2,5}"
ATTRACTOR_O2_GRID="${ATTRACTOR_O2_GRID:-}"
SUMMARY_O2="${SUMMARY_O2:-}"
ATTRACTOR_FEATURE_O2_VALUES="${ATTRACTOR_FEATURE_O2_VALUES:-}"
OVERWRITE_MODES="${OVERWRITE_MODES:-FALSE}"
OVERWRITE_PARAMETER_CONTRIBUTION="${OVERWRITE_PARAMETER_CONTRIBUTION:-FALSE}"
MODE_CONTRIBUTION_BOOTSTRAP="${MODE_CONTRIBUTION_BOOTSTRAP:-100}"
RUN_PARAMETER_CONTRIBUTION="${RUN_PARAMETER_CONTRIBUTION:-TRUE}"
RUN_UMAP="${RUN_UMAP:-TRUE}"
DRY_RUN="${DRY_RUN:-FALSE}"

for arg in "$@"; do
  case "${arg}" in
    --help|-h)
      usage
      exit 0
      ;;
    --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
    --result_root=*) RESULT_ROOT="${arg#*=}" ;;
    --invivo_input=*) INVIVO_INPUT="${arg#*=}" ;;
    --invitro_input=*) INVITRO_INPUT="${arg#*=}" ;;
    --r_module=*) R_MODULE="${arg#*=}" ;;
    --job_name=*) JOB_NAME="${arg#*=}" ;;
    --cpus=*) CPUS="${arg#*=}" ;;
    --mem=*) MEM="${arg#*=}" ;;
    --time=*) TIME_LIMIT="${arg#*=}" ;;
    --qos=*) QOS="${arg#*=}" ;;
    --partition=*) PARTITION="${arg#*=}" ;;
    --account=*) ACCOUNT="${arg#*=}" ;;
    --max_seeds=*) MAX_SEEDS="${arg#*=}" ;;
    --mode_reference_o2=*) MODE_REFERENCE_O2="${arg#*=}" ;;
    --mode_reference_o2_values=*) MODE_REFERENCE_O2_VALUES="${arg#*=}" ;;
    --attractor_o2_grid=*) ATTRACTOR_O2_GRID="${arg#*=}" ;;
    --summary_o2=*) SUMMARY_O2="${arg#*=}" ;;
    --attractor_feature_o2_values=*) ATTRACTOR_FEATURE_O2_VALUES="${arg#*=}" ;;
    --overwrite_modes=*|--force_modes=*) OVERWRITE_MODES="${arg#*=}" ;;
    --overwrite_parameter_contribution=*|--force_parameter_contribution=*|--force_extra_results=*) OVERWRITE_PARAMETER_CONTRIBUTION="${arg#*=}" ;;
    --mode_contribution_bootstrap=*) MODE_CONTRIBUTION_BOOTSTRAP="${arg#*=}" ;;
    --run_parameter_contribution=*) RUN_PARAMETER_CONTRIBUTION="${arg#*=}" ;;
    --run_umap=*) RUN_UMAP="${arg#*=}" ;;
    --dry_run=*) DRY_RUN="${arg#*=}" ;;
    *)
      echo "Unknown argument: ${arg}" >&2
      usage >&2
      exit 2
      ;;
  esac
done

RESULT_ROOT="$(resolve_hpc_path "${PROJECT_ROOT}" "${RESULT_ROOT}")"
LOG_DIR="${RESULT_ROOT}/logs"
STAMP="$(date +%Y%m%d_%H%M%S)"
SBATCH_SCRIPT="${LOG_DIR}/run_parameter_landscape_full_${STAMP}.sbatch"

mkdir -p "${LOG_DIR}"

{
  printf '%s\n' '#!/usr/bin/env bash'
  printf '#SBATCH --job-name=%s\n' "${JOB_NAME}"
  printf '%s\n' '#SBATCH --ntasks=1'
  printf '#SBATCH --cpus-per-task=%s\n' "${CPUS}"
  printf '#SBATCH --mem=%s\n' "${MEM}"
  printf '#SBATCH --time=%s\n' "${TIME_LIMIT}"
  if [[ -n "${QOS}" ]]; then
    printf '#SBATCH --qos=%s\n' "${QOS}"
  fi
  if [[ -n "${PARTITION}" ]]; then
    printf '#SBATCH --partition=%s\n' "${PARTITION}"
  fi
  if [[ -n "${ACCOUNT}" ]]; then
    printf '#SBATCH --account=%s\n' "${ACCOUNT}"
  fi
  printf '#SBATCH --output=%s/%s_%%j.out\n' "${LOG_DIR}" "${JOB_NAME}"
  printf '#SBATCH --error=%s/%s_%%j.err\n' "${LOG_DIR}" "${JOB_NAME}"
  printf '\n'
  printf '%s\n' 'set -euo pipefail'
  printf 'PROJECT_ROOT=%q\n' "${PROJECT_ROOT}"
  printf 'RESULT_ROOT=%q\n' "${RESULT_ROOT}"
  printf 'INVIVO_INPUT=%q\n' "${INVIVO_INPUT}"
  printf 'INVITRO_INPUT=%q\n' "${INVITRO_INPUT}"
  printf 'R_MODULE=%q\n' "${R_MODULE}"
  printf 'MAX_SEEDS=%q\n' "${MAX_SEEDS}"
  printf 'MODE_REFERENCE_O2=%q\n' "${MODE_REFERENCE_O2}"
  printf 'MODE_REFERENCE_O2_VALUES=%q\n' "${MODE_REFERENCE_O2_VALUES}"
  printf 'ATTRACTOR_O2_GRID=%q\n' "${ATTRACTOR_O2_GRID}"
  printf 'SUMMARY_O2=%q\n' "${SUMMARY_O2}"
  printf 'ATTRACTOR_FEATURE_O2_VALUES=%q\n' "${ATTRACTOR_FEATURE_O2_VALUES}"
  printf 'OVERWRITE_MODES=%q\n' "${OVERWRITE_MODES}"
  printf 'OVERWRITE_PARAMETER_CONTRIBUTION=%q\n' "${OVERWRITE_PARAMETER_CONTRIBUTION}"
  printf 'MODE_CONTRIBUTION_BOOTSTRAP=%q\n' "${MODE_CONTRIBUTION_BOOTSTRAP}"
  printf 'RUN_PARAMETER_CONTRIBUTION=%q\n' "${RUN_PARAMETER_CONTRIBUTION}"
  printf 'RUN_UMAP=%q\n' "${RUN_UMAP}"
  cat <<'BATCH_BODY'

SCRIPT_DIR="oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

require_file() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    echo "Missing expected output: ${path}" >&2
    return 1
  fi
}

load_r_module() {
  if [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi
  if command -v module >/dev/null 2>&1; then
    module use /app/eb/modules/all >/dev/null 2>&1 || true
  fi
  if command -v ml >/dev/null 2>&1; then
    ml "${R_MODULE}"
  elif command -v module >/dev/null 2>&1; then
    module load "${R_MODULE}"
  fi
}

run_rscript() {
  local label="$1"
  shift
  echo
  echo "[$(date '+%F %T')] ${label}"
  printf 'Command: Rscript'
  printf ' %q' "$@"
  printf '\n'
  Rscript "$@"
}

cd "${PROJECT_ROOT}"
mkdir -p "${RESULT_ROOT}/logs"
MAIN_LOG="${RESULT_ROOT}/logs/parameter_landscape_full_${SLURM_JOB_ID:-manual}.log"
exec > >(tee -a "${MAIN_LOG}") 2>&1

echo "[$(date '+%F %T')] Starting full O2 parameter-landscape analysis"
echo "Project root: ${PROJECT_ROOT}"
echo "Result root: ${RESULT_ROOT}"
echo "In vivo input: ${INVIVO_INPUT}"
echo "In vitro input: ${INVITRO_INPUT}"
echo "Threads: ${THREADS}"

if [[ ! -d "${INVIVO_INPUT}" ]]; then
  echo "Missing in vivo input directory: ${INVIVO_INPUT}" >&2
  exit 1
fi
if [[ ! -d "${INVITRO_INPUT}" ]]; then
  echo "Missing in vitro input directory: ${INVITRO_INPUT}" >&2
  exit 1
fi

load_r_module
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found after loading ${R_MODULE}." >&2
  exit 1
fi
Rscript --version

MAX_SEEDS_ARGS=()
if [[ -n "${MAX_SEEDS}" ]]; then
  MAX_SEEDS_ARGS=("--max_seeds=${MAX_SEEDS}")
fi

MODE_TABLE_ARGS=("--mode_reference_o2=${MODE_REFERENCE_O2}")
if [[ -n "${MODE_REFERENCE_O2_VALUES}" ]]; then
  MODE_TABLE_ARGS+=("--mode_reference_o2_values=${MODE_REFERENCE_O2_VALUES},${MODE_REFERENCE_O2}")
fi
if [[ -n "${ATTRACTOR_O2_GRID}" ]]; then
  MODE_TABLE_ARGS+=("--attractor_o2_grid=${ATTRACTOR_O2_GRID}")
fi
if [[ -n "${SUMMARY_O2}" ]]; then
  MODE_TABLE_ARGS+=("--summary_o2=${SUMMARY_O2}")
fi

ATTRACTOR_FEATURE_ARGS=()
if [[ -n "${ATTRACTOR_FEATURE_O2_VALUES}" ]]; then
  ATTRACTOR_FEATURE_ARGS=("--attractor_feature_o2_values=${ATTRACTOR_FEATURE_O2_VALUES}")
fi

run_rscript "Write in vivo UMAP tables and fixed-O2 mode tables" \
  "${SCRIPT_DIR}/invivo_umap_tables.R" \
  "--input_dir=${INVIVO_INPUT}" \
  "--result_root=${RESULT_ROOT}" \
  "--write_modes=TRUE" \
  "--overwrite_modes=${OVERWRITE_MODES}" \
  "--n_workers=${THREADS}" \
  "${MAX_SEEDS_ARGS[@]}" \
  "${MODE_TABLE_ARGS[@]}"

if truthy "${RUN_PARAMETER_CONTRIBUTION}"; then
  run_rscript "Estimate parameter contributions to fixed-O2 mode separation" \
    "${SCRIPT_DIR}/mode_parameter_contribution.R" \
    "--result_root=${RESULT_ROOT}" \
    "--best_csv=${RESULT_ROOT}/invivo/Tables/invivo_best_params_by_seed.csv" \
    "--mode_tables_dir=${RESULT_ROOT}/invivo/Tables/FixO2Modes" \
    "--mode_reference_o2=${MODE_REFERENCE_O2}" \
    "--mode_reference_o2_values=${MODE_REFERENCE_O2_VALUES},${MODE_REFERENCE_O2}" \
    "--n_bootstrap=${MODE_CONTRIBUTION_BOOTSTRAP}" \
    "--overwrite=${OVERWRITE_PARAMETER_CONTRIBUTION}" \
    "${MAX_SEEDS_ARGS[@]}"
  require_file "${RESULT_ROOT}/mode_parameter_contribution/mode_parameter_contribution_index.csv"
fi

if truthy "${RUN_UMAP}"; then
  run_rscript "Write in vitro UMAP tables" \
    "${SCRIPT_DIR}/invitro_umap_tables.R" \
    "--input_dir=${INVITRO_INPUT}" \
    "--result_root=${RESULT_ROOT}" \
    "${MAX_SEEDS_ARGS[@]}"

  run_rscript "Generate in vivo UMAP figures" \
    "${SCRIPT_DIR}/invivo_umap_figures.R" \
    "--result_root=${RESULT_ROOT}" \
    "--objective_seed_dir=${INVIVO_INPUT}" \
    "--mode_tables_dir=${RESULT_ROOT}/invivo/Tables/FixO2Modes" \
    "--mode_reference_o2=${MODE_REFERENCE_O2}" \
    "--shape_by_mode=TRUE" \
    "--drop_parameter_table_initial=TRUE" \
    "--umap_seed=123" \
    "--cluster_seed=123" \
    "--n_threads=${THREADS}" \
    "${ATTRACTOR_FEATURE_ARGS[@]}"

  run_rscript "Generate in vitro UMAP figures" \
    "${SCRIPT_DIR}/invitro_umap_figures.R" \
    "--result_root=${RESULT_ROOT}" \
    "--objective_seed_dir=${INVITRO_INPUT}" \
    "--drop_parameter_table_initial=TRUE" \
    "--umap_seed=123" \
    "--cluster_seed=123" \
    "--n_threads=${THREADS}"

  run_rscript "Generate pooled in vivo/in vitro UMAP figures and distance tables" \
    "${SCRIPT_DIR}/pooled_invivo_invitro_umap_figures.R" \
    "--result_root=${RESULT_ROOT}" \
    "--invivo_objective_seed_dir=${INVIVO_INPUT}" \
    "--invitro_objective_seed_dir=${INVITRO_INPUT}" \
    "--drop_parameter_table_initial=TRUE" \
    "--drop_invitro_parameter_table_initial=TRUE" \
    "--umap_seed=123" \
    "--cluster_seed=123" \
    "--n_threads=${THREADS}"

  run_rscript "Render UMAP cluster report" \
    "${SCRIPT_DIR}/umap_clusters.R" \
    "--result_root=${RESULT_ROOT}" \
    "--output_html=${RESULT_ROOT}/parameter_landscape_clustering_umap_cluster_report.html"

  require_file "${RESULT_ROOT}/parameter_landscape_clustering_umap_cluster_report.html"
fi

echo
echo "[$(date '+%F %T')] Finished full O2 parameter-landscape analysis"
echo "Main log: ${MAIN_LOG}"
echo "Result root: ${RESULT_ROOT}"
BATCH_BODY
} > "${SBATCH_SCRIPT}"

chmod +x "${SBATCH_SCRIPT}"

echo "Wrote Slurm script: ${SBATCH_SCRIPT}"
if truthy "${DRY_RUN}"; then
  echo "DRY_RUN=TRUE; not submitting."
  exit 0
fi

if ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found. Run this submitter on the HPC login node, or rerun with --dry_run=TRUE to only write the batch script." >&2
  exit 1
fi

sbatch "${SBATCH_SCRIPT}"
