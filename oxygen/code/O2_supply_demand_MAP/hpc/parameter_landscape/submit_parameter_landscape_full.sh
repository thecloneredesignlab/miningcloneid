#!/usr/bin/env bash
# Submit the O2 parameter-landscape analysis as a Slurm dependency workflow.
#
# The workflow is split into task-specific jobs:
#   A  in vivo tables + fixed-O2 modes
#   B  in vitro tables
#   C  mode parameter contribution array, one reference O2 per task
#   C2 merge mode-contribution summaries
#   D  in vivo dimensionality reductions
#   E  in vitro dimensionality reductions
#   F  pooled in vivo/in vitro dimensionality reductions
#   G  final HTML report

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash submit_parameter_landscape_full.sh [options]

Primary options:
  --project_root=/share/.../miningcloneid_soft_coupling
  --result_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/parameter_landscape
  --invivo_input=oxygen/results/fit_invivo_O2_buffering_500seed
  --invitro_input=oxygen/results/fit_invitro_O2_buffering_500seed
  --r_module=R/4.4.2-gfbf-2024a
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
  --run_report=TRUE|FALSE
  --reductions=umap,pca,tsne
  --preprocess_modes=zscore,prior_unit
  --pooled_preprocess_modes=zscore,context_prior_unit,common_prior_unit
  --run_full_tsne=TRUE|FALSE
  --dry_run=TRUE|FALSE

Default task resources:
  A invivo tables/modes: --cpus=8 --mem=96G --time=2-00:00:00
  B invitro tables:      --invitro_table_cpus=8 --invitro_table_mem=32G --invitro_table_time=4:00:00
  C mode contribution:   --mode_contribution_cpus=4 --mode_contribution_mem=24G --mode_contribution_time=4:00:00
  D invivo reductions:   --invivo_umap_cpus=20 --invivo_umap_mem=64G --invivo_umap_time=8:00:00
  E invitro reductions:  --invitro_umap_cpus=12 --invitro_umap_mem=48G --invitro_umap_time=6:00:00
  F pooled reductions:   --pooled_umap_cpus=20 --pooled_umap_mem=64G --pooled_umap_time=8:00:00
  G report:              --report_cpus=2 --report_mem=16G --report_time=2:00:00

Compatibility:
  --cpus, --mem, and --time set the A task resources. Other task resources
  can be overridden independently with the task-specific options above.
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

append_unique_csv_value() {
  local current="$1"
  local value="$2"
  value="${value//[[:space:]]/}"
  if [[ -z "${value}" ]]; then
    printf "%s" "${current}"
    return 0
  fi
  case ",${current}," in
    *",${value},"*) printf "%s" "${current}" ;;
    *) if [[ -z "${current}" ]]; then printf "%s" "${value}"; else printf "%s,%s" "${current}" "${value}"; fi ;;
  esac
}

normalize_reference_o2_values() {
  local raw="$1"
  local primary="$2"
  local out=""
  local value
  IFS=',' read -r -a values <<< "${raw},${primary}"
  for value in "${values[@]}"; do
    out="$(append_unique_csv_value "${out}" "${value}")"
  done
  if [[ -z "${out}" ]]; then
    echo "At least one reference O2 value is required." >&2
    exit 2
  fi
  printf "%s" "${out}"
}

csv_count() {
  local raw="$1"
  IFS=',' read -r -a values <<< "${raw}"
  printf "%s" "${#values[@]}"
}

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling"
DEFAULT_RESULT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/parameter_landscape"
DEFAULT_INVIVO_INPUT="oxygen/results/fit_invivo_O2_buffering_500seed"
DEFAULT_INVITRO_INPUT="oxygen/results/fit_invitro_O2_buffering_500seed"
DEFAULT_R_MODULE="R/4.4.2-gfbf-2024a"
DEFAULT_QOS="xxlarge"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
RESULT_ROOT="${RESULT_ROOT:-${DEFAULT_RESULT_ROOT}}"
INVIVO_INPUT="${INVIVO_INPUT:-${DEFAULT_INVIVO_INPUT}}"
INVITRO_INPUT="${INVITRO_INPUT:-${DEFAULT_INVITRO_INPUT}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
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
MODE_CONTRIBUTION_ARRAY_CONCURRENCY="${MODE_CONTRIBUTION_ARRAY_CONCURRENCY:-6}"
RUN_PARAMETER_CONTRIBUTION="${RUN_PARAMETER_CONTRIBUTION:-TRUE}"
RUN_UMAP="${RUN_UMAP:-TRUE}"
RUN_REPORT="${RUN_REPORT:-TRUE}"
REDUCTIONS="${REDUCTIONS:-umap,pca,tsne}"
PREPROCESS_MODES="${PREPROCESS_MODES:-zscore,prior_unit}"
POOLED_PREPROCESS_MODES="${POOLED_PREPROCESS_MODES:-zscore,context_prior_unit,common_prior_unit}"
RUN_FULL_TSNE="${RUN_FULL_TSNE:-FALSE}"
DRY_RUN="${DRY_RUN:-FALSE}"

INVIVO_TABLE_CPUS="${INVIVO_TABLE_CPUS:-${CPUS:-8}}"
INVIVO_TABLE_MEM="${INVIVO_TABLE_MEM:-${MEM:-96G}}"
INVIVO_TABLE_TIME="${INVIVO_TABLE_TIME:-${TIME_LIMIT:-2-00:00:00}}"
INVITRO_TABLE_CPUS="${INVITRO_TABLE_CPUS:-8}"
INVITRO_TABLE_MEM="${INVITRO_TABLE_MEM:-32G}"
INVITRO_TABLE_TIME="${INVITRO_TABLE_TIME:-4:00:00}"
MODE_CONTRIBUTION_CPUS="${MODE_CONTRIBUTION_CPUS:-4}"
MODE_CONTRIBUTION_MEM="${MODE_CONTRIBUTION_MEM:-24G}"
MODE_CONTRIBUTION_TIME="${MODE_CONTRIBUTION_TIME:-4:00:00}"
MODE_CONTRIBUTION_MERGE_CPUS="${MODE_CONTRIBUTION_MERGE_CPUS:-2}"
MODE_CONTRIBUTION_MERGE_MEM="${MODE_CONTRIBUTION_MERGE_MEM:-8G}"
MODE_CONTRIBUTION_MERGE_TIME="${MODE_CONTRIBUTION_MERGE_TIME:-1:00:00}"
INVIVO_UMAP_CPUS="${INVIVO_UMAP_CPUS:-20}"
INVIVO_UMAP_MEM="${INVIVO_UMAP_MEM:-64G}"
INVIVO_UMAP_TIME="${INVIVO_UMAP_TIME:-8:00:00}"
INVITRO_UMAP_CPUS="${INVITRO_UMAP_CPUS:-12}"
INVITRO_UMAP_MEM="${INVITRO_UMAP_MEM:-48G}"
INVITRO_UMAP_TIME="${INVITRO_UMAP_TIME:-6:00:00}"
POOLED_UMAP_CPUS="${POOLED_UMAP_CPUS:-20}"
POOLED_UMAP_MEM="${POOLED_UMAP_MEM:-64G}"
POOLED_UMAP_TIME="${POOLED_UMAP_TIME:-8:00:00}"
REPORT_CPUS="${REPORT_CPUS:-2}"
REPORT_MEM="${REPORT_MEM:-16G}"
REPORT_TIME="${REPORT_TIME:-2:00:00}"

for arg in "$@"; do
  case "${arg}" in
    --help|-h) usage; exit 0 ;;
    --project_root=*) PROJECT_ROOT="${arg#*=}" ;;
    --result_root=*) RESULT_ROOT="${arg#*=}" ;;
    --invivo_input=*) INVIVO_INPUT="${arg#*=}" ;;
    --invitro_input=*) INVITRO_INPUT="${arg#*=}" ;;
    --r_module=*) R_MODULE="${arg#*=}" ;;
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
    --mode_contribution_array_concurrency=*) MODE_CONTRIBUTION_ARRAY_CONCURRENCY="${arg#*=}" ;;
    --run_parameter_contribution=*) RUN_PARAMETER_CONTRIBUTION="${arg#*=}" ;;
    --run_umap=*) RUN_UMAP="${arg#*=}" ;;
    --run_report=*) RUN_REPORT="${arg#*=}" ;;
    --reductions=*) REDUCTIONS="${arg#*=}" ;;
    --preprocess_modes=*) PREPROCESS_MODES="${arg#*=}" ;;
    --pooled_preprocess_modes=*) POOLED_PREPROCESS_MODES="${arg#*=}" ;;
    --run_full_tsne=*) RUN_FULL_TSNE="${arg#*=}" ;;
    --dry_run=*) DRY_RUN="${arg#*=}" ;;
    --cpus=*) INVIVO_TABLE_CPUS="${arg#*=}" ;;
    --mem=*) INVIVO_TABLE_MEM="${arg#*=}" ;;
    --time=*) INVIVO_TABLE_TIME="${arg#*=}" ;;
    --invivo_table_cpus=*) INVIVO_TABLE_CPUS="${arg#*=}" ;;
    --invivo_table_mem=*) INVIVO_TABLE_MEM="${arg#*=}" ;;
    --invivo_table_time=*) INVIVO_TABLE_TIME="${arg#*=}" ;;
    --invitro_table_cpus=*) INVITRO_TABLE_CPUS="${arg#*=}" ;;
    --invitro_table_mem=*) INVITRO_TABLE_MEM="${arg#*=}" ;;
    --invitro_table_time=*) INVITRO_TABLE_TIME="${arg#*=}" ;;
    --mode_contribution_cpus=*) MODE_CONTRIBUTION_CPUS="${arg#*=}" ;;
    --mode_contribution_mem=*) MODE_CONTRIBUTION_MEM="${arg#*=}" ;;
    --mode_contribution_time=*) MODE_CONTRIBUTION_TIME="${arg#*=}" ;;
    --mode_contribution_merge_cpus=*) MODE_CONTRIBUTION_MERGE_CPUS="${arg#*=}" ;;
    --mode_contribution_merge_mem=*) MODE_CONTRIBUTION_MERGE_MEM="${arg#*=}" ;;
    --mode_contribution_merge_time=*) MODE_CONTRIBUTION_MERGE_TIME="${arg#*=}" ;;
    --invivo_umap_cpus=*) INVIVO_UMAP_CPUS="${arg#*=}" ;;
    --invivo_umap_mem=*) INVIVO_UMAP_MEM="${arg#*=}" ;;
    --invivo_umap_time=*) INVIVO_UMAP_TIME="${arg#*=}" ;;
    --invitro_umap_cpus=*) INVITRO_UMAP_CPUS="${arg#*=}" ;;
    --invitro_umap_mem=*) INVITRO_UMAP_MEM="${arg#*=}" ;;
    --invitro_umap_time=*) INVITRO_UMAP_TIME="${arg#*=}" ;;
    --pooled_umap_cpus=*) POOLED_UMAP_CPUS="${arg#*=}" ;;
    --pooled_umap_mem=*) POOLED_UMAP_MEM="${arg#*=}" ;;
    --pooled_umap_time=*) POOLED_UMAP_TIME="${arg#*=}" ;;
    --report_cpus=*) REPORT_CPUS="${arg#*=}" ;;
    --report_mem=*) REPORT_MEM="${arg#*=}" ;;
    --report_time=*) REPORT_TIME="${arg#*=}" ;;
    *)
      echo "Unknown argument: ${arg}" >&2
      usage >&2
      exit 2
      ;;
  esac
done

PROJECT_ROOT="$(resolve_hpc_path "${PWD}" "${PROJECT_ROOT}")"
RESULT_ROOT="$(resolve_hpc_path "${PROJECT_ROOT}" "${RESULT_ROOT}")"
MODE_REFERENCE_O2_VALUES="$(normalize_reference_o2_values "${MODE_REFERENCE_O2_VALUES}" "${MODE_REFERENCE_O2}")"
REFERENCE_O2_COUNT="$(csv_count "${MODE_REFERENCE_O2_VALUES}")"
if [[ "${MODE_CONTRIBUTION_ARRAY_CONCURRENCY}" -lt 1 ]]; then
  MODE_CONTRIBUTION_ARRAY_CONCURRENCY=1
fi

LOG_DIR="${RESULT_ROOT}/logs"
STAMP="$(date +%Y%m%d_%H%M%S)"
WORKFLOW_PREFIX="${LOG_DIR}/parameter_landscape_${STAMP}"
MANIFEST="${WORKFLOW_PREFIX}_jobs.tsv"
mkdir -p "${LOG_DIR}"

write_task_preamble() {
  local script_path="$1"
  local job_name="$2"
  local cpus="$3"
  local mem="$4"
  local time_limit="$5"
  local output_pattern="$6"
  {
    printf '%s\n' '#!/usr/bin/env bash'
    printf '#SBATCH --job-name=%s\n' "${job_name}"
    printf '%s\n' '#SBATCH --ntasks=1'
    printf '#SBATCH --cpus-per-task=%s\n' "${cpus}"
    printf '#SBATCH --mem=%s\n' "${mem}"
    printf '#SBATCH --time=%s\n' "${time_limit}"
    if [[ -n "${QOS}" ]]; then printf '#SBATCH --qos=%s\n' "${QOS}"; fi
    if [[ -n "${PARTITION}" ]]; then printf '#SBATCH --partition=%s\n' "${PARTITION}"; fi
    if [[ -n "${ACCOUNT}" ]]; then printf '#SBATCH --account=%s\n' "${ACCOUNT}"; fi
    printf '#SBATCH --output=%s/%s.out\n' "${LOG_DIR}" "${output_pattern}"
    printf '#SBATCH --error=%s/%s.err\n' "${LOG_DIR}" "${output_pattern}"
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
    printf 'RUN_REPORT=%q\n' "${RUN_REPORT}"
    printf 'REDUCTIONS=%q\n' "${REDUCTIONS}"
    printf 'PREPROCESS_MODES=%q\n' "${PREPROCESS_MODES}"
    printf 'POOLED_PREPROCESS_MODES=%q\n' "${POOLED_PREPROCESS_MODES}"
    printf 'RUN_FULL_TSNE=%q\n' "${RUN_FULL_TSNE}"
    cat <<'BATCH_PREAMBLE'

SCRIPT_DIR="oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

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

fixed_o2_slug() {
  Rscript --vanilla -e 'source("oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/parameter_landscape_utils.R"); cat(fixed_o2_o2_slug(as.numeric(commandArgs(TRUE)[[1]])))' "$1"
}

cd "${PROJECT_ROOT}"
mkdir -p "${RESULT_ROOT}/logs"
MAIN_LOG="${RESULT_ROOT}/logs/${SLURM_JOB_NAME}_${SLURM_JOB_ID:-manual}.log"
exec > >(tee -a "${MAIN_LOG}") 2>&1

echo "[$(date '+%F %T')] Starting ${SLURM_JOB_NAME}"
echo "Project root: ${PROJECT_ROOT}"
echo "Result root: ${RESULT_ROOT}"
echo "Threads: ${THREADS}"

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
  MODE_TABLE_ARGS+=("--mode_reference_o2_values=${MODE_REFERENCE_O2_VALUES}")
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
BATCH_PREAMBLE
  } > "${script_path}"
}

write_invivo_tables_script() {
  local path="$1"
  write_task_preamble "${path}" "o2pl_A_invivo_tables" "${INVIVO_TABLE_CPUS}" "${INVIVO_TABLE_MEM}" "${INVIVO_TABLE_TIME}" "o2pl_A_invivo_tables_%j"
  cat <<'BATCH_BODY' >> "${path}"

if [[ ! -d "${INVIVO_INPUT}" ]]; then
  echo "Missing in vivo input directory: ${INVIVO_INPUT}" >&2
  exit 1
fi

run_rscript "Write in vivo UMAP tables and fixed-O2 mode tables" \
  "${SCRIPT_DIR}/parameter_landscape_clustering.R" \
  "--run_parts=invivo_tables" \
  "--invivo_input=${INVIVO_INPUT}" \
  "--result_root=${RESULT_ROOT}" \
  "--overwrite_modes=${OVERWRITE_MODES}" \
  "--n_workers=${THREADS}" \
  "${MAX_SEEDS_ARGS[@]}" \
  "${MODE_TABLE_ARGS[@]}"

require_file "${RESULT_ROOT}/invivo/UMAPs/Tables/invivo_best_params_by_seed.csv"
require_file "${RESULT_ROOT}/FixO2Modes/fixed_o2_attractor_mode_by_seed_o2.tsv"
echo "[$(date '+%F %T')] Completed ${SLURM_JOB_NAME}"
BATCH_BODY
}

write_invitro_tables_script() {
  local path="$1"
  write_task_preamble "${path}" "o2pl_B_invitro_tables" "${INVITRO_TABLE_CPUS}" "${INVITRO_TABLE_MEM}" "${INVITRO_TABLE_TIME}" "o2pl_B_invitro_tables_%j"
  cat <<'BATCH_BODY' >> "${path}"

if [[ ! -d "${INVITRO_INPUT}" ]]; then
  echo "Missing in vitro input directory: ${INVITRO_INPUT}" >&2
  exit 1
fi

run_rscript "Write in vitro UMAP tables" \
  "${SCRIPT_DIR}/parameter_landscape_clustering.R" \
  "--run_parts=invitro_tables" \
  "--invitro_input=${INVITRO_INPUT}" \
  "--result_root=${RESULT_ROOT}" \
  "${MAX_SEEDS_ARGS[@]}"

require_file "${RESULT_ROOT}/invitro/UMAPs/Tables/invitro_best_params_by_seed.csv"
echo "[$(date '+%F %T')] Completed ${SLURM_JOB_NAME}"
BATCH_BODY
}

write_mode_contribution_array_script() {
  local path="$1"
  write_task_preamble "${path}" "o2pl_C_mode_contrib" "${MODE_CONTRIBUTION_CPUS}" "${MODE_CONTRIBUTION_MEM}" "${MODE_CONTRIBUTION_TIME}" "o2pl_C_mode_contrib_%A_%a"
  cat <<'BATCH_BODY' >> "${path}"

IFS=',' read -r -a REFERENCE_O2_ARRAY <<< "${MODE_REFERENCE_O2_VALUES}"
TASK_INDEX="${SLURM_ARRAY_TASK_ID:-1}"
if [[ "${TASK_INDEX}" -lt 1 || "${TASK_INDEX}" -gt "${#REFERENCE_O2_ARRAY[@]}" ]]; then
  echo "Invalid SLURM_ARRAY_TASK_ID=${TASK_INDEX}; available reference O2 count=${#REFERENCE_O2_ARRAY[@]}" >&2
  exit 1
fi
TASK_O2="${REFERENCE_O2_ARRAY[$((TASK_INDEX - 1))]}"
TASK_O2="${TASK_O2//[[:space:]]/}"

run_rscript "Estimate parameter contributions to fixed-O2 mode separation at O2=${TASK_O2}" \
  "${SCRIPT_DIR}/parameter_landscape_clustering.R" \
  "--run_parts=mode_contribution_reference" \
  "--result_root=${RESULT_ROOT}" \
  "--best_csv=${RESULT_ROOT}/invivo/UMAPs/Tables/invivo_best_params_by_seed.csv" \
  "--mode_tables_dir=${RESULT_ROOT}/FixO2Modes" \
  "--mode_reference_o2=${TASK_O2}" \
  "--n_bootstrap=${MODE_CONTRIBUTION_BOOTSTRAP}" \
  "--overwrite_parameter_contribution=${OVERWRITE_PARAMETER_CONTRIBUTION}" \
  "${MAX_SEEDS_ARGS[@]}"

TASK_SLUG="$(fixed_o2_slug "${TASK_O2}")"
require_file "${RESULT_ROOT}/mode_parameter_contribution/reference_o2_${TASK_SLUG}/mode_parameter_feature_importance.csv"
require_file "${RESULT_ROOT}/mode_parameter_contribution/reference_o2_${TASK_SLUG}/mode_parameter_contribution_summary.tsv"
echo "[$(date '+%F %T')] Completed ${SLURM_JOB_NAME} O2=${TASK_O2}"
BATCH_BODY
}

write_mode_contribution_merge_script() {
  local path="$1"
  write_task_preamble "${path}" "o2pl_C2_mode_merge" "${MODE_CONTRIBUTION_MERGE_CPUS}" "${MODE_CONTRIBUTION_MERGE_MEM}" "${MODE_CONTRIBUTION_MERGE_TIME}" "o2pl_C2_mode_merge_%j"
  cat <<'BATCH_BODY' >> "${path}"

run_rscript "Merge mode parameter contribution outputs across reference O2 values" \
  "${SCRIPT_DIR}/parameter_landscape_clustering.R" \
  "--run_parts=mode_contribution_merge" \
  "--result_root=${RESULT_ROOT}" \
  "--mode_reference_o2=${MODE_REFERENCE_O2}" \
  "--mode_reference_o2_values=${MODE_REFERENCE_O2_VALUES}"

require_file "${RESULT_ROOT}/mode_parameter_contribution/mode_parameter_contribution_index.csv"
require_file "${RESULT_ROOT}/mode_parameter_contribution/mode_parameter_top_features_across_reference_o2.csv"
echo "[$(date '+%F %T')] Completed ${SLURM_JOB_NAME}"
BATCH_BODY
}

write_invivo_umap_script() {
  local path="$1"
  write_task_preamble "${path}" "o2pl_D_invivo_reduce" "${INVIVO_UMAP_CPUS}" "${INVIVO_UMAP_MEM}" "${INVIVO_UMAP_TIME}" "o2pl_D_invivo_reduce_%j"
  cat <<'BATCH_BODY' >> "${path}"

run_rscript "Generate in vivo parameter UMAP/PCA/t-SNE figures" \
  "${SCRIPT_DIR}/parameter_landscape_clustering.R" \
  "--run_parts=invivo_reductions" \
  "--result_root=${RESULT_ROOT}" \
  "--invivo_input=${INVIVO_INPUT}" \
  "--reductions=${REDUCTIONS}" \
  "--preprocess_modes=${PREPROCESS_MODES}" \
  "--run_full_tsne=${RUN_FULL_TSNE}" \
  "--mode_tables_dir=${RESULT_ROOT}/FixO2Modes" \
  "--mode_reference_o2=${MODE_REFERENCE_O2}" \
  "--umap_seed=123" \
  "--cluster_seed=123" \
  "--n_threads=${THREADS}" \
  "${ATTRACTOR_FEATURE_ARGS[@]}"

echo "[$(date '+%F %T')] Completed ${SLURM_JOB_NAME}"
BATCH_BODY
}

write_invitro_umap_script() {
  local path="$1"
  write_task_preamble "${path}" "o2pl_E_invitro_reduce" "${INVITRO_UMAP_CPUS}" "${INVITRO_UMAP_MEM}" "${INVITRO_UMAP_TIME}" "o2pl_E_invitro_reduce_%j"
  cat <<'BATCH_BODY' >> "${path}"

run_rscript "Generate in vitro parameter UMAP/PCA/t-SNE figures" \
  "${SCRIPT_DIR}/parameter_landscape_clustering.R" \
  "--run_parts=invitro_reductions" \
  "--result_root=${RESULT_ROOT}" \
  "--invitro_input=${INVITRO_INPUT}" \
  "--reductions=${REDUCTIONS}" \
  "--preprocess_modes=${PREPROCESS_MODES}" \
  "--run_full_tsne=${RUN_FULL_TSNE}" \
  "--umap_seed=123" \
  "--cluster_seed=123" \
  "--n_threads=${THREADS}"

echo "[$(date '+%F %T')] Completed ${SLURM_JOB_NAME}"
BATCH_BODY
}

write_pooled_umap_script() {
  local path="$1"
  write_task_preamble "${path}" "o2pl_F_pooled_reduce" "${POOLED_UMAP_CPUS}" "${POOLED_UMAP_MEM}" "${POOLED_UMAP_TIME}" "o2pl_F_pooled_reduce_%j"
  cat <<'BATCH_BODY' >> "${path}"

run_rscript "Generate pooled in vivo/in vitro UMAP/PCA/t-SNE figures and distance tables" \
  "${SCRIPT_DIR}/parameter_landscape_clustering.R" \
  "--run_parts=pooled_reductions" \
  "--result_root=${RESULT_ROOT}" \
  "--invivo_input=${INVIVO_INPUT}" \
  "--invitro_input=${INVITRO_INPUT}" \
  "--reductions=${REDUCTIONS}" \
  "--pooled_preprocess_modes=${POOLED_PREPROCESS_MODES}" \
  "--run_full_tsne=${RUN_FULL_TSNE}" \
  "--umap_seed=123" \
  "--cluster_seed=123" \
  "--n_threads=${THREADS}"

echo "[$(date '+%F %T')] Completed ${SLURM_JOB_NAME}"
BATCH_BODY
}

write_report_script() {
  local path="$1"
  write_task_preamble "${path}" "o2pl_G_report" "${REPORT_CPUS}" "${REPORT_MEM}" "${REPORT_TIME}" "o2pl_G_report_%j"
  cat <<'BATCH_BODY' >> "${path}"

run_rscript "Render parameter landscape reports" \
  "${SCRIPT_DIR}/parameter_landscape_clustering.R" \
  "--run_parts=reports" \
  "--result_root=${RESULT_ROOT}" \
  "--run_umap_cluster_report=${RUN_UMAP}" \
  "--run_mode_contribution_report=${RUN_PARAMETER_CONTRIBUTION}" \
  "--umap_report_html=${RESULT_ROOT}/parameter_landscape_clustering_umap_cluster_report.html" \
  "--mode_contribution_report_html=${RESULT_ROOT}/mode_parameter_contribution/mode_parameter_contribution_report.html" \
  "--embed_assets=TRUE"

if truthy "${RUN_UMAP}"; then
  require_file "${RESULT_ROOT}/parameter_landscape_clustering_umap_cluster_report.html"
fi
if truthy "${RUN_PARAMETER_CONTRIBUTION}"; then
  require_file "${RESULT_ROOT}/mode_parameter_contribution/mode_parameter_contribution_report.html"
fi
echo "[$(date '+%F %T')] Completed ${SLURM_JOB_NAME}"
BATCH_BODY
}

submit_job() {
  local task="$1"
  local script_path="$2"
  local dependency="${3:-}"
  local array_spec="${4:-}"
  local cmd=(sbatch --parsable)
  if [[ -n "${dependency}" ]]; then cmd+=(--dependency="${dependency}"); fi
  if [[ -n "${array_spec}" ]]; then cmd+=(--array="${array_spec}"); fi
  cmd+=("${script_path}")
  if truthy "${DRY_RUN}"; then
    printf 'DRY_RUN submit %s:' "${task}" >&2
    printf ' %q' "${cmd[@]}" >&2
    printf '\n' >&2
    printf "DRYRUN_%s" "${task}"
    return 0
  fi
  local output
  output="$("${cmd[@]}")"
  output="${output%%;*}"
  printf "%s" "${output}"
}

record_job() {
  local task="$1"
  local job_id="$2"
  local dependency="$3"
  local script_path="$4"
  printf '%s\t%s\t%s\t%s\n' "${task}" "${job_id}" "${dependency}" "${script_path}" >> "${MANIFEST}"
}

SCRIPT_A="${WORKFLOW_PREFIX}_A_invivo_tables_modes.sbatch"
SCRIPT_B="${WORKFLOW_PREFIX}_B_invitro_tables.sbatch"
SCRIPT_C="${WORKFLOW_PREFIX}_C_mode_contribution_array.sbatch"
SCRIPT_C2="${WORKFLOW_PREFIX}_C2_mode_contribution_merge.sbatch"
SCRIPT_D="${WORKFLOW_PREFIX}_D_invivo_umap.sbatch"
SCRIPT_E="${WORKFLOW_PREFIX}_E_invitro_umap.sbatch"
SCRIPT_F="${WORKFLOW_PREFIX}_F_pooled_umap.sbatch"
SCRIPT_G="${WORKFLOW_PREFIX}_G_report.sbatch"

write_invivo_tables_script "${SCRIPT_A}"
write_invitro_tables_script "${SCRIPT_B}"
write_mode_contribution_array_script "${SCRIPT_C}"
write_mode_contribution_merge_script "${SCRIPT_C2}"
write_invivo_umap_script "${SCRIPT_D}"
write_invitro_umap_script "${SCRIPT_E}"
write_pooled_umap_script "${SCRIPT_F}"
write_report_script "${SCRIPT_G}"

chmod +x "${WORKFLOW_PREFIX}"_*.sbatch
printf 'task\tjob_id\tdependency\tscript\n' > "${MANIFEST}"

echo "Wrote workflow scripts under: ${LOG_DIR}"
echo "Workflow manifest: ${MANIFEST}"
echo "Reference O2 values: ${MODE_REFERENCE_O2_VALUES}"

if ! truthy "${DRY_RUN}" && ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not found. Run this submitter on the HPC login node, or rerun with --dry_run=TRUE." >&2
  exit 1
fi

JID_A="$(submit_job "A_invivo_tables_modes" "${SCRIPT_A}")"
record_job "A_invivo_tables_modes" "${JID_A}" "" "${SCRIPT_A}"
echo "Submitted A_invivo_tables_modes: ${JID_A}"

JID_B=""
JID_C=""
JID_C2=""
JID_D=""
JID_E=""
JID_F=""
JID_G=""

if truthy "${RUN_UMAP}"; then
  JID_B="$(submit_job "B_invitro_tables" "${SCRIPT_B}")"
  record_job "B_invitro_tables" "${JID_B}" "" "${SCRIPT_B}"
  echo "Submitted B_invitro_tables: ${JID_B}"
fi

if truthy "${RUN_PARAMETER_CONTRIBUTION}"; then
  ARRAY_SPEC="1-${REFERENCE_O2_COUNT}%${MODE_CONTRIBUTION_ARRAY_CONCURRENCY}"
  DEP_C="afterok:${JID_A}"
  JID_C="$(submit_job "C_mode_contribution_array" "${SCRIPT_C}" "${DEP_C}" "${ARRAY_SPEC}")"
  record_job "C_mode_contribution_array" "${JID_C}" "${DEP_C};array=${ARRAY_SPEC}" "${SCRIPT_C}"
  echo "Submitted C_mode_contribution_array: ${JID_C}"

  DEP_C2="afterok:${JID_C}"
  JID_C2="$(submit_job "C2_mode_contribution_merge" "${SCRIPT_C2}" "${DEP_C2}")"
  record_job "C2_mode_contribution_merge" "${JID_C2}" "${DEP_C2}" "${SCRIPT_C2}"
  echo "Submitted C2_mode_contribution_merge: ${JID_C2}"
fi

if truthy "${RUN_UMAP}"; then
  DEP_D="afterok:${JID_A}"
  JID_D="$(submit_job "D_invivo_umap" "${SCRIPT_D}" "${DEP_D}")"
  record_job "D_invivo_umap" "${JID_D}" "${DEP_D}" "${SCRIPT_D}"
  echo "Submitted D_invivo_umap: ${JID_D}"

  DEP_E="afterok:${JID_B}"
  JID_E="$(submit_job "E_invitro_umap" "${SCRIPT_E}" "${DEP_E}")"
  record_job "E_invitro_umap" "${JID_E}" "${DEP_E}" "${SCRIPT_E}"
  echo "Submitted E_invitro_umap: ${JID_E}"

  DEP_F="afterok:${JID_A}:${JID_B}"
  JID_F="$(submit_job "F_pooled_umap" "${SCRIPT_F}" "${DEP_F}")"
  record_job "F_pooled_umap" "${JID_F}" "${DEP_F}" "${SCRIPT_F}"
  echo "Submitted F_pooled_umap: ${JID_F}"
fi

if truthy "${RUN_REPORT}"; then
  REPORT_DEPS=()
  if [[ -n "${JID_C2}" ]]; then REPORT_DEPS+=("${JID_C2}"); fi
  if [[ -n "${JID_D}" ]]; then REPORT_DEPS+=("${JID_D}"); fi
  if [[ -n "${JID_E}" ]]; then REPORT_DEPS+=("${JID_E}"); fi
  if [[ -n "${JID_F}" ]]; then REPORT_DEPS+=("${JID_F}"); fi
  if [[ "${#REPORT_DEPS[@]}" -gt 0 ]]; then
  DEP_G="afterok:$(IFS=:; echo "${REPORT_DEPS[*]}")"
  JID_G="$(submit_job "G_report" "${SCRIPT_G}" "${DEP_G}")"
  record_job "G_report" "${JID_G}" "${DEP_G}" "${SCRIPT_G}"
  echo "Submitted G_report: ${JID_G}"
  else
    echo "Skipping G_report: no report-producing upstream tasks were enabled."
  fi
else
  echo "Skipping G_report because RUN_REPORT=${RUN_REPORT}"
fi

echo
echo "Parameter landscape dependency workflow submitted."
echo "Manifest: ${MANIFEST}"
echo "Result root: ${RESULT_ROOT}"
