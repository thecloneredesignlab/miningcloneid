#!/usr/bin/env bash
# HPC wrapper for warm_up_joint_fitting_results_extra.R.

set -euo pipefail

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

usage() {
  cat <<'EOF'
Usage:
  bash run_warm_up_joint_fitting_results_extra_hpc.sh [options] [-- extra R args]

Defaults are set for the current HPC joint multi-warmup result.

Options:
  --project_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling
  --input_root=/path/to/fit_joint_multi_warmup_result
  --output_root=/path/to/analysis_output
  --script=/path/to/warm_up_joint_fitting_results_extra.R
  --stage=prepare|embedding|curve|summary|all
  --embedding_param_set=invivo|shared|invivo,shared
  --reductions=pca,umap,tsne
  --run_full_tsne=TRUE|FALSE
  --prepare_workers=64
  --n_workers=64
  --umap_threads=32
  --overwrite=TRUE|FALSE
  --r_module=R/4.4
  --log_dir=/path/to/log_dir
  --log_path=/path/to/exact_log_file
  --background=TRUE|FALSE
  --dry_run=TRUE|FALSE
  --help

Examples:
  # Run the full embedding refresh in the foreground.
  bash run_warm_up_joint_fitting_results_extra_hpc.sh

  # Run detached from SSH. The PID and log path are printed immediately.
  bash run_warm_up_joint_fitting_results_extra_hpc.sh --background=TRUE

  # Only regenerate the missing full in vivo PCA embedding.
  bash run_warm_up_joint_fitting_results_extra_hpc.sh \
    --embedding_param_set=invivo \
    --reductions=pca \
    --run_full_tsne=FALSE
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

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling"
DEFAULT_RESULT_NAME="fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_200seed_20260708_173447"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
INPUT_ROOT=""
OUTPUT_ROOT=""
ANALYSIS_SCRIPT=""
STAGE="embedding"
EMBEDDING_PARAM_SET="invivo,shared"
REDUCTIONS="pca,umap,tsne"
RUN_FULL_TSNE="TRUE"
PREPARE_WORKERS="64"
N_WORKERS="64"
UMAP_THREADS="32"
OVERWRITE="TRUE"
R_MODULE="R/4.4"
LOG_DIR=""
LOG_PATH=""
BACKGROUND="FALSE"
DRY_RUN="FALSE"
EXTRA_R_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      usage
      exit 0
      ;;
    --)
      shift
      EXTRA_R_ARGS+=("$@")
      break
      ;;
    --project_root=*) PROJECT_ROOT="${1#*=}" ;;
    --input_root=*) INPUT_ROOT="${1#*=}" ;;
    --output_root=*) OUTPUT_ROOT="${1#*=}" ;;
    --script=*) ANALYSIS_SCRIPT="${1#*=}" ;;
    --stage=*) STAGE="${1#*=}" ;;
    --embedding_param_set=*|--embedding_param_sets=*) EMBEDDING_PARAM_SET="${1#*=}" ;;
    --reductions=*) REDUCTIONS="${1#*=}" ;;
    --run_full_tsne=*) RUN_FULL_TSNE="${1#*=}" ;;
    --prepare_workers=*) PREPARE_WORKERS="${1#*=}" ;;
    --n_workers=*) N_WORKERS="${1#*=}" ;;
    --umap_threads=*) UMAP_THREADS="${1#*=}" ;;
    --overwrite=*) OVERWRITE="${1#*=}" ;;
    --r_module=*) R_MODULE="${1#*=}" ;;
    --log_dir=*) LOG_DIR="${1#*=}" ;;
    --log_path=*) LOG_PATH="${1#*=}" ;;
    --background=*) BACKGROUND="${1#*=}" ;;
    --dry_run=*) DRY_RUN="${1#*=}" ;;
    --*)
      EXTRA_R_ARGS+=("$1")
      ;;
    *)
      echo "Unknown positional argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

if [[ -z "${INPUT_ROOT}" ]]; then
  INPUT_ROOT="${PROJECT_ROOT}/oxygen/results/${DEFAULT_RESULT_NAME}"
fi
if [[ -z "${OUTPUT_ROOT}" ]]; then
  OUTPUT_ROOT="${PROJECT_ROOT}/oxygen/results/analysis/warm_up_joint_fitting_results_extra/$(basename "${INPUT_ROOT}")"
fi
if [[ -z "${ANALYSIS_SCRIPT}" ]]; then
  ANALYSIS_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/warm_up_joint_fitting_results_extra/warm_up_joint_fitting_results_extra.R"
fi
if [[ -z "${LOG_DIR}" ]]; then
  LOG_DIR="${OUTPUT_ROOT}_hpc_logs"
fi
if [[ -z "${LOG_PATH}" ]]; then
  stage_slug="${STAGE//,/__}"
  stage_slug="${stage_slug//[^A-Za-z0-9_.-]/_}"
  LOG_PATH="${LOG_DIR}/warm_up_joint_fitting_results_extra_${stage_slug}_$(date +%Y%m%d_%H%M%S).log"
fi

CMD=(
  Rscript "${ANALYSIS_SCRIPT}"
  "--stage=${STAGE}"
  "--input_root=${INPUT_ROOT}"
  "--output_root=${OUTPUT_ROOT}"
  "--embedding_param_set=${EMBEDDING_PARAM_SET}"
  "--reductions=${REDUCTIONS}"
  "--run_full_tsne=${RUN_FULL_TSNE}"
  "--prepare_workers=${PREPARE_WORKERS}"
  "--n_workers=${N_WORKERS}"
  "--umap_threads=${UMAP_THREADS}"
  "--overwrite=${OVERWRITE}"
)
CMD+=("${EXTRA_R_ARGS[@]}")

run_analysis() {
  echo "HOST=$(hostname)"
  echo "START=$(date)"
  echo "PROJECT_ROOT=${PROJECT_ROOT}"
  echo "INPUT_ROOT=${INPUT_ROOT}"
  echo "OUTPUT_ROOT=${OUTPUT_ROOT}"
  echo "SCRIPT=${ANALYSIS_SCRIPT}"
  echo "STAGE=${STAGE}"
  echo "EMBEDDING_PARAM_SET=${EMBEDDING_PARAM_SET}"
  echo "REDUCTIONS=${REDUCTIONS}"
  echo "RUN_FULL_TSNE=${RUN_FULL_TSNE}"
  echo "COMMIT=$(git -C "${PROJECT_ROOT}" rev-parse HEAD 2>/dev/null || echo NA)"
  echo "Command: $(shell_join "${CMD[@]}")"

  cd "${PROJECT_ROOT}"
  load_r_module
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript not found after loading ${R_MODULE}." >&2
    exit 1
  fi
  Rscript --version
  "${CMD[@]}"
  echo "DONE=$(date)"
}

mkdir -p "$(dirname "${LOG_PATH}")"

echo "Log path: ${LOG_PATH}"
echo "Command: $(shell_join "${CMD[@]}")"

if truthy "${DRY_RUN}"; then
  exit 0
fi

if [[ ! -d "${PROJECT_ROOT}" ]]; then
  echo "Missing project root: ${PROJECT_ROOT}" >&2
  exit 1
fi
if [[ ! -f "${ANALYSIS_SCRIPT}" ]]; then
  echo "Missing analysis script: ${ANALYSIS_SCRIPT}" >&2
  exit 1
fi
if [[ ! -d "${INPUT_ROOT}" ]]; then
  echo "Missing input root: ${INPUT_ROOT}" >&2
  exit 1
fi

if truthy "${BACKGROUND}"; then
  child_args=(
    "--project_root=${PROJECT_ROOT}"
    "--input_root=${INPUT_ROOT}"
    "--output_root=${OUTPUT_ROOT}"
    "--script=${ANALYSIS_SCRIPT}"
    "--stage=${STAGE}"
    "--embedding_param_set=${EMBEDDING_PARAM_SET}"
    "--reductions=${REDUCTIONS}"
    "--run_full_tsne=${RUN_FULL_TSNE}"
    "--prepare_workers=${PREPARE_WORKERS}"
    "--n_workers=${N_WORKERS}"
    "--umap_threads=${UMAP_THREADS}"
    "--overwrite=${OVERWRITE}"
    "--r_module=${R_MODULE}"
    "--log_path=${LOG_PATH}"
    "--background=FALSE"
  )
  child_args+=("${EXTRA_R_ARGS[@]}")
  nohup env O2_WARMUP_JOINT_EXTRA_NO_TEE=1 bash "$0" "${child_args[@]}" > "${LOG_PATH}" 2>&1 </dev/null &
  pid=$!
  disown "${pid}" 2>/dev/null || true
  echo "PID=${pid}"
  echo "LOG=${LOG_PATH}"
  exit 0
fi

if [[ "${O2_WARMUP_JOINT_EXTRA_NO_TEE:-0}" == "1" ]]; then
  run_analysis
else
  run_analysis 2>&1 | tee "${LOG_PATH}"
fi
