#!/usr/bin/env bash

set -euo pipefail

WORKSPACE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_ROOT="${WORKSPACE_ROOT}/Code/Figures"

usage() {
  cat <<'EOF'
Usage:
  manager.sh \
    --invitro-result-root=PATH \
    --invivo-result-root=PATH \
    --joint-result-root=PATH \
    --gemcitabine-data-root=PATH \
    --ltee-data-root=PATH \
    [--n-core=N] \
    [--recompute-fixed-o2] \
    [--recompute-invivo-tsne] \
    [--rebuild-figure6-grid] \
    [--model-dependent-only] \
    [--resume-after-figure5f-de] \
    [--fresh-run]

Required scientific-input paths:
  --invitro-result-root=PATH
      Directory containing the separate in-vitro fit results.
  --invivo-result-root=PATH
      Directory containing the separate in-vivo fit results.
  --joint-result-root=PATH
      Directory containing the joint multi-warmup fit results.
  --gemcitabine-data-root=PATH
      Directory containing the in-vivo Gemcitabine input data.
  --ltee-data-root=PATH
      Directory containing the InVitroData_LTEE input tables.

Optional controls:
  --n-core=N
      Positive integer worker count. Default: 8.
  --recompute-fixed-o2
      Recompute the fixed-O2 analysis instead of reusing staged results.
  --recompute-invivo-tsne
      Recompute the in-vivo t-SNE landscape instead of reusing staged results.
  --rebuild-figure6-grid
      Rebuild all Figure 6 multi-seed response caches.
  --model-dependent-only
      Preserve existing Figure 1-3 outputs and regenerate Figure 4-6 plus
      their parent-indexed supplementary figures.
  --resume-after-figure5f-de
      Resume a validated interrupted run after the Figure 5F DE initial-
      population product has been written and checked.
  --fresh-run
      Move prior generated data/figures into an audit archive before running,
      then regenerate from empty output directories.
  -h, --help
      Print this help text and exit.

Outputs:
  The manager regenerates Figure 1-6 and parent-indexed supplementary figures, publishes figure
  artifacts, validates scientific inputs against the fixed Code/config MD5
  baseline, hashes generated intermediates, verifies published-figure MD5
  identity, and generates the embedded manuscript HTML report. The manager
  does not compile the TeX manuscript or create/replace its PDF.

Example:
  bash /absolute/path/to/workspace/manager.sh \
    --invitro-result-root=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_invitro_unified_500seed_r442_exact_20260825_032031 \
    --invivo-result-root=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_invivo_unified_500seed_r442_exact_20260825_032031 \
    --joint-result-root=/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_joint_invivo_clusters_global_invitro_best_500seed_r442_exact_20260826_033633 \
    --gemcitabine-data-root=/path/to/data/InVivoData_Gemcitabine \
    --ltee-data-root=/path/to/data/InVitroData_LTEE \
    --n-core=8
EOF
}

for argument in "$@"; do
  case "${argument}" in
    -h|--help)
      usage
      exit 0
      ;;
  esac
done

if [[ ! -d "${CODE_ROOT}" ]]; then
  echo "Portable workspace is missing Code/Figures: ${CODE_ROOT}" >&2
  exit 2
fi

export FIGURE_WORKSPACE_ROOT="${WORKSPACE_ROOT}"
export FIGURE_MODEL_CODE_ROOT="${FIGURE_MODEL_CODE_ROOT:-/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/code/O2_supply_demand_MAP}"
export OMP_NUM_THREADS=1
export KMP_USE_SHM=0
export RCPP_PARALLEL_NUM_THREADS=1

REPO_ROOT="${WORKSPACE_ROOT}"
AUDIT_ROOT="${WORKSPACE_ROOT}/audit"
MANUSCRIPT_ROOT="${WORKSPACE_ROOT}/manuscript"
REPORT_TOOL="${WORKSPACE_ROOT}/Code/vendor/O2_supply_demand_MAP_4491cd/report/text_to_html_report.py"
MODEL_CODE_ROOT="${FIGURE_MODEL_CODE_ROOT}"

N_CORE=8
RECOMPUTE_FIXED_O2=FALSE
RECOMPUTE_INVIVO_TSNE=FALSE
REBUILD_FIGURE6_GRID=FALSE
MODEL_DEPENDENT_ONLY=FALSE
RESUME_AFTER_FIGURE5F_DE=FALSE
FRESH_RUN=FALSE
INVITRO_RESULT_ROOT=""
INVIVO_RESULT_ROOT=""
JOINT_RESULT_ROOT=""
GEMCITABINE_DATA_ROOT=""
LTEE_DATA_ROOT=""

for argument in "$@"; do
  case "${argument}" in
    --invitro-result-root=*)
      INVITRO_RESULT_ROOT="${argument#*=}"
      ;;
    --invivo-result-root=*)
      INVIVO_RESULT_ROOT="${argument#*=}"
      ;;
    --joint-result-root=*)
      JOINT_RESULT_ROOT="${argument#*=}"
      ;;
    --gemcitabine-data-root=*)
      GEMCITABINE_DATA_ROOT="${argument#*=}"
      ;;
    --ltee-data-root=*)
      LTEE_DATA_ROOT="${argument#*=}"
      ;;
    --n-core=*)
      N_CORE="${argument#*=}"
      ;;
    --recompute-fixed-o2)
      RECOMPUTE_FIXED_O2=TRUE
      ;;
    --recompute-invivo-tsne)
      RECOMPUTE_INVIVO_TSNE=TRUE
      ;;
    --rebuild-figure6-grid)
      REBUILD_FIGURE6_GRID=TRUE
      ;;
    --model-dependent-only)
      MODEL_DEPENDENT_ONLY=TRUE
      ;;
    --resume-after-figure5f-de)
      RESUME_AFTER_FIGURE5F_DE=TRUE
      ;;
    --fresh-run)
      FRESH_RUN=TRUE
      ;;
    *)
      echo "Unknown option: ${argument}" >&2
      echo "Run 'manager.sh --help' for usage." >&2
      exit 2
      ;;
  esac
done

missing_options=()
[[ -n "${INVITRO_RESULT_ROOT}" ]] ||
  missing_options+=("--invitro-result-root")
[[ -n "${INVIVO_RESULT_ROOT}" ]] ||
  missing_options+=("--invivo-result-root")
[[ -n "${JOINT_RESULT_ROOT}" ]] ||
  missing_options+=("--joint-result-root")
[[ -n "${GEMCITABINE_DATA_ROOT}" ]] ||
  missing_options+=("--gemcitabine-data-root")
[[ -n "${LTEE_DATA_ROOT}" ]] ||
  missing_options+=("--ltee-data-root")

if (( ${#missing_options[@]} > 0 )); then
  echo "Missing required option(s): ${missing_options[*]}" >&2
  echo "Run 'manager.sh --help' for usage." >&2
  exit 2
fi

if ! [[ "${N_CORE}" =~ ^[1-9][0-9]*$ ]]; then
  echo "--n-core must be a positive integer." >&2
  exit 2
fi

normalize_input_directory() {
  local option_name="$1"
  local input_path="$2"
  if [[ ! -d "${input_path}" ]]; then
    echo "${option_name} is not an existing directory: ${input_path}" >&2
    return 1
  fi
  (
    cd "${input_path}"
    pwd -P
  )
}

INVITRO_RESULT_ROOT="$(
  normalize_input_directory \
    "--invitro-result-root" "${INVITRO_RESULT_ROOT}"
)"
INVIVO_RESULT_ROOT="$(
  normalize_input_directory \
    "--invivo-result-root" "${INVIVO_RESULT_ROOT}"
)"
JOINT_RESULT_ROOT="$(
  normalize_input_directory \
    "--joint-result-root" "${JOINT_RESULT_ROOT}"
)"
GEMCITABINE_DATA_ROOT="$(
  normalize_input_directory \
    "--gemcitabine-data-root" "${GEMCITABINE_DATA_ROOT}"
)"
LTEE_DATA_ROOT="$(
  normalize_input_directory \
    "--ltee-data-root" "${LTEE_DATA_ROOT}"
)"

RUNTIME_PATH_ARGS=(
  "--invitro-result-root=${INVITRO_RESULT_ROOT}"
  "--invivo-result-root=${INVIVO_RESULT_ROOT}"
  "--joint-result-root=${JOINT_RESULT_ROOT}"
  "--gemcitabine-data-root=${GEMCITABINE_DATA_ROOT}"
  "--ltee-data-root=${LTEE_DATA_ROOT}"
)

export HYPOXIA_REPO_ROOT="${REPO_ROOT}"
export FIGURE_INVITRO_RESULT_ROOT="${INVITRO_RESULT_ROOT}"
export FIGURE_INVIVO_RESULT_ROOT="${INVIVO_RESULT_ROOT}"
export FIGURE_JOINT_RESULT_ROOT="${JOINT_RESULT_ROOT}"
export FIGURE_GEMCITABINE_DATA_ROOT="${GEMCITABINE_DATA_ROOT}"
export FIGURE_LTEE_DATA_ROOT="${LTEE_DATA_ROOT}"

if [[ ! -f "${REPORT_TOOL}" ]]; then
  echo "Missing manuscript HTML report tool: ${REPORT_TOOL}" >&2
  exit 2
fi

if [[ ! -f "${MODEL_CODE_ROOT}/model/model_O2_supply_demand_MAP.R" ]] ||
   [[ ! -f "${MODEL_CODE_ROOT}/model/model_O2_supply_demand_MAP.cpp" ]]; then
  echo "Missing required external model implementation: ${MODEL_CODE_ROOT}" >&2
  exit 2
fi

Rscript \
  "${CODE_ROOT}/util/runtime/workspace_paths.R" \
  "${RUNTIME_PATH_ARGS[@]}"

mkdir -p \
  "${AUDIT_ROOT}/logs" \
  "${AUDIT_ROOT}/md5" \
  "${AUDIT_ROOT}/reports" \
  "${MANUSCRIPT_ROOT}/Figures"

RUN_ID="$(date '+%Y%m%d_%H%M%S')"
RUN_LOG="${AUDIT_ROOT}/logs/manager_${RUN_ID}.log"
LATEST_LOG="${AUDIT_ROOT}/logs/manager_latest.log"

if [[ "${FRESH_RUN}" == "TRUE" ]]; then
  FRESH_ARCHIVE="${AUDIT_ROOT}/pre_fresh_run_${RUN_ID}"
  mkdir -p "${FRESH_ARCHIVE}"
  archive_generated_path() {
    local source_path="$1"
    local archive_name="$2"
    if [[ -e "${source_path}" ]]; then
      mv "${source_path}" "${FRESH_ARCHIVE}/${archive_name}"
    fi
  }
  archive_generated_path "${WORKSPACE_ROOT}/data/Figures" "data_Figures"
  archive_generated_path "${WORKSPACE_ROOT}/Figures" "Figures"
  archive_generated_path "${MANUSCRIPT_ROOT}/Figures" "manuscript_Figures"
  mkdir -p \
    "${WORKSPACE_ROOT}/data/Figures" \
    "${WORKSPACE_ROOT}/Figures" \
    "${MANUSCRIPT_ROOT}/Figures"
fi

{
  echo "figure workspace manager start: $(date -Iseconds)"
  echo "workspace_root=${WORKSPACE_ROOT}"
  echo "repository_root=${REPO_ROOT}"
  echo "model_code_root=${MODEL_CODE_ROOT}"
  echo "invitro_result_root=${INVITRO_RESULT_ROOT}"
  echo "invivo_result_root=${INVIVO_RESULT_ROOT}"
  echo "joint_result_root=${JOINT_RESULT_ROOT}"
  echo "gemcitabine_data_root=${GEMCITABINE_DATA_ROOT}"
  echo "ltee_data_root=${LTEE_DATA_ROOT}"
  echo "n_core=${N_CORE}"
  echo "recompute_fixed_o2=${RECOMPUTE_FIXED_O2}"
  echo "recompute_invivo_tsne=${RECOMPUTE_INVIVO_TSNE}"
  echo "rebuild_figure6_grid=${REBUILD_FIGURE6_GRID}"
  echo "model_dependent_only=${MODEL_DEPENDENT_ONLY}"
  echo "resume_after_figure5f_de=${RESUME_AFTER_FIGURE5F_DE}"
  echo "fresh_run=${FRESH_RUN}"
  if [[ "${FRESH_RUN}" == "TRUE" ]]; then
    echo "fresh_run_archive=${FRESH_ARCHIVE}"
  fi

  python3 \
    "${CODE_ROOT}/util/workflow/input_validator.py" \
    --phase=preflight

  Rscript "${CODE_ROOT}/run_all_figures.R" \
    "${RUNTIME_PATH_ARGS[@]}" \
    "--n-core=${N_CORE}" \
    "--recompute-fixed-o2=${RECOMPUTE_FIXED_O2}" \
    "--recompute-invivo-tsne=${RECOMPUTE_INVIVO_TSNE}" \
    "--rebuild-figure6-grid=${REBUILD_FIGURE6_GRID}" \
    "--model-dependent-only=${MODEL_DEPENDENT_ONLY}" \
    "--resume-after-figure5f-de=${RESUME_AFTER_FIGURE5F_DE}"

  find "${WORKSPACE_ROOT}" -type f \
    \( -name '.DS_Store' -o -name '*.pyc' \) -delete
  find "${WORKSPACE_ROOT}" -type d -name '__pycache__' -empty -delete

  python3 "${CODE_ROOT}/util/workflow/artifact_publisher.py"
  python3 "${CODE_ROOT}/util/workflow/artifact_md5_validator.py"
  find "${WORKSPACE_ROOT}" -type f \
    \( -name '.DS_Store' -o -name '*.pyc' \) -delete
  find "${WORKSPACE_ROOT}" -type d -name '__pycache__' -empty -delete

  python3 "${REPORT_TOOL}" \
    --input "${MANUSCRIPT_ROOT}/ltee_hypoxia_model.tex" \
    --output "${MANUSCRIPT_ROOT}/ltee_hypoxia_model.html" \
    --asset-mode embed \
    --figure-placement source \
    --hide-source

  find "${WORKSPACE_ROOT}" -type f \
    \( -name '.DS_Store' -o -name '*.pyc' \) -delete
  find "${WORKSPACE_ROOT}" -type d -name '__pycache__' -empty -delete

  python3 \
    "${CODE_ROOT}/util/workflow/input_validator.py" \
    --phase=postflight
  echo "figure workspace manager complete: $(date -Iseconds)"
} 2>&1 | tee "${RUN_LOG}"

cp "${RUN_LOG}" "${LATEST_LOG}"
