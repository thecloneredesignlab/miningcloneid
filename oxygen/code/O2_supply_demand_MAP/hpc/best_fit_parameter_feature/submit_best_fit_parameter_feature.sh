#!/usr/bin/env bash
# Submit the integrated best-fit-parameter feature workflows to Slurm.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash submit_best_fit_parameter_feature.sh --workflow=all [options]
  bash submit_best_fit_parameter_feature.sh --workflow=fixed_o2 [options]
  bash submit_best_fit_parameter_feature.sh --workflow=parameter_landscape [options]
  bash submit_best_fit_parameter_feature.sh --workflow=dense_grid_monotonicity [options]

Workflow options:
  --workflow=NAME                 all, fixed_o2, parameter_landscape, dense_grid_monotonicity.
  --parameter_parts=PARTS         parameter-landscape parts: clustering, mode_contribution,
                                  dominant_ploidy_contribution, or all.
  --dense_parts=PARTS             Dense-grid run_parts: monotonicity, initial_ploidy, or all.

Path defaults:
  --project_root=DIR              Repo checkout root on HPC.
  --dry_run=TRUE|FALSE            Print child submit commands without submitting.

Other options are forwarded to each selected child submitter. For workflow-specific
options, the existing child submitter option names are still supported.
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

append_workflow() {
  local raw="${1:-all}"
  local part
  local parts
  IFS=',' read -r -a parts <<< "${raw}"
  for part in "${parts[@]}"; do
    part="$(echo "${part}" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d '[:space:]')"
    case "${part}" in
      all|full)
        WORKFLOWS+=("fixed_o2" "parameter_landscape" "dense_grid_monotonicity")
        ;;
      fixed_o2|fixo2|fixed)
        WORKFLOWS+=("fixed_o2")
        ;;
      parameter_landscape|parameter_landscape_clustering|landscape)
        WORKFLOWS+=("parameter_landscape")
        ;;
      dense_grid_monotonicity|dense_grid|dense_grid_monotonicity_classification|monotonicity)
        WORKFLOWS+=("dense_grid_monotonicity")
        ;;
      "")
        ;;
      *)
        echo "Unknown workflow: ${part}" >&2
        exit 2
        ;;
    esac
  done
}

parse_parameter_parts() {
  local raw="${1:-}"
  local part
  local parts
  local seen_clustering="FALSE"
  local seen_mode="FALSE"
  local seen_dominant="FALSE"

  PARAMETER_RUN_UMAP=""
  PARAMETER_RUN_CONTRIBUTION=""
  PARAMETER_CONTRIBUTION_TARGET=""

  [[ -z "${raw}" ]] && return 0
  IFS=',' read -r -a parts <<< "${raw}"
  for part in "${parts[@]}"; do
    part="$(echo "${part}" | tr '[:upper:]' '[:lower:]' | tr '-' '_' | tr -d '[:space:]')"
    case "${part}" in
      all|full)
        seen_clustering="TRUE"
        seen_mode="TRUE"
        seen_dominant="TRUE"
        ;;
      mode|mode_contribution|discrete)
        seen_mode="TRUE"
        ;;
      dominant_ploidy|dominant_ploidy_contribution|continuous)
        seen_dominant="TRUE"
        ;;
      clustering|reduction|reductions|tables)
        seen_clustering="TRUE"
        ;;
      "")
        ;;
      *)
        echo "Unknown parameter_parts value for HPC submission: ${part}" >&2
        exit 2
        ;;
    esac
  done

  if [[ "${seen_clustering}" == "TRUE" ]]; then
    PARAMETER_RUN_UMAP="TRUE"
  else
    PARAMETER_RUN_UMAP="FALSE"
  fi

  if [[ "${seen_mode}" == "TRUE" && "${seen_dominant}" == "TRUE" ]]; then
    PARAMETER_RUN_CONTRIBUTION="TRUE"
    PARAMETER_CONTRIBUTION_TARGET="all"
  elif [[ "${seen_dominant}" == "TRUE" ]]; then
    PARAMETER_RUN_CONTRIBUTION="TRUE"
    PARAMETER_CONTRIBUTION_TARGET="dominant_ploidy"
  elif [[ "${seen_mode}" == "TRUE" ]]; then
    PARAMETER_RUN_CONTRIBUTION="TRUE"
    PARAMETER_CONTRIBUTION_TARGET="mode"
  else
    PARAMETER_RUN_CONTRIBUTION="FALSE"
  fi

  if [[ "${PARAMETER_RUN_UMAP}" != "TRUE" && "${PARAMETER_RUN_CONTRIBUTION}" != "TRUE" ]]; then
    echo "parameter_parts did not select any runnable parameter-landscape step: ${raw}" >&2
    exit 2
  fi
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HPC_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORKFLOW_ROOT="$(cd "${HPC_ROOT}/.." && pwd)"
DEFAULT_PROJECT_ROOT="$(cd "${WORKFLOW_ROOT}/../../.." && pwd)"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
WORKFLOW_RAW="${WORKFLOW:-all}"
PARAMETER_PARTS="${PARAMETER_PARTS:-}"
DENSE_PARTS="${DENSE_PARTS:-}"
RUN_PARTS_VALUE="${RUN_PARTS:-}"
DRY_RUN="${DRY_RUN:-FALSE}"
FORWARD_ARGS=()

for arg in "$@"; do
  case "${arg}" in
    --help|-h)
      usage
      exit 0
      ;;
    --workflow=*) WORKFLOW_RAW="${arg#*=}" ;;
    --parameter_parts=*|--parameter-parts=*) PARAMETER_PARTS="${arg#*=}" ;;
    --dense_parts=*|--dense-parts=*) DENSE_PARTS="${arg#*=}" ;;
    --run_parts=*|--run-parts=*) RUN_PARTS_VALUE="${arg#*=}" ;;
    --project_root=*) PROJECT_ROOT="${arg#*=}"; FORWARD_ARGS+=("${arg}") ;;
    --dry_run=*) DRY_RUN="${arg#*=}"; FORWARD_ARGS+=("${arg}") ;;
    *) FORWARD_ARGS+=("${arg}") ;;
  esac
done

WORKFLOWS=()
append_workflow "${WORKFLOW_RAW}"

FIXED_SUBMIT="${HPC_ROOT}/fix_o2_simulation/submit_fix_o2_simulation_array.sh"
PARAMETER_SUBMIT="${HPC_ROOT}/parameter_landscape/submit_parameter_landscape_full.sh"
DENSE_SUBMIT="${HPC_ROOT}/dense_grid_monotonicity_classification/submit_dense_grid_monotonicity_classification.sh"

for workflow in "${WORKFLOWS[@]}"; do
  case "${workflow}" in
    fixed_o2)
      cmd=(
        bash "${FIXED_SUBMIT}"
        "--project_root=${PROJECT_ROOT}"
        "--run_parts=${RUN_PARTS_VALUE:-analysis}"
        "--analysis_out_dir=oxygen/results/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo_500seed"
        "--fixo2_script=oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo.R"
        "${FORWARD_ARGS[@]}"
      )
      ;;
    parameter_landscape)
      parameter_args=()
      if [[ -n "${PARAMETER_PARTS}" ]]; then
        parse_parameter_parts "${PARAMETER_PARTS}"
        parameter_args+=(
          "--run_umap=${PARAMETER_RUN_UMAP}"
          "--run_parameter_contribution=${PARAMETER_RUN_CONTRIBUTION}"
        )
        if [[ "${PARAMETER_RUN_CONTRIBUTION}" == "TRUE" ]]; then
          parameter_args+=("--mode_contribution_target=${PARAMETER_CONTRIBUTION_TARGET}")
        fi
      fi
      cmd=(
        bash "${PARAMETER_SUBMIT}"
        "--project_root=${PROJECT_ROOT}"
        "--result_root=oxygen/results/analysis/best_fit_parameter_feature/02_parameter_landscape_clustering/parameter_landscape"
        "${parameter_args[@]}"
        "${FORWARD_ARGS[@]}"
      )
      ;;
    dense_grid_monotonicity)
      cmd=(
        bash "${DENSE_SUBMIT}"
        "--project_root=${PROJECT_ROOT}"
        "--result_root=oxygen/results/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/monotonicity_classification"
        "--array_backend=oxygen/code/O2_supply_demand_MAP/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/dense_grid_monotonicity_array_backend.R"
        "${FORWARD_ARGS[@]}"
      )
      if [[ -n "${DENSE_PARTS}" ]]; then
        cmd+=("--run_parts=${DENSE_PARTS}")
      fi
      ;;
    *)
      echo "Internal error: unsupported workflow ${workflow}" >&2
      exit 2
      ;;
  esac
  echo "[$(date '+%F %T')] ${workflow}"
  echo "Command: $(shell_join "${cmd[@]}")"
  "${cmd[@]}"
done
