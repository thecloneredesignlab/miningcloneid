#!/bin/bash
# Run extra_results.R for one O2 result directory when the summary is missing.

set -euo pipefail

O2SD_EARLY_PROJECT_ROOT="${PROJECT_ROOT:-}"
for arg in "$@"; do
  case "${arg}" in
    --project_root=*) O2SD_EARLY_PROJECT_ROOT="${arg#*=}"; break ;;
  esac
done
if [[ -z "${O2SD_DOCKER_HPC_ROOT:-}" && -n "${O2SD_EARLY_PROJECT_ROOT}" && -d "${O2SD_EARLY_PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc" ]]; then
  O2SD_DOCKER_HPC_ROOT="${O2SD_EARLY_PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/Docker/hpc"
fi
O2SD_DOCKER_HPC_ROOT="${O2SD_DOCKER_HPC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck source=../util/o2_supply_demand_map_apptainer_runtime.sh
source "${O2SD_DOCKER_HPC_ROOT}/util/o2_supply_demand_map_apptainer_runtime.sh"

if [[ -n "${O2SD_EARLY_PROJECT_ROOT}" && -f "${O2SD_EARLY_PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_shell_utils.sh" ]]; then
  O2SD_SHELL_UTILS="${O2SD_EARLY_PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_shell_utils.sh"
else
  O2SD_SHELL_UTILS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../util" && pwd)/o2_supply_demand_map_shell_utils.sh"
fi
# shellcheck source=../../../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"

if [[ "${O2_POSTPROCESS_LOGIN_SHELL:-0}" != "1" ]]; then
  export O2_POSTPROCESS_LOGIN_SHELL=1
  exec bash -lc "$(shell_join bash "$0" "$@")"
fi

usage() {
  cat <<'EOF'
Usage:
  bash postprocess_extra_results.sh --run_dir=/path/to/run_dir [options]

Options:
  --project_root=/path/to/Rescomposite
  --extra_results_script=/path/to/extra_results.R
  --r_module=R/4.4              Compatibility option; ignored by Docker/hpc.
  --force_extra_results=TRUE|FALSE
  --dry_run=TRUE|FALSE
  --help

Behavior:
  If run_dir/extra_results/seed_summary.tsv already exists, the script exits
  without rerunning extra_results.R unless --force_extra_results=TRUE.
  Under Docker/hpc, Rscript is executed from the configured Apptainer SIF.
EOF
}

parse_args() {
  for arg in "$@"; do
    case "${arg}" in
      --help|-h)
        usage
        exit 0
        ;;
      --run_dir=*)
        RUN_DIR="${arg#*=}"
        ;;
      --project_root=*)
        PROJECT_ROOT="${arg#*=}"
        ;;
      --extra_results_script=*)
        EXTRA_RESULTS_SCRIPT="${arg#*=}"
        ;;
      --r_module=*)
        R_MODULE="${arg#*=}"
        ;;
      --force_extra_results=*)
        FORCE_EXTRA_RESULTS="${arg#*=}"
        ;;
      --dry_run=*)
        DRY_RUN="${arg#*=}"
        ;;
      *)
        echo "Unknown argument: ${arg}" >&2
        usage >&2
        exit 2
        ;;
    esac
  done
}

DEFAULT_PROJECT_ROOT="/share/lab_crd/lab_crd/taoli/Project/Rescomposite"
DEFAULT_R_MODULE="R/4.4"
DEFAULT_FORCE_EXTRA_RESULTS="FALSE"
DEFAULT_DRY_RUN="FALSE"

PROJECT_ROOT="${PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
R_MODULE="${R_MODULE:-${DEFAULT_R_MODULE}}"
FORCE_EXTRA_RESULTS="${FORCE_EXTRA_RESULTS:-${DEFAULT_FORCE_EXTRA_RESULTS}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"
RUN_DIR="${RUN_DIR:-}"
EXTRA_RESULTS_SCRIPT="${EXTRA_RESULTS_SCRIPT:-}"

parse_args "$@"

if [[ -z "${RUN_DIR}" ]]; then
  echo "--run_dir is required." >&2
  usage >&2
  exit 2
fi

PROJECT_ROOT="$(cd "${PROJECT_ROOT}" && pwd)"
if [[ ! -d "${RUN_DIR}" ]]; then
  echo "Missing run directory: ${RUN_DIR}" >&2
  exit 1
fi
RUN_DIR="$(cd "${RUN_DIR}" && pwd)"
if [[ -z "${EXTRA_RESULTS_SCRIPT}" ]]; then
  EXTRA_RESULTS_SCRIPT="${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/analysis/fit_results/extra_results.R"
fi
EXTRA_RESULTS_SCRIPT="$(cd "$(dirname "${EXTRA_RESULTS_SCRIPT}")" && pwd)/$(basename "${EXTRA_RESULTS_SCRIPT}")"

if [[ ! -f "${EXTRA_RESULTS_SCRIPT}" ]]; then
  echo "Missing extra_results script: ${EXTRA_RESULTS_SCRIPT}" >&2
  exit 1
fi

SUMMARY_PATH="${RUN_DIR}/extra_results/seed_summary.tsv"
LOG_PATH="${RUN_DIR}/extra_results_run.log"

if [[ -f "${SUMMARY_PATH}" ]] && ! truthy "${FORCE_EXTRA_RESULTS}"; then
  echo "extra_results already exists; skipping: ${SUMMARY_PATH}"
  exit 0
fi

cmd=(Rscript "${EXTRA_RESULTS_SCRIPT}" "--run_dir=${RUN_DIR}")
if truthy "${DRY_RUN}"; then
  echo "DRY_RUN=TRUE; not running extra_results.R"
  printf "Command:"
  printf " %q" "${cmd[@]}"
  printf "\n"
  exit 0
fi

load_r_module
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Container-backed Rscript not found; SIF: ${O2SD_CONTAINER_IMAGE}." >&2
  exit 1
fi

cd "${PROJECT_ROOT}"
echo "Running extra_results.R for ${RUN_DIR}"
echo "Log: ${LOG_PATH}"
"${cmd[@]}" > "${LOG_PATH}" 2>&1

if [[ ! -f "${SUMMARY_PATH}" ]]; then
  echo "extra_results.R finished but did not create ${SUMMARY_PATH}" >&2
  tail -40 "${LOG_PATH}" >&2 || true
  exit 1
fi

echo "extra_results complete: ${SUMMARY_PATH}"
