#!/usr/bin/env bash

set -euo pipefail

LOCAL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=docker_runtime.sh
source "${LOCAL_ROOT}/docker_runtime.sh"

usage() {
  cat <<'EOF'
Usage:
  bash run_full_analysis_docker.sh \
    [--invivo_fit_dir=/absolute/path/to/seed] \
    [--invitro_fit_dir=/absolute/path/to/seed] \
    [--joint_fit_dir=/absolute/path/to/seed] \
    [--run_best_fit_parameter_feature=TRUE|FALSE] \
    [--joint_result_root=/absolute/path/to/multi_warmup_result] \
    [--joint_output_root=/absolute/path/to/output] \
    [--dry_run=TRUE|FALSE]

At least one fit directory or --run_best_fit_parameter_feature=TRUE is
required. Each supplied fit seed runs the canonical simulation -> analysis ->
visualization -> report pipeline. The best-fit feature runner uses
--workflow=all. Joint coupling analysis is added when both joint paths are
supplied.
EOF
}

INVIVO_FIT_DIR=""
INVITRO_FIT_DIR=""
JOINT_FIT_DIR=""
RUN_BEST_FIT_PARAMETER_FEATURE="TRUE"
JOINT_RESULT_ROOT=""
JOINT_OUTPUT_ROOT=""
DRY_RUN="FALSE"

for arg in "$@"; do
  case "${arg}" in
    --help|-h)
      usage
      exit 0
      ;;
    --invivo_fit_dir=*) INVIVO_FIT_DIR="${arg#*=}" ;;
    --invitro_fit_dir=*) INVITRO_FIT_DIR="${arg#*=}" ;;
    --joint_fit_dir=*) JOINT_FIT_DIR="${arg#*=}" ;;
    --run_best_fit_parameter_feature=*) RUN_BEST_FIT_PARAMETER_FEATURE="${arg#*=}" ;;
    --joint_result_root=*) JOINT_RESULT_ROOT="${arg#*=}" ;;
    --joint_output_root=*) JOINT_OUTPUT_ROOT="${arg#*=}" ;;
    --dry_run=*) DRY_RUN="${arg#*=}" ;;
    *)
      echo "Unknown argument: ${arg}" >&2
      usage >&2
      exit 2
      ;;
  esac
done

as_true() {
  case "$(printf '%s' "$1" | tr '[:lower:]' '[:upper:]')" in
    TRUE|T|YES|Y|1) return 0 ;;
    *) return 1 ;;
  esac
}

if [[ -z "${INVIVO_FIT_DIR}" && -z "${INVITRO_FIT_DIR}" \
  && -z "${JOINT_FIT_DIR}" ]] && ! as_true "${RUN_BEST_FIT_PARAMETER_FEATURE}"; then
  echo "No analysis stage was selected." >&2
  usage >&2
  exit 2
fi

run_postfit() {
  local fit_dir="$1"
  local scope="$2"
  o2sd_docker_exec Rscript \
    "${O2SD_WORKFLOW_ROOT}/runner/run_postfit_pipeline.R" \
    "--fit_dir=${fit_dir}" \
    "--scope=${scope}" \
    "--dry_run=${DRY_RUN}"
}

[[ -z "${INVIVO_FIT_DIR}" ]] || run_postfit "${INVIVO_FIT_DIR}" invivo
[[ -z "${INVITRO_FIT_DIR}" ]] || run_postfit "${INVITRO_FIT_DIR}" invitro
[[ -z "${JOINT_FIT_DIR}" ]] || run_postfit "${JOINT_FIT_DIR}" joint

if as_true "${RUN_BEST_FIT_PARAMETER_FEATURE}"; then
  o2sd_docker_exec Rscript \
    "${O2SD_WORKFLOW_ROOT}/runner/best_fit_parameter_feature/runner.R" \
    --workflow=all \
    "--dry_run=${DRY_RUN}"
fi

if [[ -n "${JOINT_RESULT_ROOT}" || -n "${JOINT_OUTPUT_ROOT}" ]]; then
  if [[ -z "${JOINT_RESULT_ROOT}" || -z "${JOINT_OUTPUT_ROOT}" ]]; then
    echo "--joint_result_root and --joint_output_root must be supplied together." >&2
    exit 2
  fi
  o2sd_docker_exec Rscript \
    "${O2SD_WORKFLOW_ROOT}/runner/joint_coupling/run_joint_coupling_pipeline.R" \
    "--result_root=${JOINT_RESULT_ROOT}" \
    "--output_root=${JOINT_OUTPUT_ROOT}" \
    "--dry_run=${DRY_RUN}"
fi
