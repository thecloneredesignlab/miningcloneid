#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS_ROOT="$(cd "${WORKFLOW_ROOT}/../../results" && pwd)"

UNIFIED_RUNNER="${SCRIPT_DIR}/run_fit_model_O2G_supply_demand_MAP.sh"
STAMP="$(date +%Y%m%d_%H%M%S)"
PIPELINE_ROOT="${RESULTS_ROOT}/fit_O2G_two_stage_${STAMP}"

stage1_args=()
stage2_args=()
stage1_out_dir=""
stage1_glucose=""
stage2_config=""

has_arg_prefix() {
  local needle="$1"
  shift
  local arg
  for arg in "$@"; do
    if [[ "${arg}" == "${needle}"* ]]; then
      return 0
    fi
  done
  return 1
}

for arg in "$@"; do
  case "${arg}" in
    --stage1_out_dir=*)
      stage1_out_dir="${arg#*=}"
      ;;
    --stage1_parameter_table=*)
      stage1_args+=("--parameter_table=${arg#*=}")
      ;;
    --stage1_x_data=*)
      stage1_args+=("--x_data=${arg#*=}")
      ;;
    --stage1_growth_data=*)
      stage1_args+=("--growth_data=${arg#*=}")
      ;;
    --stage1_seed=*)
      stage1_args+=("--seed=${arg#*=}")
      ;;
    --stage1_itermax=*)
      stage1_args+=("--itermax=${arg#*=}")
      ;;
    --stage1_NP=*)
      stage1_args+=("--NP=${arg#*=}")
      ;;
    --stage1_dt=*)
      stage1_args+=("--dt=${arg#*=}")
      ;;
    --stage1_pop_growth_factor=*)
      stage1_args+=("--pop_growth_factor=${arg#*=}")
      ;;
    --stage1_report_passages=*)
      stage1_args+=("--report_passages=${arg#*=}")
      ;;
    --stage1_glucose=*)
      stage1_glucose="${arg#*=}"
      stage1_args+=("--glucose=${stage1_glucose}")
      ;;
    --stage1_*=*)
      echo "Unsupported stage1 option: ${arg}" >&2
      exit 1
      ;;
    --config=*)
      stage2_config="${arg#*=}"
      stage2_args+=("${arg}")
      ;;
    --glucose=*)
      if [[ -z "${stage1_glucose}" ]]; then
        stage1_glucose="${arg#*=}"
        stage1_args+=("--glucose=${stage1_glucose}")
      fi
      stage2_args+=("${arg}")
      ;;
    --parameter_table=*|--parameters=*)
      echo "Do not pass ${arg%%=*} to stage2 in the two-stage runner." >&2
      echo "Use --stage1_parameter_table=... for the stage1 input table; stage2 always receives the stage1 locked table." >&2
      exit 1
      ;;
    *)
      stage2_args+=("${arg}")
      ;;
  esac
done

if [[ -z "${stage1_out_dir}" ]]; then
  stage1_out_dir="${PIPELINE_ROOT}/stage1_invitro"
fi
mkdir -p "${stage1_out_dir}"

if [[ -z "${stage1_glucose}" ]]; then
  cfg_probe="${stage2_config:-${WORKFLOW_ROOT}/../../config/O2G_supply_demand.yaml}"
  if [[ -f "${cfg_probe}" ]]; then
    probed_value="$(awk 'BEGIN{FS=":"} /^[[:space:]]*glucose[[:space:]]*:/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}' "${cfg_probe}")"
    if [[ -n "${probed_value}" ]]; then
      stage1_args+=("--glucose=${probed_value}")
    fi
  fi
fi

if ! has_arg_prefix "--out_root=" "${stage2_args[@]}"; then
  stage2_args+=("--out_root=${PIPELINE_ROOT}/stage2_invivo")
fi
if ! has_arg_prefix "--run_prefix=" "${stage2_args[@]}"; then
  stage2_args+=("--run_prefix=fit_model_O2G_supply_demand_MAP")
fi
if ! has_arg_prefix "--append_run_prefix_timestamp=" "${stage2_args[@]}"; then
  stage2_args+=("--append_run_prefix_timestamp=FALSE")
fi

echo "[O2G two-stage] Stage 1 output: ${stage1_out_dir}"
echo "[O2G two-stage] Running in vitro anoxia fit"
bash "${UNIFIED_RUNNER}" \
  "--fit_invitro" \
  "--out_dir=${stage1_out_dir}" \
  "${stage1_args[@]}"

LOCKED_TABLE="${stage1_out_dir}/parameter_table.invitro_locked.csv"
if [[ ! -f "${LOCKED_TABLE}" ]]; then
  echo "Missing locked parameter table from stage1: ${LOCKED_TABLE}" >&2
  exit 1
fi

echo "[O2G two-stage] Stage 2 locked parameter table: ${LOCKED_TABLE}"
echo "[O2G two-stage] Running in vivo MAP fit"
bash "${UNIFIED_RUNNER}" \
  "--fit_invivo" \
  "--mode=run" \
  "--parameter_table=${LOCKED_TABLE}" \
  "${stage2_args[@]}"
