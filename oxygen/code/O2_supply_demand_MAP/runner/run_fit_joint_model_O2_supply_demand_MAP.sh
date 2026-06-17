#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UNIFIED_RUNNER="${SCRIPT_DIR}/run_fit_model_O2_supply_demand_MAP.sh"

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

if [[ -z "${O2SD_RUN_COMMAND:-}" ]]; then
  export O2SD_RUN_COMMAND="$(shell_join bash "${BASH_SOURCE[0]}" "$@")"
fi
export O2SD_ENTRYPOINT_SCRIPT="${BASH_SOURCE[0]}"

has_mode=false
for arg in "$@"; do
  case "${arg}" in
    --mode|--mode=*)
      has_mode=true
      break
      ;;
  esac
done

if "${has_mode}"; then
  exec bash "${UNIFIED_RUNNER}" --fit_joint "$@"
else
  exec bash "${UNIFIED_RUNNER}" --fit_joint --mode=run "$@"
fi
