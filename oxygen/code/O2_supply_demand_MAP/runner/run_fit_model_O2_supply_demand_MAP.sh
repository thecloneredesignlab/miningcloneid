#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
FIT_SCRIPT="${WORKFLOW_ROOT}/optimizer/fit_model_O2_supply_demand_MAP.R"

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
if [[ -z "${O2SD_ENTRYPOINT_SCRIPT:-}" ]]; then
  export O2SD_ENTRYPOINT_SCRIPT="${BASH_SOURCE[0]}"
fi

exec Rscript "${FIT_SCRIPT}" "$@"
