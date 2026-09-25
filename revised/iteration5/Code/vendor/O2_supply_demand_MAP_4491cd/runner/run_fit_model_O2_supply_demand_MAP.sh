#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
FIT_SCRIPT="${WORKFLOW_ROOT}/optimizer/fit_model_O2_supply_demand_MAP.R"
O2SD_SHELL_UTILS="${WORKFLOW_ROOT}/util/o2_supply_demand_map_shell_utils.sh"
# shellcheck source=../util/o2_supply_demand_map_shell_utils.sh
source "${O2SD_SHELL_UTILS}"

if [[ -z "${O2SD_RUN_COMMAND:-}" ]]; then
  export O2SD_RUN_COMMAND="$(shell_join bash "${BASH_SOURCE[0]}" "$@")"
fi
if [[ -z "${O2SD_ENTRYPOINT_SCRIPT:-}" ]]; then
  export O2SD_ENTRYPOINT_SCRIPT="${BASH_SOURCE[0]}"
fi

exec Rscript "${FIT_SCRIPT}" "$@"
