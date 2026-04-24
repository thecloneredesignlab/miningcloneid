#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
FIT_SCRIPT="${WORKFLOW_ROOT}/optimizer/fit_invivo_model_O2_supply_demand_MAP.R"

exec Rscript "${FIT_SCRIPT}" "--mode=run" "$@"
