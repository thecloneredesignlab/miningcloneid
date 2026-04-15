#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
FIT_SCRIPT="${WORKFLOW_ROOT}/optimizer/fit_model_O2G_supply_demand_MAP.R"

exec Rscript "${FIT_SCRIPT}" "$@"
