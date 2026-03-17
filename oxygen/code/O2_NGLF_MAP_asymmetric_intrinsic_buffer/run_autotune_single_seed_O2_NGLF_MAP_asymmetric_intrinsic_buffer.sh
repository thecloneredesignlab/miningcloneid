#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/autotune_single_seed_O2_NGLF_MAP_asymmetric_intrinsic_buffer.R"

exec Rscript "${R_SCRIPT}" "$@"
