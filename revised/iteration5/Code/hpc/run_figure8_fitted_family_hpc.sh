#!/usr/bin/env bash
# Thin compatibility routing; old positional arguments are rejected.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
exec bash "${script_dir}/figure78_legacy_hpc.sh" "run_figure8_fitted_family_hpc.sh" "$@"
