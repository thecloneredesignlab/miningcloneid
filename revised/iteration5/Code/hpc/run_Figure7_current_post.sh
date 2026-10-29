#!/usr/bin/env bash
# Numerical collection/acceptance only. No publication, rendering, or submission.
set -euo pipefail
workspace="${1:?Absolute workspace}"
data_root="${2:?Absolute prepared data root}"
[[ "$workspace" == /* && "$data_root" == /* ]] || exit 2
export FIGURE78_WORKSPACE="$workspace" FIGURE78_FIGURE=7 FIGURE78_DATA_ROOT="$data_root"
source "${workspace}/Code/config/figure78_hpc_paths.env"
python3 "${workspace}/Code/Figures/figure78_entry_tools.py" check-prepared \
  "--workspace=$workspace" --figure=7 "--data-root=$data_root"
for stage in aggregate validate; do
  "$FIGURE78_RSCRIPT" "${workspace}/Code/Figures/run_Figure7.R" \
    "--stage=$stage" --profile=full "--workspace=$workspace" "--data-root=$data_root"
done
echo "FIGURE7_NUMERICAL_COLLECTION_VALIDATED; publication and rendering remain separately gated"
