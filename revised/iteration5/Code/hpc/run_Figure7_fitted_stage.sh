#!/usr/bin/env bash
set -euo pipefail
workspace="$1";campaign="$2";stage="$3"
source "${workspace}/Code/config/figure78_hpc_paths.env"
unset FIGURE78_DATA_ROOT FIGURE78_FIGURE
case "$stage" in
  prepare|smoke|acceptance|aggregate|validate) ;;
  collect) for part in aggregate validate; do
    "$FIGURE78_RSCRIPT" "${workspace}/Code/Figures/run_Figure7_fitted_p_misseg.R" \
      "--stage=$part" "--workspace=$workspace" "--campaign-root=$campaign"
    done
    exit 0 ;;
  *) echo "Unsupported explicit stage" >&2;exit 2 ;;
esac
exec "$FIGURE78_RSCRIPT" "${workspace}/Code/Figures/run_Figure7_fitted_p_misseg.R" \
  "--stage=$stage" "--workspace=$workspace" "--campaign-root=$campaign"
