#!/usr/bin/env bash
# Current host/container entry. Stages are explicit; this file never submits jobs.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "${script_dir}/../config/figure78_hpc_paths.env"
figure=;stage=;args=()
for arg in "$@"; do
  case "$arg" in
    --figure=7|--figure=8) figure="${arg#*=}" ;;
    --stage=*) stage="${arg#*=}" ;;
    --*=*) args+=("$arg") ;;
    *) echo "Use explicit --figure=7|8 --stage=...; legacy positional arguments are not accepted" >&2; exit 2 ;;
  esac
done
[[ "$figure" == 7 || "$figure" == 8 ]] && [[ -n "$stage" ]] || { echo "Specify --figure and --stage" >&2; exit 2; }
case "$stage" in
  static) exec python3 "${script_dir}/../Figures/figure78_entry_tools.py" static "--workspace=$FIGURE78_WORKSPACE" "${args[@]}" ;;
  plot)
    if [[ "$figure" == 7 && " ${args[*]} " != *" --analysis="* ]]; then args+=(--analysis=all); fi
    exec "$FIGURE78_RSCRIPT" "${script_dir}/../Figures/draw_Figure${figure}.R" "${args[@]}" ;;
  retained-validate) exec "$FIGURE78_RSCRIPT" "${script_dir}/../Figures/validate_Figure78.R" "--figure=$figure" "${args[@]}" ;;
  plan|prepare|simulate|aggregate|validate|all)
    exec "$FIGURE78_RSCRIPT" "${script_dir}/../Figures/run_Figure${figure}.R" "--stage=$stage" "${args[@]}" ;;
  *) echo "Unknown current stage: $stage" >&2; exit 2 ;;
esac
