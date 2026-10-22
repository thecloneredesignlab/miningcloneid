#!/usr/bin/env bash
# Future separately authorized submission only. Never executed by static checks.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "${script_dir}/../config/figure78_hpc_paths.env"
figure=;data_root=;submit=no
for arg in "$@"; do
  case "$arg" in
    --figure=7|--figure=8) figure="${arg#*=}" ;;
    --data-root=*) data_root="${arg#*=}" ;;
    --submit=yes) submit=yes ;;
    *) echo "Use --figure=7|8 --data-root=ABSOLUTE_PREPARED_ROOT --submit=yes; old positional/task arguments are forbidden" >&2; exit 2 ;;
  esac
done
[[ "$submit" == yes && -n "$figure" && "$data_root" == /* ]] || { echo "Submission requires explicit --submit=yes, figure and prepared data root" >&2; exit 2; }
count="$(python3 "${script_dir}/../Figures/figure78_entry_tools.py" check-prepared \
  "--workspace=$FIGURE78_WORKSPACE" "--figure=$figure" "--data-root=$data_root")"
[[ "$count" =~ ^[1-9][0-9]*$ ]] || exit 2
export FIGURE78_FIGURE="$figure" FIGURE78_DATA_ROOT="$data_root"
exec sbatch --parsable --partition=red --qos=xxlarge --cpus-per-task=1 --mem=8G --time=12:00:00 \
  --array="1-${count}" --export=ALL --job-name="figure${figure}_current" \
  "${script_dir}/run_current_Figure78_array.sub" "$FIGURE78_WORKSPACE"
