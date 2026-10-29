#!/usr/bin/env bash
# Authorized Figure7 current campaign only: no Figure8/historical generation.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "${script_dir}/../config/figure78_hpc_paths.env"
data_root=;submit=no;task_ids=;postprocess=yes
for arg in "$@"; do
  case "$arg" in
    --data-root=*) data_root="${arg#*=}" ;;
    --submit=yes) submit=yes ;;
    --task-ids=*) task_ids="${arg#*=}" ;;
    --postprocess=no) postprocess=no ;;
    *) echo "Use --data-root=ABSOLUTE_PREPARED_ROOT --submit=yes [--task-ids=IDS] [--postprocess=no]" >&2; exit 2 ;;
  esac
done
[[ "$submit" == yes && "$data_root" == /* ]] || exit 2
count="$(python3 "${script_dir}/../Figures/figure78_entry_tools.py" check-prepared \
  "--workspace=$FIGURE78_WORKSPACE" --figure=7 "--data-root=$data_root")"
[[ "$count" =~ ^[1-9][0-9]*$ ]] || exit 2
array="1-${count}"
if [[ -n "$task_ids" ]]; then
  [[ "$task_ids" =~ ^[0-9]+(,[0-9]+)*$ ]] || exit 2
  IFS=, read -ra figure7_ids <<< "$task_ids"
  for id in "${figure7_ids[@]}"; do (( id >= 1 && id <= count )) || exit 2; done
  array="$task_ids"
fi
export FIGURE78_FIGURE=7 FIGURE78_DATA_ROOT="$data_root"
mkdir -p "${data_root}/Figure7/logs"
# Deliberately NO %N concurrency cap. Scheduling belongs to Slurm/xxlarge.
job="$(sbatch --parsable --partition=red --qos=xxlarge --cpus-per-task=1 --mem=16G --time=12:00:00 \
  --array="$array" --export=ALL --job-name=figure7_native_last5 \
  --output="${data_root}/Figure7/logs/array_%A_%a.out" \
  "${script_dir}/run_current_Figure78_array.sub" "$FIGURE78_WORKSPACE")"
job="${job%%;*}"
[[ "$job" =~ ^[0-9]+$ ]] || exit 2
echo "ARRAY_JOB_ID=$job"
if [[ "$postprocess" == yes ]]; then
  post="$(sbatch --parsable --partition=red --qos=xxlarge --cpus-per-task=1 --mem=32G --time=12:00:00 \
    --dependency="afterok:${job}" --export=ALL --job-name=figure7_native_validate \
    --output="${data_root}/Figure7/logs/post_%j.out" \
    "${script_dir}/run_Figure7_current_post.sh" "$FIGURE78_WORKSPACE" "$data_root")"
  echo "POSTPROCESS_JOB_ID=${post%%;*}"
fi
