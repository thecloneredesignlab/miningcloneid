#!/usr/bin/env bash

set -euo pipefail

ITERATION_ROOT="/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration4"
REPO_ROOT="/share/lab_crd/taoli/Project/HypoxiaLTEEFigures"
RUN_ID="supp6_17_18_full_20260914_2033"
MAX_CONCURRENT=256
QOS="xxlarge"
ARRAY_TIME="12:00:00"
ARRAY_MEM="8G"
FINALIZE_MEM="128G"

for argument in "$@"; do
  case "${argument}" in
    --run-id=*) RUN_ID="${argument#*=}" ;;
    --max-concurrent=*) MAX_CONCURRENT="${argument#*=}" ;;
    --qos=*) QOS="${argument#*=}" ;;
    --array-time=*) ARRAY_TIME="${argument#*=}" ;;
    --array-mem=*) ARRAY_MEM="${argument#*=}" ;;
    --finalize-mem=*) FINALIZE_MEM="${argument#*=}" ;;
    -h|--help)
      echo "Usage: $0 [--run-id=ID] [--max-concurrent=N] [--qos=xxlarge]"
      exit 0 ;;
    *) echo "Unknown option: ${argument}" >&2; exit 2 ;;
  esac
done

[[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$ ]] || {
  echo "Invalid run id: ${RUN_ID}" >&2; exit 2
}
[[ "${MAX_CONCURRENT}" =~ ^[1-9][0-9]*$ ]] || {
  echo "--max-concurrent must be a positive integer." >&2; exit 2
}
command -v sbatch >/dev/null 2>&1 || { echo "sbatch is required." >&2; exit 2; }

FIGURE6_DATA_ROOT="${ITERATION_ROOT}/data/Figures/Figure6/fixed_pmisseg_v1"
RUN_ROOT="${FIGURE6_DATA_ROOT}/net_growth_full_range_q10_runs/${RUN_ID}"
MANIFEST="${RUN_ROOT}/net_growth_task_manifest.tsv"
STATUS_PATH="${ITERATION_ROOT}/audit/hpc_supp_figure6_17_18/${RUN_ID}/status.tsv"
WORKER="${ITERATION_ROOT}/Code/hpc/run_supp_figure6_17_18_array.sub"
FINALIZER="${ITERATION_ROOT}/Code/hpc/finalize_supp_figure6_17_18_array.sub"
for path in "${MANIFEST}" "${WORKER}" "${FINALIZER}"; do
  [[ -r "${path}" ]] || { echo "Missing required file: ${path}" >&2; exit 2; }
done

if [[ -d "${ITERATION_ROOT}/audit/locks/supp_figure6_17_18.lock" ]]; then
  echo "The monolithic Figure 6 runner still owns its lock; refusing overlap." >&2
  exit 2
fi

CAMPAIGN_ID="array_resume_$(date '+%Y%m%d_%H%M%S')"
CAMPAIGN_ROOT="${ITERATION_ROOT}/audit/hpc_supp_figure6_17_18_array/${RUN_ID}/${CAMPAIGN_ID}"
mkdir -p "${CAMPAIGN_ROOT}/logs"
TASK_FILE="${CAMPAIGN_ROOT}/missing_tasks.tsv"
TEMP_TASK_FILE="${TASK_FILE}.tmp.$$"
printf 'array_index\ttask_id\tcache_path\tcheckpoint_path\n' > "${TEMP_TASK_FILE}"
index=0
while IFS=$'\t' read -r task_id model_context propagation_mode pair_label \
    p_misseg o2_chunk_index o2_index_start o2_index_end endpoint_indices \
    n_unique_endpoint represented_optimizer_endpoint cache_path; do
  [[ "${task_id}" == "task_id" ]] && continue
  [[ "${task_id}" =~ ^NG[0-9]{4}$ ]] || {
    echo "Malformed manifest task id: ${task_id}" >&2; exit 2
  }
  if [[ ! -f "${cache_path}" ]]; then
    index=$((index + 1))
    printf '%s\t%s\t%s\t%s\n' "${index}" "${task_id}" "${cache_path}" \
      "${cache_path}.checkpoint.rds" >> "${TEMP_TASK_FILE}"
  fi
done < "${MANIFEST}"
mv "${TEMP_TASK_FILE}" "${TASK_FILE}"

total_tasks="$(( $(wc -l < "${MANIFEST}") - 1 ))"
missing_tasks="${index}"
completed_tasks="$((total_tasks - missing_tasks))"
[[ "${total_tasks}" -eq 4020 ]] || {
  echo "Expected 4020 manifest tasks; found ${total_tasks}." >&2; exit 2
}
cp "${STATUS_PATH}" "${CAMPAIGN_ROOT}/status_before_array.tsv"
sha256sum "${MANIFEST}" "${TASK_FILE}" > "${CAMPAIGN_ROOT}/input_sha256.tsv"
{
  printf 'key\tvalue\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'campaign_id\t%s\n' "${CAMPAIGN_ID}"
  printf 'git_head\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse HEAD)"
  printf 'total_tasks\t%s\n' "${total_tasks}"
  printf 'completed_tasks_at_submission\t%s\n' "${completed_tasks}"
  printf 'missing_tasks_at_submission\t%s\n' "${missing_tasks}"
  printf 'qos\t%s\n' "${QOS}"
  printf 'max_concurrent\t%s\n' "${MAX_CONCURRENT}"
  printf 'array_time\t%s\n' "${ARRAY_TIME}"
  printf 'array_mem\t%s\n' "${ARRAY_MEM}"
  printf 'submitted_at\t%s\n' "$(date -Iseconds)"
} > "${CAMPAIGN_ROOT}/campaign.tsv"

export_list="ALL,RUN_ID=${RUN_ID},CAMPAIGN_ROOT=${CAMPAIGN_ROOT}"
prewarm_id="$(sbatch --parsable --qos="${QOS}" --job-name=s6ng_prewarm \
  --cpus-per-task=1 --mem=16G --time=01:00:00 \
  --output="${CAMPAIGN_ROOT}/logs/prewarm_%j.log" \
  --export="${export_list},ARRAY_MODE=prewarm" "${WORKER}")"

array_id=""
dependency="afterok:${prewarm_id}"
if (( missing_tasks > 0 )); then
  array_id="$(sbatch --parsable --qos="${QOS}" --job-name=s6ng_tasks \
    --dependency="${dependency}" --array="1-${missing_tasks}%${MAX_CONCURRENT}" \
    --cpus-per-task=1 --mem="${ARRAY_MEM}" --time="${ARRAY_TIME}" \
    --output="${CAMPAIGN_ROOT}/logs/task_%A_%a.log" \
    --export="${export_list},ARRAY_MODE=task,TASK_FILE=${TASK_FILE}" "${WORKER}")"
  dependency="afterany:${array_id}"
fi

finalize_id="$(sbatch --parsable --qos="${QOS}" --job-name=s6ng_finalize \
  --dependency="${dependency}" --cpus-per-task=1 --mem="${FINALIZE_MEM}" \
  --time=12:00:00 --output="${CAMPAIGN_ROOT}/logs/finalize_%j.log" \
  --export="${export_list}" "${FINALIZER}")"

{
  printf 'job_role\tjob_id\tdependency\n'
  printf 'prewarm\t%s\t\n' "${prewarm_id}"
  if [[ -n "${array_id}" ]]; then
    printf 'array\t%s\tafterok:%s\n' "${array_id}" "${prewarm_id}"
  fi
  printf 'finalize\t%s\t%s\n' "${finalize_id}" "${dependency}"
} > "${CAMPAIGN_ROOT}/jobs.tsv"

{
  printf 'run_id\tstatus\texit_code\tstage\thost\tn_core\tgit_head\tupdated_at\n'
  printf '%s\tRUNNING\t0\tSLURM_ARRAY_RESUME\tRED\t1\t%s\t%s\n' \
    "${RUN_ID}" "$(git -C "${REPO_ROOT}" rev-parse HEAD)" "$(date -Iseconds)"
} > "${STATUS_PATH}"

echo "campaign_root=${CAMPAIGN_ROOT}"
echo "completed_tasks=${completed_tasks}"
echo "missing_tasks=${missing_tasks}"
echo "prewarm_job_id=${prewarm_id}"
echo "array_job_id=${array_id}"
echo "finalize_job_id=${finalize_id}"
