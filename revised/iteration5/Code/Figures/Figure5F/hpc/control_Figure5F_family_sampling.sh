#!/usr/bin/env bash
set -euo pipefail

iteration_root="${FIGURE5F_ITERATION_ROOT:?FIGURE5F_ITERATION_ROOT is required}"
result_root="${FIGURE5F_RESULT_ROOT:?FIGURE5F_RESULT_ROOT is required}"
family="${FIGURE5F_FAMILY:?FIGURE5F_FAMILY is required}"
current_iter="${FIGURE5F_SAMPLE_ITER:?FIGURE5F_SAMPLE_ITER is required}"
burnin="${FIGURE5F_SAMPLE_BURNIN:-3000}"
chunk_iter="${FIGURE5F_CHUNK_ITER:-10000}"
max_iter="${FIGURE5F_MAX_ITER:-200000}"
runner="${iteration_root}/Code/Figures/Figure5F/hpc/run_Figure5F_family_stage.sh"
controller="${iteration_root}/Code/Figures/Figure5F/hpc/control_Figure5F_family_sampling.sh"
output_dir="${iteration_root}/data/Figures/Figure5/generalized_posterior_family_conditioned"
checkpoint_dir="${iteration_root}/tmp/figure5f_generalized_posterior_family_conditioned"
log_dir="${iteration_root}/tmp/figure5f_family_conditioned_hpc_logs"
readiness="${output_dir}/figure5f_${family}_readiness.tsv"
state_file="${checkpoint_dir}/family_${family}_controller_state.tsv"
complete_file="${checkpoint_dir}/family_${family}_complete.tsv"

case "${family}" in C01|C02|C03) ;; *) echo "Invalid family: ${family}" >&2; exit 64 ;; esac
mkdir -p "${checkpoint_dir}" "${log_dir}"
if [[ ! -f "${readiness}" ]]; then
  echo "Missing family readiness file: ${readiness}" >&2
  exit 2
fi

passed_now=1
if ! awk -F '\t' -v expected_iter="${current_iter}" -v expected_family="${family}" '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "passed") p = i
      if ($i == "n_iter") n = i
      if ($i == "family") f = i
    }
    next
  }
  p == 0 || n == 0 || f == 0 || $p != "TRUE" ||
    $n != expected_iter || $f != expected_family { failed = 1 }
  END { exit failed }
' "${readiness}"; then
  passed_now=0
fi

previous_iter=0
previous_streak=0
if [[ -f "${state_file}" ]]; then
  previous_iter="$(awk -F '\t' 'NR == 2 {print $2}' "${state_file}")"
  previous_streak="$(awk -F '\t' 'NR == 2 {print $4}' "${state_file}")"
fi
if (( passed_now == 1 )); then
  if (( previous_iter == current_iter )); then
    streak="${previous_streak}"
  else
    streak="$((previous_streak + 1))"
  fi
else
  streak=0
fi

state_tmp="${state_file}.tmp_${SLURM_JOB_ID:-$$}"
printf 'family\tn_iter\tpassed\tconsecutive_passes\tupdated_at\n' > "${state_tmp}"
printf '%s\t%s\t%s\t%s\t%s\n' \
  "${family}" "${current_iter}" "${passed_now}" "${streak}" "$(date -Is)" >> "${state_tmp}"
mv "${state_tmp}" "${state_file}"

if (( streak >= 2 )); then
  complete_tmp="${complete_file}.tmp_${SLURM_JOB_ID:-$$}"
  printf 'family\tn_iter\tconsecutive_passes\tcompleted_at\n' > "${complete_tmp}"
  printf '%s\t%s\t%s\t%s\n' "${family}" "${current_iter}" "${streak}" "$(date -Is)" >> "${complete_tmp}"
  mv "${complete_tmp}" "${complete_file}"
  echo "${family} passed its frozen gates twice consecutively at n_iter=${current_iter}; this family is stopped."

  all_complete=1
  for candidate in C01 C02 C03; do
    if [[ ! -f "${checkpoint_dir}/family_${candidate}_complete.tsv" ]]; then
      all_complete=0
    fi
  done
  if (( all_complete == 1 )); then
    combine_lock="${checkpoint_dir}/combine_submit.lock"
    if mkdir "${combine_lock}" 2>/dev/null; then
      combine_job="$({ sbatch --parsable \
        --job-name=f5f_family_combine \
        --account=tao.li2 \
        --partition=red \
        --qos=small \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem=128G \
        --time=2-00:00:00 \
        --dependency="afterok:${SLURM_JOB_ID}" \
        --output="${log_dir}/combine_%j.out" \
        --error="${log_dir}/combine_%j.err" \
        --export=ALL \
        "${runner}" combine; } | cut -d';' -f1)"
      printf 'job_id\tsubmitted_at\n%s\t%s\n' "${combine_job}" "$(date -Is)" > \
        "${checkpoint_dir}/combine_submitted.tsv"
      echo "All families complete; submitted combine job ${combine_job}."
    else
      echo "All families complete; combine submission is already locked/submitted."
    fi
  fi
  exit 0
fi

if (( current_iter >= max_iter )); then
  halted="${checkpoint_dir}/family_${family}_halted_at_max.tsv"
  printf 'family\tn_iter\tconsecutive_passes\treadiness\thalted_at\n' > "${halted}"
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "${family}" "${current_iter}" "${streak}" "${readiness}" "$(date -Is)" >> "${halted}"
  echo "${family} reached frozen max n_iter=${max_iter} without two consecutive passes." >&2
  exit 2
fi

next_iter="$((current_iter + chunk_iter))"
if (( next_iter > max_iter )); then next_iter="${max_iter}"; fi
submission_file="${checkpoint_dir}/family_${family}_submitted_${next_iter}.tsv"
if [[ -f "${submission_file}" ]]; then
  echo "The n_iter=${next_iter} extension for ${family} was already submitted: ${submission_file}"
  exit 0
fi

common_export="ALL,FIGURE5F_ITERATION_ROOT=${iteration_root},FIGURE5F_RESULT_ROOT=${result_root},FIGURE5F_R_MODULE=R/4.4,FIGURE5F_FAMILY=${family},FIGURE5F_SAMPLE_ITER=${next_iter},FIGURE5F_SAMPLE_BURNIN=${burnin},FIGURE5F_CHUNK_ITER=${chunk_iter},FIGURE5F_MAX_ITER=${max_iter}"
sample_job="$({ sbatch --parsable \
  --job-name="f5f_${family}_sample" \
  --account=tao.li2 \
  --partition=red \
  --qos=small \
  --ntasks=2 \
  --cpus-per-task=5 \
  --mem=128G \
  --time=2-00:00:00 \
  --dependency="afterok:${SLURM_JOB_ID}" \
  --output="${log_dir}/${family}_sample_%j.out" \
  --error="${log_dir}/${family}_sample_%j.err" \
  --export="${common_export}" \
  "${runner}" sample; } | cut -d';' -f1)"

aggregate_job="$({ sbatch --parsable \
  --job-name="f5f_${family}_aggregate" \
  --account=tao.li2 \
  --partition=red \
  --qos=small \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem=64G \
  --time=2-00:00:00 \
  --dependency="afterok:${sample_job}" \
  --output="${log_dir}/${family}_aggregate_%j.out" \
  --error="${log_dir}/${family}_aggregate_%j.err" \
  --export="${common_export}" \
  "${runner}" aggregate; } | cut -d';' -f1)"

controller_job="$({ sbatch --parsable \
  --job-name="f5f_${family}_controller" \
  --account=tao.li2 \
  --partition=red \
  --qos=small \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem=2G \
  --time=2-00:00:00 \
  --dependency="afterok:${aggregate_job}" \
  --output="${log_dir}/${family}_controller_%j.out" \
  --error="${log_dir}/${family}_controller_%j.err" \
  --export="${common_export}" \
  "${controller}"; } | cut -d';' -f1)"

printf 'family\tn_iter\tsample_job\taggregate_job\tcontroller_job\tsubmitted_at\n' > "${submission_file}"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${family}" "${next_iter}" "${sample_job}" "${aggregate_job}" "${controller_job}" "$(date -Is)" >> "${submission_file}"
echo "Extended only ${family} to n_iter=${next_iter}: ${sample_job} -> ${aggregate_job} -> ${controller_job}"
