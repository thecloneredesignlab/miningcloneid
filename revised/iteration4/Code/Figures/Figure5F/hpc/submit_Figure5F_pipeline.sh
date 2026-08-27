#!/usr/bin/env bash
set -euo pipefail

iteration_root="${FIGURE5F_ITERATION_ROOT:-/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration2}"
result_root="${FIGURE5F_RESULT_ROOT:-/share/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540}"
runner="${iteration_root}/Code/Figures/Figure5F/hpc/run_Figure5F_stage.sh"
log_dir="${iteration_root}/tmp/figure5f_hpc_logs"

mkdir -p "${log_dir}"
for path in "${iteration_root}" "${result_root}" "${runner}"; do
  if [[ ! -e "${path}" ]]; then
    echo "Missing required submission input: ${path}" >&2
    exit 2
  fi
done

export FIGURE5F_ITERATION_ROOT="${iteration_root}"
export FIGURE5F_RESULT_ROOT="${result_root}"
export FIGURE5F_R_MODULE="${FIGURE5F_R_MODULE:-R/4.4}"

prepare_job="$({ sbatch --parsable \
  --job-name=f5f_prepare \
  --account=tao.li2 \
  --partition=red \
  --qos=xxlarge \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem=16G \
  --time=12:00:00 \
  --output="${log_dir}/prepare_%j.out" \
  --error="${log_dir}/prepare_%j.err" \
  --export=ALL \
  "${runner}" prepare; } | cut -d';' -f1)"

benchmark_job="$({ sbatch --parsable \
  --job-name=f5f_benchmark \
  --account=tao.li2 \
  --partition=red \
  --qos=xxlarge \
  --ntasks=1 \
  --cpus-per-task=30 \
  --mem=192G \
  --time=12:00:00 \
  --dependency="afterok:${prepare_job}" \
  --output="${log_dir}/benchmark_%j.out" \
  --error="${log_dir}/benchmark_%j.err" \
  --export=ALL \
  "${runner}" benchmark; } | cut -d';' -f1)"

pilot_ladders_job="$({ sbatch --parsable \
  --job-name=f5f_pilot_ladders \
  --account=tao.li2 \
  --partition=red \
  --qos=xxlarge \
  --ntasks=6 \
  --cpus-per-task=5 \
  --mem=192G \
  --time=12:00:00 \
  --dependency="afterok:${benchmark_job}" \
  --output="${log_dir}/pilot_ladders_%j.out" \
  --error="${log_dir}/pilot_ladders_%j.err" \
  --export=ALL \
  "${runner}" pilot_ladders; } | cut -d';' -f1)"

pilot_aggregate_job="$({ sbatch --parsable \
  --job-name=f5f_pilot_aggregate \
  --account=tao.li2 \
  --partition=red \
  --qos=xxlarge \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem=32G \
  --time=12:00:00 \
  --dependency="afterok:${pilot_ladders_job}" \
  --output="${log_dir}/pilot_aggregate_%j.out" \
  --error="${log_dir}/pilot_aggregate_%j.err" \
  --export=ALL \
  "${runner}" pilot_aggregate; } | cut -d';' -f1)"

calibrate_job="$({ sbatch --parsable \
  --job-name=f5f_calibrate \
  --account=tao.li2 \
  --partition=red \
  --qos=xxlarge \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem=8G \
  --time=12:00:00 \
  --dependency="afterok:${pilot_aggregate_job}" \
  --output="${log_dir}/calibrate_%j.out" \
  --error="${log_dir}/calibrate_%j.err" \
  --export=ALL \
  "${runner}" calibrate; } | cut -d';' -f1)"

sample_ladders_job="$({ sbatch --parsable \
  --job-name=f5f_sample_ladders \
  --account=tao.li2 \
  --partition=red \
  --qos=xxlarge \
  --ntasks=6 \
  --cpus-per-task=5 \
  --mem=192G \
  --time=12:00:00 \
  --dependency="afterok:${calibrate_job}" \
  --output="${log_dir}/sample_ladders_%j.out" \
  --error="${log_dir}/sample_ladders_%j.err" \
  --export=ALL \
  "${runner}" sample_ladders; } | cut -d';' -f1)"

sample_aggregate_job="$({ sbatch --parsable \
  --job-name=f5f_sample_aggregate \
  --account=tao.li2 \
  --partition=red \
  --qos=xxlarge \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem=32G \
  --time=12:00:00 \
  --dependency="afterok:${sample_ladders_job}" \
  --output="${log_dir}/sample_aggregate_%j.out" \
  --error="${log_dir}/sample_aggregate_%j.err" \
  --export=ALL \
  "${runner}" sample_aggregate; } | cut -d';' -f1)"

printf 'stage\tjob_id\tdependency\n'
printf 'prepare\t%s\tnone\n' "${prepare_job}"
printf 'benchmark\t%s\tafterok:%s\n' "${benchmark_job}" "${prepare_job}"
printf 'pilot_ladders\t%s\tafterok:%s\n' "${pilot_ladders_job}" "${benchmark_job}"
printf 'pilot_aggregate\t%s\tafterok:%s\n' "${pilot_aggregate_job}" "${pilot_ladders_job}"
printf 'calibrate\t%s\tafterok:%s\n' "${calibrate_job}" "${pilot_aggregate_job}"
printf 'sample_ladders\t%s\tafterok:%s\n' "${sample_ladders_job}" "${calibrate_job}"
printf 'sample_aggregate\t%s\tafterok:%s\n' "${sample_aggregate_job}" "${sample_ladders_job}"
