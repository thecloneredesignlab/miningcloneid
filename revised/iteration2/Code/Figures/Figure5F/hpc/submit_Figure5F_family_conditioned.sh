#!/usr/bin/env bash
set -euo pipefail

iteration_root="${FIGURE5F_ITERATION_ROOT:-/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration2}"
result_root="${FIGURE5F_RESULT_ROOT:-/share/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540}"
initial_iter="${FIGURE5F_INITIAL_ITER:-20000}"
burnin="${FIGURE5F_SAMPLE_BURNIN:-3000}"
chunk_iter="${FIGURE5F_CHUNK_ITER:-10000}"
max_iter="${FIGURE5F_MAX_ITER:-200000}"
runner="${iteration_root}/Code/Figures/Figure5F/hpc/run_Figure5F_family_stage.sh"
controller="${iteration_root}/Code/Figures/Figure5F/hpc/control_Figure5F_family_sampling.sh"
log_dir="${iteration_root}/tmp/figure5f_family_conditioned_hpc_logs"

mkdir -p "${log_dir}"
for path in "${iteration_root}" "${result_root}" "${runner}" "${controller}"; do
  if [[ ! -e "${path}" ]]; then
    echo "Missing required family-conditioned input: ${path}" >&2
    exit 2
  fi
done

export FIGURE5F_ITERATION_ROOT="${iteration_root}"
export FIGURE5F_RESULT_ROOT="${result_root}"
export FIGURE5F_R_MODULE="R/4.4"

prepare_job="$({ sbatch --parsable \
  --job-name=f5f_family_prepare \
  --account=tao.li2 \
  --partition=red \
  --qos=small \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem=64G \
  --time=2-00:00:00 \
  --output="${log_dir}/prepare_%j.out" \
  --error="${log_dir}/prepare_%j.err" \
  --export=ALL \
  "${runner}" prepare; } | cut -d';' -f1)"

printf 'family\tstage\tjob_id\tdependency\tn_iter\n'
printf 'ALL\tprepare\t%s\tnone\tNA\n' "${prepare_job}"

for family in C01 C02 C03; do
  common_export="ALL,FIGURE5F_ITERATION_ROOT=${iteration_root},FIGURE5F_RESULT_ROOT=${result_root},FIGURE5F_R_MODULE=R/4.4,FIGURE5F_FAMILY=${family},FIGURE5F_SAMPLE_ITER=${initial_iter},FIGURE5F_SAMPLE_BURNIN=${burnin},FIGURE5F_CHUNK_ITER=${chunk_iter},FIGURE5F_MAX_ITER=${max_iter}"
  sample_job="$({ sbatch --parsable \
    --job-name="f5f_${family}_sample" \
    --account=tao.li2 \
    --partition=red \
    --qos=small \
    --ntasks=2 \
    --cpus-per-task=5 \
    --mem=128G \
    --time=2-00:00:00 \
    --dependency="afterok:${prepare_job}" \
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

  printf '%s\tsample\t%s\tafterok:%s\t%s\n' "${family}" "${sample_job}" "${prepare_job}" "${initial_iter}"
  printf '%s\taggregate\t%s\tafterok:%s\t%s\n' "${family}" "${aggregate_job}" "${sample_job}" "${initial_iter}"
  printf '%s\tcontroller\t%s\tafterok:%s\t%s\n' "${family}" "${controller_job}" "${aggregate_job}" "${initial_iter}"
done
