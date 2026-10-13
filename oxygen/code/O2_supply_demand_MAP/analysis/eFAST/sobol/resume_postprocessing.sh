#!/usr/bin/env bash
set -euo pipefail
repo=${SOBOL_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
out=${SOBOL_OUT_DIR:-$repo/oxygen/results/eFAST/sobol}
sif=${SOBOL_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol"
exec 8>"$out/slurm/submit.lock"
flock -n 8 || exit 1
apptainer exec --cleanenv "$sif" env OPENBLAS_NUM_THREADS=1 python3 \
  "$code/prepare_postprocessing_resume.py" --out-dir "$out"
common=(--parsable --qos=xxlarge --time=12:00:00)
exports="ALL,SOBOL_REPO_ROOT=$repo,SOBOL_OUT_DIR=$out,SOBOL_SIF=$sif"
record() {
  printf '%s\t%s\t%s\txxlarge\t12:00:00\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$out/slurm/submission_jobs.tsv"
}
finalize=$(sbatch "${common[@]}" --cpus-per-task=1 --mem=8G \
  --job-name=sobol_pilot_indices --output="$out/slurm/logs/indices_%j.out" \
  --error="$out/slurm/logs/indices_%j.err" --export="$exports,SOBOL_ROLE=finalize" "$code/run_array.sbatch")
finalize=${finalize%%;*}; record finalize_retry "$finalize" '' 1 8G ''
response=$(sbatch "${common[@]}" --array=1-54 --cpus-per-task=1 --mem=4G --dependency="afterok:$finalize" \
  --job-name=sobol_pilot_response --output="$out/slurm/logs/response_%A_%a.out" \
  --error="$out/slurm/logs/response_%A_%a.err" --export="$exports,SOBOL_ROLE=response" "$code/run_array.sbatch")
response=${response%%;*}; record response_retry "$response" 1-54 1 4G "afterok:$finalize"
finish=$(sbatch "${common[@]}" --cpus-per-task=1 --mem=8G --dependency="afterany:$response" \
  --job-name=sobol_pilot_report --output="$out/slurm/logs/report_%j.out" \
  --error="$out/slurm/logs/report_%j.err" --export="$exports,SOBOL_ROLE=finish" "$code/run_array.sbatch")
finish=${finish%%;*}; record finish_retry "$finish" '' 1 8G "afterany:$response"
tail -3 "$out/slurm/submission_jobs.tsv"
