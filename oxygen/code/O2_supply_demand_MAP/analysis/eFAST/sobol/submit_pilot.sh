#!/usr/bin/env bash
set -euo pipefail
repo=${SOBOL_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
fit=${SOBOL_FIT_ROOT:-/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253}
old=${SOBOL_EFAST_ROOT:-$repo/oxygen/results/eFAST}
figure=${SOBOL_FIGURE4_DIR:-/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration4/data/Figures/Figure4}
sif=${SOBOL_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
out=${SOBOL_OUT_DIR:-$repo/oxygen/results/eFAST/sobol}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol"
mkdir -p "$out/slurm/logs"
exec 8>"$out/slurm/submit.lock"
flock -n 8 || { echo 'Another pilot submission is in progress' >&2; exit 1; }
[[ ! -e "$out/slurm/submission_jobs.tsv" ]] || { echo 'Pilot already submitted; inspect recorded jobs before another attempt' >&2; exit 1; }
apptainer exec --cleanenv "$sif" env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 "$code/run.py" prepare --out-dir "$out" --fit-root "$fit" --efast-root "$old" \
  --figure4-dir "$figure" --sif "$sif"
common=(--parsable --qos=xxlarge --time=12:00:00)
exports="ALL,SOBOL_REPO_ROOT=$repo,SOBOL_OUT_DIR=$out,SOBOL_SIF=$sif"
printf 'role\tjob_id\tarray_spec\tqos\ttime\tcpus\tmem\tdependency\tsubmitted_utc\n' > "$out/slurm/submission_jobs.tsv"
record() {
  printf '%s\t%s\t%s\txxlarge\t12:00:00\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$out/slurm/submission_jobs.tsv"
}
validation=$(sbatch "${common[@]}" --cpus-per-task=1 --mem=4G \
  --job-name=sobol_pilot_validate --output="$out/slurm/logs/validate_%j.out" \
  --error="$out/slurm/logs/validate_%j.err" --export="$exports,SOBOL_ROLE=validate" "$code/run_array.sbatch")
validation=${validation%%;*}; record validate "$validation" '' 1 4G ''
evaluation=$(sbatch "${common[@]}" --array=1-1350 --cpus-per-task=1 --mem=4G \
  --dependency="afterok:$validation" \
  --job-name=sobol_pilot_eval --output="$out/slurm/logs/evaluate_%A_%a.out" \
  --error="$out/slurm/logs/evaluate_%A_%a.err" --export="$exports,SOBOL_ROLE=evaluate" "$code/run_array.sbatch")
evaluation=${evaluation%%;*}; record evaluate "$evaluation" 1-1350 1 4G "afterok:$validation"
summary=$(sbatch "${common[@]}" --array=1-9 --cpus-per-task=1 --mem=8G --dependency="afterany:$evaluation" \
  --job-name=sobol_pilot_summary --output="$out/slurm/logs/summary_%A_%a.out" \
  --error="$out/slurm/logs/summary_%A_%a.err" --export="$exports,SOBOL_ROLE=summarize" "$code/run_array.sbatch")
summary=${summary%%;*}; record summarize "$summary" 1-9 1 8G "afterany:$evaluation"
finalize=$(sbatch "${common[@]}" --cpus-per-task=1 --mem=8G --dependency="afterany:$summary" \
  --job-name=sobol_pilot_indices --output="$out/slurm/logs/indices_%j.out" \
  --error="$out/slurm/logs/indices_%j.err" --export="$exports,SOBOL_ROLE=finalize" "$code/run_array.sbatch")
finalize=${finalize%%;*}; record finalize "$finalize" '' 1 8G "afterany:$summary"
response=$(sbatch "${common[@]}" --array=1-54 --cpus-per-task=1 --mem=4G --dependency="afterok:$finalize" \
  --job-name=sobol_pilot_response --output="$out/slurm/logs/response_%A_%a.out" \
  --error="$out/slurm/logs/response_%A_%a.err" --export="$exports,SOBOL_ROLE=response" "$code/run_array.sbatch")
response=${response%%;*}; record response "$response" 1-54 1 4G "afterok:$finalize"
finish=$(sbatch "${common[@]}" --cpus-per-task=1 --mem=8G --dependency="afterany:$response" \
  --job-name=sobol_pilot_report --output="$out/slurm/logs/report_%j.out" \
  --error="$out/slurm/logs/report_%j.err" --export="$exports,SOBOL_ROLE=finish" "$code/run_array.sbatch")
finish=${finish%%;*}; record finish "$finish" '' 1 8G "afterany:$response"
cat "$out/slurm/submission_jobs.tsv"
