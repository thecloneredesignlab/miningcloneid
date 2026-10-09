#!/usr/bin/env bash
set -euo pipefail
repo=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
fit=${EFAST_FIT_ROOT:-/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253}
figure=${EFAST_FIGURE4_DIR:-/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration4/data/Figures/Figure4}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
out=${EFAST_OUT_DIR:-$repo/oxygen/results/eFAST/neighborhood10pct}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
mkdir -p "$out/slurm/logs"
exec 8>"$out/slurm/submit.lock"
flock -n 8 || { echo "Another Slurm submission is in progress" >&2; exit 1; }
[[ ! -s "$out/slurm/submission_jobs.tsv" ]] || { echo "Jobs already submitted; use the recorded task IDs to retry only failed tasks" >&2; exit 1; }
maximum=$(scontrol show config | awk '/^[[:space:]]*MaxArraySize[[:space:]]*=/ {print $3}')
[[ "$maximum" =~ ^[0-9]+$ && "$maximum" -gt 35000 ]] || { echo "This submission requires MaxArraySize > 35000; split arrays for this cluster" >&2; exit 1; }
apptainer exec --cleanenv "$sif" env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 "$code/neighborhood_slurm.py" prepare --out-dir "$out" \
  --fit-root "$fit" --figure4-dir "$figure" --sif "$sif"
common=(--parsable --qos=xxlarge --time=12:00:00)
exports="ALL,EFAST_REPO_ROOT=$repo,EFAST_OUT_DIR=$out,EFAST_SIF=$sif"
printf 'role\tjob_id\tarray_spec\tqos\ttime\tcpus\tmem\tdependency\tsubmitted_utc\n' > "$out/slurm/submission_jobs.tsv"
record() {
  printf '%s\t%s\t%s\txxlarge\t12:00:00\t%s\t%s\t%s\t%s\n' \
    "$1" "$2" "$3" "$4" "$5" "$6" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$out/slurm/submission_jobs.tsv"
}
trajectory=$(sbatch "${common[@]}" --array=1-35000 --cpus-per-task=1 --mem=4G \
  --job-name=efast_trajectory --output="$out/slurm/logs/trajectory_%A_%a.out" \
  --error="$out/slurm/logs/trajectory_%A_%a.err" \
  --export="$exports,EFAST_ARRAY_ROLE=trajectory" "$code/run_neighborhood_array.sbatch")
trajectory=${trajectory%%;*}
record trajectory "$trajectory" 1-35000 1 4G ''
summary=$(sbatch "${common[@]}" --array=1-500 --cpus-per-task=1 --mem=8G \
  --dependency="afterany:$trajectory" --job-name=efast_seed_summary \
  --output="$out/slurm/logs/seed_summary_%A_%a.out" --error="$out/slurm/logs/seed_summary_%A_%a.err" \
  --export="$exports,EFAST_ARRAY_ROLE=seed_summary" "$code/run_neighborhood_array.sbatch")
summary=${summary%%;*}
record seed_summary "$summary" 1-500 1 8G "afterany:$trajectory"
final=$(sbatch "${common[@]}" --cpus-per-task=4 --mem=32G --dependency="afterany:$summary" \
  --job-name=efast_finalize --output="$out/slurm/logs/finalize_%j.out" \
  --error="$out/slurm/logs/finalize_%j.err" \
  --export="$exports,EFAST_ARRAY_ROLE=finalize" "$code/run_neighborhood_array.sbatch")
final=${final%%;*}
record finalize "$final" '' 4 32G "afterany:$summary"
cat "$out/slurm/submission_jobs.tsv"
