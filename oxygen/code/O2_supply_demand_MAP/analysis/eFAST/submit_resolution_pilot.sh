#!/usr/bin/env bash
set -euo pipefail
repo=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
baseline=${EFAST_BASELINE_DIR:-$repo/oxygen/results/eFAST/neighborhood10pct}
out=${EFAST_OUT_DIR:-$baseline/pilot/N1025_N2049}
fit=${EFAST_FIT_ROOT:-/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253}
figure=${EFAST_FIGURE4_DIR:-/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration4/data/Figures/Figure4}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
mode=${1:-submit}
[[ "$mode" == prepare || "$mode" == smoke || "$mode" == submit ]] || exit 2
mkdir -p "$out/slurm/logs"
exec 8>"$out/slurm/submit.lock"
flock -n 8 || { echo 'Another pilot submission is in progress' >&2; exit 1; }
common=(--parsable --qos=xxlarge --time=12:00:00)
exports="ALL,EFAST_REPO_ROOT=$repo,EFAST_OUT_DIR=$out,EFAST_BASELINE_DIR=$baseline,EFAST_SIF=$sif"
if [[ "$mode" == prepare ]]; then
  apptainer exec --cleanenv "$sif" env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    python3 "$code/resolution_pilot.py" prepare --out-dir "$out" --baseline-dir "$baseline" \
    --fit-root "$fit" --figure4-dir "$figure" --sif "$sif"
  exit
fi
[[ -s "$out/slurm/submission_plan.json" ]] || { echo 'Run prepare first' >&2; exit 1; }
record() {
  printf '%s\t%s\t%s\txxlarge\t12:00:00\t%s\t%s\t%s\t%s\n' \
    "$1" "$2" "$3" "$4" "$5" "$6" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$out/slurm/submission_jobs.tsv"
}
if [[ ! -f "$out/slurm/submission_jobs.tsv" ]]; then
  printf 'role\tjob_id\tarray_spec\tqos\ttime\tcpus\tmem\tdependency\tsubmitted_utc\n' > "$out/slurm/submission_jobs.tsv"
fi
if [[ "$mode" == smoke ]]; then
  # A prior smoke job cannot be silently duplicated.
  [[ ! -s "$out/slurm/smoke_job_id.txt" ]] || { echo 'Smoke already submitted; inspect its logs/receipt' >&2; exit 1; }
  job=$(sbatch "${common[@]}" --cpus-per-task=1 --mem=4G --job-name=efast_N_smoke \
    --output="$out/slurm/logs/smoke_%j.out" --error="$out/slurm/logs/smoke_%j.err" \
    --export="$exports,EFAST_ARRAY_ROLE=smoke" "$code/run_resolution_pilot.sbatch")
  job=${job%%;*}
  printf '%s\n' "$job" > "$out/slurm/smoke_job_id.txt"
  record smoke "$job" '' 1 4G ''
  echo "$job"
  exit
fi
[[ ! -s "$out/slurm/trajectory_job_id.txt" ]] || { echo 'Pilot trajectory array already submitted' >&2; exit 1; }
# The plan hash ties the passed smoke to exactly the submitted scientific code.
apptainer exec --cleanenv "$sif" python3 - "$out" <<'PY'
import hashlib,json,pathlib,sys
root=pathlib.Path(sys.argv[1]); receipt=json.loads((root/'smoke/smoke_receipt.json').read_text())
assert receipt['status']=='passed', 'Require passed compute-node smoke'
assert receipt['submission_plan_sha256']==hashlib.sha256((root/'slurm/submission_plan.json').read_bytes()).hexdigest(), 'Smoke plan changed'
PY
trajectory=$(sbatch "${common[@]}" --array=1-1680 --cpus-per-task=1 --mem=4G \
  --job-name=efast_N_trajectory --output="$out/slurm/logs/trajectory_%A_%a.out" \
  --error="$out/slurm/logs/trajectory_%A_%a.err" \
  --export="$exports,EFAST_ARRAY_ROLE=trajectory" "$code/run_resolution_pilot.sbatch")
trajectory=${trajectory%%;*}
printf '%s\n' "$trajectory" > "$out/slurm/trajectory_job_id.txt"
record trajectory "$trajectory" 1-1680 1 4G ''
summary=$(sbatch "${common[@]}" --array=1-24 --cpus-per-task=1 --mem=8G \
  --dependency="afterany:$trajectory" --job-name=efast_N_summary \
  --output="$out/slurm/logs/summary_%A_%a.out" --error="$out/slurm/logs/summary_%A_%a.err" \
  --export="$exports,EFAST_ARRAY_ROLE=seed_summary" "$code/run_resolution_pilot.sbatch")
summary=${summary%%;*}
record seed_summary "$summary" 1-24 1 8G "afterany:$trajectory"
final=$(sbatch "${common[@]}" --cpus-per-task=4 --mem=32G --dependency="afterany:$summary" \
  --job-name=efast_N_compare --output="$out/slurm/logs/finalize_%j.out" \
  --error="$out/slurm/logs/finalize_%j.err" \
  --export="$exports,EFAST_ARRAY_ROLE=finalize" "$code/run_resolution_pilot.sbatch")
final=${final%%;*}
record finalize "$final" '' 4 32G "afterany:$summary"
cat "$out/slurm/submission_jobs.tsv"
