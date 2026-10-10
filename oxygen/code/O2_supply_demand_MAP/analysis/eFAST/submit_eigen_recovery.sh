#!/usr/bin/env bash
set -euo pipefail
repo=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
out=${EFAST_OUT_DIR:-$repo/oxygen/results/eFAST/neighborhood10pct}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
mode=${1:-submit}
[[ "$mode" == submit || "$mode" == --retry ]] || exit 2
exec 8>"$out/slurm/eigen_recovery_submit.lock"
flock -n 8 || { echo "Another recovery submission is in progress" >&2; exit 1; }
jobs="$out/slurm/numerical_recovery_jobs.tsv"
if [[ -s "$jobs" ]]; then
  [[ "$mode" == --retry ]] || { echo "Recovery jobs already recorded; inspect before retry" >&2; exit 1; }
  known=$(awk -F '\t' 'NR>1 {print $2}' "$jobs" | paste -sd, -)
  active=$(squeue -h -j "$known" -o '%i' 2>/dev/null || true)
  [[ -z "$active" ]] || { echo "Recovery jobs are still active" >&2; exit 1; }
fi
python3 - "$out" <<'PY'
import hashlib,json,pathlib,sys
root=pathlib.Path(sys.argv[1]); previous=root/'slurm/numerical_recovery_plan.json'
if previous.exists():
    archive=root/'slurm'/('numerical_recovery_plan_before_'+hashlib.sha256(previous.read_bytes()).hexdigest()[:12]+'.json')
    if not archive.exists(): archive.write_bytes(previous.read_bytes())
PY
apptainer exec --cleanenv "$sif" env OPENBLAS_NUM_THREADS=1 python3 "$code/neighborhood_slurm.py" prepare \
  --out-dir "$out" --fit-root /share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253 \
  --figure4-dir /share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration4/data/Figures/Figure4 \
  --sif "$sif" --replace-plan > "$out/slurm/eigen_recovery/parent_prepare.log"
apptainer exec --cleanenv "$sif" env OPENBLAS_NUM_THREADS=1 python3 "$code/recover_eigen.py" prepare --out-dir "$out" > "$out/slurm/eigen_recovery/prepare.log"
arrays=$(python3 - "$out/slurm/numerical_recovery_plan.json" <<'PY'
import json,sys
p=json.load(open(sys.argv[1])); print(p['array_spec']); print(p['summary_array_spec'])
PY
)
repair_spec=${arrays%%$'\n'*}; summary_spec=${arrays#*$'\n'}
[[ -n "$repair_spec" && -n "$summary_spec" ]] || { echo "No pending repair or summary array" >&2; exit 1; }
[[ -s "$jobs" ]] || printf 'role\tjob_id\tarray_spec\tqos\ttime\tcpus\tmem\tdependency\tsubmitted_utc\n' > "$jobs"
common=(--parsable --qos=xxlarge --time=12:00:00)
exports="ALL,EFAST_REPO_ROOT=$repo,EFAST_OUT_DIR=$out,EFAST_SIF=$sif"
record() {
  printf '%s\t%s\t%s\txxlarge\t12:00:00\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$jobs"
}
repair=$(sbatch "${common[@]}" --array="$repair_spec" --cpus-per-task=1 --mem=4G \
  --job-name=efast_point_recovery --output="$out/slurm/logs/point_recovery_%A_%a.out" \
  --error="$out/slurm/logs/point_recovery_%A_%a.err" --export="$exports,EFAST_RECOVERY_ROLE=repair" "$code/run_eigen_recovery.sbatch")
repair=${repair%%;*}; record repair "$repair" "$repair_spec" 1 4G ''
summary=$(sbatch "${common[@]}" --array="$summary_spec" --cpus-per-task=1 --mem=8G \
  --dependency="afterany:$repair" --job-name=efast_recovered_summary \
  --output="$out/slurm/logs/recovered_summary_%A_%a.out" --error="$out/slurm/logs/recovered_summary_%A_%a.err" \
  --export="$exports,EFAST_RECOVERY_ROLE=summary" "$code/run_eigen_recovery.sbatch")
summary=${summary%%;*}; record seed_summary "$summary" "$summary_spec" 1 8G "afterany:$repair"
final=$(sbatch "${common[@]}" --cpus-per-task=4 --mem=32G --dependency="afterany:$summary" \
  --job-name=efast_recovered_final --output="$out/slurm/logs/recovered_final_%j.out" \
  --error="$out/slurm/logs/recovered_final_%j.err" --export="$exports,EFAST_RECOVERY_ROLE=finalize" "$code/run_eigen_recovery.sbatch")
final=${final%%;*}; record finalize "$final" '' 4 32G "afterany:$summary"
cat "$jobs"
