#!/usr/bin/env bash
set -euo pipefail
repo=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
baseline=${EFAST_BASELINE_DIR:-$repo/oxygen/results/eFAST/neighborhood10pct}
pilot=${EFAST_PILOT_ROOT:-$baseline/pilot/N8193_M4_and_N4097_M8}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
mode=${1:-submit}
[[ "$mode" == prepare || "$mode" == submit ]] || exit 2
mkdir -p "$pilot"
exec 8>"$pilot/submit.lock"
flock -n 8 || { echo 'Another N/M submission is in progress' >&2; exit 1; }
common=(--parsable --qos=xxlarge --time=12:00:00)
for config in N8193_M4 N4097_M8; do
  out="$pilot/$config"
  mkdir -p "$out/slurm/logs"
  if [[ "$mode" == prepare ]]; then
    apptainer exec --cleanenv "$sif" env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
      python3 "$code/resolution_nm.py" prepare --out-dir "$out" --baseline-dir "$baseline" --config "$config"
    continue
  fi
  apptainer exec --cleanenv "$sif" env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python3 - "$out" "$sif" "$code" <<'PY'
import pathlib,sys
sys.path.insert(0,sys.argv[3])
import resolution_nm as run
plan=run.context(pathlib.Path(sys.argv[1]),verify_sources=True)
assert pathlib.Path(sys.argv[2]).resolve()==pathlib.Path(plan['sif'])
assert plan['qos']=='xxlarge' and plan['time']=='12:00:00'
assert plan['specified_compute_node'] is None and plan['array_concurrency_limit'] is None
PY
  exports="ALL,EFAST_REPO_ROOT=$repo,EFAST_OUT_DIR=$out,EFAST_BASELINE_DIR=$baseline,EFAST_SIF=$sif"
  record="$out/slurm/submission_jobs.tsv"
  if [[ ! -f "$record" ]]; then
    printf 'role\tjob_id\tarray_spec\tqos\ttime\tcpus\tmem\tdependency\tsubmitted_utc\n' > "$record"
  fi
  submit_role() {
    local role=$1 array=$2 cpus=$3 mem=$4 dependency=$5 job
    local marker="$out/slurm/${role}_job_id.txt"
    if [[ -s "$marker" ]]; then
      cat "$marker"
      return
    fi
    local args=("${common[@]}" --cpus-per-task="$cpus" --mem="$mem" --job-name="efast_${config}_${role}"
      --output="$out/slurm/logs/${role}_%A_%a.out" --error="$out/slurm/logs/${role}_%A_%a.err"
      --export="$exports,EFAST_ARRAY_ROLE=$role")
    [[ -z "$array" ]] || args+=(--array="$array")
    [[ -z "$dependency" ]] || args+=(--dependency="$dependency")
    job=$(sbatch "${args[@]}" "$code/run_resolution_nm.sbatch")
    job=${job%%;*}
    [[ "$job" =~ ^[0-9]+$ ]] || { echo "Invalid sbatch output: $job" >&2; return 1; }
    printf '%s\n' "$job" > "$marker"
    printf '%s\t%s\t%s\txxlarge\t12:00:00\t%s\t%s\t%s\t%s\n' \
      "$role" "$job" "$array" "$cpus" "$mem" "$dependency" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$record"
    printf '%s\n' "$job"
  }
  smoke=$(submit_role smoke '' 1 4G '')
  trajectory=$(submit_role trajectory 1-140 1 4G "afterok:$smoke")
  summary=$(submit_role seed_summary 1-2 1 8G "afterany:$trajectory")
  finalize=$(submit_role finalize '' 4 32G "afterany:$summary")
  echo "$config"
  cat "$record"
done
