#!/usr/bin/env bash
set -euo pipefail
repo=${SOBOL_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
fit=${SOBOL_FIT_ROOT:-/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253}
old=${SOBOL_EFAST_ROOT:-$repo/oxygen/results/eFAST}
figure=${SOBOL_FIGURE4_DIR:-/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration4/data/Figures/Figure4}
sif=${SOBOL_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
out=${SOBOL_OUT_DIR:-$repo/oxygen/results/eFAST/sobol/full_N4096}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol"
mkdir -p "$out/slurm/logs"
exec 8>"$out/slurm/submit.lock"
flock -n 8 || { echo 'Another submission is in progress' >&2; exit 1; }
[[ ! -e "$out/slurm/submission_jobs.tsv" ]] || { echo 'Already submitted; inspect recorded jobs' >&2; exit 1; }
snapshot="$out/slurm/code_snapshot"
[[ ! -d "$snapshot" ]] || { echo 'Snapshot exists; inspect interrupted preparation' >&2; exit 1; }
mkdir -p "$snapshot"
cd "$repo"
revision=$(git rev-parse HEAD)
git archive --format=tar "$revision" oxygen/code/O2_supply_demand_MAP | tar -xf - -C "$snapshot"
mkdir -p "$snapshot/oxygen/code/O2_supply_demand_MAP/model/.rcpp_cache_o2_supply_demand_map"
apptainer exec --cleanenv "$sif" env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 "$code/full_run.py" prepare --out-dir "$out" --fit-root "$fit" --efast-root "$old" \
  --figure4-dir "$figure" --sif "$sif" --snapshot-root "$snapshot" --git-commit "$revision"
python3 - "$snapshot" <<'PY'
import pathlib,sys
p=pathlib.Path(sys.argv[1])
for f in p.rglob('*'):
    f.chmod(0o555 if f.is_dir() else 0o444)
p.chmod(0o555)
PY
common=(--parsable --qos=xxlarge --time=12:00:00)
exports="ALL,SOBOL_OUT_DIR=$out,SOBOL_SIF=$sif"
printf 'role\tjob_id\tarray_spec\tqos\ttime\tcpus\tmem\tdependency\tsubmitted_utc\n' > "$out/slurm/submission_jobs.tsv"
record() {
  printf '%s\t%s\t%s\txxlarge\t12:00:00\t1\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$out/slurm/submission_jobs.tsv"
}
submit() {
  local role=$1 array=$2 memory=$3 dependency=$4
  local options=(--cpus-per-task=1 --mem="$memory" --job-name="sobol_full_$role" --output="$out/slurm/logs/${role}_%A_%a.out" --error="$out/slurm/logs/${role}_%A_%a.err" --export="$exports,SOBOL_ROLE=$role")
  [[ -z "$array" ]] || options+=(--array="$array")
  [[ -z "$dependency" ]] || options+=(--dependency="$dependency")
  local job
  job=$(sbatch "${common[@]}" "${options[@]}" "$snapshot/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/run_full_array.sbatch")
  job=${job%%;*}; record "$role" "$job" "$array" "$memory" "$dependency"; printf '%s' "$job"
}
bootstrap=$(submit bootstrap '' 4G '')
validation=$(submit validate 1-501 4G "afterok:$bootstrap")
gate=$(submit gate '' 4G "afterany:$validation")
evaluation=$(submit evaluate 1-75150 4G "afterok:$gate")
summary=$(submit summarize 1-501 16G "afterany:$evaluation")
aggregate=$(submit aggregate '' 16G "afterany:$summary")
response=$(submit response 1-54 4G "afterok:$aggregate")
finish=$(submit finish '' 8G "afterany:$response")
cat "$out/slurm/submission_jobs.tsv"
