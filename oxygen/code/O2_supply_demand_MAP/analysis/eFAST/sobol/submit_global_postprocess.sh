#!/usr/bin/env bash
set -euo pipefail
repo=${SOBOL_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
source_root=${GLOBAL_SOURCE_ROOT:-$repo/oxygen/results/eFAST/sobol/full_N4096}
out=${GLOBAL_OUT_DIR:-$source_root/global/standalone_N4096}
sif=${GLOBAL_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
mkdir -p "$out/slurm/logs"
exec 8>"$out/slurm/submit.lock"
flock -n 8 || { echo 'Another global submission is in progress' >&2; exit 1; }
[[ ! -e "$out/slurm/submission_jobs.tsv" ]] || { echo 'Global job already recorded; inspect rather than resubmit' >&2; exit 1; }
snapshot="$out/slurm/code_snapshot"
[[ ! -e "$snapshot" ]] || { echo 'Snapshot already exists; inspect interrupted preparation' >&2; exit 1; }
mkdir "$snapshot"
cd "$repo"
revision=$(git rev-parse HEAD)
git archive --format=tar "$revision" oxygen/code/O2_supply_demand_MAP | tar -xf - -C "$snapshot"
code="$snapshot/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol"
fit=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["fit_root"])' "$source_root/slurm/submission_plan.json")
apptainer exec --cleanenv --bind "$source_root:$source_root:ro,$out:$out:rw,$fit:$fit:ro" "$sif" \
  env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
  python3 "$code/global_postprocess.py" prepare --source-root "$source_root" --out-dir "$out" \
  --snapshot-root "$snapshot" --git-commit "$revision"
python3 - "$snapshot" <<'PY'
import pathlib,sys
p=pathlib.Path(sys.argv[1])
for f in p.rglob('*'):
    f.chmod(0o555 if f.is_dir() else 0o444)
p.chmod(0o555)
PY
squeue -h -j 20839424,20839425,20839426,20839427 -o '%i|%T|%E' > "$out/slurm/existing_pipeline_dependencies.txt"
job=$(sbatch --parsable --qos=xxlarge --time=12:00:00 --cpus-per-task=1 --mem=16G \
  --job-name=sobol_global_postprocess --output="$out/slurm/logs/global_%j.out" \
  --error="$out/slurm/logs/global_%j.err" \
  --export="ALL,GLOBAL_SOURCE_ROOT=$source_root,GLOBAL_OUT_DIR=$out,GLOBAL_SIF=$sif" \
  "$code/run_global_postprocess.sbatch")
job=${job%%;*}
printf 'role\tjob_id\tqos\ttime\tcpus\tmem\tdependency\tsubmitted_utc\n' > "$out/slurm/submission_jobs.tsv"
printf 'global_postprocess\t%s\txxlarge\t12:00:00\t1\t16G\t\t%s\n' "$job" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$out/slurm/submission_jobs.tsv"
cat "$out/slurm/submission_jobs.tsv"
