#!/usr/bin/env bash
set -euo pipefail
root=$1
role=$2
n=$3
snapshot="$root/slurm/code_snapshot"
scratch=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/efast_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
trap 'exit 143' TERM INT
source "$snapshot/retry_global_startup.sh"
readarray -t settings < <(python3 - "$root/plan.json" <<'PY'
import json,sys
p=json.load(open(sys.argv[1]));old=json.load(open(p['original_plan']))
print(p['sif']['path']);print(old['code_snapshot_root']);print(p['source_initial_root']);print(p['canonical_helpers'])
PY
)
sif=${settings[0]};frozen=${settings[1]};source_root=${settings[2]};helpers=${settings[3]}
mkdir "$scratch/model_cache" "$scratch/perron_cache"
expected=$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["archive_sha256"])' "$source_root/slurm/runtime_cache_manifest.json")
actual=$(sha256sum "$source_root/slurm/rcpp_template.tar.gz");actual=${actual%% *}
[[ "$actual" = "$expected" ]] || { echo 'Runtime template changed' >&2;exit 2; }
tar -xzf "$source_root/slurm/rcpp_template.tar.gz" -C "$scratch/model_cache"
model_cache="$frozen/oxygen/code/O2_supply_demand_MAP/model/.rcpp_cache_o2_supply_demand_map"
echo "eFAST role=$role N=$n task=${SLURM_ARRAY_TASK_ID:-1} host=$(hostname -s) start=$(date -Is)"
launch apptainer exec --cleanenv --bind "$scratch:$scratch,$scratch/model_cache:$model_cache" "$sif" \
  env R_HOME= R_ENVIRON_USER=/dev/null R_LIBS_USER=/opt/R/4.4.2/lib64/R/library \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR="$scratch/matplotlib" \
  MININGCLONEID_RCPP_REBUILD=FALSE SOBOL_PERRON_CACHE="$scratch/perron_cache" \
  EFAST_CANONICAL_HELPERS="$helpers" SOBOL_SCRATCH="$scratch" SLURM_JOB_ID="$SLURM_JOB_ID" \
  python3 -u "$snapshot/workflow.py" "$role" --root "$root" --n "$n" --task "${SLURM_ARRAY_TASK_ID:-1}"
echo "eFAST payload completed $(date -Is)"
