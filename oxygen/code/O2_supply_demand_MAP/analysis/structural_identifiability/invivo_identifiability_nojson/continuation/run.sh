#!/usr/bin/env bash
set -euo pipefail
study="${1:?study}"; mode="${2:?mode}"; task="${3:?task}"; release="${4:?release}"
allocation_start="$(date +%s)"
project="${study%/oxygen/results/structural_identifiability/*}"
image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
model="${study}/source/numerical/model"
runtime="${study}/runtime/${task}"
mkdir -p "${runtime}/home" "${runtime}/tmp" "${runtime}/rcpp_cache"
cache_source="${project}/oxygen/results/structural_identifiability/20261111_invivo_identifiability_202133Z/cache_sources/validated01"
if [[ -d "${study}/cache_sources/validated01" ]]; then cache_source="${study}/cache_sources/validated01"; fi
cp -a "${cache_source}/." "${runtime}/rcpp_cache/"
python3 - "${runtime}/rcpp_cache" "${model}/.rcpp_cache_o2_supply_demand_map" <<'PY'
from pathlib import Path
import re,sys
for p in Path(sys.argv[1]).rglob('*.cpp.R'):
    text=p.read_text()
    text=re.sub(r"/[^'\"\n]*?/\.rcpp_cache_o2_supply_demand_map/",sys.argv[2]+'/',text)
    p.write_text(text)
PY
actual_hash="$(sha256sum "${image}")"; actual_hash="${actual_hash%% *}"
test "${actual_hash}" = dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372
printf '%s\n' "${actual_hash}" > "${runtime}/image_sha256.txt"
trap 'rc=$?; printf "%s\n" "$rc" > "'"${runtime}"'/launch_exit_code.txt"' EXIT
/usr/bin/time -v apptainer exec --cleanenv --home "${runtime}/home" --bind /share:/share \
 --bind "${study}/source:${study}/source:ro" --bind "${study}/releases/${release}:${study}/releases/${release}:ro" \
 --bind "${study}/manifests:${study}/manifests:ro" --bind "${study}/reference:${study}/reference:ro" \
 --bind "${study}/prior:${study}/prior:ro" --bind "${study}/dispatch/${task}:${study}/dispatch/${task}:ro" \
 --bind "${runtime}/tmp:/tmp" --bind "${runtime}/rcpp_cache:${model}/.rcpp_cache_o2_supply_demand_map" \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
 --env RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript --env TMPDIR=/tmp \
 --env IVI_ALLOCATION_START="${allocation_start}" --env IVI_EVALUATION_MODE=fresh_fork --env IVI_IMAGE_SHA256="${actual_hash}" \
 --env SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-}" --env SLURM_JOB_ID="${SLURM_JOB_ID:-}" \
 --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null \
 "${image}" bash "${project}/scripts/agentRrunner.sh" "${study}/releases/${release}/run.R" "${study}" "${mode}" "${task}"
