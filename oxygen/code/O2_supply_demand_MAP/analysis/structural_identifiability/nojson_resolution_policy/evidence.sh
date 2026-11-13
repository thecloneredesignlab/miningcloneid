#!/usr/bin/env bash
set -euo pipefail
project="${1:?project}"
study="${2:?study}"
release="${3:?release}"
task="${4:?task}"
image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
wp0="${project}/oxygen/results/structural_identifiability/20261109_seed144_wp0"
parent="${project}/oxygen/results/structural_identifiability/20261109_culture_v2_pilot"
model="${wp0}/source/code/oxygen/code/O2_supply_demand_MAP/model"
runtime="${study}/runtime/${task}_policy"
mkdir -p "${runtime}/rcpp_cache" "${runtime}/home"
cp -a "${parent}/runtime/validation/rcpp_cache/." "${runtime}/rcpp_cache/"
actual_hash="$(sha256sum "${image}")"
actual_hash="${actual_hash%% *}"
if [[ "${actual_hash}" != dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372 ]]; then
  printf '%s\n' 'Image fingerprint mismatch' >&2
  exit 1
fi
printf '%s\n' "${actual_hash}" > "${runtime}/image_sha256.txt"
trap 'rc=$?; printf "%s\n" "$rc" > "'"${runtime}"'/launch_exit_code.txt"' EXIT
/usr/bin/time -v apptainer exec --cleanenv \
  --home "${runtime}/home" --bind /share:/share \
  --bind "${wp0}/source:${wp0}/source:ro" --bind "${parent}:${parent}:ro" \
  --bind "${study}/implementation:${study}/implementation:ro" --bind "${release}:${release}:ro" \
  --bind "${runtime}/rcpp_cache:${model}/.rcpp_cache_o2_supply_demand_map" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --env O2PROF_IMAGE_SHA256="${actual_hash}" --env RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript \
  --env SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-}" --env SLURM_JOB_ID="${SLURM_JOB_ID:-}" \
  --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null \
  "${image}" bash "${project}/scripts/agentRrunner.sh" "${release}/evidence.R" "${project}" "${study}" "${release}" "${task}"
