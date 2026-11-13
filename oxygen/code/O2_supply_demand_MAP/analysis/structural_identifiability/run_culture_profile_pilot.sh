#!/usr/bin/env bash
set -euo pipefail
project_root="${1:?project root}"
study="${2:?study root}"
mode="${3:?validate or fit}"
task_id="${4:-}"
image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
model="${project_root}/oxygen/results/structural_identifiability/20261109_seed144_wp0/source/code/oxygen/code/O2_supply_demand_MAP/model"
wp0="${project_root}/oxygen/results/structural_identifiability/20261109_seed144_wp0"
code="${study}/implementation/oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability"
case "${mode}" in
  validate) entry=validate_culture_profile_pilot.R; runtime="${study}/runtime/validation"; extra=() ;;
  fit) test -n "${task_id}"; entry=run_culture_profile_pilot.R; runtime="${study}/runtime/${task_id}"; extra=("${task_id}") ;;
  *) exit 64 ;;
esac
mkdir -p "${runtime}/rcpp_cache" "${runtime}/home"
if [[ "${mode}" == fit && -d "${study}/runtime/validation/rcpp_cache" ]]; then
  cp -a "${study}/runtime/validation/rcpp_cache/." "${runtime}/rcpp_cache/"
fi
actual_hash="$(sha256sum "${image}")"
actual_hash="${actual_hash%% *}"
test "${actual_hash}" = dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372
printf '%s\n' "${actual_hash}" > "${runtime}/image_sha256.txt"
trap 'exit_code=$?; printf "%s\n" "$exit_code" > "'"${runtime}"'/launch_exit_code.txt"' EXIT
/usr/bin/time -v apptainer exec --cleanenv \
  --home "${runtime}/home" --bind /share:/share \
  --bind "${wp0}/source:${wp0}/source:ro" \
  --bind "${study}/source:${study}/source:ro" \
  --bind "${study}/implementation:${study}/implementation:ro" \
  --bind "${runtime}/rcpp_cache:${model}/.rcpp_cache_o2_supply_demand_map" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --env O2SD_CONTAINER_IMAGE="${image}" --env O2PROF_IMAGE_SHA256="${actual_hash}" \
  --env SLURM_JOB_ID="${SLURM_JOB_ID:-}" --env RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript \
  --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null \
  "${image}" bash "${project_root}/scripts/agentRrunner.sh" "${code}/${entry}" \
  "${project_root}" "${study}" "${extra[@]}"
