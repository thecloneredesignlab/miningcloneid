#!/usr/bin/env bash
set -euo pipefail
project_root="${1:?project root}"
study="${2:?study root}"
image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
wp0="${project_root}/oxygen/results/structural_identifiability/20261109_seed144_wp0"
model="${wp0}/source/code/oxygen/code/O2_supply_demand_MAP/model"
runtime="${study}/runtime/final_report"
code="${study}/scripts"

# Final outputs require closed fitting tasks and unchanged accepted payloads.
python3 "${code}/audit_culture_pilot_artifacts.py" "${study}" --require-terminal
mkdir -p "${runtime}/rcpp_cache" "${runtime}/home"
cp -a "${study}/runtime/validation/rcpp_cache/." "${runtime}/rcpp_cache/"
actual_hash="$(sha256sum "${image}")"
actual_hash="${actual_hash%% *}"
test "${actual_hash}" = dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372
printf '%s\n' "${actual_hash}" > "${runtime}/image_sha256.txt"
trap 'exit_code=$?; printf "%s\n" "$exit_code" > "'"${runtime}"'/launch_exit_code.txt"' EXIT
apptainer exec --cleanenv --bind /share:/share \
  --env RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null \
  "${image}" bash "${project_root}/scripts/agentRrunner.sh" "${code}/test_culture_offspring_identity.R" "${code}/culture_profile_diagnostics.R"
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
  "${image}" bash "${project_root}/scripts/agentRrunner.sh" "${code}/report_culture_profile_pilot.R" \
  "${project_root}" "${study}"
python3 "${code}/audit_culture_pilot_artifacts.py" "${study}" --accounting --accounting-only
python3 "${code}/summarize_culture_profile_findings.py" "${project_root}" "${study}"
python3 "${code}/summarize_culture_profile_pilot.py" "${project_root}" "${study}"
