#!/usr/bin/env bash
set -euo pipefail
study="${1:?study}"; task="${2:?task}"
release_name="$(awk -F '\t' '$1=="release" {print $2}' "${study}/dispatch/${task}/spec.tsv")"
release="${study}/releases/${release_name}"
project=/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability
image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
expected=dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372
wp0="${project}/oxygen/results/structural_identifiability/20261109_seed144_wp0"
archive="${project}/oxygen/results/structural_identifiability/20261112_culture_v2_seven_point_review_nojson"
runtime="${study}/runtime/${task}"
mkdir -p "${runtime}/home" "${runtime}/tmp"
actual="$(sha256sum "${image}")"; actual="${actual%% *}"; [[ "${actual}" == "${expected}" ]]
printf '%s\n' "${actual}" > "${runtime}/image_sha256.txt"
trap 'rc=$?; printf "%s\n" "$rc" > "'"${runtime}"'/launch_exit_code.txt"' EXIT
python3 -B "${release}/workflow.py" verify "${study}" "${release}"
mode="$(awk -F '\t' '$1=="mode" {print $2}' "${study}/dispatch/${task}/spec.tsv")"
case "${mode}" in
  prepare|replay|replay_final|parallel|engineering_check|local_api_check) script=engineering.R;;
  screen_*|global_*|globalcheck_*|local_*|check_*|positive_*) script=fit.R;;
  report|report_resume01) script=report.R;;
  *) exit 64;;
esac
apptainer exec --cleanenv --home "${runtime}/home" \
  --bind /share:/share --bind "${wp0}:${wp0}:ro" --bind "${archive}:${archive}:ro" \
  --bind "${project}/oxygen/results/synthetic_identifiability/20261113_culture_conditional_d0_seed144_02:${project}/oxygen/results/synthetic_identifiability/20261113_culture_conditional_d0_seed144_02:ro" \
  --bind "${study}/releases:${study}/releases:ro" --bind "${study}/inputs:${study}/inputs:ro" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --env PYTHONDONTWRITEBYTECODE=1 --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null \
  --env TMPDIR="${runtime}/tmp" --env D0_RUNTIME="${runtime}" --env D0_RELEASE="${release}" \
  --env D0_IMAGE_SHA256="${actual}" --env RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript \
  --env SLURM_JOB_ID="${SLURM_JOB_ID}" --env SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK}" \
  "${image}" bash "${project}/scripts/agentRrunner.sh" "${release}/${script}" "${study}" "${task}" "${mode}"
python3 -B "${release}/workflow.py" verify "${study}" "${release}"
