#!/usr/bin/env bash
set -euo pipefail
project="${1:?project}"; study="${2:?study}"; phase="${3:?phase}"
image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
expected=dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372
model="${project}/oxygen/results/structural_identifiability/20261109_seed144_wp0/source/code/oxygen/code/O2_supply_demand_MAP/model"
wp0="${project}/oxygen/results/structural_identifiability/20261109_seed144_wp0"
archive="${project}/oxygen/results/structural_identifiability/20261112_culture_v2_seven_point_review_nojson"
actual="$(sha256sum "${image}")"; actual="${actual%% *}"; [[ "${actual}" == "${expected}" ]]
[[ "${SLURM_CPUS_PER_TASK}" == 1 ]]
mkdir -p "${study}/results/${phase}"
printf '%s\n' "${actual}" > "${study}/results/${phase}/image_sha256.txt"
trap 'rc=$?; printf "%s\n" "$rc" > "'"${study}"'/results/'"${phase}"'/launch_exit_code.txt"' EXIT
inside() {
  local task="$1"; shift
  local runtime="${study}/runtime/${task}"
  mkdir -p "${runtime}/home" "${runtime}/rcpp_cache" "${runtime}/tmp"
  apptainer exec --cleanenv --home "${runtime}/home" \
    --bind /share:/share --bind "${wp0}/source:${wp0}/source:ro" \
    --bind "${archive}:${archive}:ro" --bind "${study}/inputs:${study}/inputs:ro" \
    --bind "${study}/implementation:${study}/implementation:ro" \
    --bind "${runtime}/rcpp_cache:${model}/.rcpp_cache_o2_supply_demand_map" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
    --env PYTHONDONTWRITEBYTECODE=1 --env TMPDIR="${runtime}/tmp" \
    --env RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript \
    --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null \
    --env SIF_VALIDATION_IMAGE_SHA256="${actual}" \
    --env SLURM_JOB_ID="${SLURM_JOB_ID}" --env SLURM_CPUS_PER_TASK=1 \
    "${image}" "$@"
}
case "${phase}" in
metadata)
  inside metadata_check python3 "${study}/implementation/check.py" verify "${project}" "${study}" metadata
  inside metadata_r bash "${project}/scripts/agentRrunner.sh" "${study}/implementation/metadata.R" "${project}" "${study}"
  ;;
replay)
  test -f "${study}/results/metadata/COMPLETE"
  inside replay_before python3 "${study}/implementation/check.py" verify "${project}" "${study}" replay_before
  for case_id in production_seed144 extended_free induction_zero survival_zero; do
    for repeat_id in 1 2; do
      task="${case_id}_rep${repeat_id}"
      test ! -d "${study}/results/full_forward_slots/${task}"
      inside "${task}" bash "${project}/scripts/agentRrunner.sh" "${study}/implementation/replay.R" "${project}" "${study}" "${case_id}" "${repeat_id}"
    done
  done
  inside replay_after python3 "${study}/implementation/check.py" verify "${project}" "${study}" replay_after
  ;;
summary)
  inside independent_compare bash "${project}/scripts/agentRrunner.sh" "${study}/implementation/compare.R" "${project}" "${study}"
  inside summary python3 "${study}/implementation/check.py" summary "${project}" "${study}" summary
  ;;
*) printf '%s\n' 'Unsupported validation mode' >&2; exit 64;;
esac
after="$(sha256sum "${image}")"; after="${after%% *}"; [[ "${after}" == "${expected}" ]]
