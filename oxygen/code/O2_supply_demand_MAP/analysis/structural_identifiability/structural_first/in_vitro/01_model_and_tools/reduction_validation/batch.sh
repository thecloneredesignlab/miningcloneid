#!/usr/bin/env bash
set -euo pipefail
study="${1:?}"; release="${study}/releases/${2:?}"; task="${3:?}"; phase="${4:?}"; case_id="${5:?}"
project=/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability
image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
scratch="$(mktemp -d "/tmp/rr96_${SLURM_JOB_ID}_XXXXXX")"
mkdir -p "${study}/runtime/${task}" "${scratch}/home"
trap 'rc=$?; printf "%s\n" "$rc" > "'"${study}"'/runtime/'"${task}"'/launch_exit_code.txt"' EXIT
python3 -B "${release}/workflow.py" "${study}" verify --release "${2}"
actual="$(sha256sum "${image}")";actual="${actual%% *}"
[[ "${actual}" == dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372 ]]
printf '%s\n' "${actual}" > "${study}/runtime/${task}/image_sha256.txt"
apptainer exec --cleanenv --home "${scratch}/home" --bind /share:/share --bind "${scratch}:${scratch}" \
 --bind "${project}/oxygen/results/structural_identifiability/20261109_seed144_wp0:${project}/oxygen/results/structural_identifiability/20261109_seed144_wp0:ro" \
 --bind "${project}/oxygen/results/structural_identifiability/20261112_culture_v2_seven_point_review_nojson:${project}/oxygen/results/structural_identifiability/20261112_culture_v2_seven_point_review_nojson:ro" \
 --bind "${project}/oxygen/results/structural_identifiability/20261113_culture_symbolic_01:${project}/oxygen/results/structural_identifiability/20261113_culture_symbolic_01:ro" \
 --bind "${project}/oxygen/results/synthetic_identifiability:${project}/oxygen/results/synthetic_identifiability:ro" \
 --bind "${study}/releases:${study}/releases:ro" --bind "${study}/inputs:${study}/inputs:ro" \
 --env D0_RUNTIME="${study}/runtime/${task}" --env D0_IMAGE_SHA256="${actual}" --env RR_RELEASE="${release}" --env TMPDIR="${scratch}" \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
 --env PYTHONDONTWRITEBYTECODE=1 --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null \
 --env SLURM_JOB_ID="${SLURM_JOB_ID}" --env RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript \
 "${image}" bash "${project}/scripts/agentRrunner.sh" "${release}/$([[ "${phase}" == report ]] && printf report || printf science).R" "${study}" "${task}" "${phase}" "${case_id}"
if [[ "${phase}" == report ]]; then python3 -B "${release}/collect.py" "${study}"; fi
python3 -B "${release}/workflow.py" "${study}" verify --release "${2}"
