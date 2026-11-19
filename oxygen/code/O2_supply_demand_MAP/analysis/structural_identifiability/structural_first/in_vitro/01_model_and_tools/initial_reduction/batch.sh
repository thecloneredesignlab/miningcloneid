#!/usr/bin/env bash
set -euo pipefail
study="${1:?}"; release="${study}/releases/${2:?}"; task="${3:?}"; mode="${4:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 PYTHONDONTWRITEBYTECODE=1
scratch="$(mktemp -d "/tmp/si4_${SLURM_JOB_ID}_XXXXXX")"
export JULIA_DEPOT_PATH="${scratch}/depot"
export TMPDIR="${scratch}"
mkdir -p "${study}/runtime/${task}" "${study}/results/${task}"
trap 'rc=$?; printf "%s\n" "$rc" > "'"${study}"'/runtime/'"${task}"'/launch_exit_code.txt"' EXIT
python3 -B "${release}/workflow.py" "${study}" verify --release "${2}"
if [[ "${mode}" != tool_setup ]]; then python3 -B "${release}/verify_inputs.py" "${study}"; fi
case "${mode}" in
 tool_setup) python3 -B "${release}/install_tool.py" "${study}"; "${study}/toolchain/julia-1.10.13/bin/julia" --startup-file=no --history-file=no --project="${study}/toolchain/environment" "${release}/tool_setup.jl" "${study}/results/${task}"; python3 -B "${release}/freeze_tool.py" "${study}" "${scratch}/depot" "${task}";;
 tool_postfreeze|julia_global)
  python3 -B "${release}/unpack_tool.py" "${study}" "${scratch}"
  if [[ "${mode}" == tool_postfreeze ]]; then
   "${study}/toolchain/julia-1.10.13/bin/julia" --startup-file=no --history-file=no --project="${study}/toolchain/environment" "${release}/tool_smoke.jl" "${study}/results/${task}"
  else python3 -B "${release}/run_julia_case.py" "${study}" "${release}" "${task}"; fi;;
 reduction)
  image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
  actual="$(sha256sum "${image}")";actual="${actual%% *}";[[ "${actual}" == dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372 ]]
  printf '%s\n' "${actual}" > "${study}/runtime/${task}/image_sha256.txt"
  mkdir -p "${scratch}/home"
  project=/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability
  apptainer exec --cleanenv --home "${scratch}/home" --bind /share:/share --bind "${scratch}:${scratch}" \
   --bind "${project}/oxygen/results/structural_identifiability/20261109_seed144_wp0:${project}/oxygen/results/structural_identifiability/20261109_seed144_wp0:ro" \
   --bind "${project}/oxygen/results/structural_identifiability/20261112_culture_v2_seven_point_review_nojson:${project}/oxygen/results/structural_identifiability/20261112_culture_v2_seven_point_review_nojson:ro" \
   --bind "${study}/releases:${study}/releases:ro" --bind "${study}/inputs:${study}/inputs:ro" \
   --env D0_RUNTIME="${study}/runtime/${task}" --env D0_IMAGE_SHA256="${actual}" --env SI4_RELEASE="${release}" --env TMPDIR="${scratch}" \
   --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
   --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null --env RSCRIPT_BIN=/usr/local/bin/o2-hpc-exact-rscript \
   "${image}" bash "${project}/scripts/agentRrunner.sh" "${release}/reduction.R" "${study}" "${task}";;
 exact_analysis) python3 -B "${release}/certify.py" "${study}" "${task}";python3 -B "${release}/verify_certificates.py" "${study}/results/${task}";;
 *) echo "Unregistered Step04 mode" >&2; exit 64;;
esac
if [[ "${mode}" != tool_setup ]]; then python3 -B "${release}/verify_inputs.py" "${study}"; fi
python3 -B "${release}/workflow.py" "${study}" verify --release "${2}"
