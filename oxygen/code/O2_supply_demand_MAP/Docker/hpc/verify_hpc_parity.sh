#!/usr/bin/env bash

set -euo pipefail

DOCKER_HPC_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${DOCKER_HPC_ROOT}/../.." && pwd)"
MODULE_HPC_ROOT="${WORKFLOW_ROOT}/hpc"
MAPPING="${DOCKER_HPC_ROOT}/hpc_script_mapping.tsv"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

[[ -f "${MAPPING}" ]] || fail "Missing mapping: ${MAPPING}"

tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/o2sd-hpc-parity.XXXXXX")"
trap 'rm -rf "${tmp_dir}"' EXIT

tail -n +2 "${MAPPING}" | cut -f1 | sort > "${tmp_dir}/mapped-module.txt"
tail -n +2 "${MAPPING}" | cut -f2 | sort > "${tmp_dir}/mapped-container.txt"
find "${MODULE_HPC_ROOT}" -type f ! -name '.DS_Store' \
  | sed "s#^${MODULE_HPC_ROOT}/##" | sort > "${tmp_dir}/actual-module.txt"

diff -u "${tmp_dir}/actual-module.txt" "${tmp_dir}/mapped-module.txt" \
  || fail "hpc_script_mapping.tsv does not cover the complete hpc/ inventory."

resource_pattern='(#SBATCH|--(account|partition|qos|time|mem|mem-per-cpu|cpus-per-task|ntasks|nodes|array|dependency|output|error|job-name)(=|[[:space:]]))'

while IFS=$'\t' read -r module_path container_path; do
  [[ "${module_path}" == "module_path" ]] && continue

  source_path="${MODULE_HPC_ROOT}/${module_path}"
  target_path="${DOCKER_HPC_ROOT}/${container_path}"
  [[ -f "${source_path}" ]] || fail "Missing module source: ${source_path}"
  [[ -f "${target_path}" ]] || fail "Missing container counterpart: ${target_path}"

  case "${target_path}" in
    *.sh|*.sub)
      bash -n "${target_path}"
      grep -q 'o2_supply_demand_map_apptainer_runtime.sh' "${target_path}" \
        || fail "Missing Apptainer runtime initialization: ${target_path}"
      ;;
  esac

  grep -E "${resource_pattern}" "${source_path}" \
    | sed 's#/Docker/hpc#/hpc#g' > "${tmp_dir}/source-resources.txt" || true
  grep -E "${resource_pattern}" "${target_path}" \
    | sed 's#/Docker/hpc#/hpc#g' > "${tmp_dir}/target-resources.txt" || true
  diff -u "${tmp_dir}/source-resources.txt" "${tmp_dir}/target-resources.txt" \
    || fail "Slurm resource configuration changed in ${container_path}."
done < "${MAPPING}"

if grep -R -n -E '^[[:space:]]*(ml|module)[[:space:]]' \
  --include='*.sh' --include='*.sub' "${DOCKER_HPC_ROOT}"; then
  fail "Docker/hpc contains an active host module command."
fi

echo "HPC parity verification passed: 28 mapped files; Slurm resources unchanged."
