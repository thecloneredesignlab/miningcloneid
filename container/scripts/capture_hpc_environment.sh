#!/usr/bin/env bash
# Capture the environment seen by a Slurm compute node for the soft_coupling O2 workflow.

#SBATCH --job-name=o2_env_audit
#SBATCH --qos=small
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

if [[ "$#" -lt 3 ]]; then
  echo "Usage: capture_hpc_environment.sh AUDIT_ROOT STATIC_MANIFEST_DIR CODE_SNAPSHOT_DIR" >&2
  exit 2
fi

AUDIT_ROOT="$1"
STATIC_MANIFEST_DIR="$2"
CODE_SNAPSHOT_DIR="$3"
SCRIPT_DIR="${AUDIT_ROOT}/scripts"
SNAPSHOT_DIR="${AUDIT_ROOT}/hpc_snapshot"
SYSTEM_DIR="${SNAPSHOT_DIR}/system"
MODULE_DIR="${SNAPSHOT_DIR}/modules"
R_DIR="${SNAPSHOT_DIR}/r"
PYTHON_DIR="${SNAPSHOT_DIR}/python"
LIB_DIR="${SNAPSHOT_DIR}/shared_libraries"
LOG_DIR="${SNAPSHOT_DIR}/logs"

mkdir -p "${SYSTEM_DIR}" "${MODULE_DIR}" "${R_DIR}" "${PYTHON_DIR}" "${LIB_DIR}" "${LOG_DIR}"

run_capture() {
  local output="$1"
  shift
  {
    echo "COMMAND: $*"
    "$@"
  } >"${output}" 2>&1 || true
}

capture_shell() {
  local output="$1"
  shift
  {
    echo "COMMAND: $*"
    "$@"
  } >"${output}" 2>&1 || true
}

{
  printf "capture_utc\t%s\n" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  printf "hostname\t%s\n" "$(hostname -f 2>/dev/null || hostname)"
  printf "slurm_job_id\t%s\n" "${SLURM_JOB_ID:-NA}"
  printf "slurm_job_name\t%s\n" "${SLURM_JOB_NAME:-NA}"
  printf "slurm_node_list\t%s\n" "${SLURM_JOB_NODELIST:-NA}"
  printf "slurm_partition\t%s\n" "${SLURM_JOB_PARTITION:-NA}"
  printf "slurm_cpus_per_task\t%s\n" "${SLURM_CPUS_PER_TASK:-NA}"
  printf "slurm_mem_per_node\t%s\n" "${SLURM_MEM_PER_NODE:-NA}"
  printf "requested_r_module\tR/4.4\n"
} >"${SNAPSHOT_DIR}/capture_metadata.tsv"

cp /etc/os-release "${SYSTEM_DIR}/os-release.txt" 2>/dev/null || true
cp /etc/redhat-release "${SYSTEM_DIR}/redhat-release.txt" 2>/dev/null || true
run_capture "${SYSTEM_DIR}/uname.txt" uname -a
run_capture "${SYSTEM_DIR}/architecture.txt" uname -m
run_capture "${SYSTEM_DIR}/lscpu.txt" lscpu
run_capture "${SYSTEM_DIR}/glibc.txt" getconf GNU_LIBC_VERSION
run_capture "${SYSTEM_DIR}/ldd-version.txt" ldd --version
run_capture "${SYSTEM_DIR}/ulimit.txt" bash -lc "ulimit -a"
run_capture "${SYSTEM_DIR}/filesystem.txt" df -Th

if [[ -f /etc/profile.d/modules.sh ]]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi
if command -v module >/dev/null 2>&1; then
  module use /app/eb/modules/all >/dev/null 2>&1 || true
fi

capture_shell "${MODULE_DIR}/module-avail-R.txt" bash -lc "module -t avail R 2>&1 || true"
capture_shell "${MODULE_DIR}/module-avail-Python.txt" bash -lc "module -t avail Python 2>&1 || true"
capture_shell "${MODULE_DIR}/module-spider-R-4.4.txt" bash -lc "module spider R/4.4 2>&1 || true"

if command -v ml >/dev/null 2>&1; then
  ml R/4.4
else
  module load R/4.4
fi

capture_shell "${MODULE_DIR}/module-list-after-R-4.4.txt" bash -lc "module list 2>&1 || true"
capture_shell "${MODULE_DIR}/module-show-R-4.4.txt" bash -lc "module show R/4.4 2>&1 || true"
capture_shell "${MODULE_DIR}/module-show-loaded-R.txt" bash -lc "module show \"\${LMOD_FAMILY_R_VERSION:-R/4.4}\" 2>&1 || true"

{
  printf "variable\tvalue\n"
  for name in EBROOTR EBVERSIONR LMOD_FAMILY_R LMOD_FAMILY_R_VERSION R_HOME R_LIBS R_LIBS_USER MODULEPATH LD_LIBRARY_PATH PATH OMP_NUM_THREADS; do
    value="${!name:-}"
    printf "%s\t%s\n" "${name}" "${value//$'\t'/ }"
  done
} >"${MODULE_DIR}/selected-environment.tsv"

run_capture "${R_DIR}/which-R.txt" bash -lc "type -a R Rscript"
run_capture "${R_DIR}/R-version-command.txt" R --version
run_capture "${R_DIR}/Rscript-version-command.txt" Rscript --version

{
  printf "key\tvalue\n"
  for key in CC CXX CXX11 CXX14 CXX17 FC F77 CPPFLAGS CFLAGS CXXFLAGS FFLAGS FCFLAGS LDFLAGS BLAS_LIBS LAPACK_LIBS SHLIB_LD SHLIB_CXXLD; do
    value="$(R CMD config "${key}" 2>&1 || true)"
    printf "%s\t%s\n" "${key}" "${value//$'\t'/ }"
  done
} >"${R_DIR}/r-cmd-config.tsv"

Rscript "${SCRIPT_DIR}/capture_r_environment.R" \
  "${R_DIR}" \
  "${STATIC_MANIFEST_DIR}/r-direct-packages.tsv" \
  "${STATIC_MANIFEST_DIR}/r-package-usage-by-file.tsv" \
  >"${LOG_DIR}/capture-r.log" 2>&1

{
  printf "tool\tpath\tversion_first_line\n"
  for tool in gcc g++ gfortran make cmake pkg-config pandoc quarto magick convert identify git curl openssl tar gzip bzip2 xz java; do
    path="$(command -v "${tool}" 2>/dev/null || true)"
    if [[ -n "${path}" ]]; then
      version="$("${tool}" --version 2>&1 | head -n 1 || true)"
      printf "%s\t%s\t%s\n" "${tool}" "${path}" "${version//$'\t'/ }"
    else
      printf "%s\tNA\tNA\n" "${tool}"
    fi
  done
} >"${SYSTEM_DIR}/external-tools.tsv"

if command -v rpm >/dev/null 2>&1; then
  rpm -qa --qf '%{NAME}\t%{EPOCHNUM}\t%{VERSION}\t%{RELEASE}\t%{ARCH}\n' \
    | LC_ALL=C sort >"${SYSTEM_DIR}/system-packages-all.tsv" || true
else
  printf "RPM command unavailable\n" >"${SYSTEM_DIR}/system-packages-all.tsv"
fi

{
  printf "kind\tpath\tldd_output\n"
  while IFS= read -r lib_path; do
    [[ -d "${lib_path}" ]] || continue
    while IFS= read -r so_path; do
      while IFS= read -r ldd_line; do
        printf "R_package\t%s\t%s\n" "${so_path}" "${ldd_line//$'\t'/ }"
      done < <(ldd "${so_path}" 2>&1 || true)
    done < <(find "${lib_path}" -type f \( -name '*.so' -o -name '*.dylib' \) 2>/dev/null)
  done < <(Rscript --vanilla -e 'writeLines(.libPaths())')
} >"${LIB_DIR}/r-package-shared-libraries-ldd.tsv"

PYTHON_BIN="$(command -v python3 || true)"
{
  printf "name\tpath\tversion\n"
  for name in python python3 pip pip3; do
    path="$(command -v "${name}" 2>/dev/null || true)"
    version=""
    if [[ -n "${path}" ]]; then
      version="$("${path}" --version 2>&1 | head -n 1 || true)"
    fi
    printf "%s\t%s\t%s\n" "${name}" "${path:-NA}" "${version:-NA}"
  done
} >"${PYTHON_DIR}/python-command-resolution.tsv"

if [[ -n "${PYTHON_BIN}" ]]; then
  "${PYTHON_BIN}" "${SCRIPT_DIR}/capture_python_environment.py" \
    "${PYTHON_DIR}" \
    "${STATIC_MANIFEST_DIR}/python-direct-modules.tsv" \
    >"${LOG_DIR}/capture-python.log" 2>&1 || true

  {
    printf "file\tstatus\tmessage\n"
    while IFS= read -r py_file; do
      message="$("${PYTHON_BIN}" -m py_compile "${py_file}" 2>&1 || true)"
      if "${PYTHON_BIN}" -m py_compile "${py_file}" >/dev/null 2>&1; then
        status="ok"
      else
        status="failed"
      fi
      printf "%s\t%s\t%s\n" "${py_file#${CODE_SNAPSHOT_DIR}/}" "${status}" "${message//$'\t'/ }"
    done < <(find "${CODE_SNAPSHOT_DIR}" -type f -name '*.py' | LC_ALL=C sort)
  } >"${PYTHON_DIR}/python-script-compile-status.tsv"

  {
    printf "kind\tpath\tldd_output\n"
    while IFS= read -r so_path; do
      while IFS= read -r ldd_line; do
        printf "Python_extension\t%s\t%s\n" "${so_path}" "${ldd_line//$'\t'/ }"
      done < <(ldd "${so_path}" 2>&1 || true)
    done < <("${PYTHON_BIN}" - <<'PY'
import os
import site
paths = []
try:
    paths.extend(site.getsitepackages())
except Exception:
    pass
try:
    paths.append(site.getusersitepackages())
except Exception:
    pass
seen = set()
for root in paths:
    if not root or root in seen or not os.path.isdir(root):
        continue
    seen.add(root)
    for current, _dirs, files in os.walk(root):
        for name in files:
            if name.endswith((".so", ".dylib")):
                print(os.path.join(current, name))
PY
    )
  } >"${LIB_DIR}/python-extension-shared-libraries-ldd.tsv"
else
  printf "python3 not found after loading R/4.4\n" >"${LOG_DIR}/capture-python.log"
fi

awk -F '\t' '
  NR > 1 {
    line = $3
    if (index(line, " => ") > 0) {
      split(line, arrow, " => ")
      split(arrow[2], fields, " ")
      if (substr(fields[1], 1, 1) == "/") print fields[1]
    } else {
      split(line, fields, " ")
      if (substr(fields[1], 1, 1) == "/") print fields[1]
    }
  }
' "${LIB_DIR}/r-package-shared-libraries-ldd.tsv" \
  "${LIB_DIR}/python-extension-shared-libraries-ldd.tsv" \
  | LC_ALL=C sort -u >"${LIB_DIR}/resolved-shared-library-paths.txt"

{
  printf "library_path\trpm_owner\n"
  while IFS= read -r library_path; do
    owner="NA"
    if command -v rpm >/dev/null 2>&1; then
      owner="$(rpm -qf "${library_path}" 2>&1 || true)"
    fi
    printf "%s\t%s\n" "${library_path}" "${owner//$'\t'/ }"
  done <"${LIB_DIR}/resolved-shared-library-paths.txt"
} >"${LIB_DIR}/shared-library-rpm-owners.tsv"

run_capture "${LIB_DIR}/R-binary-ldd.txt" ldd "$(command -v R)"
if [[ -n "${PYTHON_BIN}" ]]; then
  run_capture "${LIB_DIR}/python-binary-ldd.txt" ldd "${PYTHON_BIN}"
fi

find "${SNAPSHOT_DIR}" -type f ! -name SHA256SUMS -print0 \
  | LC_ALL=C sort -z \
  | xargs -0 sha256sum >"${SNAPSHOT_DIR}/SHA256SUMS"

printf "HPC environment capture completed: %s\n" "${SNAPSHOT_DIR}"
