#!/usr/bin/env bash

# Shared shell helpers for local runners and HPC launchers. Sourcing this file
# defines functions only; it does not submit jobs, load modules, or run a stage.

truthy() {
  case "${1:-FALSE}" in
    TRUE|true|True|1|yes|YES|y|Y|on|ON) return 0 ;;
    *) return 1 ;;
  esac
}

is_null_value() {
  local val
  val="$(echo "${1:-}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')"
  [[ -z "${val}" || "${val}" == "null" || "${val}" == "none" || "${val}" == "na" ]]
}

normalize_fitting_mode() {
  local val
  val="$(echo "${1:-}" | tr '[:upper:]' '[:lower:]')"
  val="${val// /}"
  val="${val//-/}"
  val="${val//_/}"
  case "${val}" in
    invivo) echo "invivo" ;;
    invitro) echo "invitro" ;;
    joint) echo "joint" ;;
    *) echo "" ;;
  esac
}

require_positive_int() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]] || (( value <= 0 )); then
    echo "${name} must be a positive integer, got: ${value}" >&2
    exit 2
  fi
}

require_nonnegative_int() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]]; then
    echo "${name} must be a non-negative integer, got: ${value}" >&2
    exit 2
  fi
}

log_msg() {
  local msg="$1"
  local stamp
  stamp="$(date '+%Y-%m-%d %H:%M:%S')"
  printf '[%s] %s\n' "${stamp}" "${msg}" | tee -a "${PROGRESS_LOG}"
}

print_command() {
  local label="$1"
  shift
  if [[ -n "${PROGRESS_LOG:-}" ]]; then
    printf "%s:" "${label}" | tee -a "${PROGRESS_LOG}"
    printf " %q" "$@" | tee -a "${PROGRESS_LOG}"
    printf "\n" | tee -a "${PROGRESS_LOG}"
  else
    printf "%s:" "${label}"
    printf " %q" "$@"
    printf "\n"
  fi
}

run_or_print() {
  local label="$1"
  shift
  print_command "${label}" "$@"
  if ! truthy "${DRY_RUN}"; then
    "$@"
  fi
}

shell_join() {
  local out=""
  local token
  local quoted
  for token in "$@"; do
    printf -v quoted "%q" "${token}"
    out+="${quoted} "
  done
  printf "%s" "${out% }"
}

load_r_module() {
  if [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi
  if command -v module >/dev/null 2>&1; then
    module use /app/eb/modules/all >/dev/null 2>&1 || true
  fi
  if [[ -z "${R_MODULE:-}" ]]; then
    return 0
  fi
  if truthy "${O2SD_R_MODULE_OPTIONAL:-FALSE}"; then
    if command -v ml >/dev/null 2>&1; then
      ml "${R_MODULE}" >/dev/null 2>&1 || true
    elif command -v module >/dev/null 2>&1; then
      module load "${R_MODULE}" >/dev/null 2>&1 || true
    fi
  elif command -v ml >/dev/null 2>&1; then
    ml "${R_MODULE}"
  elif command -v module >/dev/null 2>&1; then
    module load "${R_MODULE}"
  fi
  return 0
}

derive_joint_warmup_seed_label() {
  local invivo_label
  local invitro_label
  invivo_label="$(sanitize_label "$(basename "${INVIVO_BEST_SEED_DIR}")")"
  invitro_label="$(sanitize_label "$(basename "${INVITRO_BEST_SEED_DIR}")")"
  printf "invivo_%s__invitro_%s" "${invivo_label}" "${invitro_label}"
}

label_joint_run_prefix() {
  if truthy "${JOINT_WARMUP_ENABLE}" && ! is_null_value "${JOINT_WARMUP_SEED_LABEL}"; then
    if [[ "${JOINT_RUN_PREFIX}" != *"${JOINT_WARMUP_SEED_LABEL}"* ]]; then
      JOINT_RUN_PREFIX="${JOINT_RUN_PREFIX}__${JOINT_WARMUP_SEED_LABEL}"
    fi
  fi
}

resolve_existing_dir() {
  local label="$1"
  local path="$2"
  if [[ "${path}" != /* && -d "${PROJECT_ROOT}/${path}" ]]; then
    path="${PROJECT_ROOT}/${path}"
  fi
  if [[ ! -d "${path}" ]]; then
    echo "Missing ${label}: ${path}" >&2
    exit 1
  fi
  (cd "${path}" && pwd)
}

first_line() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    echo "Missing file: ${path}" >&2
    exit 1
  fi
  local line
  line="$(head -n 1 "${path}")"
  line="$(printf "%s" "${line}" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  if [[ -z "${line}" ]]; then
    echo "Empty first line in ${path}" >&2
    exit 1
  fi
  printf "%s" "${line}"
}

o2sd_prov_shell_join() {
  shell_join "$@"
}

o2sd_prov_cell() {
  local val="${1-}"
  val="${val//$'\t'/ }"
  val="${val//$'\r'/ }"
  val="${val//$'\n'/ }"
  printf "%s" "${val}"
}

o2sd_prov_init() {
  local run_dir="$1"
  mkdir -p "${run_dir}"
  local path="${run_dir}/run_provenance.tsv"
  if [[ ! -f "${path}" ]]; then
    printf "section\tkey\tvalue\n" > "${path}"
  fi
}

o2sd_prov_append() {
  local run_dir="$1"
  local section="$2"
  local key="$3"
  local value="${4-}"
  o2sd_prov_init "${run_dir}"
  printf "%s\t%s\t%s\n" \
    "$(o2sd_prov_cell "${section}")" \
    "$(o2sd_prov_cell "${key}")" \
    "$(o2sd_prov_cell "${value}")" >> "${run_dir}/run_provenance.tsv"
}

o2sd_prov_write_many() {
  local run_dir="$1"
  shift
  while (( $# >= 3 )); do
    o2sd_prov_append "${run_dir}" "$1" "$2" "$3"
    shift 3
  done
}

o2sd_prov_write_text() {
  local path="$1"
  shift
  mkdir -p "$(dirname "${path}")"
  printf "%s\n" "$*" > "${path}"
}

o2sd_prov_write_standard() {
  local run_dir="$1"
  local entrypoint="${2-}"
  local command_text="${3-}"
  o2sd_prov_write_many "${run_dir}" \
    execution timestamp "$(date '+%Y-%m-%d %H:%M:%S %Z')" \
    execution hostname "$(hostname 2>/dev/null || printf NA)" \
    execution user "${USER:-NA}" \
    execution cwd "$(pwd)" \
    execution entrypoint "${entrypoint}" \
    execution run_dir "${run_dir}"
  if [[ -n "${command_text}" ]]; then
    o2sd_prov_write_text "${run_dir}/run_command.txt" "${command_text}"
    o2sd_prov_append "${run_dir}" execution run_command_file "${run_dir}/run_command.txt"
  fi
  if command -v git >/dev/null 2>&1; then
    local git_root=""
    git_root="$(git -C "${run_dir}" rev-parse --show-toplevel 2>/dev/null || true)"
    if [[ -z "${git_root}" && -n "${PROJECT_ROOT:-}" ]]; then
      git_root="$(git -C "${PROJECT_ROOT}" rev-parse --show-toplevel 2>/dev/null || true)"
    fi
    if [[ -n "${git_root}" ]]; then
      o2sd_prov_write_many "${run_dir}" \
        git root "${git_root}" \
        git branch "$(git -C "${git_root}" rev-parse --abbrev-ref HEAD 2>/dev/null || printf NA)" \
        git commit "$(git -C "${git_root}" rev-parse HEAD 2>/dev/null || printf NA)" \
        git dirty_status "$(if git -C "${git_root}" diff --quiet --ignore-submodules -- 2>/dev/null && git -C "${git_root}" diff --cached --quiet --ignore-submodules -- 2>/dev/null; then printf clean; else printf dirty; fi)"
    fi
  fi
}

o2sd_prov_record_sbatch() {
  local run_dir="$1"
  local command_text="$2"
  local job_id="${3-NA}"
  o2sd_prov_write_text "${run_dir}/sbatch_command.txt" "${command_text}"
  o2sd_prov_write_many "${run_dir}" \
    slurm sbatch_command_file "${run_dir}/sbatch_command.txt" \
    slurm job_id "${job_id}"
}
