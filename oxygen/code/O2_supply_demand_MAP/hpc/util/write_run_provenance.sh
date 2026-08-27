#!/bin/bash
# Shared helpers for recording how an O2 fit/report run was launched.

o2sd_prov_shell_join() {
  local out=""
  local token
  local quoted
  for token in "$@"; do
    printf -v quoted "%q" "${token}"
    out+="${quoted} "
  done
  printf "%s" "${out% }"
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
