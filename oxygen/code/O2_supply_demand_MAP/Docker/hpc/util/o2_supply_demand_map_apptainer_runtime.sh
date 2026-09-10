#!/usr/bin/env bash

# Shared runtime for the Docker-image-backed Slurm scripts.
#
# Slurm submission and task orchestration stay on the HPC host. R and Python
# commands are intercepted through Docker/hpc/bin and executed in the verified
# Apptainer SIF.

O2SD_DOCKER_HPC_ROOT="${O2SD_DOCKER_HPC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
O2SD_WORKFLOW_ROOT="${O2SD_WORKFLOW_ROOT:-$(cd "${O2SD_DOCKER_HPC_ROOT}/../.." && pwd)}"
O2SD_CONTAINER_IMAGE="${O2SD_CONTAINER_IMAGE:-/share/lab_crd/lab_crd/taoli/Docker/o2_supply_demand_map_r44.sif}"
O2SD_CONTAINER_R_LIBRARY="${O2SD_CONTAINER_R_LIBRARY:-/opt/R/4.4.2/lib64/R/library}"
O2SD_CONTAINER_RCPP_CACHE="${O2SD_CONTAINER_RCPP_CACHE:-/tmp/o2sd-rcpp-cache-${UID:-$(id -u)}}"
O2SD_CONTAINER_RUNTIME_ACTIVE=TRUE

case ":${PATH}:" in
  *":${O2SD_DOCKER_HPC_ROOT}/bin:"*) ;;
  *) PATH="${O2SD_DOCKER_HPC_ROOT}/bin:${PATH}" ;;
esac

export O2SD_DOCKER_HPC_ROOT
export O2SD_WORKFLOW_ROOT
export O2SD_CONTAINER_IMAGE
export O2SD_CONTAINER_R_LIBRARY
export O2SD_CONTAINER_RCPP_CACHE
export O2SD_CONTAINER_RUNTIME_ACTIVE
export PATH

o2sd_container_prepare() {
  if ! command -v apptainer >/dev/null 2>&1; then
    echo "apptainer is required for Docker/hpc execution." >&2
    return 127
  fi
  if [[ ! -r "${O2SD_CONTAINER_IMAGE}" ]]; then
    echo "Container SIF is not readable: ${O2SD_CONTAINER_IMAGE}" >&2
    return 2
  fi
}

o2sd_container_ignore_r_module() {
  local requested="${1:-}"
  if [[ -n "${requested}" && "${O2SD_CONTAINER_MODULE_NOTICE_SHOWN:-FALSE}" != "TRUE" ]]; then
    echo "Container runtime active; ignoring host R module request: ${requested}" >&2
    O2SD_CONTAINER_MODULE_NOTICE_SHOWN=TRUE
    export O2SD_CONTAINER_MODULE_NOTICE_SHOWN
  fi
  o2sd_container_prepare
}

o2sd_apptainer_exec() {
  o2sd_container_prepare

  local container_home="/tmp/o2sd-container-home-${UID:-$(id -u)}"
  mkdir -p "${container_home}/cache"
  mkdir -p "${O2SD_CONTAINER_RCPP_CACHE}"

  local -a cmd=(
    apptainer exec
    --cleanenv
    --home "${container_home}"
    --env "XDG_CACHE_HOME=${container_home}/cache"
    --env "R_LIBS_USER=${O2SD_CONTAINER_R_LIBRARY}"
    --env "R_PROFILE_USER=/dev/null"
    --env "R_ENVIRON_USER=/dev/null"
    --env "PYTHONUSERBASE=${container_home}/.local"
    --env "PYTHONPATH="
    --env "O2SD_CONTAINER_RUNTIME_ACTIVE=TRUE"
    --env "O2SD_CONTAINER_IMAGE=${O2SD_CONTAINER_IMAGE}"
  )

  local name
  while IFS= read -r name; do
    case "${name}" in
      SLURM_*|O2SD_*|O2_*|O2PL_*|MININGCLONEID_*|FIXO2_*|AGREEMENT_*|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|MKL_NUM_THREADS|VECLIB_MAXIMUM_THREADS)
        cmd+=(--env "${name}=${!name}")
        ;;
    esac
  done < <(compgen -e)

  if [[ -d /share ]]; then
    cmd+=(--bind /share:/share)
  fi
  if [[ -n "${PROJECT_ROOT:-}" && -d "${PROJECT_ROOT}" ]]; then
    cmd+=(--bind "${PROJECT_ROOT}:${PROJECT_ROOT}" --pwd "${PROJECT_ROOT}")
  elif [[ -d "${PWD}" ]]; then
    cmd+=(--bind "${PWD}:${PWD}" --pwd "${PWD}")
  fi
  cmd+=(
    --bind
    "${O2SD_CONTAINER_RCPP_CACHE}:${O2SD_WORKFLOW_ROOT}/model/.rcpp_cache_o2_supply_demand_map"
  )

  if [[ -n "${O2SD_CONTAINER_BINDS:-}" ]]; then
    local bind_path
    local old_ifs="${IFS}"
    IFS=","
    for bind_path in ${O2SD_CONTAINER_BINDS}; do
      [[ -n "${bind_path}" ]] && cmd+=(--bind "${bind_path}")
    done
    IFS="${old_ifs}"
  fi

  cmd+=("${O2SD_CONTAINER_IMAGE}")
  cmd+=("$@")
  "${cmd[@]}"
}
