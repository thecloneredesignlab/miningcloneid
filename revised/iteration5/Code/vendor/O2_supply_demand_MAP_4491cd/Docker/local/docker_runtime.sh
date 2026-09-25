#!/usr/bin/env bash

# Shared Docker runtime for local fitting and analysis wrappers.

O2SD_DOCKER_LOCAL_ROOT="${O2SD_DOCKER_LOCAL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
O2SD_WORKFLOW_ROOT="${O2SD_WORKFLOW_ROOT:-$(cd "${O2SD_DOCKER_LOCAL_ROOT}/../.." && pwd)}"
O2SD_REPO_ROOT="${O2SD_REPO_ROOT:-$(cd "${O2SD_WORKFLOW_ROOT}/../../.." && pwd)}"
O2SD_DOCKER_IMAGE="${O2SD_DOCKER_IMAGE:-zafiro/o2_supply_demand_map:r44}"
O2SD_DOCKER_PLATFORM="${O2SD_DOCKER_PLATFORM:-linux/amd64}"
O2SD_DOCKER_RCPP_CACHE="${O2SD_DOCKER_RCPP_CACHE:-${TMPDIR:-/tmp}/o2sd-docker-rcpp-cache-$(id -u)}"

export O2SD_DOCKER_LOCAL_ROOT
export O2SD_WORKFLOW_ROOT
export O2SD_REPO_ROOT
export O2SD_DOCKER_IMAGE
export O2SD_DOCKER_PLATFORM
export O2SD_DOCKER_RCPP_CACHE

o2sd_docker_prepare() {
  if ! command -v docker >/dev/null 2>&1; then
    echo "Docker is required for Docker/local execution." >&2
    return 127
  fi
  if ! docker image inspect "${O2SD_DOCKER_IMAGE}" >/dev/null 2>&1; then
    echo "Docker image is not available locally: ${O2SD_DOCKER_IMAGE}" >&2
    echo "Pull it with: docker pull ${O2SD_DOCKER_IMAGE}" >&2
    return 2
  fi
}

o2sd_docker_exec() {
  o2sd_docker_prepare

  local uid_value
  local gid_value
  uid_value="$(id -u)"
  gid_value="$(id -g)"
  mkdir -p "${O2SD_DOCKER_RCPP_CACHE}"

  local -a command=(
    docker run
    --rm
    --init
    --platform "${O2SD_DOCKER_PLATFORM}"
    --user "${uid_value}:${gid_value}"
    --volume "${O2SD_REPO_ROOT}:${O2SD_REPO_ROOT}"
    --volume "${O2SD_DOCKER_RCPP_CACHE}:${O2SD_WORKFLOW_ROOT}/model/.rcpp_cache_o2_supply_demand_map"
    --workdir "${O2SD_REPO_ROOT}"
    --env "HOME=/tmp/o2sd-container-home"
    --env "XDG_CACHE_HOME=/tmp/o2sd-container-home/cache"
    --env "R_PROFILE_USER=/dev/null"
    --env "R_ENVIRON_USER=/dev/null"
    --env "PYTHONUSERBASE=/tmp/o2sd-container-home/.local"
    --env "PYTHONPATH="
    --env "O2SD_CONTAINER_RUNTIME_ACTIVE=TRUE"
  )

  local name
  while IFS= read -r name; do
    case "${name}" in
      O2SD_*|O2_*|O2PL_*|MININGCLONEID_*|FIXO2_*|AGREEMENT_*|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|MKL_NUM_THREADS|VECLIB_MAXIMUM_THREADS)
        command+=(--env "${name}=${!name}")
        ;;
    esac
  done < <(compgen -e)

  if [[ -n "${O2SD_DOCKER_BINDS:-}" ]]; then
    local bind_spec
    local old_ifs="${IFS}"
    IFS=","
    for bind_spec in ${O2SD_DOCKER_BINDS}; do
      [[ -n "${bind_spec}" ]] && command+=(--volume "${bind_spec}")
    done
    IFS="${old_ifs}"
  fi

  command+=("${O2SD_DOCKER_IMAGE}")
  command+=("$@")
  "${command[@]}"
}
