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
O2SD_CONTAINER_BITMAP_TYPE="${O2SD_CONTAINER_BITMAP_TYPE:-cairo}"
O2SD_CONTAINER_R_PROFILE="${O2SD_CONTAINER_R_PROFILE:-${O2SD_DOCKER_HPC_ROOT}/util/o2_supply_demand_map_container.Rprofile}"
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
export O2SD_CONTAINER_BITMAP_TYPE
export O2SD_CONTAINER_R_PROFILE
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
  if [[ ! -r "${O2SD_CONTAINER_R_PROFILE}" ]]; then
    echo "Container R profile is not readable: ${O2SD_CONTAINER_R_PROFILE}" >&2
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
    --env "R_PROFILE_USER=${O2SD_CONTAINER_R_PROFILE}"
    --env "R_ENVIRON_USER=/dev/null"
    --env "R_BITMAP_TYPE=${O2SD_CONTAINER_BITMAP_TYPE}"
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

o2sd_container_r_sanity_check() {
  local max_attempts="${1:-3}"
  local retry_delay_seconds="${2:-2}"
  if ! [[ "${max_attempts}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Container R sanity-check attempts must be a positive integer: ${max_attempts}" >&2
    return 2
  fi
  if ! [[ "${retry_delay_seconds}" =~ ^[0-9]+$ ]]; then
    echo "Container R sanity-check delay must be a non-negative integer: ${retry_delay_seconds}" >&2
    return 2
  fi

  local log_path
  log_path="$(mktemp "${TMPDIR:-/tmp}/o2sd-r-sanity.XXXXXX")"
  local attempt
  for ((attempt = 1; attempt <= max_attempts; attempt++)); do
    if Rscript -e '
      expected <- Sys.getenv("R_BITMAP_TYPE", unset = "")
      actual <- getOption("bitmapType")
      if (nzchar(expected) && !identical(actual, expected)) {
        stop("bitmapType mismatch: expected ", expected, ", got ", actual)
      }
      png_path <- tempfile(fileext = ".png")
      on.exit(unlink(png_path), add = TRUE)
      grDevices::png(png_path, width = 64, height = 64)
      graphics::par(mar = rep(0, 4))
      graphics::plot.new()
      grDevices::dev.off()
      if (!file.exists(png_path) || file.info(png_path)$size <= 0) {
        stop("container PNG sanity output was not created")
      }
      cat("Container R sanity check OK; bitmapType=", actual, "\n", sep = "")
    ' >"${log_path}" 2>&1; then
      rm -f "${log_path}"
      return 0
    fi
    echo "Container R sanity check attempt ${attempt}/${max_attempts} failed." >&2
    if (( attempt < max_attempts && retry_delay_seconds > 0 )); then
      sleep "${retry_delay_seconds}"
    fi
  done

  echo "Container R sanity check failed after ${max_attempts} attempts." >&2
  echo "SIF: ${O2SD_CONTAINER_IMAGE}" >&2
  echo "Last Rscript output:" >&2
  sed -n '1,160p' "${log_path}" >&2
  rm -f "${log_path}"
  return 1
}
