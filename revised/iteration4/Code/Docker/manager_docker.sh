#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ITERATION_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
REPO_ROOT="$(cd "${ITERATION_ROOT}/../.." && pwd -P)"
MANAGER="${ITERATION_ROOT}/manager.sh"
FONT_FILE="${REPO_ROOT}/data/fonts/Arial Bold.ttf"
EXPECTED_IMAGE_ID="sha256:32c49db0ad27a0b5832b601ba96e2b72bfc1e2f1ccbf34687f8f596f1f7cdcd5"

usage() {
  cat <<'EOF'
Usage:
  manager_docker.sh \
    --docker-image=IMAGE \
    --invitro-result-root=PATH \
    --invivo-result-root=PATH \
    --joint-result-root=PATH \
    --gemcitabine-data-root=PATH \
    --ltee-data-root=PATH \
    [--n-core=N] \
    [--recompute-fixed-o2] \
    [--recompute-invivo-tsne] \
    [--rebuild-figure6-grid]

Docker image:
  --docker-image=IMAGE
      Locally available Docker image used to run manager.sh. The image must
      resolve to:
        sha256:32c49db0ad27a0b5832b601ba96e2b72bfc1e2f1ccbf34687f8f596f1f7cdcd5

      Known equivalent local tags:
        o2_supply_demand_map:r44
        zafiro/o2_supply_demand_map:r44

      The wrapper does not automatically pull or rebuild an image.

Required scientific-input paths:
  --invitro-result-root=PATH
      Directory containing the separate in-vitro fit results.
  --invivo-result-root=PATH
      Directory containing the separate in-vivo fit results.
  --joint-result-root=PATH
      Directory containing the joint multi-warmup fit results.
  --gemcitabine-data-root=PATH
      Directory containing the in-vivo Gemcitabine input data.
  --ltee-data-root=PATH
      Directory containing the InVitroData_LTEE input tables.

Optional controls:
  --n-core=N
      Positive integer worker count. Default: 8.
  --recompute-fixed-o2
      Recompute the fixed-O2 analysis instead of reusing staged results.
  --recompute-invivo-tsne
      Recompute the in-vivo t-SNE landscape instead of reusing staged results.
  --rebuild-figure6-grid
      Rebuild all Figure 6 multi-seed response caches.
  -h, --help, --hlep
      Print this help text and exit. --hlep is accepted as a compatibility
      alias for the standard --help spelling.

Container behavior:
  The repository is mounted read-write at its existing absolute path. Input
  roots outside the repository are mounted read-only at their existing
  absolute paths. The container runs as the current host UID:GID, without
  network access, and is removed after the run. The wrapper requires a
  linux/amd64 image and checks for Rscript, python3, and ImageMagick's magick
  command, plus the R packages used directly by the figure workflow. It forces
  R's bitmapType to cairo and smoke-tests headless PNG output before starting
  manager.sh.

Outputs:
  manager.sh regenerates Figure 1-6 and parent-indexed supplementary figures, publishes figure
  artifacts, validates scientific inputs against the fixed Code/config MD5
  baseline, hashes intermediates, verifies published-figure MD5 identity, and
  generates the embedded manuscript HTML report. It does not compile the TeX
  manuscript or create/replace PDF.

Example:
  bash /absolute/path/to/workspace/Code/Docker/manager_docker.sh \
    --docker-image=o2_supply_demand_map:r44 \
    --invitro-result-root=/path/to/fit_invitro_O2_buffering_500seed \
    --invivo-result-root=/path/to/fit_invivo_O2_buffering_500seed \
    --joint-result-root=/path/to/fit_joint_multi_warmup_results \
    --gemcitabine-data-root=/path/to/data/InVivoData_Gemcitabine \
    --ltee-data-root=/path/to/data/InVitroData_LTEE \
    --n-core=8
EOF
}

normalize_input_directory() {
  local option_name="$1"
  local input_path="$2"
  if [[ ! -d "${input_path}" ]]; then
    echo "${option_name} is not an existing directory: ${input_path}" >&2
    return 1
  fi
  (
    cd "${input_path}"
    pwd -P
  )
}

DOCKER_IMAGE=""
INVITRO_RESULT_ROOT=""
INVIVO_RESULT_ROOT=""
JOINT_RESULT_ROOT=""
GEMCITABINE_DATA_ROOT=""
LTEE_DATA_ROOT=""
N_CORE=8
RECOMPUTE_FIXED_O2=FALSE
RECOMPUTE_INVIVO_TSNE=FALSE
REBUILD_FIGURE6_GRID=FALSE

for argument in "$@"; do
  case "${argument}" in
    -h|--help|--hlep)
      usage
      exit 0
      ;;
    --docker-image=*)
      DOCKER_IMAGE="${argument#*=}"
      ;;
    --invitro-result-root=*)
      INVITRO_RESULT_ROOT="${argument#*=}"
      ;;
    --invivo-result-root=*)
      INVIVO_RESULT_ROOT="${argument#*=}"
      ;;
    --joint-result-root=*)
      JOINT_RESULT_ROOT="${argument#*=}"
      ;;
    --gemcitabine-data-root=*)
      GEMCITABINE_DATA_ROOT="${argument#*=}"
      ;;
    --ltee-data-root=*)
      LTEE_DATA_ROOT="${argument#*=}"
      ;;
    --n-core=*)
      N_CORE="${argument#*=}"
      ;;
    --recompute-fixed-o2)
      RECOMPUTE_FIXED_O2=TRUE
      ;;
    --recompute-invivo-tsne)
      RECOMPUTE_INVIVO_TSNE=TRUE
      ;;
    --rebuild-figure6-grid)
      REBUILD_FIGURE6_GRID=TRUE
      ;;
    *)
      echo "Unknown option: ${argument}" >&2
      echo "Run 'manager_docker.sh --help' for usage." >&2
      exit 2
      ;;
  esac
done

missing_options=()
[[ -n "${DOCKER_IMAGE}" ]] ||
  missing_options+=("--docker-image")
[[ -n "${INVITRO_RESULT_ROOT}" ]] ||
  missing_options+=("--invitro-result-root")
[[ -n "${INVIVO_RESULT_ROOT}" ]] ||
  missing_options+=("--invivo-result-root")
[[ -n "${JOINT_RESULT_ROOT}" ]] ||
  missing_options+=("--joint-result-root")
[[ -n "${GEMCITABINE_DATA_ROOT}" ]] ||
  missing_options+=("--gemcitabine-data-root")
[[ -n "${LTEE_DATA_ROOT}" ]] ||
  missing_options+=("--ltee-data-root")

if (( ${#missing_options[@]} > 0 )); then
  echo "Missing required option(s): ${missing_options[*]}" >&2
  echo "Run 'manager_docker.sh --help' for usage." >&2
  exit 2
fi

if ! [[ "${N_CORE}" =~ ^[1-9][0-9]*$ ]]; then
  echo "--n-core must be a positive integer." >&2
  exit 2
fi

if [[ ! -f "${MANAGER}" ]]; then
  echo "Could not find manager.sh: ${MANAGER}" >&2
  exit 2
fi

if [[ ! -f "${FONT_FILE}" ]]; then
  echo "Could not find repository font: ${FONT_FILE}" >&2
  exit 2
fi

INVITRO_RESULT_ROOT="$(
  normalize_input_directory \
    "--invitro-result-root" "${INVITRO_RESULT_ROOT}"
)"
INVIVO_RESULT_ROOT="$(
  normalize_input_directory \
    "--invivo-result-root" "${INVIVO_RESULT_ROOT}"
)"
JOINT_RESULT_ROOT="$(
  normalize_input_directory \
    "--joint-result-root" "${JOINT_RESULT_ROOT}"
)"
GEMCITABINE_DATA_ROOT="$(
  normalize_input_directory \
    "--gemcitabine-data-root" "${GEMCITABINE_DATA_ROOT}"
)"
LTEE_DATA_ROOT="$(
  normalize_input_directory \
    "--ltee-data-root" "${LTEE_DATA_ROOT}"
)"

if ! command -v docker >/dev/null 2>&1; then
  echo "Docker is not available on PATH." >&2
  exit 2
fi

if ! docker info >/dev/null 2>&1; then
  echo "The Docker daemon is not available." >&2
  exit 2
fi

if ! RESOLVED_IMAGE_ID="$(
  docker image inspect --format '{{.Id}}' "${DOCKER_IMAGE}" 2>/dev/null
)"; then
  echo "Docker image is not available locally: ${DOCKER_IMAGE}" >&2
  echo "The wrapper does not pull images automatically." >&2
  exit 2
fi

if [[ "${RESOLVED_IMAGE_ID}" != "${EXPECTED_IMAGE_ID}" ]]; then
  echo "Docker image ID mismatch for ${DOCKER_IMAGE}." >&2
  echo "Expected: ${EXPECTED_IMAGE_ID}" >&2
  echo "Resolved: ${RESOLVED_IMAGE_ID}" >&2
  exit 2
fi

RESOLVED_PLATFORM="$(
  docker image inspect --format '{{.Os}}/{{.Architecture}}' "${DOCKER_IMAGE}"
)"
if [[ "${RESOLVED_PLATFORM}" != "linux/amd64" ]]; then
  echo "Docker image platform mismatch for ${DOCKER_IMAGE}." >&2
  echo "Expected: linux/amd64" >&2
  echo "Resolved: ${RESOLVED_PLATFORM}" >&2
  exit 2
fi

MANAGER_ARGS=(
  "--invitro-result-root=${INVITRO_RESULT_ROOT}"
  "--invivo-result-root=${INVIVO_RESULT_ROOT}"
  "--joint-result-root=${JOINT_RESULT_ROOT}"
  "--gemcitabine-data-root=${GEMCITABINE_DATA_ROOT}"
  "--ltee-data-root=${LTEE_DATA_ROOT}"
  "--n-core=${N_CORE}"
)
if [[ "${RECOMPUTE_FIXED_O2}" == "TRUE" ]]; then
  MANAGER_ARGS+=("--recompute-fixed-o2")
fi
if [[ "${RECOMPUTE_INVIVO_TSNE}" == "TRUE" ]]; then
  MANAGER_ARGS+=("--recompute-invivo-tsne")
fi
if [[ "${REBUILD_FIGURE6_GRID}" == "TRUE" ]]; then
  MANAGER_ARGS+=("--rebuild-figure6-grid")
fi

DOCKER_ARGS=(
  run
  --rm
  --init
  --platform linux/amd64
  --user "$(id -u):$(id -g)"
  --network none
  --workdir "${ITERATION_ROOT}"
  --env HOME=/tmp
  --env TMPDIR=/tmp
  --env XDG_CACHE_HOME=/tmp/.cache
  --mount "type=bind,source=${REPO_ROOT},target=${REPO_ROOT}"
)

MOUNTED_INPUT_ROOTS=()
append_read_only_input_mount() {
  local input_root="$1"
  local mounted_root

  if [[ "${input_root}" == "${REPO_ROOT}" ||
        "${input_root}" == "${REPO_ROOT}/"* ]]; then
    return 0
  fi

  for mounted_root in "${MOUNTED_INPUT_ROOTS[@]}"; do
    if [[ "${input_root}" == "${mounted_root}" ]]; then
      return 0
    fi
  done

  DOCKER_ARGS+=(
    --mount "type=bind,source=${input_root},target=${input_root},readonly"
  )
  MOUNTED_INPUT_ROOTS+=("${input_root}")
}

append_read_only_input_mount "${INVITRO_RESULT_ROOT}"
append_read_only_input_mount "${INVIVO_RESULT_ROOT}"
append_read_only_input_mount "${JOINT_RESULT_ROOT}"
append_read_only_input_mount "${GEMCITABINE_DATA_ROOT}"
append_read_only_input_mount "${LTEE_DATA_ROOT}"

CONTAINER_COMMAND='
set -euo pipefail
r_profile_user="${TMPDIR:-/tmp}/figure-workspace-Rprofile"
printf "%s\n" "options(bitmapType = \"cairo\")" > "${r_profile_user}"
export R_PROFILE_USER="${r_profile_user}"
preflight_failed=0
for required_command in Rscript python3 magick; do
  if ! command -v "${required_command}" >/dev/null 2>&1; then
    echo "Container is missing required command: ${required_command}" >&2
    preflight_failed=1
  fi
done
if command -v Rscript >/dev/null 2>&1; then
  missing_r_packages="$(
    Rscript --vanilla -e "packages <- c(\"Matrix\", \"Rcpp\", \"Rtsne\", \"cluster\", \"cowplot\", \"data.table\", \"dplyr\", \"ggnewscale\", \"ggplot2\", \"ggrepel\", \"magick\", \"patchwork\", \"readxl\", \"scales\", \"shadowtext\", \"svglite\", \"tidyr\"); missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]; cat(missing, sep = \"\\n\")"
  )"
  if [[ -n "${missing_r_packages}" ]]; then
    while IFS= read -r missing_r_package; do
      [[ -n "${missing_r_package}" ]] &&
        echo "Container is missing required R package: ${missing_r_package}" >&2
    done <<<"${missing_r_packages}"
    preflight_failed=1
  fi
  if ! Rscript -e "stopifnot(isTRUE(capabilities(\"cairo\")), identical(getOption(\"bitmapType\"), \"cairo\")); output <- tempfile(fileext = \".png\"); grDevices::png(output, width = 64, height = 64, type = getOption(\"bitmapType\")); graphics::par(mar = rep(0, 4)); graphics::plot.new(); grDevices::dev.off(); stopifnot(file.exists(output), as.numeric(file.info(output)[[\"size\"]]) > 0); unlink(output)" >/dev/null; then
    echo "Container failed the headless Cairo PNG smoke test." >&2
    preflight_failed=1
  fi
fi
if (( preflight_failed != 0 )); then
  exit 2
fi
manager_path="$1"
shift
exec bash "${manager_path}" "$@"
'

echo "Docker image: ${DOCKER_IMAGE}"
echo "Docker image ID: ${RESOLVED_IMAGE_ID}"
echo "Docker platform: ${RESOLVED_PLATFORM}"
echo "Repository mount (read-write): ${REPO_ROOT}"

exec docker "${DOCKER_ARGS[@]}" \
  "${DOCKER_IMAGE}" \
  bash -c "${CONTAINER_COMMAND}" bash \
  "${MANAGER}" \
  "${MANAGER_ARGS[@]}"
