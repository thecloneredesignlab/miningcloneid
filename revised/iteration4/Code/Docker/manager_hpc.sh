#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ITERATION_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
REPO_ROOT="$(cd "${ITERATION_ROOT}/../.." && pwd -P)"
MANAGER="${ITERATION_ROOT}/manager.sh"
FONT_FILE="${REPO_ROOT}/data/fonts/Arial Bold.ttf"
DEFAULT_SIF_IMAGE="/share/lab_crd/lab_crd/taoli/Docker/o2_supply_demand_map_r44.sif"

usage() {
  cat <<'EOF'
Usage:
  manager_hpc.sh \
    [--sif-image=PATH] \
    --invitro-result-root=PATH \
    --invivo-result-root=PATH \
    --joint-result-root=PATH \
    --gemcitabine-data-root=PATH \
    --ltee-data-root=PATH \
    [--n-core=N] \
    [--recompute-fixed-o2] \
    [--recompute-invivo-tsne] \
    [--rebuild-figure6-grid]

HPC container image:
  --sif-image=PATH
      Apptainer/Singularity SIF file used to run manager.sh.

      Default:
        /share/lab_crd/lab_crd/taoli/Docker/o2_supply_demand_map_r44.sif

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
  The wrapper uses apptainer when available and otherwise tries singularity.
  It starts a clean, contained environment. The repository is bound
  read-write at its existing absolute path, while input roots outside the
  repository are bound read-only at their existing absolute paths. It checks
  for Rscript, python3, ImageMagick's magick command, and the R packages used
  directly by the figure workflow. It forces R's bitmapType to cairo and
  smoke-tests headless PNG output before starting manager.sh.

HPC scheduling:
  This script directly runs the container. It does not call sbatch and does
  not choose CPU, memory, partition, or time limits. Invoke it from an
  interactive allocation or from a separately managed Slurm submission
  script when compute-node execution is required.

Outputs:
  manager.sh regenerates Figure 1-6 and parent-indexed supplementary figures, publishes figure
  artifacts, validates scientific inputs against the fixed Code/config MD5
  baseline, hashes intermediates, verifies published-figure MD5 identity, and
  generates the embedded manuscript HTML report. It does not compile the TeX
  manuscript or create/replace PDF.

Example:
  bash /absolute/path/to/workspace/Code/Docker/manager_hpc.sh \
    --sif-image=/share/lab_crd/lab_crd/taoli/Docker/o2_supply_demand_map_r44.sif \
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

SIF_IMAGE="${DEFAULT_SIF_IMAGE}"
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
    --sif-image=*)
      SIF_IMAGE="${argument#*=}"
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
      echo "Run 'manager_hpc.sh --help' for usage." >&2
      exit 2
      ;;
  esac
done

missing_options=()
[[ -n "${SIF_IMAGE}" ]] ||
  missing_options+=("--sif-image")
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
  echo "Run 'manager_hpc.sh --help' for usage." >&2
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

if [[ ! -r "${SIF_IMAGE}" || ! -f "${SIF_IMAGE}" ]]; then
  echo "SIF image is not a readable file: ${SIF_IMAGE}" >&2
  exit 2
fi
SIF_IMAGE_DIRECTORY="$(cd "$(dirname "${SIF_IMAGE}")" && pwd -P)"
SIF_IMAGE="${SIF_IMAGE_DIRECTORY}/$(basename "${SIF_IMAGE}")"

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

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v apptainer)"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v singularity)"
else
  echo "Neither apptainer nor singularity is available on PATH." >&2
  exit 2
fi

if command -v sha256sum >/dev/null 2>&1; then
  SIF_SHA256="$(sha256sum "${SIF_IMAGE}")"
  SIF_SHA256="${SIF_SHA256%% *}"
elif command -v shasum >/dev/null 2>&1; then
  SIF_SHA256="$(shasum -a 256 "${SIF_IMAGE}")"
  SIF_SHA256="${SIF_SHA256%% *}"
else
  echo "Neither sha256sum nor shasum is available to verify the SIF file." >&2
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

CONTAINER_ARGS=(
  exec
  --cleanenv
  --containall
  --pwd "${ITERATION_ROOT}"
  --env HOME=/tmp
  --env TMPDIR=/tmp
  --env XDG_CACHE_HOME=/tmp/.cache
  --bind "${REPO_ROOT}:${REPO_ROOT}:rw"
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

  CONTAINER_ARGS+=(--bind "${input_root}:${input_root}:ro")
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

echo "Container runtime: ${CONTAINER_RUNTIME}"
echo "SIF image: ${SIF_IMAGE}"
echo "SIF SHA-256: ${SIF_SHA256}"
echo "Repository bind (read-write): ${REPO_ROOT}"

exec "${CONTAINER_RUNTIME}" "${CONTAINER_ARGS[@]}" \
  "${SIF_IMAGE}" \
  bash -c "${CONTAINER_COMMAND}" bash \
  "${MANAGER}" \
  "${MANAGER_ARGS[@]}"
