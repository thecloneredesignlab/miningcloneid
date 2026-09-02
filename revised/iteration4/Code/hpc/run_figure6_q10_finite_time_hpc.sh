#!/usr/bin/env bash

# Fresh Figure 6 q10 finite-time computation, diagnostics, and headless drawing
# on the allocated RED compute node hpctpa3pc0028.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
ITERATION_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"
REPO_ROOT="$(cd -- "${ITERATION_ROOT}/../.." && pwd -P)"
CODE_ROOT="${ITERATION_ROOT}/Code/Figures"

EXPECTED_REPO_ROOT="/share/lab_crd/taoli/Project/HypoxiaLTEEFigures"
EXPECTED_NODE="hpctpa3pc0028"
SIF_IMAGE="/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif"
MODEL_CODE_ROOT="/share/lab_crd/taoli/Project/soft_couping_org/oxygen/code/O2_supply_demand_MAP"
RESULTS_ROOT="/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results"
INVIVO_RESULT_ROOT="${RESULTS_ROOT}/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253"
INVITRO_RESULT_ROOT="${RESULTS_ROOT}/fit_invitro_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253"
JOINT_RESULT_ROOT="${RESULTS_ROOT}/fit_joint_unified_global_invitro_500seed_all_xxlarge_r442_exact_20260828_145253"

N_CORE=63
RUN_ID="$(date '+%Y%m%d_%H%M%S')"
PREFLIGHT_ONLY=FALSE
SMOKE_ONLY=FALSE
RED_EASYBUILD_ROOT="/app/eb"
RED_FLEXIBLAS_LIB="${RED_EASYBUILD_ROOT}/software/FlexiBLAS/3.4.4-GCC-13.3.0/lib64"
RED_OPENBLAS_LIB="${RED_EASYBUILD_ROOT}/software/OpenBLAS/0.3.27-GCC-13.3.0/lib"
RED_GCC_LIB="${RED_EASYBUILD_ROOT}/software/GCCcore/13.3.0/lib64"
RED_BINUTILS_LIB="${RED_EASYBUILD_ROOT}/software/binutils/2.42-GCCcore-13.3.0/lib"
CONTAINER_LD_LIBRARY_PATH="${RED_FLEXIBLAS_LIB}:${RED_OPENBLAS_LIB}:${RED_GCC_LIB}:${RED_BINUTILS_LIB}:/opt/rh/gcc-toolset-13/root/usr/lib64"

usage() {
  printf '%s\n' \
    'Usage: run_figure6_q10_finite_time_hpc.sh [options]' \
    '' \
    'Options:' \
    '  --n-core=N       Worker count; coordinator plus workers must fit 64 CPUs.' \
    '                   Default: 63.' \
    '  --run-id=ID      Fresh output namespace. Default: current timestamp.' \
    '  --preflight-only Validate node, paths, container, packages, and inputs.' \
    '  --smoke-only     Run the reduced finite-time calculation without drawing.' \
    '  -h, --help       Show this help.'
}

for argument in "$@"; do
  case "${argument}" in
    --n-core=*) N_CORE="${argument#*=}" ;;
    --run-id=*) RUN_ID="${argument#*=}" ;;
    --preflight-only) PREFLIGHT_ONLY=TRUE ;;
    --smoke-only) SMOKE_ONLY=TRUE ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: ${argument}" >&2; usage >&2; exit 2 ;;
  esac
done

if ! [[ "${N_CORE}" =~ ^[1-9][0-9]*$ ]] || (( N_CORE > 63 )); then
  echo "--n-core must be an integer from 1 through 63." >&2
  exit 2
fi
if ! [[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$ ]]; then
  echo "--run-id contains unsupported characters or is too long." >&2
  exit 2
fi
if [[ "${PREFLIGHT_ONLY}" == TRUE && "${SMOKE_ONLY}" == TRUE ]]; then
  echo "--preflight-only and --smoke-only are mutually exclusive." >&2
  exit 2
fi

HOST_SHORT="$(hostname -s)"
if [[ "${HOST_SHORT}" != "${EXPECTED_NODE}" ]]; then
  echo "This runner must execute on ${EXPECTED_NODE}; observed ${HOST_SHORT}." >&2
  exit 2
fi
if [[ "${REPO_ROOT}" != "${EXPECTED_REPO_ROOT}" ]]; then
  echo "Unexpected HPC repository root: ${REPO_ROOT}" >&2
  exit 2
fi

require_file() {
  local path="$1"
  local role="$2"
  if [[ ! -f "${path}" || ! -r "${path}" ]]; then
    echo "Missing readable ${role}: ${path}" >&2
    exit 2
  fi
}

require_directory() {
  local path="$1"
  local role="$2"
  if [[ ! -d "${path}" || ! -r "${path}" ]]; then
    echo "Missing readable ${role}: ${path}" >&2
    exit 2
  fi
}

require_file "${SIF_IMAGE}" "SIF image"
require_directory "${MODEL_CODE_ROOT}" "external model-code root"
require_directory "${INVIVO_RESULT_ROOT}" "in-vivo fit-result root"
require_directory "${INVITRO_RESULT_ROOT}" "in-vitro fit-result root"
require_directory "${JOINT_RESULT_ROOT}" "joint fit-result root"
require_directory "${RED_EASYBUILD_ROOT}" "RED EasyBuild root"
require_file "${RED_FLEXIBLAS_LIB}/libflexiblas.so.3" "RED FlexiBLAS runtime"
require_directory "${RED_OPENBLAS_LIB}" "RED OpenBLAS runtime"
require_directory "${RED_GCC_LIB}" "RED GCC runtime"
require_directory "${RED_BINUTILS_LIB}" "RED binutils runtime"

MODEL_SOURCE_FILES=(
  "${MODEL_CODE_ROOT}/model/model_O2_supply_demand_MAP.R"
  "${MODEL_CODE_ROOT}/model/model_O2_supply_demand_MAP.cpp"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_shared.R"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_common_semantics.R"
  "${MODEL_CODE_ROOT}/simulation/o2/fixed_o2/run_fixed_o2_simulation.R"
  "${MODEL_CODE_ROOT}/simulation/o2/fixed_o2/fixed_o2_numerical_producers.R"
)
for path in "${MODEL_SOURCE_FILES[@]}"; do
  require_file "${path}" "external model source"
done

FIGURE6_SOURCE_FILES=(
  "${ITERATION_ROOT}/data/Figures/Figure6/joint_multiseed_surface_summary.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure6/joint_multiseed_surface_summary_invitro.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure6/joint_multiseed_trajectory_summary.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure6/joint_multiseed_trajectory_summary_invitro.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure6/figure6_inverse_response_summary.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure6/figure6_invitro_inverse_response_summary.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure6/figure6d_fixed_p_curve_family.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure6/figure6_invitro_fixed_p_curve_family.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure5/figure5_frozen_inputs/selected_results.tsv"
  "${ITERATION_ROOT}/data/Figures/Supp_Figure5_1/soft_coupling_master_long.tsv"
)
for path in "${FIGURE6_SOURCE_FILES[@]}"; do
  require_file "${path}" "Figure 6 source input"
done

while IFS=$'\t' read -r record_type warmup_label rest; do
  if [[ "${record_type}" != "joint_pair_best" ]]; then
    continue
  fi
  require_file \
    "${ITERATION_ROOT}/data/Figures/Figure5/figure5_frozen_inputs/winners/${warmup_label}/fit_config.rds" \
    "joint-pair fit config"
  require_file \
    "${ITERATION_ROOT}/data/Figures/Figure5/figure5_frozen_inputs/winners/${warmup_label}/best_params.tsv" \
    "joint-pair best parameters"
  require_file \
    "${ITERATION_ROOT}/data/Figures/Supp_Figure5_2/joint/${warmup_label}/seed_objective_simple.tsv" \
    "joint-pair 500-seed objective table"
done < <(tail -n +2 "${ITERATION_ROOT}/data/Figures/Figure5/figure5_frozen_inputs/selected_results.tsv")

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v apptainer)"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v singularity)"
else
  echo "Neither apptainer nor singularity is available on ${EXPECTED_NODE}." >&2
  exit 2
fi

AVAILABLE_CPU="$(getconf _NPROCESSORS_ONLN)"
if (( AVAILABLE_CPU < N_CORE + 1 )); then
  echo "Need at least $((N_CORE + 1)) online CPUs; observed ${AVAILABLE_CPU}." >&2
  exit 2
fi
MEMORY_KB="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
if [[ -z "${MEMORY_KB}" ]] || (( MEMORY_KB < 500000000 )); then
  echo "Expected approximately 512 GB RAM; observed ${MEMORY_KB:-unknown} kB." >&2
  exit 2
fi

AUDIT_ROOT="${ITERATION_ROOT}/audit/hpc_figure6_finite_time/${RUN_ID}"
LOG_ROOT="${ITERATION_ROOT}/audit/logs"
LOCK_ROOT="${ITERATION_ROOT}/audit/locks"
LOCK_DIR="${LOCK_ROOT}/figure6_q10_finite_time.lock"
RUN_LOG="${LOG_ROOT}/figure6_q10_finite_time_${RUN_ID}.log"
LATEST_LOG="${LOG_ROOT}/figure6_q10_finite_time_latest.log"
STATUS_PATH="${AUDIT_ROOT}/status.tsv"
INPUT_MD5_PATH="${AUDIT_ROOT}/input_md5.tsv"
OUTPUT_SHA256_PATH="${AUDIT_ROOT}/output_sha256.tsv"
TASK_TMP_DIR=""
RUN_STATUS="INITIALIZING"
LOCK_ACQUIRED=FALSE

mkdir -p "${AUDIT_ROOT}" "${LOG_ROOT}" "${LOCK_ROOT}"
if ! mkdir "${LOCK_DIR}" 2>/dev/null; then
  echo "Another Figure 6 finite-time run owns the lock: ${LOCK_DIR}" >&2
  exit 2
fi
LOCK_ACQUIRED=TRUE
printf 'host=%s\npid=%s\nrun_id=%s\n' \
  "${HOST_SHORT}" "$$" "${RUN_ID}" > "${LOCK_DIR}/owner"

write_status() {
  local status="$1"
  local exit_code="$2"
  local stage="$3"
  {
    printf 'run_id\tstatus\texit_code\tstage\thost\tn_core\tgit_head\tupdated_at\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${RUN_ID}" "${status}" "${exit_code}" "${stage}" "${HOST_SHORT}" \
      "${N_CORE}" "$(git -C "${REPO_ROOT}" rev-parse HEAD)" "$(date -Iseconds)"
  } > "${STATUS_PATH}"
}

cleanup() {
  local exit_code="$?"
  set +e
  if [[ "${RUN_STATUS}" != "COMPLETE" ]]; then
    write_status "FAILED" "${exit_code}" "${RUN_STATUS}"
  fi
  [[ -f "${RUN_LOG}" ]] && cp "${RUN_LOG}" "${LATEST_LOG}"
  if [[ -n "${TASK_TMP_DIR}" && -d "${TASK_TMP_DIR}" ]]; then
    rm -rf -- "${TASK_TMP_DIR}"
  fi
  if [[ "${LOCK_ACQUIRED}" == TRUE ]]; then
    rm -f -- "${LOCK_DIR}/owner"
    rmdir "${LOCK_DIR}" 2>/dev/null || true
  fi
  trap - EXIT
  exit "${exit_code}"
}
trap cleanup EXIT

exec > >(tee -a "${RUN_LOG}") 2>&1

echo "Figure 6 q10 finite-time run start: $(date -Iseconds)"
echo "run_id=${RUN_ID}"
echo "host=${HOST_SHORT}"
echo "n_core=${N_CORE}"
echo "memory_kb=${MEMORY_KB}"
echo "git_head=$(git -C "${REPO_ROOT}" rev-parse HEAD)"
echo "model_code_root=${MODEL_CODE_ROOT}"
echo "joint_result_root=${JOINT_RESULT_ROOT}"

RUN_STATUS="PREFLIGHT"
write_status "RUNNING" "0" "${RUN_STATUS}"

{
  printf 'md5\tsize_bytes\tpath\n'
  for path in "${MODEL_SOURCE_FILES[@]}" "${FIGURE6_SOURCE_FILES[@]}"; do
    printf '%s\t%s\t%s\n' \
      "$(md5sum "${path}" | awk '{print $1}')" \
      "$(stat -c '%s' "${path}")" "${path}"
  done
} > "${INPUT_MD5_PATH}"

SIF_SHA256="$(sha256sum "${SIF_IMAGE}" | awk '{print $1}')"
echo "sif_sha256=${SIF_SHA256}"

TASK_TMP_DIR="$(mktemp -d /tmp/figure6-q10-finite-time.XXXXXX)"
mkdir -p \
  "${TASK_TMP_DIR}/home" \
  "${TASK_TMP_DIR}/cache" \
  "${TASK_TMP_DIR}/model_rcpp_cache"
printf '%s\n' \
  'options(bitmapType = "cairo", device = "png", warn = 1)' \
  > "${TASK_TMP_DIR}/Rprofile"

CONTAINER_ARGS=(
  exec
  --cleanenv
  --containall
  --pwd "${ITERATION_ROOT}"
  --home "${TASK_TMP_DIR}/home"
  --env "TMPDIR=${TASK_TMP_DIR}"
  --env "XDG_CACHE_HOME=${TASK_TMP_DIR}/cache"
  --env "R_PROFILE_USER=${TASK_TMP_DIR}/Rprofile"
  --env "R_HOME=/opt/R/4.4.2/lib64/R"
  --env "LD_LIBRARY_PATH=${CONTAINER_LD_LIBRARY_PATH}"
  --env "DISPLAY="
  --env "QT_QPA_PLATFORM=offscreen"
  --env "FIGURE_WORKSPACE_ROOT=${ITERATION_ROOT}"
  --env "FIGURE_MODEL_CODE_ROOT=${MODEL_CODE_ROOT}"
  --env "FIGURE_INVIVO_RESULT_ROOT=${INVIVO_RESULT_ROOT}"
  --env "FIGURE_INVITRO_RESULT_ROOT=${INVITRO_RESULT_ROOT}"
  --env "FIGURE_JOINT_RESULT_ROOT=${JOINT_RESULT_ROOT}"
  --env "FIGURE6_FINITE_TIME_RUN_ID=${RUN_ID}"
  --env "FIGURE6_FINITE_TIME_FUTURE_PLAN=multicore"
  --env "OMP_NUM_THREADS=1"
  --env "OPENBLAS_NUM_THREADS=1"
  --env "MKL_NUM_THREADS=1"
  --env "VECLIB_MAXIMUM_THREADS=1"
  --env "RCPP_PARALLEL_NUM_THREADS=1"
  --env "KMP_USE_SHM=0"
  --env "KMP_INIT_AT_FORK=FALSE"
  --bind "${RED_EASYBUILD_ROOT}:${RED_EASYBUILD_ROOT}:ro"
  --bind "${ITERATION_ROOT}:${ITERATION_ROOT}:rw"
  --bind "${MODEL_CODE_ROOT}:${MODEL_CODE_ROOT}:ro"
  --bind "${TASK_TMP_DIR}/model_rcpp_cache:${MODEL_CODE_ROOT}/model/.rcpp_cache_o2_supply_demand_map:rw"
  --bind "${RESULTS_ROOT}:${RESULTS_ROOT}:ro"
  --bind "${TASK_TMP_DIR}:${TASK_TMP_DIR}:rw"
)

container_command() {
  "${CONTAINER_RUNTIME}" "${CONTAINER_ARGS[@]}" "${SIF_IMAGE}" "$@"
}

container_command Rscript --vanilla -e '
packages <- c(
  "Matrix", "Rcpp", "data.table", "future", "future.apply", "ggplot2",
  "patchwork", "magick", "isoband", "scales"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
workspace <- normalizePath(Sys.getenv("FIGURE_WORKSPACE_ROOT"), mustWork = TRUE)
model_root <- normalizePath(Sys.getenv("FIGURE_MODEL_CODE_ROOT"), mustWork = TRUE)
source(file.path(workspace, "Code", "Figures", "util", "analysis", "figure6_robustness.R"))
source(file.path(workspace, "Code", "Figures", "util", "analysis", "figure6_context_extension.R"))
source(file.path(workspace, "Code", "Figures", "util", "analysis", "figure6_finite_time_q10.R"))
source(file.path(workspace, "Code", "Figures", "util", "analysis", "figure6_finite_time_plots.R"))
paths <- f6r_paths(workspace)
stopifnot(identical(paths$oxygen_code, model_root))
f6r_load_response_engine(paths)
stopifnot(
  future::supportsMulticore(),
  exists("fixo2_fixed_matrix", envir = globalenv()),
  exists("response_force_effective_p_misseg", envir = globalenv()),
  exists("f6ft_data", envir = globalenv()),
  exists("f6ft_draw_main", envir = globalenv())
)
headless_png <- file.path(tempdir(), "figure6_headless_cairo_smoke.png")
grDevices::png(
  headless_png, width = 480, height = 480, type = "cairo", bg = "white"
)
graphics::par(mar = rep(0, 4))
graphics::plot.new()
grDevices::dev.off()
stopifnot(file.exists(headless_png), file.info(headless_png)$size > 0)
unlink(headless_png)
cat("figure6_finite_time_container_preflight_ok\n")
'

if [[ "${PREFLIGHT_ONLY}" == TRUE ]]; then
  RUN_STATUS="COMPLETE"
  write_status "COMPLETE" "0" "PREFLIGHT_ONLY"
  echo "Figure 6 q10 finite-time preflight complete: $(date -Iseconds)"
  exit 0
fi

if [[ "${SMOKE_ONLY}" == TRUE ]]; then
  RUN_STATUS="SMOKE"
  write_status "RUNNING" "0" "${RUN_STATUS}"
  container_command Rscript --vanilla \
    "${CODE_ROOT}/data_Figure6_finite_time_q10.R" \
    "--n-core=${N_CORE}" "--run-id=${RUN_ID}_smoke" \
    --smoke=TRUE --publish-current=FALSE --compute-diagnostics=FALSE
  RUN_STATUS="COMPLETE"
  write_status "COMPLETE" "0" "SMOKE_ONLY"
  echo "Figure 6 q10 finite-time smoke test complete: $(date -Iseconds)"
  exit 0
fi

RUN_STATUS="FINITE_TIME_DATA"
write_status "RUNNING" "0" "${RUN_STATUS}"
container_command Rscript --vanilla \
  "${CODE_ROOT}/data_Figure6_finite_time_q10.R" \
  "--n-core=${N_CORE}" "--run-id=${RUN_ID}" \
  --smoke=FALSE --publish-current=TRUE --compute-diagnostics=TRUE

RUN_STATUS="DRAW_FIGURE6"
write_status "RUNNING" "0" "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/draw_Figure6.R"

RUN_STATUS="DRAW_SUPP_FIGURE6_5"
write_status "RUNNING" "0" "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/draw_Supp_Figure6_5.R"

RUN_STATUS="DRAW_SUPP_FIGURE6_6"
write_status "RUNNING" "0" "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/draw_Supp_Figure6_6.R"

RUN_STATUS="DRAW_SUPP_FIGURE6_7"
write_status "RUNNING" "0" "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/draw_Supp_Figure6_7.R"

RUN_STATUS="VALIDATE"
write_status "RUNNING" "0" "${RUN_STATUS}"

check_validation() {
  local path="$1"
  require_file "${path}" "validation table"
  awk -F '\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) if ($i == "passed") passed_col = i
      next
    }
    passed_col > 0 && toupper($passed_col) != "TRUE" {failed = 1}
    END {if (NR < 2 || passed_col == 0 || failed) exit 1}
  ' "${path}" || {
    echo "Validation table failed: ${path}" >&2
    exit 1
  }
}

RUN_ROOT="${ITERATION_ROOT}/data/Figures/Figure6/finite_time_q10_runs/${RUN_ID}"
check_validation "${RUN_ROOT}/q10_endpoint_manifest_validation.tsv"
check_validation "${RUN_ROOT}/finite_time_task_validation.tsv"
check_validation "${RUN_ROOT}/finite_time_panel_validation.tsv"
check_validation "${RUN_ROOT}/finite_time_diagnostic_endpoint_validation.tsv"
check_validation "${RUN_ROOT}/finite_time_euler_step_sensitivity.tsv"
check_validation "${RUN_ROOT}/figure6_finite_time_render_validation.tsv"
check_validation "${ITERATION_ROOT}/data/Figures/Supp_Figure6_5/supp_fig6-5_render_validation.tsv"
check_validation "${ITERATION_ROOT}/data/Figures/Supp_Figure6_6/supp_fig6-6_render_validation.tsv"
check_validation "${ITERATION_ROOT}/data/Figures/Supp_Figure6_7/supp_fig6-7_render_validation.tsv"

OUTPUT_FILES=(
  "${ITERATION_ROOT}/Figures/assembled_fig6.png"
  "${ITERATION_ROOT}/Figures/assembled_fig6.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig6-5_finite_time_eigen_expm_vs_euler.png"
  "${ITERATION_ROOT}/Figures/supp_fig6-5_finite_time_eigen_expm_vs_euler.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig6-6_eigen_expm_agreement.png"
  "${ITERATION_ROOT}/Figures/supp_fig6-6_eigen_expm_agreement.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig6-7_finite_time_vs_steady_attractor.png"
  "${ITERATION_ROOT}/Figures/supp_fig6-7_finite_time_vs_steady_attractor.pdf"
)
for path in "${OUTPUT_FILES[@]}"; do
  require_file "${path}" "rendered Figure 6 output"
done

{
  printf 'sha256\tsize_bytes\tpath\n'
  for path in "${OUTPUT_FILES[@]}"; do
    printf '%s\t%s\t%s\n' \
      "$(sha256sum "${path}" | awk '{print $1}')" \
      "$(stat -c '%s' "${path}")" "${path}"
  done
} > "${OUTPUT_SHA256_PATH}"

RUN_STATUS="COMPLETE"
write_status "COMPLETE" "0" "${RUN_STATUS}"
echo "Figure 6 q10 finite-time run complete: $(date -Iseconds)"
echo "status_path=${STATUS_PATH}"
echo "output_sha256=${OUTPUT_SHA256_PATH}"
