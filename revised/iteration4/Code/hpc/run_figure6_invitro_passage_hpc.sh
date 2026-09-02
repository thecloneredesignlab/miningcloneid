#!/usr/bin/env bash

# Compute passage-constrained Figure 6E-F and render Figure 6 plus the archived
# no-passaging Supplementary Figure 6-8 on hpctpa3pc0028.

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
JOINT_RESULT_ROOT="${RESULTS_ROOT}/fit_joint_unified_global_invitro_500seed_all_xxlarge_r442_exact_20260828_145253"

N_CORE=63
RUN_ID="figure6_passage_$(date '+%Y%m%d_%H%M%S')"
PREFLIGHT_ONLY=FALSE
SMOKE_ONLY=FALSE
RED_EASYBUILD_ROOT="/app/eb"
RED_FLEXIBLAS_LIB="${RED_EASYBUILD_ROOT}/software/FlexiBLAS/3.4.4-GCC-13.3.0/lib64"
RED_OPENBLAS_LIB="${RED_EASYBUILD_ROOT}/software/OpenBLAS/0.3.27-GCC-13.3.0/lib"
RED_GCC_LIB="${RED_EASYBUILD_ROOT}/software/GCCcore/13.3.0/lib64"
RED_BINUTILS_LIB="${RED_EASYBUILD_ROOT}/software/binutils/2.42-GCCcore-13.3.0/lib"
CONTAINER_LD_LIBRARY_PATH="${RED_FLEXIBLAS_LIB}:${RED_OPENBLAS_LIB}:${RED_GCC_LIB}:${RED_BINUTILS_LIB}:/opt/rh/gcc-toolset-13/root/usr/lib64"

for argument in "$@"; do
  case "${argument}" in
    --n-core=*) N_CORE="${argument#*=}" ;;
    --run-id=*) RUN_ID="${argument#*=}" ;;
    --preflight-only) PREFLIGHT_ONLY=TRUE ;;
    --smoke-only) SMOKE_ONLY=TRUE ;;
    -h|--help)
      echo "Usage: $0 [--n-core=1..63] [--run-id=ID] [--preflight-only|--smoke-only]"
      exit 0 ;;
    *) echo "Unknown option: ${argument}" >&2; exit 2 ;;
  esac
done

if ! [[ "${N_CORE}" =~ ^[1-9][0-9]*$ ]] || (( N_CORE > 63 )); then
  echo "--n-core must be an integer from 1 through 63." >&2
  exit 2
fi
if ! [[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$ ]]; then
  echo "Invalid --run-id." >&2
  exit 2
fi
if [[ "${PREFLIGHT_ONLY}" == TRUE && "${SMOKE_ONLY}" == TRUE ]]; then
  echo "--preflight-only and --smoke-only are mutually exclusive." >&2
  exit 2
fi
if [[ "$(hostname -s)" != "${EXPECTED_NODE}" ]]; then
  echo "This runner must execute on ${EXPECTED_NODE}." >&2
  exit 2
fi
if [[ "${REPO_ROOT}" != "${EXPECTED_REPO_ROOT}" ]]; then
  echo "Unexpected repository root: ${REPO_ROOT}" >&2
  exit 2
fi

require_file() {
  [[ -f "$1" && -r "$1" ]] || { echo "Missing readable $2: $1" >&2; exit 2; }
}
require_dir() {
  [[ -d "$1" && -r "$1" ]] || { echo "Missing readable $2: $1" >&2; exit 2; }
}

require_file "${SIF_IMAGE}" "SIF image"
require_dir "${MODEL_CODE_ROOT}" "external model-code root"
require_dir "${JOINT_RESULT_ROOT}" "joint-result root"
require_dir "${RED_EASYBUILD_ROOT}" "RED EasyBuild root"
require_file "${RED_FLEXIBLAS_LIB}/libflexiblas.so.3" "FlexiBLAS runtime"
require_dir "${RED_OPENBLAS_LIB}" "OpenBLAS runtime"
require_dir "${RED_GCC_LIB}" "GCC runtime"
require_dir "${RED_BINUTILS_LIB}" "binutils runtime"

MODEL_FILES=(
  "${MODEL_CODE_ROOT}/model/model_O2_supply_demand_MAP.R"
  "${MODEL_CODE_ROOT}/model/model_O2_supply_demand_MAP.cpp"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_shared.R"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_common_semantics.R"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_invitro_lineage_utils.R"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_invitro_lineage_simulation_utils.R"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_invitro_objective_utils.R"
  "${MODEL_CODE_ROOT}/simulation/fix_o2_simulation.R"
  "${MODEL_CODE_ROOT}/simulation/fix_o2_simulation_shared_utils.R"
  "${MODEL_CODE_ROOT}/simulation/o2/fixed_o2/run_fixed_o2_simulation.R"
)
for path in "${MODEL_FILES[@]}"; do require_file "${path}" "external model source"; done

BASE_POINTER="${ITERATION_ROOT}/data/Figures/Figure6/finite_time_q10_current.tsv"
require_file "${BASE_POINTER}" "continuous Figure 6 run pointer"
BASE_RELATIVE_RUN="$(awk -F '\t' 'NR==2 {print $2}' "${BASE_POINTER}")"
BASE_RUN="${ITERATION_ROOT}/data/Figures/Figure6/${BASE_RELATIVE_RUN}"
require_file "${BASE_RUN}/q10_unique_endpoint_manifest.tsv" "q10 endpoint manifest"
require_file "${BASE_RUN}/q10_optimizer_seed_manifest.tsv" "q10 seed manifest"
for letter in c d e f; do
  require_file "${BASE_RUN}/finite_time_panel_${letter}.rds" "continuous finite-time panel ${letter}"
done

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v apptainer)"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v singularity)"
else
  echo "No apptainer/singularity runtime found." >&2
  exit 2
fi
AVAILABLE_CPU="$(getconf _NPROCESSORS_ONLN)"
MEMORY_KB="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
if (( AVAILABLE_CPU < N_CORE + 1 )); then
  echo "Need $((N_CORE + 1)) CPUs; observed ${AVAILABLE_CPU}." >&2
  exit 2
fi
if [[ -z "${MEMORY_KB}" ]] || (( MEMORY_KB < 500000000 )); then
  echo "Expected approximately 512 GB RAM; observed ${MEMORY_KB:-unknown} kB." >&2
  exit 2
fi

AUDIT_ROOT="${ITERATION_ROOT}/audit/hpc_figure6_invitro_passage/${RUN_ID}"
LOG_ROOT="${ITERATION_ROOT}/audit/logs"
LOCK_ROOT="${ITERATION_ROOT}/audit/locks"
LOCK_DIR="${LOCK_ROOT}/figure6_invitro_passage.lock"
RUN_LOG="${LOG_ROOT}/figure6_invitro_passage_${RUN_ID}.log"
LATEST_LOG="${LOG_ROOT}/figure6_invitro_passage_latest.log"
STATUS_PATH="${AUDIT_ROOT}/status.tsv"
INPUT_MD5_PATH="${AUDIT_ROOT}/input_md5.tsv"
OUTPUT_SHA256_PATH="${AUDIT_ROOT}/output_sha256.tsv"
TASK_TMP_DIR=""
RUN_STATUS="INITIALIZING"
LOCK_ACQUIRED=FALSE

mkdir -p "${AUDIT_ROOT}" "${LOG_ROOT}" "${LOCK_ROOT}"
if ! mkdir "${LOCK_DIR}" 2>/dev/null; then
  echo "Another Figure 6 passage run owns ${LOCK_DIR}." >&2
  exit 2
fi
LOCK_ACQUIRED=TRUE
printf 'host=%s\npid=%s\nrun_id=%s\n' "$(hostname -s)" "$$" "${RUN_ID}" > "${LOCK_DIR}/owner"

write_status() {
  {
    printf 'run_id\tstatus\texit_code\tstage\thost\tn_core\tgit_head\tupdated_at\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${RUN_ID}" "$1" "$2" "$3" "$(hostname -s)" "${N_CORE}" \
      "$(git -C "${REPO_ROOT}" rev-parse HEAD)" "$(date -Iseconds)"
  } > "${STATUS_PATH}"
}

cleanup() {
  exit_code="$?"
  set +e
  if [[ "${RUN_STATUS}" != "COMPLETE" ]]; then write_status FAILED "${exit_code}" "${RUN_STATUS}"; fi
  [[ -f "${RUN_LOG}" ]] && cp "${RUN_LOG}" "${LATEST_LOG}"
  [[ -n "${TASK_TMP_DIR}" && -d "${TASK_TMP_DIR}" ]] && rm -rf -- "${TASK_TMP_DIR}"
  if [[ "${LOCK_ACQUIRED}" == TRUE ]]; then
    rm -f -- "${LOCK_DIR}/owner"
    rmdir "${LOCK_DIR}" 2>/dev/null || true
  fi
  trap - EXIT
  exit "${exit_code}"
}
trap cleanup EXIT
exec > >(tee -a "${RUN_LOG}") 2>&1

echo "Figure 6 passage run start: $(date -Iseconds)"
echo "run_id=${RUN_ID} host=$(hostname -s) n_core=${N_CORE} memory_kb=${MEMORY_KB}"
echo "git_head=$(git -C "${REPO_ROOT}" rev-parse HEAD)"
echo "model_code_root=${MODEL_CODE_ROOT}"
echo "joint_result_root=${JOINT_RESULT_ROOT}"
RUN_STATUS="PREFLIGHT"
write_status RUNNING 0 "${RUN_STATUS}"

{
  printf 'md5\tsize_bytes\tpath\n'
  for path in "${MODEL_FILES[@]}" \
    "${BASE_RUN}/q10_unique_endpoint_manifest.tsv" \
    "${BASE_RUN}/q10_optimizer_seed_manifest.tsv"; do
    printf '%s\t%s\t%s\n' "$(md5sum "${path}" | awk '{print $1}')" "$(stat -c '%s' "${path}")" "${path}"
  done
} > "${INPUT_MD5_PATH}"

TASK_TMP_DIR="$(mktemp -d /tmp/figure6-invitro-passage.XXXXXX)"
mkdir -p "${TASK_TMP_DIR}/home" "${TASK_TMP_DIR}/cache" "${TASK_TMP_DIR}/model_rcpp_cache"
printf '%s\n' 'options(bitmapType = "cairo", device = "png", warn = 1)' > "${TASK_TMP_DIR}/Rprofile"

CONTAINER_ARGS=(
  exec --cleanenv --containall --pwd "${ITERATION_ROOT}"
  --home "${TASK_TMP_DIR}/home"
  --env "TMPDIR=${TASK_TMP_DIR}"
  --env "XDG_CACHE_HOME=${TASK_TMP_DIR}/cache"
  --env "R_PROFILE_USER=${TASK_TMP_DIR}/Rprofile"
  --env "R_HOME=/opt/R/4.4.2/lib64/R"
  --env "LD_LIBRARY_PATH=${CONTAINER_LD_LIBRARY_PATH}"
  --env "DISPLAY=" --env "QT_QPA_PLATFORM=offscreen"
  --env "FIGURE_WORKSPACE_ROOT=${ITERATION_ROOT}"
  --env "FIGURE_MODEL_CODE_ROOT=${MODEL_CODE_ROOT}"
  --env "FIGURE_JOINT_RESULT_ROOT=${JOINT_RESULT_ROOT}"
  --env "FIGURE6_FINITE_TIME_FUTURE_PLAN=multicore"
  --env "OMP_NUM_THREADS=1" --env "OPENBLAS_NUM_THREADS=1"
  --env "MKL_NUM_THREADS=1" --env "VECLIB_MAXIMUM_THREADS=1"
  --env "RCPP_PARALLEL_NUM_THREADS=1" --env "KMP_USE_SHM=0"
  --env "KMP_INIT_AT_FORK=FALSE"
  --bind "${RED_EASYBUILD_ROOT}:${RED_EASYBUILD_ROOT}:ro"
  --bind "${ITERATION_ROOT}:${ITERATION_ROOT}:rw"
  --bind "${MODEL_CODE_ROOT}:${MODEL_CODE_ROOT}:ro"
  --bind "${TASK_TMP_DIR}/model_rcpp_cache:${MODEL_CODE_ROOT}/model/.rcpp_cache_o2_supply_demand_map:rw"
  --bind "${RESULTS_ROOT}:${RESULTS_ROOT}:ro"
  --bind "${TASK_TMP_DIR}:${TASK_TMP_DIR}:rw"
)
container_command() { "${CONTAINER_RUNTIME}" "${CONTAINER_ARGS[@]}" "${SIF_IMAGE}" "$@"; }

container_command Rscript --vanilla -e '
packages <- c("Matrix", "Rcpp", "data.table", "future", "future.apply", "ggplot2", "patchwork", "magick", "scales")
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
workspace <- normalizePath(Sys.getenv("FIGURE_WORKSPACE_ROOT"), mustWork = TRUE)
source(file.path(workspace, "Code/Figures/util/analysis/figure6_robustness.R"))
source(file.path(workspace, "Code/Figures/util/analysis/figure6_context_extension.R"))
source(file.path(workspace, "Code/Figures/util/analysis/figure6_finite_time_q10.R"))
source(file.path(workspace, "Code/Figures/util/analysis/figure6_invitro_passage_q10.R"))
source(file.path(workspace, "Code/Figures/util/analysis/figure6_finite_time_plots.R"))
paths <- f6r_paths(workspace)
stopifnot(identical(paths$oxygen_code, normalizePath(Sys.getenv("FIGURE_MODEL_CODE_ROOT"), mustWork = TRUE)))
f6r_load_response_engine(paths)
stopifnot(exists("fixo2_fixed_matrix"), exists("f6p_data"), exists("f6ft_draw_supp6_8"))
png_path <- file.path(tempdir(), "headless.png")
png(png_path, width = 200, height = 200, type = "cairo"); plot.new(); dev.off()
stopifnot(file.exists(png_path), file.info(png_path)$size > 0)
cat("figure6_invitro_passage_preflight_ok\n")
'

if [[ "${PREFLIGHT_ONLY}" == TRUE ]]; then
  RUN_STATUS="COMPLETE"; write_status COMPLETE 0 PREFLIGHT_ONLY; exit 0
fi

if [[ "${SMOKE_ONLY}" == TRUE ]]; then
  RUN_STATUS="SMOKE"; write_status RUNNING 0 "${RUN_STATUS}"
  container_command Rscript --vanilla "${CODE_ROOT}/data_Figure6_invitro_passage_q10.R" \
    "--n-core=${N_CORE}" "--run-id=${RUN_ID}_smoke" --smoke=TRUE --publish-current=FALSE
  RUN_STATUS="COMPLETE"; write_status COMPLETE 0 SMOKE_ONLY; exit 0
fi

RUN_STATUS="PASSAGE_DATA"; write_status RUNNING 0 "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/data_Figure6_invitro_passage_q10.R" \
  "--n-core=${N_CORE}" "--run-id=${RUN_ID}" --smoke=FALSE --publish-current=TRUE

RUN_STATUS="DRAW_SUPP_6_8"; write_status RUNNING 0 "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/draw_Supp_Figure6_8.R"

RUN_STATUS="DRAW_FIGURE6"; write_status RUNNING 0 "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/draw_Figure6.R"

RUN_STATUS="VALIDATE"; write_status RUNNING 0 "${RUN_STATUS}"
RUN_ROOT="${ITERATION_ROOT}/data/Figures/Figure6/invitro_passage_q10_runs/${RUN_ID}"
check_validation() {
  require_file "$1" "validation table"
  awk -F '\t' '
    NR == 1 {for (i=1; i<=NF; i++) if ($i=="passed") c=i; next}
    c>0 && toupper($c)!="TRUE" {bad=1}
    END {if (NR<2 || c==0 || bad) exit 1}
  ' "$1" || { echo "Validation failed: $1" >&2; exit 1; }
}
check_validation "${RUN_ROOT}/passage_task_validation.tsv"
check_validation "${RUN_ROOT}/passage_panel_validation.tsv"
check_validation "${RUN_ROOT}/figure6_finite_time_render_validation.tsv"
check_validation "${ITERATION_ROOT}/data/Figures/Supp_Figure6_8/supp_fig6-8_render_validation.tsv"
require_file "${RUN_ROOT}/passage_schedule.tsv" "passage schedule"
require_file "${RUN_ROOT}/passage_time_coverage.tsv" "passage coverage"
require_file "${RUN_ROOT}/passage_vs_continuous_comparison.tsv" "passage comparison"

OUTPUT_FILES=(
  "${ITERATION_ROOT}/Figures/assembled_fig6.png"
  "${ITERATION_ROOT}/Figures/assembled_fig6.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig6-8_continuous_invitro_no_passaging.png"
  "${ITERATION_ROOT}/Figures/supp_fig6-8_continuous_invitro_no_passaging.pdf"
)
for path in "${OUTPUT_FILES[@]}"; do require_file "${path}" "rendered output"; done
{
  printf 'sha256\tsize_bytes\tpath\n'
  for path in "${OUTPUT_FILES[@]}"; do
    printf '%s\t%s\t%s\n' "$(sha256sum "${path}" | awk '{print $1}')" "$(stat -c '%s' "${path}")" "${path}"
  done
} > "${OUTPUT_SHA256_PATH}"

RUN_STATUS="COMPLETE"
write_status COMPLETE 0 COMPLETE
echo "Figure 6 passage run complete: $(date -Iseconds)"
echo "status_path=${STATUS_PATH}"
echo "output_sha256=${OUTPUT_SHA256_PATH}"
