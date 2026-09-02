#!/usr/bin/env bash

# Fresh, headless Figure 7 fixed-p_misseg computation and rendering on
# hpctpa3pc0028. Figure 6 code, data, and rendered assets are guarded by a
# byte-level SHA256 inventory before and after the run.

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
RUN_ID="figure7_pmisseg_$(date '+%Y%m%d_%H%M%S')"
PREFLIGHT_ONLY=FALSE
DATA_ONLY=FALSE
DRAW_ONLY=FALSE
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
    --data-only) DATA_ONLY=TRUE ;;
    --draw-only) DRAW_ONLY=TRUE ;;
    -h|--help)
      printf '%s\n' \
        "Usage: $0 [--n-core=1..63] [--run-id=ID]" \
        "          [--preflight-only|--data-only|--draw-only]"
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
N_SPECIAL_MODE=0
[[ "${PREFLIGHT_ONLY}" == TRUE ]] && N_SPECIAL_MODE=$((N_SPECIAL_MODE + 1))
[[ "${DATA_ONLY}" == TRUE ]] && N_SPECIAL_MODE=$((N_SPECIAL_MODE + 1))
[[ "${DRAW_ONLY}" == TRUE ]] && N_SPECIAL_MODE=$((N_SPECIAL_MODE + 1))
if (( N_SPECIAL_MODE > 1 )); then
  echo "Special modes are mutually exclusive." >&2
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
require_dir "${INVIVO_RESULT_ROOT}" "in-vivo result root"
require_dir "${INVITRO_RESULT_ROOT}" "in-vitro result root"
require_dir "${JOINT_RESULT_ROOT}" "joint result root"
require_file "${ITERATION_ROOT}/provenance/figure6_frozen_before_figure7.tsv" "Figure 6 freeze baseline"
require_file "${CODE_ROOT}/util/analysis/figure7_freeze_guard.R" "Figure 6 freeze guard"
require_dir "${RED_EASYBUILD_ROOT}" "RED EasyBuild root"
require_file "${RED_FLEXIBLAS_LIB}/libflexiblas.so.3" "FlexiBLAS runtime"

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

CODE_FILES=(
  "${CODE_ROOT}/data_Figure7.R"
  "${CODE_ROOT}/data_Figure7_finite_time_q10.R"
  "${CODE_ROOT}/data_Figure7_invitro_passage_q10.R"
  "${CODE_ROOT}/data_Supp_Figure7_1.R"
  "${CODE_ROOT}/data_Supp_Figure7_2.R"
  "${CODE_ROOT}/data_Supp_Figure7_3.R"
  "${CODE_ROOT}/data_Supp_Figure7_4.R"
  "${CODE_ROOT}/data_Supp_Figure7_11_12.R"
  "${CODE_ROOT}/draw_Figure7.R"
  "${CODE_ROOT}/util/analysis/figure7_robustness.R"
  "${CODE_ROOT}/util/analysis/figure7_context_extension.R"
  "${CODE_ROOT}/util/analysis/figure7_finite_time_q10.R"
  "${CODE_ROOT}/util/analysis/figure7_invitro_passage_q10.R"
  "${CODE_ROOT}/util/analysis/figure7_extended_time_o2.R"
  "${CODE_ROOT}/util/analysis/figure7_finite_time_plots.R"
)
for index in $(seq 1 12); do CODE_FILES+=("${CODE_ROOT}/draw_Supp_Figure7_${index}.R"); done
for path in "${CODE_FILES[@]}"; do require_file "${path}" "Figure 7 source"; done

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

AUDIT_ROOT="${ITERATION_ROOT}/audit/hpc_figure7_pmisseg/${RUN_ID}"
LOG_ROOT="${ITERATION_ROOT}/audit/logs"
LOCK_ROOT="${ITERATION_ROOT}/audit/locks"
LOCK_DIR="${LOCK_ROOT}/figure7_pmisseg.lock"
RUN_LOG="${LOG_ROOT}/figure7_pmisseg_${RUN_ID}.log"
LATEST_LOG="${LOG_ROOT}/figure7_pmisseg_latest.log"
STATUS_PATH="${AUDIT_ROOT}/status.tsv"
INPUT_SHA256_PATH="${AUDIT_ROOT}/input_sha256.tsv"
OUTPUT_SHA256_PATH="${AUDIT_ROOT}/output_sha256.tsv"
PDF_FONT_VALIDATION_PATH="${AUDIT_ROOT}/pdf_font_validation.tsv"
FIGURE6_GUARD_PATH="${AUDIT_ROOT}/figure6_frozen_after_figure7_check.tsv"
TASK_TMP_DIR=""
RUN_STATUS="INITIALIZING"
LOCK_ACQUIRED=FALSE

mkdir -p "${AUDIT_ROOT}" "${LOG_ROOT}" "${LOCK_ROOT}"
if ! mkdir "${LOCK_DIR}" 2>/dev/null; then
  echo "Another Figure 7 run owns ${LOCK_DIR}." >&2
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
  if [[ "${RUN_STATUS}" != "COMPLETE" ]]; then
    write_status FAILED "${exit_code}" "${RUN_STATUS}"
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

echo "Figure 7 fixed-p_misseg run start: $(date -Iseconds)"
echo "run_id=${RUN_ID} host=$(hostname -s) n_core=${N_CORE} memory_kb=${MEMORY_KB}"
echo "model_code_root=${MODEL_CODE_ROOT}"
echo "joint_result_root=${JOINT_RESULT_ROOT}"
RUN_STATUS="PREFLIGHT"
write_status RUNNING 0 "${RUN_STATUS}"

{
  printf 'sha256\tsize_bytes\tpath\n'
  for path in "${MODEL_FILES[@]}" "${CODE_FILES[@]}"; do
    printf '%s\t%s\t%s\n' \
      "$(sha256sum "${path}" | awk '{print $1}')" \
      "$(stat -c '%s' "${path}")" "${path}"
  done
} > "${INPUT_SHA256_PATH}"

TASK_TMP_DIR="$(mktemp -d /tmp/figure7-pmisseg.XXXXXX)"
mkdir -p "${TASK_TMP_DIR}/home" "${TASK_TMP_DIR}/cache" \
  "${TASK_TMP_DIR}/fontconfig" "${TASK_TMP_DIR}/model_rcpp_cache"
printf '%s\n' 'options(bitmapType = "cairo", device = "png", warn = 1)' > "${TASK_TMP_DIR}/Rprofile"

CONTAINER_ARGS=(
  exec --cleanenv --containall --pwd "${ITERATION_ROOT}"
  --home "${TASK_TMP_DIR}/home"
  --env "TMPDIR=${TASK_TMP_DIR}"
  --env "XDG_CACHE_HOME=${TASK_TMP_DIR}/cache"
  --env "R_PROFILE_USER=${TASK_TMP_DIR}/Rprofile"
  --env "R_HOME=/opt/R/4.4.2/lib64/R"
  --env "LD_LIBRARY_PATH=${CONTAINER_LD_LIBRARY_PATH}"
  --env "DISPLAY=" --env "QT_QPA_PLATFORM=offscreen" --env "MPLBACKEND=Agg"
  --env "FIGURE_WORKSPACE_ROOT=${ITERATION_ROOT}"
  --env "FIGURE_MODEL_CODE_ROOT=${MODEL_CODE_ROOT}"
  --env "FIGURE_INVIVO_RESULT_ROOT=${INVIVO_RESULT_ROOT}"
  --env "FIGURE_INVITRO_RESULT_ROOT=${INVITRO_RESULT_ROOT}"
  --env "FIGURE_JOINT_RESULT_ROOT=${JOINT_RESULT_ROOT}"
  --env "FIGURE7_AUDIT_RUN_ID=${RUN_ID}"
  --env "FIGURE7_FINITE_TIME_FUTURE_PLAN=multicore"
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
packages <- c(
  "Matrix", "Rcpp", "RcppEigen", "cluster", "data.table", "dplyr",
  "future", "future.apply", "ggplot2", "matrixStats", "tidyr",
  "patchwork", "magick", "isoband", "scales"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
workspace <- normalizePath(Sys.getenv("FIGURE_WORKSPACE_ROOT"), mustWork = TRUE)
model_root <- normalizePath(Sys.getenv("FIGURE_MODEL_CODE_ROOT"), mustWork = TRUE)
code_root <- file.path(workspace, "Code", "Figures")
files <- c(
  list.files(code_root, pattern = "Figure7.*[.]R$", full.names = TRUE),
  list.files(file.path(code_root, "util", "analysis"),
             pattern = "figure7.*[.]R$", full.names = TRUE),
  file.path(code_root, "util", "analysis", "si_figure7_eigenmodes.R")
)
invisible(lapply(unique(files), parse))
source(file.path(code_root, "util", "analysis", "figure7_robustness.R"))
paths <- f7r_paths(workspace)
stopifnot(identical(paths$oxygen_code, model_root))
f7r_load_response_engine(paths)
stopifnot(
  exists("fixo2_fixed_matrix", envir = globalenv()),
  exists("fixo2_dominant_attractor_one", envir = globalenv()),
  exists(".pmisseg_of_O2", envir = globalenv()),
  exists("figure7_force_p_misseg", envir = globalenv())
)
png_path <- file.path(tempdir(), "figure7_headless.png")
grDevices::png(png_path, width = 240, height = 240, type = "cairo", bg = "white")
graphics::plot.new(); grDevices::dev.off()
stopifnot(file.exists(png_path), file.info(png_path)$size > 0)
cat("figure7_container_preflight_ok\n")
'

container_command Rscript --vanilla \
  "${CODE_ROOT}/util/analysis/figure7_freeze_guard.R" verify

if [[ "${PREFLIGHT_ONLY}" == TRUE ]]; then
  RUN_STATUS="COMPLETE"; write_status COMPLETE 0 PREFLIGHT_ONLY; exit 0
fi

if [[ "${DRAW_ONLY}" != TRUE ]]; then
  RUN_STATUS="DATA_FIGURE7_AB"
  write_status RUNNING 0 "${RUN_STATUS}"
  container_command Rscript --vanilla "${CODE_ROOT}/data_Figure7.R" \
    "--n-core=${N_CORE}" --rebuild=TRUE --n-resample=100

  RUN_STATUS="DATA_SUPP_FIGURE7_1_TO_4"
  write_status RUNNING 0 "${RUN_STATUS}"
  container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure7_1.R" \
    "--n-core=${N_CORE}" --rebuild=TRUE
  container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure7_2.R"
  container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure7_3.R" \
    "--n-core=${N_CORE}" --rebuild=FALSE
  container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure7_4.R" \
    "--n-core=${N_CORE}" --rebuild=TRUE

  RUN_STATUS="DATA_FIGURE7_CD_AND_DIAGNOSTICS"
  write_status RUNNING 0 "${RUN_STATUS}"
  container_command Rscript --vanilla "${CODE_ROOT}/data_Figure7_finite_time_q10.R" \
    "--n-core=${N_CORE}" "--run-id=${RUN_ID}_finite" \
    --smoke=FALSE --publish-current=TRUE --compute-diagnostics=TRUE

  RUN_STATUS="DATA_FIGURE7_EF_PASSAGING"
  write_status RUNNING 0 "${RUN_STATUS}"
  container_command Rscript --vanilla "${CODE_ROOT}/data_Figure7_invitro_passage_q10.R" \
    "--n-core=${N_CORE}" "--run-id=${RUN_ID}_passage" \
    --smoke=FALSE --publish-current=TRUE

  RUN_STATUS="DATA_SUPP_FIGURE7_11_12"
  write_status RUNNING 0 "${RUN_STATUS}"
  container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure7_11_12.R" \
    --mode=all "--n-core=${N_CORE}" "--run-id=${RUN_ID}_extended" \
    --smoke=FALSE --publish-current=TRUE
fi

if [[ "${DATA_ONLY}" == TRUE ]]; then
  container_command Rscript --vanilla \
    "${CODE_ROOT}/util/analysis/figure7_freeze_guard.R" verify
  RUN_STATUS="COMPLETE"; write_status COMPLETE 0 DATA_ONLY; exit 0
fi

RUN_STATUS="DRAW_FIGURE7_AND_SUPPLEMENTS"
write_status RUNNING 0 "${RUN_STATUS}"
for script in \
  draw_Supp_Figure7_1.R draw_Supp_Figure7_2.R \
  draw_Supp_Figure7_3.R draw_Supp_Figure7_4.R \
  draw_Supp_Figure7_5.R draw_Supp_Figure7_6.R \
  draw_Supp_Figure7_7.R draw_Supp_Figure7_8.R \
  draw_Supp_Figure7_9.R draw_Supp_Figure7_10.R \
  draw_Supp_Figure7_11.R draw_Supp_Figure7_12.R \
  draw_Figure7.R; do
  container_command Rscript --vanilla "${CODE_ROOT}/${script}"
done

RUN_STATUS="VALIDATE"
write_status RUNNING 0 "${RUN_STATUS}"

check_validation() {
  local path="$1"
  require_file "${path}" "validation table"
  awk -F '\t' '
    NR == 1 {for (i=1; i<=NF; i++) if ($i=="passed") c=i; next}
    c>0 && toupper($c)!="TRUE" {bad=1}
    END {if (NR<2 || c==0 || bad) exit 1}
  ' "${path}" || { echo "Validation failed: ${path}" >&2; exit 1; }
}

VALIDATION_COUNT=0
while IFS= read -r path; do
  if head -n 1 "${path}" | tr '\t' '\n' | grep -Fx passed >/dev/null; then
    check_validation "${path}"
    VALIDATION_COUNT=$((VALIDATION_COUNT + 1))
  fi
done < <(find \
  "${ITERATION_ROOT}/data/Figures/Figure7" \
  "${ITERATION_ROOT}/data/Figures"/Supp_Figure7_* \
  -type f -name '*validation.tsv' -print | sort)
if (( VALIDATION_COUNT < 20 )); then
  echo "Expected at least 20 Figure 7 validation tables; observed ${VALIDATION_COUNT}." >&2
  exit 1
fi

mapfile -t OUTPUT_FILES < <(find "${ITERATION_ROOT}/Figures" -maxdepth 1 -type f \
  \( -name 'assembled_fig7.png' -o -name 'assembled_fig7.pdf' \
     -o -name 'supp_fig7-*.png' -o -name 'supp_fig7-*.pdf' \) -print | sort)
if (( ${#OUTPUT_FILES[@]} != 26 )); then
  echo "Expected 26 Figure 7 PNG/PDF outputs; observed ${#OUTPUT_FILES[@]}." >&2
  exit 1
fi
for path in "${OUTPUT_FILES[@]}"; do require_file "${path}" "Figure 7 output"; done

require_file "$(command -v pdffonts)" "pdffonts executable"
require_file "$(command -v pdftotext)" "pdftotext executable"
{
  printf 'path\ttype3_detected\ttext_bytes\tpassed\n'
  for path in "${OUTPUT_FILES[@]}"; do
    [[ "${path}" == *.pdf ]] || continue
    type3_detected=FALSE
    if pdffonts "${path}" | grep -F 'Type 3' >/dev/null; then
      type3_detected=TRUE
    fi
    text_bytes="$(pdftotext -layout "${path}" - | wc -c | tr -d ' ')"
    passed=FALSE
    if [[ "${type3_detected}" == FALSE ]] && (( text_bytes > 50 )); then passed=TRUE; fi
    printf '%s\t%s\t%s\t%s\n' "${path}" "${type3_detected}" "${text_bytes}" "${passed}"
    [[ "${passed}" == TRUE ]] || { echo "PDF text/font validation failed: ${path}" >&2; exit 1; }
  done
} > "${PDF_FONT_VALIDATION_PATH}"

{
  printf 'sha256\tsize_bytes\tpath\n'
  for path in "${OUTPUT_FILES[@]}"; do
    printf '%s\t%s\t%s\n' \
      "$(sha256sum "${path}" | awk '{print $1}')" \
      "$(stat -c '%s' "${path}")" "${path}"
  done
} > "${OUTPUT_SHA256_PATH}"

container_command Rscript --vanilla -e '
workspace <- normalizePath(Sys.getenv("FIGURE_WORKSPACE_ROOT"), mustWork = TRUE)
source(file.path(workspace, "Code/Figures/util/analysis/figure7_freeze_guard.R"))
f7g_verify(
  workspace_root = workspace,
  output = file.path(
    workspace, "audit/hpc_figure7_pmisseg", Sys.getenv("FIGURE7_AUDIT_RUN_ID"),
    "figure6_frozen_after_figure7_check.tsv"
  )
)
'

RUN_STATUS="COMPLETE"
write_status COMPLETE 0 COMPLETE
echo "Figure 7 fixed-p_misseg run complete: $(date -Iseconds)"
echo "status_path=${STATUS_PATH}"
echo "output_sha256=${OUTPUT_SHA256_PATH}"
