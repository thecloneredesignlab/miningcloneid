#!/usr/bin/env bash

# Compute Supplementary Figure 6-14 net-growth data and render Supplementary
# Figures 6-14 through 6-16 headlessly on hpctpa3pc0028. Supplementary Figures
# 6-15/16 are slices of the existing full-range Figure 6 cache.

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
RED_EASYBUILD_ROOT="/app/eb"
RED_FLEXIBLAS_LIB="${RED_EASYBUILD_ROOT}/software/FlexiBLAS/3.4.4-GCC-13.3.0/lib64"
RED_OPENBLAS_LIB="${RED_EASYBUILD_ROOT}/software/OpenBLAS/0.3.27-GCC-13.3.0/lib"
RED_GCC_LIB="${RED_EASYBUILD_ROOT}/software/GCCcore/13.3.0/lib64"
RED_BINUTILS_LIB="${RED_EASYBUILD_ROOT}/software/binutils/2.42-GCCcore-13.3.0/lib"
CONTAINER_LD_LIBRARY_PATH="${RED_FLEXIBLAS_LIB}:${RED_OPENBLAS_LIB}:${RED_GCC_LIB}:${RED_BINUTILS_LIB}:/opt/rh/gcc-toolset-13/root/usr/lib64"

N_CORE=60
O2_CHUNK_SIZE=1
RUN_ID="supp6_14_net_growth_$(date '+%Y%m%d_%H%M%S')"
PREFLIGHT_ONLY=FALSE
DRAW_ONLY=FALSE
FINALIZE_EXISTING=FALSE

for argument in "$@"; do
  case "${argument}" in
    --n-core=*) N_CORE="${argument#*=}" ;;
    --o2-chunk-size=*) O2_CHUNK_SIZE="${argument#*=}" ;;
    --run-id=*) RUN_ID="${argument#*=}" ;;
    --preflight-only) PREFLIGHT_ONLY=TRUE ;;
    --draw-only) DRAW_ONLY=TRUE ;;
    --finalize-existing) FINALIZE_EXISTING=TRUE ;;
    -h|--help)
      printf '%s\n' \
        "Usage: $0 [--n-core=1..63] [--o2-chunk-size=N] [--run-id=ID]" \
        "          [--preflight-only|--draw-only] [--finalize-existing]"
      exit 0 ;;
    *) echo "Unknown option: ${argument}" >&2; exit 2 ;;
  esac
done

if ! [[ "${N_CORE}" =~ ^[1-9][0-9]*$ ]] || (( N_CORE > 63 )); then
  echo "--n-core must be an integer from 1 through 63." >&2; exit 2
fi
if ! [[ "${O2_CHUNK_SIZE}" =~ ^[1-9][0-9]*$ ]]; then
  echo "--o2-chunk-size must be a positive integer." >&2; exit 2
fi
if ! [[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$ ]]; then
  echo "Invalid --run-id." >&2; exit 2
fi
if [[ "${PREFLIGHT_ONLY}" == TRUE && "${DRAW_ONLY}" == TRUE ]]; then
  echo "--preflight-only and --draw-only are mutually exclusive." >&2; exit 2
fi
if [[ "${FINALIZE_EXISTING}" == TRUE && "${DRAW_ONLY}" != TRUE ]]; then
  echo "--finalize-existing requires --draw-only." >&2; exit 2
fi
[[ "$(hostname -s)" == "${EXPECTED_NODE}" ]] || {
  echo "This runner must execute on ${EXPECTED_NODE}." >&2; exit 2;
}
[[ "${REPO_ROOT}" == "${EXPECTED_REPO_ROOT}" ]] || {
  echo "Unexpected repository root: ${REPO_ROOT}" >&2; exit 2;
}

require_file() { [[ -f "$1" && -r "$1" ]] || { echo "Missing readable $2: $1" >&2; exit 2; }; }
require_dir() { [[ -d "$1" && -r "$1" ]] || { echo "Missing readable $2: $1" >&2; exit 2; }; }
require_command() { command -v "$1" >/dev/null 2>&1 || { echo "Missing $2: $1" >&2; exit 2; }; }

require_file "${SIF_IMAGE}" "SIF image"
require_dir "${MODEL_CODE_ROOT}" "external model-code root"
require_dir "${INVIVO_RESULT_ROOT}" "in-vivo result root"
require_dir "${INVITRO_RESULT_ROOT}" "in-vitro result root"
require_dir "${JOINT_RESULT_ROOT}" "joint result root"
require_file "${RED_FLEXIBLAS_LIB}/libflexiblas.so.3" "FlexiBLAS runtime"
require_file "${CODE_ROOT}/data_Supp_Figure6_14.R" "net-growth data entry point"
require_file "${CODE_ROOT}/finalize_Supp_Figure6_14.R" "net-growth finalizer"
require_file "${CODE_ROOT}/util/analysis/figure6_net_growth_q10.R" "net-growth implementation"
require_file "${CODE_ROOT}/util/analysis/figure6_net_growth_propagator.cpp" "net-growth C++ propagator"
require_file "${CODE_ROOT}/util/analysis/figure6_supplementary_b_layout.R" "B-only layout"
for index in 14 15 16; do
  require_file "${CODE_ROOT}/draw_Supp_Figure6_${index}.R" "supplement drawing entry point"
done
SOURCE_FIXED_ROOT="${ITERATION_ROOT}/data/Figures/Figure7/fixed_pmisseg_v1"
SOURCE_POINTER="${SOURCE_FIXED_ROOT}/finite_time_full_q10_current.tsv"
require_file "${SOURCE_POINTER}" "full-range source pointer"
SOURCE_RELATIVE_RUN="$(awk -F '\t' 'NR==2 {print $2}' "${SOURCE_POINTER}")"
SOURCE_RUN_ROOT="${SOURCE_FIXED_ROOT}/${SOURCE_RELATIVE_RUN}"
require_dir "${SOURCE_RUN_ROOT}" "full-range source run"
require_command pdftotext "PDF text validation"
require_command pdffonts "PDF font validation"

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v apptainer)"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v singularity)"
else
  echo "No apptainer/singularity runtime found." >&2; exit 2
fi

AVAILABLE_CPU="$(getconf _NPROCESSORS_ONLN)"
MEMORY_KB="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
(( AVAILABLE_CPU >= N_CORE + 1 )) || {
  echo "Need $((N_CORE + 1)) CPUs; observed ${AVAILABLE_CPU}." >&2; exit 2;
}
[[ -n "${MEMORY_KB}" ]] && (( MEMORY_KB >= 500000000 )) || {
  echo "Expected approximately 512 GB RAM; observed ${MEMORY_KB:-unknown} kB." >&2; exit 2;
}

AUDIT_ROOT="${ITERATION_ROOT}/audit/hpc_supp_figure6_14_16/${RUN_ID}"
LOG_ROOT="${ITERATION_ROOT}/audit/logs"
LOCK_ROOT="${ITERATION_ROOT}/audit/locks"
LOCK_DIR="${LOCK_ROOT}/supp_figure6_14_16.lock"
RUN_LOG="${LOG_ROOT}/supp_figure6_14_16_${RUN_ID}.log"
STATUS_PATH="${AUDIT_ROOT}/status.tsv"
OUTPUT_SHA256="${AUDIT_ROOT}/output_sha256.tsv"
TASK_TMP_DIR=""
RUN_STATUS="INITIALIZING"

mkdir -p "${AUDIT_ROOT}" "${LOG_ROOT}" "${LOCK_ROOT}"
mkdir "${LOCK_DIR}" 2>/dev/null || {
  echo "Another Supplementary Figure 6-14/16 run owns ${LOCK_DIR}." >&2; exit 2;
}
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
  if [[ -n "${TASK_TMP_DIR}" && -d "${TASK_TMP_DIR}" ]]; then
    rm -rf -- "${TASK_TMP_DIR}"
  fi
  rm -f -- "${LOCK_DIR}/owner"
  rmdir "${LOCK_DIR}" 2>/dev/null || true
  trap - EXIT
  exit "${exit_code}"
}
trap cleanup EXIT
exec > >(tee -a "${RUN_LOG}") 2>&1

echo "Supplementary Figure 6-14/16 run start: $(date -Iseconds)"
echo "run_id=${RUN_ID} host=$(hostname -s) n_core=${N_CORE}"
echo "model_code_root=${MODEL_CODE_ROOT}"
RUN_STATUS="PREFLIGHT"
write_status RUNNING 0 "${RUN_STATUS}"
sha256sum "${SIF_IMAGE}" > "${AUDIT_ROOT}/container_image.sha256"

mkdir -p "${ITERATION_ROOT}/audit/tmp"
TASK_TMP_DIR="$(mktemp -d "${ITERATION_ROOT}/audit/tmp/supp6-14-16.XXXXXX")"
mkdir -p "${TASK_TMP_DIR}/home" "${TASK_TMP_DIR}/cache" \
  "${TASK_TMP_DIR}/model_rcpp_cache"
printf '%s\n' 'options(bitmapType = "cairo", device = "png", warn = 1)' > "${TASK_TMP_DIR}/Rprofile"

CONTAINER_ARGS=(
  exec --cleanenv --containall --pwd "${ITERATION_ROOT}"
  --home "${TASK_TMP_DIR}/home"
  --env "TMPDIR=/tmp" --env "TMP=/tmp" --env "TEMP=/tmp"
  --env "XDG_CACHE_HOME=${TASK_TMP_DIR}/cache"
  --env "FONTCONFIG_PATH=/etc/fonts"
  --env "R_PROFILE_USER=${TASK_TMP_DIR}/Rprofile"
  --env "R_MAKEVARS_USER=${SCRIPT_DIR}/figure6_compiler_makevars"
  --env "R_HOME=/opt/R/4.4.2/lib64/R"
  --env "LD_LIBRARY_PATH=${CONTAINER_LD_LIBRARY_PATH}"
  --env "DISPLAY=" --env "QT_QPA_PLATFORM=offscreen" --env "MPLBACKEND=Agg"
  --env "FIGURE_WORKSPACE_ROOT=${ITERATION_ROOT}"
  --env "FIGURE_MODEL_CODE_ROOT=${MODEL_CODE_ROOT}"
  --env "FIGURE_INVIVO_RESULT_ROOT=${INVIVO_RESULT_ROOT}"
  --env "FIGURE_INVITRO_RESULT_ROOT=${INVITRO_RESULT_ROOT}"
  --env "FIGURE_JOINT_RESULT_ROOT=${JOINT_RESULT_ROOT}"
  --env "FIGURE6_FULL_RANGE_SOURCE_RUN_ROOT=${SOURCE_RUN_ROOT}"
  --env "FIGURE6_FINITE_TIME_FUTURE_PLAN=multicore"
  --env "OMP_NUM_THREADS=1" --env "OPENBLAS_NUM_THREADS=1"
  --env "MKL_NUM_THREADS=1" --env "VECLIB_MAXIMUM_THREADS=1"
  --env "RCPP_PARALLEL_NUM_THREADS=1" --env "KMP_USE_SHM=0"
  --env "KMP_INIT_AT_FORK=FALSE"
  --bind "${ITERATION_ROOT}:${ITERATION_ROOT}:rw"
  --bind "${RED_EASYBUILD_ROOT}:${RED_EASYBUILD_ROOT}:ro"
  --bind "${MODEL_CODE_ROOT}:${MODEL_CODE_ROOT}:ro"
  --bind "${TASK_TMP_DIR}/model_rcpp_cache:${MODEL_CODE_ROOT}/model/.rcpp_cache_o2_supply_demand_map:rw"
  --bind "${RESULTS_ROOT}:${RESULTS_ROOT}:ro"
  --bind "${TASK_TMP_DIR}:${TASK_TMP_DIR}:rw"
  --bind "${TASK_TMP_DIR}:/tmp:rw"
  --bind "${TASK_TMP_DIR}:/var/tmp:rw"
)
container_command() {
  "${CONTAINER_RUNTIME}" "${CONTAINER_ARGS[@]}" "${SIF_IMAGE}" "$@"
}

container_command Rscript --vanilla -e '
packages <- c("Matrix", "Rcpp", "future", "future.apply", "ggplot2", "scales")
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
workspace <- normalizePath(Sys.getenv("FIGURE_WORKSPACE_ROOT"), mustWork = TRUE)
model_root <- normalizePath(Sys.getenv("FIGURE_MODEL_CODE_ROOT"), mustWork = TRUE)
code_root <- file.path(workspace, "Code", "Figures")
files <- c(
  "util/analysis/figure6_net_growth_q10.R",
  "util/analysis/figure6_supplementary_b_layout.R",
  "data_Supp_Figure6_14.R", "draw_Supp_Figure6_14.R",
  "draw_Supp_Figure6_15.R", "draw_Supp_Figure6_16.R"
)
invisible(lapply(file.path(code_root, files), parse))
source(file.path(code_root, "util", "analysis", "figure6_robustness.R"))
paths <- f6r_paths(workspace)
stopifnot(identical(paths$oxygen_code, model_root))
png_path <- file.path(tempdir(), "supp6_headless.png")
grDevices::png(png_path, width = 240, height = 240, type = "cairo", bg = "white")
graphics::plot.new(); grDevices::dev.off()
stopifnot(file.exists(png_path), file.info(png_path)$size > 0)
cat("supp_figure6_14_16_preflight_ok\n")
'

if [[ "${PREFLIGHT_ONLY}" == TRUE ]]; then
  RUN_STATUS="COMPLETE"; write_status COMPLETE 0 PREFLIGHT_ONLY; exit 0
fi

if [[ "${FINALIZE_EXISTING}" == TRUE ]]; then
  RUN_STATUS="FINALIZE_EXISTING"
  write_status RUNNING 0 "${RUN_STATUS}"
  container_command Rscript --vanilla "${CODE_ROOT}/finalize_Supp_Figure6_14.R" \
    "--run-id=${RUN_ID}"
fi

if [[ "${DRAW_ONLY}" != TRUE ]]; then
  RUN_STATUS="COMPUTE_NET_GROWTH"
  write_status RUNNING 0 "${RUN_STATUS}"
  container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure6_14.R" \
    "--n-core=${N_CORE}" "--o2-chunk-size=${O2_CHUNK_SIZE}" \
    "--run-id=${RUN_ID}" --smoke=FALSE --publish-current=TRUE
fi

RUN_STATUS="HEADLESS_RENDER"
write_status RUNNING 0 "${RUN_STATUS}"
for script in draw_Supp_Figure6_14.R draw_Supp_Figure6_15.R draw_Supp_Figure6_16.R; do
  container_command Rscript --vanilla "${CODE_ROOT}/${script}"
done

RUN_STATUS="VALIDATE"
write_status RUNNING 0 "${RUN_STATUS}"
OUTPUTS=(
  "${ITERATION_ROOT}/Figures/supp_fig6-14_population_net_live_growth_rate.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig6-14_population_net_live_growth_rate.png"
  "${ITERATION_ROOT}/Figures/supp_fig6-15_finite_time_ploidy_o2_0_5_day500.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig6-15_finite_time_ploidy_o2_0_5_day500.png"
  "${ITERATION_ROOT}/Figures/supp_fig6-16_finite_time_ploidy_o2_0_5_day1000.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig6-16_finite_time_ploidy_o2_0_5_day1000.png"
)
VALIDATIONS=(
  "${ITERATION_ROOT}/data/Figures/Supp_Figure6_14/supp_fig6-14_population_net_live_growth_rate_render_validation.tsv"
  "${ITERATION_ROOT}/data/Figures/Supp_Figure6_15/supp_fig6-15_finite_time_ploidy_o2_0_5_day500_render_validation.tsv"
  "${ITERATION_ROOT}/data/Figures/Supp_Figure6_16/supp_fig6-16_finite_time_ploidy_o2_0_5_day1000_render_validation.tsv"
)
for path in "${VALIDATIONS[@]}"; do
  require_file "${path}" "render validation"
  awk -F '\t' '
    NR == 1 {for (i=1; i<=NF; i++) if ($i=="passed") c=i; next}
    c>0 && toupper($c)!="TRUE" {bad=1}
    END {if (NR<2 || c==0 || bad) exit 1}
  ' "${path}" || { echo "Validation failed: ${path}" >&2; exit 1; }
done

{
  printf 'sha256\tsize_bytes\tpath\n'
  for path in "${OUTPUTS[@]}"; do
    require_file "${path}" "supplement output"
    printf '%s\t%s\t%s\n' "$(sha256sum "${path}" | awk '{print $1}')" \
      "$(stat -c '%s' "${path}")" "${path}"
  done
} > "${OUTPUT_SHA256}"

{
  printf 'pdf\ttext_words_present\tfont_count\tno_type3\tpassed\n'
  for pdf in "${OUTPUTS[@]}"; do
    [[ "${pdf}" == *.pdf ]] || continue
    text_path="${TASK_TMP_DIR}/$(basename "${pdf}").txt"
    pdftotext -layout "${pdf}" "${text_path}"
    words_present=FALSE
    grep -Eiq 'oxygen|ploidy|growth|misseg' "${text_path}" && words_present=TRUE
    font_count="$(pdffonts "${pdf}" | awk 'NR>2 {n++} END {print n+0}')"
    no_type3=TRUE
    pdffonts "${pdf}" | awk 'NR>2 && $2=="Type 3" {bad=1} END {exit bad}' || no_type3=FALSE
    passed=FALSE
    [[ "${words_present}" == TRUE && "${font_count}" -gt 0 && "${no_type3}" == TRUE ]] && passed=TRUE
    printf '%s\t%s\t%s\t%s\t%s\n' "${pdf}" "${words_present}" \
      "${font_count}" "${no_type3}" "${passed}"
    [[ "${passed}" == TRUE ]] || exit 1
  done
} > "${AUDIT_ROOT}/pdf_text_font_validation.tsv"

RUN_STATUS="COMPLETE"
write_status COMPLETE 0 COMPLETE
echo "Supplementary Figure 6-14/16 run complete: $(date -Iseconds)"
echo "output_manifest=${OUTPUT_SHA256}"
