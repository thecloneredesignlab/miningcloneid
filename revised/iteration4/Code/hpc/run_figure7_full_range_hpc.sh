#!/usr/bin/env bash

# One fresh full-range Figure 7C computation followed by headless rendering of
# Figure 7 and Supplementary Figures 7-8 through 7-11 on hpctpa3pc0028.

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
O2_CHUNK_SIZE=10
RUN_ID="figure7_full_range_$(date '+%Y%m%d_%H%M%S')"
PREFLIGHT_ONLY=FALSE
DRAW_ONLY=FALSE
COMPUTE_ONLY=FALSE
REUSE_CONTINUOUS_RUN=""

for argument in "$@"; do
  case "${argument}" in
    --n-core=*) N_CORE="${argument#*=}" ;;
    --o2-chunk-size=*) O2_CHUNK_SIZE="${argument#*=}" ;;
    --run-id=*) RUN_ID="${argument#*=}" ;;
    --preflight-only) PREFLIGHT_ONLY=TRUE ;;
    --draw-only) DRAW_ONLY=TRUE ;;
    --compute-only) COMPUTE_ONLY=TRUE ;;
    --reuse-continuous-run=*) REUSE_CONTINUOUS_RUN="${argument#*=}" ;;
    -h|--help)
      printf '%s\n' \
        "Usage: $0 [--n-core=1..63] [--o2-chunk-size=N] [--run-id=ID]" \
        "          [--preflight-only|--draw-only]"
      exit 0 ;;
    *) echo "Unknown option: ${argument}" >&2; exit 2 ;;
  esac
done

if ! [[ "${N_CORE}" =~ ^[1-9][0-9]*$ ]] || (( N_CORE > 63 )); then
  echo "--n-core must be an integer from 1 through 63." >&2
  exit 2
fi
if ! [[ "${O2_CHUNK_SIZE}" =~ ^[1-9][0-9]*$ ]]; then
  echo "--o2-chunk-size must be a positive integer." >&2
  exit 2
fi
if ! [[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$ ]]; then
  echo "Invalid --run-id." >&2
  exit 2
fi
if [[ "${PREFLIGHT_ONLY}" == TRUE && "${DRAW_ONLY}" == TRUE ]]; then
  echo "--preflight-only and --draw-only are mutually exclusive." >&2
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
  [[ -f "$1" && -r "$1" ]] || {
    echo "Missing readable $2: $1" >&2
    exit 2
  }
}
require_dir() {
  [[ -d "$1" && -r "$1" ]] || {
    echo "Missing readable $2: $1" >&2
    exit 2
  }
}
require_command() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Missing required command $2: $1" >&2
    exit 2
  }
}

require_file "${SIF_IMAGE}" "SIF image"
require_dir "${MODEL_CODE_ROOT}" "external model-code root"
require_dir "${INVIVO_RESULT_ROOT}" "in-vivo result root"
require_dir "${INVITRO_RESULT_ROOT}" "in-vitro result root"
require_dir "${JOINT_RESULT_ROOT}" "joint result root"
require_dir "${RED_EASYBUILD_ROOT}" "RED EasyBuild root"
require_file "${RED_FLEXIBLAS_LIB}/libflexiblas.so.3" "FlexiBLAS runtime"
require_file "${CODE_ROOT}/data_Figure7_full_range_q10.R" "full-range entry point"
require_file "${CODE_ROOT}/draw_Figure7.R" "Figure 7 drawing entry point"
require_file "${CODE_ROOT}/util/analysis/figure7_full_range_q10.R" "full-range implementation"
require_file "${CODE_ROOT}/util/analysis/figure7_full_range_propagator.cpp" "C++ propagator"
for index in 8 9 10 11; do
  require_file "${CODE_ROOT}/draw_Supp_Figure7_${index}.R" "supplement drawing entry point"
done
require_command pdftotext "for PDF word-level validation"
require_command pdffonts "for PDF font validation"

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

AUDIT_ROOT="${ITERATION_ROOT}/audit/hpc_figure7_full_range/${RUN_ID}"
LOG_ROOT="${ITERATION_ROOT}/audit/logs"
LOCK_ROOT="${ITERATION_ROOT}/audit/locks"
LOCK_DIR="${LOCK_ROOT}/figure7_full_range.lock"
RUN_LOG="${LOG_ROOT}/figure7_full_range_${RUN_ID}.log"
LATEST_LOG="${LOG_ROOT}/figure7_full_range_latest.log"
STATUS_PATH="${AUDIT_ROOT}/status.tsv"
FIGURE6_BEFORE="${AUDIT_ROOT}/figure6_before.sha256"
FIGURE6_AFTER="${AUDIT_ROOT}/figure6_after.sha256"
OUTPUT_SHA256="${AUDIT_ROOT}/output_sha256.tsv"
PDF_TEXT_VALIDATION="${AUDIT_ROOT}/pdf_text_validation.tsv"
TASK_TMP_DIR=""
RUN_STATUS="INITIALIZING"

mkdir -p "${AUDIT_ROOT}" "${LOG_ROOT}" "${LOCK_ROOT}"
if ! mkdir "${LOCK_DIR}" 2>/dev/null; then
  echo "Another Figure 7 full-range run owns ${LOCK_DIR}." >&2
  exit 2
fi
printf 'host=%s\npid=%s\nrun_id=%s\n' \
  "$(hostname -s)" "$$" "${RUN_ID}" > "${LOCK_DIR}/owner"

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
  rm -f -- "${LOCK_DIR}/owner"
  rmdir "${LOCK_DIR}" 2>/dev/null || true
  trap - EXIT
  exit "${exit_code}"
}
trap cleanup EXIT
exec > >(tee -a "${RUN_LOG}") 2>&1

inventory_figure6() {
  local output="$1"
  find "${ITERATION_ROOT}/data/Figures/Figure6" "${ITERATION_ROOT}/Figures" \
    -maxdepth 3 -type f \
    \( -iname '*fig6*' -o -iname '*figure6*' \) -print0 2>/dev/null |
    sort -z | xargs -0 sha256sum > "${output}"
}

echo "Figure 7 full-range run start: $(date -Iseconds)"
echo "run_id=${RUN_ID} host=$(hostname -s) n_core=${N_CORE} memory_kb=${MEMORY_KB}"
echo "model_code_root=${MODEL_CODE_ROOT}"
echo "joint_result_root=${JOINT_RESULT_ROOT}"
RUN_STATUS="PREFLIGHT"
write_status RUNNING 0 "${RUN_STATUS}"
inventory_figure6 "${FIGURE6_BEFORE}"

TASK_TMP_DIR="$(mktemp -d /tmp/figure7-full-range.XXXXXX)"
mkdir -p "${TASK_TMP_DIR}/home" "${TASK_TMP_DIR}/cache" \
  "${TASK_TMP_DIR}/fontconfig" "${TASK_TMP_DIR}/model_rcpp_cache"
printf '%s\n' \
  'options(bitmapType = "cairo", device = "png", warn = 1)' \
  > "${TASK_TMP_DIR}/Rprofile"

CONTAINER_ARGS=(
  exec --cleanenv --containall --pwd "${ITERATION_ROOT}"
  --home "${TASK_TMP_DIR}/home"
  --env "TMPDIR=${TASK_TMP_DIR}"
  --env "XDG_CACHE_HOME=${TASK_TMP_DIR}/cache"
  --env "FONTCONFIG_PATH=/etc/fonts"
  --env "R_PROFILE_USER=${TASK_TMP_DIR}/Rprofile"
  --env "R_HOME=/opt/R/4.4.2/lib64/R"
  --env "LD_LIBRARY_PATH=${CONTAINER_LD_LIBRARY_PATH}"
  --env "DISPLAY=" --env "QT_QPA_PLATFORM=offscreen" --env "MPLBACKEND=Agg"
  --env "FIGURE_WORKSPACE_ROOT=${ITERATION_ROOT}"
  --env "FIGURE_MODEL_CODE_ROOT=${MODEL_CODE_ROOT}"
  --env "FIGURE_INVIVO_RESULT_ROOT=${INVIVO_RESULT_ROOT}"
  --env "FIGURE_INVITRO_RESULT_ROOT=${INVITRO_RESULT_ROOT}"
  --env "FIGURE_JOINT_RESULT_ROOT=${JOINT_RESULT_ROOT}"
  --env "FIGURE7_FINITE_TIME_FUTURE_PLAN=multicore"
  --env "FIGURE7_REUSE_CONTINUOUS_RUN=${REUSE_CONTINUOUS_RUN}"
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
)
container_command() {
  "${CONTAINER_RUNTIME}" "${CONTAINER_ARGS[@]}" "${SIF_IMAGE}" "$@"
}

container_command Rscript --vanilla -e '
packages <- c(
  "Matrix", "Rcpp", "data.table", "future", "future.apply", "ggplot2",
  "ggh4x", "patchwork", "magick", "isoband", "scales"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
workspace <- normalizePath(Sys.getenv("FIGURE_WORKSPACE_ROOT"), mustWork = TRUE)
model_root <- normalizePath(Sys.getenv("FIGURE_MODEL_CODE_ROOT"), mustWork = TRUE)
code_root <- file.path(workspace, "Code", "Figures")
invisible(lapply(list.files(
  file.path(code_root, "util", "analysis"), pattern = "figure7.*[.]R$",
  full.names = TRUE
), parse))
source(file.path(code_root, "util", "analysis", "figure7_robustness.R"))
source(file.path(code_root, "util", "analysis", "figure7_context_extension.R"))
source(file.path(code_root, "util", "analysis", "figure7_finite_time_q10.R"))
source(file.path(code_root, "util", "analysis", "figure7_invitro_passage_q10.R"))
source(file.path(code_root, "util", "analysis", "figure7_full_range_q10.R"))
paths <- f7r_paths(workspace)
stopifnot(identical(paths$oxygen_code, model_root))
f7r_load_response_engine(paths)
f7g_load_propagator(paths)
stopifnot(
  identical(f7g_modes("in vivo"), "continuous"),
  identical(f7g_modes("in vitro"), c("continuous", "passage")),
  identical(f7g_days(), 0:10000),
  isTRUE(all.equal(f7g_o2("in vivo"), seq(0, 5, length.out = 201L))),
  isTRUE(all.equal(f7g_o2("in vitro"), seq(0, 20, by = 0.1)))
)
png_path <- file.path(tempdir(), "figure7_headless.png")
grDevices::png(png_path, width = 240, height = 240, type = "cairo", bg = "white")
graphics::plot.new(); grDevices::dev.off()
stopifnot(file.exists(png_path), file.info(png_path)$size > 0)
cat("figure7_full_range_preflight_ok\n")
'

if [[ "${PREFLIGHT_ONLY}" == TRUE ]]; then
  RUN_STATUS="COMPLETE"
  write_status COMPLETE 0 PREFLIGHT_ONLY
  exit 0
fi

if [[ "${DRAW_ONLY}" != TRUE ]]; then
  RUN_STATUS="COMPUTE_FULL_RANGE"
  write_status RUNNING 0 "${RUN_STATUS}"
  container_command Rscript --vanilla \
    "${CODE_ROOT}/data_Figure7_full_range_q10.R" \
    "--n-core=${N_CORE}" "--o2-chunk-size=${O2_CHUNK_SIZE}" \
    "--run-id=${RUN_ID}" --smoke=FALSE --publish-current=TRUE
fi

if [[ "${COMPUTE_ONLY}" == TRUE ]]; then
  inventory_figure6 "${FIGURE6_AFTER}"
  cmp -s "${FIGURE6_BEFORE}" "${FIGURE6_AFTER}"
  RUN_STATUS="COMPLETE"
  write_status COMPLETE 0 COMPUTE_ONLY
  exit 0
fi

RUN_STATUS="HEADLESS_RENDER"
write_status RUNNING 0 "${RUN_STATUS}"
LEGACY_ROOT="${AUDIT_ROOT}/legacy_pre_reorganization"
mkdir -p "${LEGACY_ROOT}"
LEGACY_BASENAMES=(
  supp_fig7-8_continuous_invitro_no_passaging
  supp_fig7-9_invivo_finite_time_full_grid
  supp_fig7-10_invitro_passage_full_grid
  supp_fig7-11_invitro_passage_1000d_extended_o2
  supp_fig7-12_continuous_invitro_1000d_extended_o2
)
for base in "${LEGACY_BASENAMES[@]}"; do
  while IFS= read -r path; do
    [[ -n "${path}" ]] || continue
    relative="${path#${ITERATION_ROOT}/}"
    destination="${LEGACY_ROOT}/${relative}"
    mkdir -p "$(dirname "${destination}")"
    mv "${path}" "${destination}"
  done < <(find \
    "${ITERATION_ROOT}/Figures" \
    "${ITERATION_ROOT}/manuscript/Figures" \
    "${ITERATION_ROOT}/data/Figures" \
    -type f -name "${base}*" -print 2>/dev/null)
done
for script in \
  draw_Supp_Figure7_8.R draw_Supp_Figure7_9.R \
  draw_Supp_Figure7_10.R draw_Supp_Figure7_11.R draw_Figure7.R; do
  container_command Rscript --vanilla "${CODE_ROOT}/${script}"
done

RUN_STATUS="VALIDATE"
write_status RUNNING 0 "${RUN_STATUS}"
CURRENT_POINTER="${ITERATION_ROOT}/data/Figures/Figure7/fixed_pmisseg_v1/finite_time_full_q10_current.tsv"
require_file "${CURRENT_POINTER}" "full-range current pointer"
CURRENT_RUN_ROOT="$(awk -F '\t' 'NR==2 {print $2}' "${CURRENT_POINTER}")"
CURRENT_RUN_ROOT="${ITERATION_ROOT}/data/Figures/Figure7/fixed_pmisseg_v1/${CURRENT_RUN_ROOT}"
require_dir "${CURRENT_RUN_ROOT}" "full-range current run"

VALIDATIONS=(
  "${CURRENT_RUN_ROOT}/full_range_task_validation.tsv"
  "${CURRENT_RUN_ROOT}/full_range_panel_validation.tsv"
  "${CURRENT_RUN_ROOT}/passage_vs_continuous_validation.tsv"
  "${CURRENT_RUN_ROOT}/figure7_full_range_render_validation.tsv"
  "${ITERATION_ROOT}/data/Figures/Supp_Figure7_8/supp_fig7-8_inverse_response_render_validation.tsv"
  "${ITERATION_ROOT}/data/Figures/Supp_Figure7_9/supp_fig7-9_invivo_continuous_full_range_render_validation.tsv"
  "${ITERATION_ROOT}/data/Figures/Supp_Figure7_10/supp_fig7-10_invitro_continuous_full_range_render_validation.tsv"
  "${ITERATION_ROOT}/data/Figures/Supp_Figure7_11/supp_fig7-11_invitro_passage_full_range_render_validation.tsv"
)
for path in "${VALIDATIONS[@]}"; do
  require_file "${path}" "validation table"
  awk -F '\t' '
    NR == 1 {for (i=1; i<=NF; i++) if ($i=="passed") c=i; next}
    c>0 && toupper($c)!="TRUE" {bad=1}
    END {if (NR<2 || c==0 || bad) exit 1}
  ' "${path}" || { echo "Validation failed: ${path}" >&2; exit 1; }
done

if find "${CURRENT_RUN_ROOT}" -maxdepth 1 -type f \
    -name 'full_range_panel_invivo_passage_*.rds' | grep -q .; then
  echo "Forbidden in-vivo passage panel was produced." >&2
  exit 1
fi

OUTPUTS=(
  "${ITERATION_ROOT}/Figures/assembled_fig7.pdf"
  "${ITERATION_ROOT}/Figures/assembled_fig7.png"
  "${ITERATION_ROOT}/Figures/supp_fig7-8_inverse_response.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig7-8_inverse_response.png"
  "${ITERATION_ROOT}/Figures/supp_fig7-9_invivo_continuous_full_range.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig7-9_invivo_continuous_full_range.png"
  "${ITERATION_ROOT}/Figures/supp_fig7-10_invitro_continuous_full_range.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig7-10_invitro_continuous_full_range.png"
  "${ITERATION_ROOT}/Figures/supp_fig7-11_invitro_passage_full_range.pdf"
  "${ITERATION_ROOT}/Figures/supp_fig7-11_invitro_passage_full_range.png"
)
{
  printf 'sha256\tsize_bytes\tpath\n'
  for path in "${OUTPUTS[@]}"; do
    require_file "${path}" "Figure 7 output"
    printf '%s\t%s\t%s\n' \
      "$(sha256sum "${path}" | awk '{print $1}')" \
      "$(stat -c '%s' "${path}")" "${path}"
  done
} > "${OUTPUT_SHA256}"

{
  printf 'pdf\ttext_words_present\tfont_count\tpassed\n'
  for pdf in "${OUTPUTS[@]}"; do
    [[ "${pdf}" == *.pdf ]] || continue
    text_path="${TASK_TMP_DIR}/$(basename "${pdf}").txt"
    pdftotext -layout "${pdf}" "${text_path}"
    words_present=FALSE
    if grep -Eiq 'oxygen|ploidy|misseg' "${text_path}"; then
      words_present=TRUE
    fi
    font_count="$(pdffonts "${pdf}" | awk 'NR>2 {n++} END {print n+0}')"
    passed=FALSE
    if [[ "${words_present}" == TRUE && "${font_count}" -gt 0 ]]; then
      passed=TRUE
    fi
    printf '%s\t%s\t%s\t%s\n' \
      "${pdf}" "${words_present}" "${font_count}" "${passed}"
    [[ "${passed}" == TRUE ]] || exit 1
  done
} > "${PDF_TEXT_VALIDATION}"

inventory_figure6 "${FIGURE6_AFTER}"
cmp -s "${FIGURE6_BEFORE}" "${FIGURE6_AFTER}" || {
  echo "Figure 6 changed during the Figure 7 run." >&2
  diff -u "${FIGURE6_BEFORE}" "${FIGURE6_AFTER}" || true
  exit 1
}

RUN_STATUS="COMPLETE"
write_status COMPLETE 0 COMPLETE
echo "Figure 7 full-range run complete: $(date -Iseconds)"
echo "current_run_root=${CURRENT_RUN_ROOT}"
echo "output_manifest=${OUTPUT_SHA256}"
