#!/usr/bin/env bash

# Recompute day-1000 trajectory-level net-growth conditioning for Supplementary
# Figures 6-20/21, then render their PDF/PNG files on the displayless node.
# All writable paths are under revised/iteration4.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
ITERATION_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"
REPO_ROOT="$(cd -- "${ITERATION_ROOT}/../.." && pwd -P)"
CODE_ROOT="${ITERATION_ROOT}/Code/Figures"
EXPECTED_NODE="hpctpa3pc0028"
EXPECTED_REPO_ROOT="/share/lab_crd/taoli/Project/HypoxiaLTEEFigures"
MODEL_CODE_ROOT="/share/lab_crd/taoli/Project/soft_couping_org/oxygen/code/O2_supply_demand_MAP"
RESULTS_ROOT="/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results"
SIF_IMAGE="/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif"
RED_EASYBUILD_ROOT="/app/eb"
RED_FLEXIBLAS_LIB="${RED_EASYBUILD_ROOT}/software/FlexiBLAS/3.4.4-GCC-13.3.0/lib64"
RED_OPENBLAS_LIB="${RED_EASYBUILD_ROOT}/software/OpenBLAS/0.3.27-GCC-13.3.0/lib"
RED_GCC_LIB="${RED_EASYBUILD_ROOT}/software/GCCcore/13.3.0/lib64"
RED_BINUTILS_LIB="${RED_EASYBUILD_ROOT}/software/binutils/2.42-GCCcore-13.3.0/lib"
CONTAINER_LD_LIBRARY_PATH="${RED_FLEXIBLAS_LIB}:${RED_OPENBLAS_LIB}:${RED_GCC_LIB}:${RED_BINUTILS_LIB}:/opt/rh/gcc-toolset-13/root/usr/lib64"

N_CORE=56
DRAW_ONLY=FALSE
PREFLIGHT_ONLY=FALSE
RUN_ID="supp6_20_21_$(date '+%Y%m%d_%H%M%S')"
for argument in "$@"; do
  case "${argument}" in
    --n-core=*) N_CORE="${argument#*=}" ;;
    --run-id=*) RUN_ID="${argument#*=}" ;;
    --draw-only) DRAW_ONLY=TRUE ;;
    --preflight-only) PREFLIGHT_ONLY=TRUE ;;
    -h|--help)
      echo "Usage: $0 [--n-core=1..60] [--run-id=ID] [--draw-only|--preflight-only]"
      exit 0 ;;
    *) echo "Unknown argument: ${argument}" >&2; exit 2 ;;
  esac
done
[[ "${N_CORE}" =~ ^[1-9][0-9]*$ ]] && (( N_CORE <= 60 )) || {
  echo "--n-core must be 1..60" >&2; exit 2;
}
[[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$ ]] || {
  echo "Invalid run ID" >&2; exit 2;
}
[[ "$(hostname -s)" == "${EXPECTED_NODE}" ]] || {
  echo "Use ${EXPECTED_NODE}; current node is $(hostname -s)" >&2; exit 2;
}
[[ "${REPO_ROOT}" == "${EXPECTED_REPO_ROOT}" ]] || {
  echo "Unexpected repository root: ${REPO_ROOT}" >&2; exit 2;
}
require_file() { [[ -r "$1" && -f "$1" ]] || { echo "Missing $2: $1" >&2; exit 2; }; }
require_dir() { [[ -r "$1" && -d "$1" ]] || { echo "Missing $2: $1" >&2; exit 2; }; }
require_file "${SIF_IMAGE}" "container image"
require_dir "${MODEL_CODE_ROOT}" "external model code"
require_dir "${RESULTS_ROOT}" "fit results"
require_file "${RED_FLEXIBLAS_LIB}/libflexiblas.so.3" "FlexiBLAS"
for file in data_Supp_Figure6_20_21.R draw_Supp_Figure6_20.R draw_Supp_Figure6_21.R; do
  require_file "${CODE_ROOT}/${file}" "figure entry point"
done
SOURCE_FIXED_ROOT="${ITERATION_ROOT}/data/Figures/Figure7/fixed_pmisseg_v1"
SOURCE_POINTER="${SOURCE_FIXED_ROOT}/finite_time_full_q10_current.tsv"
require_file "${SOURCE_POINTER}" "Figure 6B source pointer"
SOURCE_RELATIVE_RUN="$(awk -F '\t' 'NR==2 {print $2}' "${SOURCE_POINTER}")"
SOURCE_RUN_ROOT="${SOURCE_FIXED_ROOT}/${SOURCE_RELATIVE_RUN}"
require_dir "${SOURCE_RUN_ROOT}" "Figure 6B source run"
require_file "${ITERATION_ROOT}/data/Figures/Figure6/fixed_pmisseg_v1/net_growth_full_range_q10_current.tsv" "Figure 6B growth pointer"
if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v apptainer)"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v singularity)"
else
  echo "No Apptainer/Singularity runtime" >&2; exit 2
fi

AUDIT_ROOT="${ITERATION_ROOT}/audit/hpc_supp_figure6_20_21/${RUN_ID}"
LOCK_DIR="${ITERATION_ROOT}/audit/locks/supp_figure6_20_21.lock"
mkdir -p "${AUDIT_ROOT}" "${ITERATION_ROOT}/audit/locks" "${ITERATION_ROOT}/audit/tmp"
mkdir "${LOCK_DIR}" 2>/dev/null || {
  echo "Another 6-20/21 computation holds ${LOCK_DIR}" >&2; exit 2;
}
printf 'host=%s\npid=%s\n' "$(hostname -s)" "$$" > "${LOCK_DIR}/owner"
TASK_TMP_DIR="$(mktemp -d "${ITERATION_ROOT}/audit/tmp/supp6-20-21.XXXXXX")"
mkdir -p "${TASK_TMP_DIR}/home" "${TASK_TMP_DIR}/cache" "${TASK_TMP_DIR}/model_rcpp_cache"
printf '%s\n' 'options(bitmapType = "cairo", device = "png", warn = 1)' > "${TASK_TMP_DIR}/Rprofile"
STATUS_PATH="${AUDIT_ROOT}/status.tsv"
STAGE="PREFLIGHT"
write_status() {
  printf 'run_id\tstatus\tstage\thost\tn_core\tupdated_at\n%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${RUN_ID}" "$1" "${STAGE}" "$(hostname -s)" "${N_CORE}" "$(date -Iseconds)" \
    > "${STATUS_PATH}"
}
cleanup() {
  code="$?"
  if (( code != 0 )); then write_status FAILED; fi
  if [[ -f "${LOCK_DIR}/owner" ]]; then rm -f -- "${LOCK_DIR}/owner"; fi
  rmdir "${LOCK_DIR}" 2>/dev/null || true
}
trap cleanup EXIT
write_status RUNNING

CONTAINER_ARGS=(
  exec --cleanenv --containall --pwd "${ITERATION_ROOT}"
  --home "${TASK_TMP_DIR}/home"
  --env "TMPDIR=/tmp" --env "TMP=/tmp" --env "TEMP=/tmp"
  --env "XDG_CACHE_HOME=${TASK_TMP_DIR}/cache"
  --env "R_PROFILE_USER=${TASK_TMP_DIR}/Rprofile"
  --env "R_MAKEVARS_USER=${SCRIPT_DIR}/figure6_compiler_makevars"
  --env "R_HOME=/opt/R/4.4.2/lib64/R"
  --env "LD_LIBRARY_PATH=${CONTAINER_LD_LIBRARY_PATH}"
  --env "DISPLAY=" --env "QT_QPA_PLATFORM=offscreen" --env "MPLBACKEND=Agg"
  --env "FIGURE_WORKSPACE_ROOT=${ITERATION_ROOT}"
  --env "FIGURE_MODEL_CODE_ROOT=${MODEL_CODE_ROOT}"
  --env "FIGURE_INVIVO_RESULT_ROOT=${RESULTS_ROOT}/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253"
  --env "FIGURE_INVITRO_RESULT_ROOT=${RESULTS_ROOT}/fit_invitro_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253"
  --env "FIGURE_JOINT_RESULT_ROOT=${RESULTS_ROOT}/fit_joint_unified_global_invitro_500seed_all_xxlarge_r442_exact_20260828_145253"
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
files <- c("data_Supp_Figure6_20_21.R", "draw_Supp_Figure6_20.R",
 "draw_Supp_Figure6_21.R", "util/analysis/figure6_growth_permissive.R",
 "util/analysis/figure6_growth_permissive_layout.R")
base <- file.path(Sys.getenv("FIGURE_WORKSPACE_ROOT"), "Code", "Figures")
invisible(lapply(file.path(base, files), parse))
path <- file.path(tempdir(), "headless_test.png")
grDevices::png(path, width=200, height=200, type="cairo")
graphics::plot.new(); grDevices::dev.off()
stopifnot(file.exists(path), file.info(path)$size > 0)
cat("supp_figure6_20_21_headless_preflight_ok\n")
'
if [[ "${PREFLIGHT_ONLY}" == TRUE ]]; then
  STAGE="COMPLETE"; write_status COMPLETE; exit 0
fi

if [[ "${DRAW_ONLY}" != TRUE ]]; then
  STAGE="COMPUTE"; write_status RUNNING
  container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure6_20_21.R" \
    --phase=all "--n-core=${N_CORE}"
fi
STAGE="RENDER"; write_status RUNNING
container_command Rscript --vanilla "${CODE_ROOT}/draw_Supp_Figure6_20.R"
container_command Rscript --vanilla "${CODE_ROOT}/draw_Supp_Figure6_21.R"

STAGE="VALIDATE"; write_status RUNNING
VALIDATION="${ITERATION_ROOT}/data/Figures/Supp_Figure6_20_21_day1000_positive_growth_v6/day1000_positive_growth_validation.tsv"
require_file "${VALIDATION}" "Figure 6B replay validation"
awk -F '\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="passed") c=i; next}
 c>0 && toupper($c)!="TRUE" {bad=1}
 END {if(NR!=5 || c==0 || bad) exit 1}' "${VALIDATION}"
for label in c01 c02; do
  number=20
  [[ "${label}" == c02 ]] && number=21
  stem="supp_fig6-${number}_growth_permissive_ploidy_vs_oxygen_${label}"
  pdf="${ITERATION_ROOT}/Figures/${stem}.pdf"
  png="${ITERATION_ROOT}/Figures/${stem}.png"
  require_file "${pdf}" "supplement PDF"
  require_file "${png}" "supplement PNG"
  pdftotext "${pdf}" - | grep -q "Day-1000"
  pdffonts "${pdf}" | awk 'NR>2 && $2=="Type 3" {bad=1} END {exit bad}'
done
{
  printf 'sha256\tpath\n'
  for file in "${ITERATION_ROOT}/Figures/supp_fig6-20_growth_permissive_ploidy_vs_oxygen_c01.pdf" \
              "${ITERATION_ROOT}/Figures/supp_fig6-20_growth_permissive_ploidy_vs_oxygen_c01.png" \
              "${ITERATION_ROOT}/Figures/supp_fig6-21_growth_permissive_ploidy_vs_oxygen_c02.pdf" \
              "${ITERATION_ROOT}/Figures/supp_fig6-21_growth_permissive_ploidy_vs_oxygen_c02.png" \
              "${ITERATION_ROOT}/data/Figures/Supp_Figure6_20_21_day1000_positive_growth_v6/day1000_positive_growth_curves.tsv"; do
    printf '%s\t%s\n' "$(sha256sum "${file}" | awk '{print $1}')" "${file}"
  done
} > "${AUDIT_ROOT}/output_sha256.tsv"
STAGE="COMPLETE"; write_status COMPLETE
echo "Supplementary Figure 6-20/21 complete: ${RUN_ID}"
