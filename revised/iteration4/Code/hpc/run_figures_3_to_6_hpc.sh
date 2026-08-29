#!/usr/bin/env bash

# Headless Figure 3--6 recomputation and rendering for hpctpa3pc0028.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
ITERATION_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"
REPO_ROOT="$(cd -- "${ITERATION_ROOT}/../.." && pwd -P)"
CODE_ROOT="${ITERATION_ROOT}/Code/Figures"
EXPECTED_NODE="hpctpa3pc0028"
EXPECTED_REPO_ROOT="/share/lab_crd/taoli/Project/HypoxiaLTEEFigures"
SIF_IMAGE="/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif"
MODEL_CODE_ROOT="/share/lab_crd/taoli/Project/soft_couping_org/oxygen/code/O2_supply_demand_MAP"
RESULTS_ROOT="/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results"
INVIVO_RESULT_ROOT="${RESULTS_ROOT}/fit_invivo_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253"
INVITRO_RESULT_ROOT="${RESULTS_ROOT}/fit_invitro_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253"
JOINT_RESULT_ROOT="${RESULTS_ROOT}/fit_joint_unified_global_invitro_500seed_all_xxlarge_r442_exact_20260828_145253"
GEMCITABINE_DATA_ROOT="${REPO_ROOT}/data/InVivoData_Gemcitabine"
LTEE_DATA_ROOT="${REPO_ROOT}/data/InVitroData_LTEE"
RED_EASYBUILD_ROOT="/app/eb"
RED_FLEXIBLAS_LIB="${RED_EASYBUILD_ROOT}/software/FlexiBLAS/3.4.4-GCC-13.3.0/lib64"
RED_OPENBLAS_LIB="${RED_EASYBUILD_ROOT}/software/OpenBLAS/0.3.27-GCC-13.3.0/lib"
RED_GCC_LIB="${RED_EASYBUILD_ROOT}/software/GCCcore/13.3.0/lib64"
RED_BINUTILS_LIB="${RED_EASYBUILD_ROOT}/software/binutils/2.42-GCCcore-13.3.0/lib"
CONTAINER_LD_LIBRARY_PATH="${RED_FLEXIBLAS_LIB}:${RED_OPENBLAS_LIB}:${RED_GCC_LIB}:${RED_BINUTILS_LIB}:/opt/rh/gcc-toolset-13/root/usr/lib64"

N_CORE=56
BASELINE_ONLY=FALSE
for argument in "$@"; do
  case "${argument}" in
    --n-core=*) N_CORE="${argument#*=}" ;;
    --write-baseline-only) BASELINE_ONLY=TRUE ;;
    -h|--help)
      echo "Usage: run_figures_3_to_6_hpc.sh [--n-core=N] [--write-baseline-only]"
      exit 0
      ;;
    *) echo "Unknown option: ${argument}" >&2; exit 2 ;;
  esac
done
[[ "${N_CORE}" =~ ^[1-9][0-9]*$ ]] && (( N_CORE <= 63 )) || {
  echo "--n-core must be 1..63." >&2
  exit 2
}
[[ "$(hostname -s)" == "${EXPECTED_NODE}" ]] || {
  echo "Must run on ${EXPECTED_NODE}; observed $(hostname -s)." >&2
  exit 2
}
[[ "${REPO_ROOT}" == "${EXPECTED_REPO_ROOT}" ]] || {
  echo "Unexpected repository root: ${REPO_ROOT}" >&2
  exit 2
}
for path in \
  "${SIF_IMAGE}" \
  "${MODEL_CODE_ROOT}/model/model_O2_supply_demand_MAP.R" \
  "${INVIVO_RESULT_ROOT}" "${INVITRO_RESULT_ROOT}" "${JOINT_RESULT_ROOT}" \
  "${GEMCITABINE_DATA_ROOT}" "${LTEE_DATA_ROOT}"; do
  [[ -e "${path}" && -r "${path}" ]] || {
    echo "Missing required input: ${path}" >&2
    exit 2
  }
done

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v apptainer)"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v singularity)"
else
  echo "Apptainer/Singularity is unavailable." >&2
  exit 2
fi

RUN_ID="$(date '+%Y%m%d_%H%M%S')"
RUN_ROOT="${ITERATION_ROOT}/audit/hpc_figures_3_6/${RUN_ID}"
TASK_TMP_DIR="${RUN_ROOT}/tmp"
RUN_LOG="${RUN_ROOT}/run.log"
STATUS_PATH="${RUN_ROOT}/status.tsv"
mkdir -p "${TASK_TMP_DIR}/home" "${TASK_TMP_DIR}/cache" \
  "${TASK_TMP_DIR}/model_rcpp_cache"
printf '%s\n' 'options(bitmapType = "cairo")' 'Sys.unsetenv("DISPLAY")' \
  > "${TASK_TMP_DIR}/Rprofile"
exec > >(tee -a "${RUN_LOG}") 2>&1

status() {
  printf 'run_id\tstatus\tstage\thost\tn_core\tgit_head\tupdated_at\n' > "${STATUS_PATH}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${RUN_ID}" "$1" "$2" "$(hostname -s)" "${N_CORE}" \
    "$(git -C "${REPO_ROOT}" rev-parse HEAD)" "$(date -Iseconds)" \
    >> "${STATUS_PATH}"
}

CONTAINER_ARGS=(
  exec --cleanenv --containall --pwd "${ITERATION_ROOT}"
  --home "${TASK_TMP_DIR}/home"
  --env "TMPDIR=${TASK_TMP_DIR}"
  --env "XDG_CACHE_HOME=${TASK_TMP_DIR}/cache"
  --env "R_PROFILE_USER=${TASK_TMP_DIR}/Rprofile"
  --env "R_HOME=/opt/R/4.4.2/lib64/R"
  --env "LD_LIBRARY_PATH=${CONTAINER_LD_LIBRARY_PATH}"
  --env "FIGURE_WORKSPACE_ROOT=${ITERATION_ROOT}"
  --env "FIGURE_MODEL_CODE_ROOT=${MODEL_CODE_ROOT}"
  --env "FIGURE_INVIVO_RESULT_ROOT=${INVIVO_RESULT_ROOT}"
  --env "FIGURE_INVITRO_RESULT_ROOT=${INVITRO_RESULT_ROOT}"
  --env "FIGURE_JOINT_RESULT_ROOT=${JOINT_RESULT_ROOT}"
  --env "FIGURE_GEMCITABINE_DATA_ROOT=${GEMCITABINE_DATA_ROOT}"
  --env "FIGURE_LTEE_DATA_ROOT=${LTEE_DATA_ROOT}"
  --env MININGCLONEID_RCPP_REBUILD=FALSE
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
  --env MKL_NUM_THREADS=1 --env RCPP_PARALLEL_NUM_THREADS=1
  --env KMP_USE_SHM=0
  --bind "${RED_EASYBUILD_ROOT}:${RED_EASYBUILD_ROOT}:ro"
  --bind "${ITERATION_ROOT}:${ITERATION_ROOT}:rw"
  --bind "${MODEL_CODE_ROOT}:${MODEL_CODE_ROOT}:ro"
  --bind "${TASK_TMP_DIR}/model_rcpp_cache:${MODEL_CODE_ROOT}/model/.rcpp_cache_o2_supply_demand_map:rw"
  --bind "${RESULTS_ROOT}:${RESULTS_ROOT}:ro"
  --bind "${GEMCITABINE_DATA_ROOT}:${GEMCITABINE_DATA_ROOT}:ro"
  --bind "${LTEE_DATA_ROOT}:${LTEE_DATA_ROOT}:ro"
  --bind "${TASK_TMP_DIR}:${TASK_TMP_DIR}:rw"
)
container_command() {
  "${CONTAINER_RUNTIME}" "${CONTAINER_ARGS[@]}" "${SIF_IMAGE}" "$@"
}

RUNTIME_ARGS=(
  "--invitro-result-root=${INVITRO_RESULT_ROOT}"
  "--invivo-result-root=${INVIVO_RESULT_ROOT}"
  "--joint-result-root=${JOINT_RESULT_ROOT}"
  "--gemcitabine-data-root=${GEMCITABINE_DATA_ROOT}"
  "--ltee-data-root=${LTEE_DATA_ROOT}"
)

echo "run_id=${RUN_ID}"
echo "git_head=$(git -C "${REPO_ROOT}" rev-parse HEAD)"
echo "model_code_root=${MODEL_CODE_ROOT}"
echo "invivo_result_root=${INVIVO_RESULT_ROOT}"
echo "invitro_result_root=${INVITRO_RESULT_ROOT}"
echo "joint_result_root=${JOINT_RESULT_ROOT}"

status RUNNING HEADLESS_PREFLIGHT
container_command Rscript -e '
stopifnot(identical(getOption("bitmapType"), "cairo"), !nzchar(Sys.getenv("DISPLAY")))
png(file.path(Sys.getenv("TMPDIR"), "smoke.png"), type = "cairo"); plot(1:3); dev.off()
cairo_pdf(file.path(Sys.getenv("TMPDIR"), "smoke.pdf")); plot(1:3); dev.off()
stopifnot(file.info(file.path(Sys.getenv("TMPDIR"), "smoke.png"))$size > 0)
cat("headless_preflight_ok\n")
'

if [[ "${BASELINE_ONLY}" == "TRUE" ]]; then
  status RUNNING WRITE_INPUT_BASELINE
  container_command python3 "${CODE_ROOT}/util/workflow/input_validator.py" \
    --write-baseline
  status COMPLETE WRITE_INPUT_BASELINE
  echo "baseline_path=${ITERATION_ROOT}/Code/config/manifests/expected_scientific_input_md5.tsv"
  exit 0
fi

status RUNNING RUN_FIGURES_3_TO_6
if ! container_command bash "${ITERATION_ROOT}/manager.sh" \
  "${RUNTIME_ARGS[@]}" "--n-core=${N_CORE}" --first-main-figure=3 \
  --recompute-fixed-o2 --recompute-invivo-tsne --rebuild-figure6-grid; then
  status FAILED RUN_FIGURES_3_TO_6
  exit 1
fi

status RUNNING VERIFY_RENDERED_PAIRS
for basename in \
  assembled_fig3 assembled_fig4 assembled_fig5 assembled_fig6 \
  supp_fig4-1_all18_cluster_prior_violins \
  supp_fig4-2_invivo_optimizer_diagnostics \
  supp_fig5-1_joint_parameter_stability \
  supp_fig5-2_joint_fit_optimizer_diagnostics \
  supp_fig6-1_response_class_diagnostics \
  supp_fig6-2_joint_ensemble_robustness \
  supp_fig6-3_weak_gap_regime_robustness \
  supp_fig6-4_extended_invitro_o2_response; do
  for extension in png pdf; do
    path="${ITERATION_ROOT}/Figures/${basename}.${extension}"
    [[ -s "${path}" ]] || {
      status FAILED VERIFY_RENDERED_PAIRS
      echo "Missing rendered artifact: ${path}" >&2
      exit 1
    }
  done
done
status COMPLETE VERIFY_RENDERED_PAIRS
echo "status_path=${STATUS_PATH}"
echo "Figure 3-6 HPC recomputation and rendering complete."
