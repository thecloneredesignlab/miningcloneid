#!/usr/bin/env bash

# Compute only the data products required by Figure 6 and Supplementary
# Figures 6-1 through 6-4 on the allocated RED compute node. Drawing and
# manuscript rendering remain local tasks.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
ITERATION_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"
REPO_ROOT="$(cd -- "${ITERATION_ROOT}/../.." && pwd -P)"
CODE_ROOT="${ITERATION_ROOT}/Code/Figures"
AUDIT_ROOT="${ITERATION_ROOT}/audit"

EXPECTED_REPO_ROOT="/share/lab_crd/taoli/Project/HypoxiaLTEEFigures"
EXPECTED_NODE="hpctpa3pc0028"
SIF_IMAGE="/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif"
MODEL_CODE_ROOT="/share/lab_crd/taoli/Project/soft_couping_org/oxygen/code/O2_supply_demand_MAP"
RESULTS_ROOT="/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results"
INVIVO_RESULT_ROOT="${RESULTS_ROOT}/fit_invivo_unified_500seed_r442_exact_20260825_032031"
INVITRO_RESULT_ROOT="${RESULTS_ROOT}/fit_invitro_unified_500seed_r442_exact_20260825_032031"
JOINT_RESULT_ROOT="${RESULTS_ROOT}/fit_joint_invivo_clusters_global_invitro_best_500seed_r442_exact_20260826_033633"
GEMCITABINE_DATA_ROOT="${REPO_ROOT}/data/InVivoData_Gemcitabine"
LTEE_DATA_ROOT="${REPO_ROOT}/data/InVitroData_LTEE"
N_CORE=63
PREFLIGHT_ONLY=FALSE
RED_EASYBUILD_ROOT="/app/eb"
RED_FLEXIBLAS_LIB="${RED_EASYBUILD_ROOT}/software/FlexiBLAS/3.4.4-GCC-13.3.0/lib64"
RED_OPENBLAS_LIB="${RED_EASYBUILD_ROOT}/software/OpenBLAS/0.3.27-GCC-13.3.0/lib"
RED_GCC_LIB="${RED_EASYBUILD_ROOT}/software/GCCcore/13.3.0/lib64"
RED_BINUTILS_LIB="${RED_EASYBUILD_ROOT}/software/binutils/2.42-GCCcore-13.3.0/lib"
CONTAINER_LD_LIBRARY_PATH="${RED_FLEXIBLAS_LIB}:${RED_OPENBLAS_LIB}:${RED_GCC_LIB}:${RED_BINUTILS_LIB}:/opt/rh/gcc-toolset-13/root/usr/lib64"

usage() {
  printf '%s\n' \
    'Usage:' \
    '  run_figure6_data_hpc.sh [--n-core=N] [--preflight-only]' \
    '' \
    'This data-only runner is specific to the active RED allocation on' \
    'hpctpa3pc0028. It uses the r442 exact SIF and the synchronized current model' \
    'and fit-result roots under /share/lab_crd/taoli/Project.' \
    '' \
    'Options:' \
    '  --n-core=N       Independent R worker count. Default: 63.' \
    '  --preflight-only Validate node, paths, container, packages, and staged' \
    '                   Figure 3/4/5 inputs without archiving or computing.' \
    '  -h, --help       Show this help text.' \
    '' \
    'The script archives prior Figure 6 data products inside iteration4, performs' \
    'a fresh rebuild, and runs only:' \
    '  data_Figure6.R' \
    '  data_Supp_Figure6_1.R' \
    '  data_Supp_Figure6_2.R' \
    '  data_Supp_Figure6_3.R' \
    '  data_Supp_Figure6_4.R' \
    '' \
    'No draw_*.R script, figure renderer, manuscript builder, or LaTeX command is' \
    'executed.'
}

for argument in "$@"; do
  case "${argument}" in
    --n-core=*)
      N_CORE="${argument#*=}"
      ;;
    --preflight-only)
      PREFLIGHT_ONLY=TRUE
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: ${argument}" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if ! [[ "${N_CORE}" =~ ^[1-9][0-9]*$ ]] || (( N_CORE > 63 )); then
  echo "--n-core must be an integer from 1 through 63." >&2
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
require_directory "${GEMCITABINE_DATA_ROOT}" "in-vivo experimental-data root"
require_directory "${LTEE_DATA_ROOT}" "in-vitro experimental-data root"
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
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_fixed_o2_utils.R"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_fixed_o2_mode_utils.R"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_fixed_o2_format_utils.R"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_fixed_o2_validation_utils.R"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_postfit_input_utils.R"
  "${MODEL_CODE_ROOT}/util/o2_supply_demand_map_postfit_probe_utils.R"
  "${MODEL_CODE_ROOT}/simulation/fix_o2_simulation.R"
  "${MODEL_CODE_ROOT}/simulation/fix_o2_simulation_shared_utils.R"
  "${MODEL_CODE_ROOT}/simulation/o2/fixed_o2/run_fixed_o2_simulation.R"
)
for path in "${MODEL_SOURCE_FILES[@]}"; do
  require_file "${path}" "external model source"
done

WARM_START_LABELS=(
  tsne_vi_seed392_C01_vt_seed228
  tsne_vi_seed338_C02_vt_seed228
  tsne_vi_seed25_C03_vt_seed228
  tsne_vi_seed366_C04_vt_seed228
  tsne_vi_seed77_C05_vt_seed228
  tsne_vi_seed119_C06_vt_seed228
)

UPSTREAM_FILES=(
  "${ITERATION_ROOT}/data/Figures/Figure3/figure3g_parameter_endpoints_500seeds.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure4/fixed_o2_dominant_ploidy_201grid.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure4/invivo_fit_objective_ranking_500seeds.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure4/all_parameter_fitted_endpoint_values.tsv"
  "${ITERATION_ROOT}/data/Figures/Figure5/pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_coordinates.csv"
  "${ITERATION_ROOT}/data/Figures/Figure5/figure5_frozen_inputs/selected_results.tsv"
  "${ITERATION_ROOT}/data/Figures/Supp_Figure5_1/soft_coupling_master_long.tsv"
)
for label in "${WARM_START_LABELS[@]}"; do
  UPSTREAM_FILES+=(
    "${ITERATION_ROOT}/data/Figures/Figure5/figure5_frozen_inputs/winners/${label}/best_params.tsv"
    "${ITERATION_ROOT}/data/Figures/Figure5/figure5_frozen_inputs/winners/${label}/fit_config.rds"
    "${ITERATION_ROOT}/data/Figures/Supp_Figure5_2/joint/${label}/seed_objective_simple.tsv"
  )
done
for path in "${UPSTREAM_FILES[@]}"; do
  require_file "${path}" "fresh staged Figure 3/4/5 input"
done

if command -v apptainer >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v apptainer)"
elif command -v singularity >/dev/null 2>&1; then
  CONTAINER_RUNTIME="$(command -v singularity)"
else
  echo "Neither apptainer nor singularity is available on the compute node." >&2
  exit 2
fi

AVAILABLE_CPU="$(getconf _NPROCESSORS_ONLN)"
if (( AVAILABLE_CPU < N_CORE + 1 )); then
  echo "Need at least $((N_CORE + 1)) online CPUs; observed ${AVAILABLE_CPU}." >&2
  exit 2
fi

RUN_ID="$(date '+%Y%m%d_%H%M%S')"
RUN_ROOT="${AUDIT_ROOT}/hpc_figure6_data/${RUN_ID}"
LOG_ROOT="${AUDIT_ROOT}/logs"
LOCK_ROOT="${AUDIT_ROOT}/locks"
LOCK_DIR="${LOCK_ROOT}/figure6_data_hpc.lock"
RUN_LOG="${LOG_ROOT}/figure6_data_hpc_${RUN_ID}.log"
LATEST_LOG="${LOG_ROOT}/figure6_data_hpc_latest.log"
STATUS_PATH="${RUN_ROOT}/status.tsv"
UPSTREAM_MD5_PATH="${RUN_ROOT}/upstream_input_md5.tsv"
MODEL_MD5_PATH="${RUN_ROOT}/model_source_md5.tsv"
OUTPUT_SUMMARY_PATH="${RUN_ROOT}/output_summary.tsv"
TASK_TMP_DIR=""
RUN_STATUS="INITIALIZING"
LOCK_ACQUIRED=FALSE

mkdir -p "${RUN_ROOT}" "${LOG_ROOT}" "${LOCK_ROOT}"
if ! mkdir "${LOCK_DIR}" 2>/dev/null; then
  echo "Another Figure 6 HPC data run owns the lock: ${LOCK_DIR}" >&2
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
  if [[ "${LOCK_ACQUIRED}" == "TRUE" ]]; then
    rm -f -- "${LOCK_DIR}/owner"
    rmdir "${LOCK_DIR}" 2>/dev/null || true
  fi
  trap - EXIT
  exit "${exit_code}"
}
trap cleanup EXIT

exec > >(tee -a "${RUN_LOG}") 2>&1

echo "Figure 6 HPC data run start: $(date -Iseconds)"
echo "run_id=${RUN_ID}"
echo "host=${HOST_SHORT}"
echo "iteration_root=${ITERATION_ROOT}"
echo "container_runtime=${CONTAINER_RUNTIME}"
echo "sif_image=${SIF_IMAGE}"
echo "model_code_root=${MODEL_CODE_ROOT}"
echo "invivo_result_root=${INVIVO_RESULT_ROOT}"
echo "invitro_result_root=${INVITRO_RESULT_ROOT}"
echo "joint_result_root=${JOINT_RESULT_ROOT}"
echo "n_core=${N_CORE}"
echo "git_head=$(git -C "${REPO_ROOT}" rev-parse HEAD)"

RUN_STATUS="PREFLIGHT"
write_status "RUNNING" "0" "${RUN_STATUS}"

{
  printf 'md5\tsize_bytes\tpath\n'
  for path in "${UPSTREAM_FILES[@]}"; do
    printf '%s\t%s\t%s\n' \
      "$(md5sum "${path}" | awk '{print $1}')" \
      "$(stat -c '%s' "${path}")" "${path}"
  done
} > "${UPSTREAM_MD5_PATH}"

{
  printf 'md5\tsize_bytes\tpath\n'
  for path in "${MODEL_SOURCE_FILES[@]}"; do
    printf '%s\t%s\t%s\n' \
      "$(md5sum "${path}" | awk '{print $1}')" \
      "$(stat -c '%s' "${path}")" "${path}"
  done
} > "${MODEL_MD5_PATH}"

SIF_SHA256="$(sha256sum "${SIF_IMAGE}" | awk '{print $1}')"
echo "sif_sha256=${SIF_SHA256}"

TASK_TMP_DIR="$(mktemp -d /tmp/figure6-data-hpc.XXXXXX)"
mkdir -p \
  "${TASK_TMP_DIR}/home" \
  "${TASK_TMP_DIR}/cache" \
  "${TASK_TMP_DIR}/model_rcpp_cache"
printf '%s\n' 'options(bitmapType = "cairo")' > "${TASK_TMP_DIR}/Rprofile"

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
  --env "FIGURE_WORKSPACE_ROOT=${ITERATION_ROOT}"
  --env "FIGURE_MODEL_CODE_ROOT=${MODEL_CODE_ROOT}"
  --env "FIGURE_INVIVO_RESULT_ROOT=${INVIVO_RESULT_ROOT}"
  --env "FIGURE_INVITRO_RESULT_ROOT=${INVITRO_RESULT_ROOT}"
  --env "FIGURE_JOINT_RESULT_ROOT=${JOINT_RESULT_ROOT}"
  --env "FIGURE_GEMCITABINE_DATA_ROOT=${GEMCITABINE_DATA_ROOT}"
  --env "FIGURE_LTEE_DATA_ROOT=${LTEE_DATA_ROOT}"
  --env MININGCLONEID_RCPP_REBUILD=FALSE
  --env MININGCLONEID_RCPP_LOCK_TIMEOUT_SEC=600
  --env OMP_NUM_THREADS=1
  --env OPENBLAS_NUM_THREADS=1
  --env MKL_NUM_THREADS=1
  --env VECLIB_MAXIMUM_THREADS=1
  --env RCPP_PARALLEL_NUM_THREADS=1
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

container_command Rscript --vanilla -e '
packages <- c(
  "Matrix", "Rcpp", "RcppEigen", "cluster", "data.table", "dplyr",
  "future", "future.apply", "ggplot2", "matrixStats", "tidyr"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
model_root <- normalizePath(Sys.getenv("FIGURE_MODEL_CODE_ROOT"), mustWork = TRUE)
stopifnot(file.exists(file.path(model_root, "model", "model_O2_supply_demand_MAP.R")))
workspace <- normalizePath(Sys.getenv("FIGURE_WORKSPACE_ROOT"), mustWork = TRUE)
source(file.path(workspace, "Code", "Figures", "util", "analysis", "figure6_robustness.R"))
paths <- f6r_paths(workspace)
stopifnot(identical(paths$oxygen_code, model_root))
f6r_load_response_engine(paths)
stopifnot(
  exists("cpp_o2simps_build_G_for_o2_triplet", envir = globalenv()),
  exists("fixo2_dominant_attractor_one", envir = globalenv())
)
cat("container_preflight_ok\n")
'

if [[ "${PREFLIGHT_ONLY}" == "TRUE" ]]; then
  RUN_STATUS="COMPLETE"
  write_status "COMPLETE" "0" "PREFLIGHT_ONLY"
  echo "Figure 6 HPC data preflight complete: $(date -Iseconds)"
  exit 0
fi

RUN_STATUS="ARCHIVE_PRIOR_OUTPUTS"
write_status "RUNNING" "0" "${RUN_STATUS}"
FRESH_ARCHIVE="${AUDIT_ROOT}/pre_hpc_figure6_data_${RUN_ID}"
TARGET_RELATIVE_PATHS=(
  data/Figures/Figure6
  data/Figures/Supp_Figure6_1
  data/Figures/Supp_Figure6_2
  data/Figures/Supp_Figure6_3
  data/Figures/Supp_Figure6_4
  manuscript/tables/data/fixed_o2
)
for relative_path in "${TARGET_RELATIVE_PATHS[@]}"; do
  source_path="${ITERATION_ROOT}/${relative_path}"
  if [[ -e "${source_path}" ]]; then
    destination="${FRESH_ARCHIVE}/${relative_path}"
    mkdir -p "$(dirname "${destination}")"
    mv "${source_path}" "${destination}"
  fi
done
mkdir -p "${ITERATION_ROOT}/data/Figures"
echo "fresh_archive=${FRESH_ARCHIVE}"

RUN_STATUS="DATA_FIGURE6"
write_status "RUNNING" "0" "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/data_Figure6.R" \
  "--n-core=${N_CORE}" --rebuild=TRUE --n-resample=100

RUN_STATUS="DATA_SUPP_FIGURE6_1"
write_status "RUNNING" "0" "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure6_1.R" \
  "--n-core=${N_CORE}" --rebuild=TRUE

RUN_STATUS="DATA_SUPP_FIGURE6_2"
write_status "RUNNING" "0" "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure6_2.R"

RUN_STATUS="DATA_SUPP_FIGURE6_3"
write_status "RUNNING" "0" "${RUN_STATUS}"
# The q10/q20 endpoint caches were freshly created by data_Figure6.R in this
# same run. Reuse only those same-run products to avoid a duplicate rebuild.
container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure6_3.R" \
  "--n-core=${N_CORE}" --rebuild=FALSE

RUN_STATUS="DATA_SUPP_FIGURE6_4"
write_status "RUNNING" "0" "${RUN_STATUS}"
container_command Rscript --vanilla "${CODE_ROOT}/data_Supp_Figure6_4.R" \
  "--n-core=${N_CORE}" --rebuild=TRUE

RUN_STATUS="VALIDATE"
write_status "RUNNING" "0" "${RUN_STATUS}"
FIGURE6_DIR="${ITERATION_ROOT}/data/Figures/Figure6"

count_rds() {
  local directory="$1"
  find "${directory}" -type f -name '*.rds' -print | wc -l | tr -d ' '
}

INVIVO_DENSE_COUNT="$(count_rds "${FIGURE6_DIR}/figure6d_dense_endpoint_cache")"
INVITRO_DENSE_COUNT="$(count_rds "${FIGURE6_DIR}/figure6_invitro_dense_endpoint_cache")"
INVIVO_Q20_COUNT="$(count_rds "${FIGURE6_DIR}/multiseed_seed_cache")"
INVITRO_Q20_COUNT="$(count_rds "${FIGURE6_DIR}/multiseed_endpoint_cache_invitro")"
INVITRO_SEPARATE_COUNT="$(count_rds "${FIGURE6_DIR}/separate_invitro_seed_cache")"

[[ "${INVIVO_DENSE_COUNT}" == "186" ]] || {
  echo "Expected 186 in-vivo dense endpoint caches; observed ${INVIVO_DENSE_COUNT}." >&2
  exit 1
}
[[ "${INVITRO_DENSE_COUNT}" == "186" ]] || {
  echo "Expected 186 in-vitro dense endpoint caches; observed ${INVITRO_DENSE_COUNT}." >&2
  exit 1
}
[[ "${INVIVO_Q20_COUNT}" == "600" ]] || {
  echo "Expected 600 in-vivo q20 seed caches; observed ${INVIVO_Q20_COUNT}." >&2
  exit 1
}
[[ "${INVITRO_Q20_COUNT}" == "372" ]] || {
  echo "Expected 372 unique in-vitro q20 endpoint caches; observed ${INVITRO_Q20_COUNT}." >&2
  exit 1
}
[[ "${INVITRO_SEPARATE_COUNT}" == "500" ]] || {
  echo "Expected 500 separate in-vitro seed caches; observed ${INVITRO_SEPARATE_COUNT}." >&2
  exit 1
}

check_validation() {
  local path="$1"
  require_file "${path}" "validation table"
  awk -F '\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) if ($i == "passed") passed_col = i
      next
    }
    passed_col > 0 && toupper($passed_col) != "TRUE" { failed = 1 }
    END {
      if (NR < 2 || passed_col == 0 || failed) exit 1
    }
  ' "${path}" || {
    echo "Validation table contains a failed or missing check: ${path}" >&2
    exit 1
  }
}

check_validation "${FIGURE6_DIR}/response_class_validation.tsv"
check_validation "${FIGURE6_DIR}/response_class_invitro_validation.tsv"
check_validation "${FIGURE6_DIR}/figure6d_dense_validation.tsv"
check_validation "${FIGURE6_DIR}/figure6_invitro_dense_validation.tsv"
check_validation "${ITERATION_ROOT}/data/Figures/Supp_Figure6_3/supp_figure6-3_data_validation.tsv"
check_validation "${ITERATION_ROOT}/data/Figures/Supp_Figure6_3/supp_figure6-3_context_validation.tsv"
check_validation "${ITERATION_ROOT}/data/Figures/Supp_Figure6_4/supp_figure6-4_validation.tsv"

{
  printf 'metric\tvalue\n'
  printf 'invivo_dense_endpoint_cache_count\t%s\n' "${INVIVO_DENSE_COUNT}"
  printf 'invitro_dense_endpoint_cache_count\t%s\n' "${INVITRO_DENSE_COUNT}"
  printf 'invivo_q20_seed_cache_count\t%s\n' "${INVIVO_Q20_COUNT}"
  printf 'invitro_q20_unique_endpoint_cache_count\t%s\n' "${INVITRO_Q20_COUNT}"
  printf 'separate_invitro_seed_cache_count\t%s\n' "${INVITRO_SEPARATE_COUNT}"
  printf 'sif_sha256\t%s\n' "${SIF_SHA256}"
  printf 'model_code_root\t%s\n' "${MODEL_CODE_ROOT}"
  printf 'git_head\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse HEAD)"
} > "${OUTPUT_SUMMARY_PATH}"

RUN_STATUS="COMPLETE"
write_status "COMPLETE" "0" "${RUN_STATUS}"
echo "Figure 6 HPC data run complete: $(date -Iseconds)"
echo "status_path=${STATUS_PATH}"
echo "output_summary=${OUTPUT_SUMMARY_PATH}"
