#!/usr/bin/env bash

# Rebuild Figure 4B and Figure 4B lite on hpctpa3pc0028 from staged
# iteration5 inputs. All generated files and logs remain inside iteration5.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
ITERATION_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"
REPO_ROOT="$(cd -- "${ITERATION_ROOT}/../.." && pwd -P)"
EXPECTED_NODE="hpctpa3pc0028"
EXPECTED_ITERATION_ROOT="/share/lab_crd/taoli/Project/HypoxiaLTEEFigures/revised/iteration5"
SIF_IMAGE="/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif"
RED_EASYBUILD_ROOT="/app/eb"
RED_FLEXIBLAS_LIB="${RED_EASYBUILD_ROOT}/software/FlexiBLAS/3.4.4-GCC-13.3.0/lib64"
RED_OPENBLAS_LIB="${RED_EASYBUILD_ROOT}/software/OpenBLAS/0.3.27-GCC-13.3.0/lib"
RED_GCC_LIB="${RED_EASYBUILD_ROOT}/software/GCCcore/13.3.0/lib64"
RED_BINUTILS_LIB="${RED_EASYBUILD_ROOT}/software/binutils/2.42-GCCcore-13.3.0/lib"
RED_LIBICONV_LIB="${RED_EASYBUILD_ROOT}/software/libiconv/1.17-GCCcore-13.3.0/lib"
RED_LIBICONV_SO="${RED_LIBICONV_LIB}/libiconv.so.2"
RED_ICU_LIB="${RED_EASYBUILD_ROOT}/software/ICU/75.1-GCCcore-13.3.0/lib"
CONTAINER_LD_LIBRARY_PATH="${RED_FLEXIBLAS_LIB}:${RED_OPENBLAS_LIB}:${RED_GCC_LIB}:${RED_BINUTILS_LIB}:${RED_LIBICONV_LIB}:${RED_ICU_LIB}:/opt/rh/gcc-toolset-13/root/usr/lib64"

[[ "$(hostname -s)" == "${EXPECTED_NODE}" ]] || {
  echo "Must run on ${EXPECTED_NODE}; observed $(hostname -s)." >&2
  exit 2
}
[[ "${ITERATION_ROOT}" == "${EXPECTED_ITERATION_ROOT}" ]] || {
  echo "Unexpected iteration root: ${ITERATION_ROOT}" >&2
  exit 2
}
for path in "${SIF_IMAGE}" "${RED_LIBICONV_SO}" \
  "${RED_ICU_LIB}/libicui18n.so.75"; do
  [[ -r "${path}" ]] || {
    echo "Missing required runtime input: ${path}" >&2
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

DATA_DIR="${ITERATION_ROOT}/data/Figures/Figure4"
PANEL_DIR="${DATA_DIR}/panels"
DELIVERABLE_DIR="${ITERATION_ROOT}/Figures"
ASSOCIATION_SCRIPT="${ITERATION_ROOT}/Code/Figures/util/analysis/figure4_continuous_ploidy_association.R"
LANDSCAPE_SCRIPT="${ITERATION_ROOT}/Code/Figures/util/analysis/parameter_landscape.R"
for path in "${DATA_DIR}" "${ASSOCIATION_SCRIPT}" "${LANDSCAPE_SCRIPT}"; do
  [[ -e "${path}" ]] || {
    echo "Missing iteration5 input: ${path}" >&2
    exit 2
  }
done

RUN_ID="$(date '+%Y%m%d_%H%M%S')"
RUN_ROOT="${ITERATION_ROOT}/audit/hpc_figure4b_peak_o2/${RUN_ID}"
TASK_TMP_DIR="${RUN_ROOT}/tmp"
RUN_LOG="${RUN_ROOT}/run.log"
STATUS_PATH="${RUN_ROOT}/status.tsv"
CHECKSUM_PATH="${RUN_ROOT}/output_sha256.tsv"
mkdir -p "${TASK_TMP_DIR}/home" "${TASK_TMP_DIR}/cache" \
  "${PANEL_DIR}" "${DELIVERABLE_DIR}"
printf '%s\n' 'options(bitmapType = "cairo")' 'Sys.unsetenv("DISPLAY")' \
  > "${TASK_TMP_DIR}/Rprofile"
exec > >(tee -a "${RUN_LOG}") 2>&1

CODE_COMMIT="$(git -C "${REPO_ROOT}" rev-parse origin/HypoxiaLTEEFigures)"
CURRENT_STAGE="START"
status() {
  printf 'run_id\tstatus\tstage\thost\tcode_commit\tupdated_at\n' > "${STATUS_PATH}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${RUN_ID}" "$1" "$2" "$(hostname -s)" "${CODE_COMMIT}" \
    "$(date -Iseconds)" >> "${STATUS_PATH}"
}
on_error() {
  local exit_code=$?
  status FAILED "${CURRENT_STAGE}"
  exit "${exit_code}"
}
trap on_error ERR

CONTAINER_ARGS=(
  exec --cleanenv --containall --pwd "${ITERATION_ROOT}"
  --home "${TASK_TMP_DIR}/home"
  --env TMPDIR=/tmp --env TMP=/tmp --env TEMP=/tmp
  --env "XDG_CACHE_HOME=${TASK_TMP_DIR}/cache"
  --env "R_PROFILE_USER=${TASK_TMP_DIR}/Rprofile"
  --env "R_HOME=/opt/R/4.4.2/lib64/R"
  --env "LD_LIBRARY_PATH=${CONTAINER_LD_LIBRARY_PATH}"
  --env "LD_PRELOAD=${RED_LIBICONV_SO}"
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1
  --env MKL_NUM_THREADS=1 --env RCPP_PARALLEL_NUM_THREADS=1
  --env KMP_USE_SHM=0
  --env "ANALYSIS_DATA_DIR=${DATA_DIR}"
  --env "PLOT_OUTPUT_DIR=${PANEL_DIR}"
  --env "DELIVERABLE_OUTPUT_DIR=${DELIVERABLE_DIR}"
  --bind "${RED_EASYBUILD_ROOT}:${RED_EASYBUILD_ROOT}:ro"
  --bind "${ITERATION_ROOT}:${ITERATION_ROOT}:rw"
  --bind "${TASK_TMP_DIR}:${TASK_TMP_DIR}:rw"
  --bind "${TASK_TMP_DIR}:/tmp:rw"
  --bind "${TASK_TMP_DIR}:/var/tmp:rw"
)
container_command() {
  "${CONTAINER_RUNTIME}" "${CONTAINER_ARGS[@]}" "${SIF_IMAGE}" "$@"
}

echo "run_id=${RUN_ID}"
echo "code_commit=${CODE_COMMIT}"
echo "iteration_root=${ITERATION_ROOT}"
echo "container=${SIF_IMAGE}"

CURRENT_STAGE="PARSE_R_SOURCES"
status RUNNING "${CURRENT_STAGE}"
container_command Rscript -e \
  "invisible(parse(file='${ASSOCIATION_SCRIPT}')); invisible(parse(file='${LANDSCAPE_SCRIPT}')); cat('parse_ok\\n')"

CURRENT_STAGE="DERIVE_PEAK_O2_RANKING"
status RUNNING "${CURRENT_STAGE}"
container_command Rscript "${ASSOCIATION_SCRIPT}" \
  "--data-dir=${DATA_DIR}"

CURRENT_STAGE="RENDER_FIGURE4B_LINEAR_AND_LOGX"
status RUNNING "${CURRENT_STAGE}"
container_command Rscript "${LANDSCAPE_SCRIPT}"

CURRENT_STAGE="VALIDATE_RANKING_AND_OUTPUTS"
status RUNNING "${CURRENT_STAGE}"
container_command Rscript -e '
suppressPackageStartupMessages(library(data.table))
data_dir <- Sys.getenv("ANALYSIS_DATA_DIR")
ranking <- fread(file.path(data_dir, "continuous_ploidy_parameter_ranking.tsv"))
ranking <- ranking[order(display_order)]
expected_levels <- c("High O2", "Medium O2", "Low O2")
expected_group <- fcase(
  ranking$O2_at_max_abs >= 0 & ranking$O2_at_max_abs < 1.5, "Low O2",
  ranking$O2_at_max_abs >= 1.5 & ranking$O2_at_max_abs < 3.5, "Medium O2",
  ranking$O2_at_max_abs >= 3.5 & ranking$O2_at_max_abs <= 5, "High O2",
  default = NA_character_
)
stopifnot(
  nrow(ranking) == 18L,
  identical(ranking$display_order, seq_len(18L)),
  identical(ranking$peak_o2_group, expected_group),
  !anyNA(ranking$peak_o2_group),
  !any(diff(ranking$peak_o2_group_order) < 0)
)
counts <- ranking[, .N, by = peak_o2_group]
stopifnot(identical(
  counts[match(expected_levels, peak_o2_group), N],
  c(3L, 6L, 9L)
))
within_group_fail <- ranking[, any(diff(max_abs_rho) > 1e-12),
                             by = peak_o2_group_order]$V1
expected_parameter_order <- ranking[
  order(peak_o2_group_order, -max_abs_rho, parameter_order),
  parameter
]
stopifnot(
  !any(within_group_fail),
  identical(ranking$parameter, expected_parameter_order)
)
validation <- fread(file.path(data_dir, "parameter_landscape_layout_validation.tsv"))
stopifnot(
  validation[metric == "parameter_sort_secondary", value] ==
    "descending maximum absolute Spearman rho within peak O2 group",
  validation[metric == "parameter_sort_tertiary", value] ==
    "configured parameter order for exact max-|rho| ties",
  validation[metric == "row_annotation_field", value] == "peak_o2_group",
  validation[metric == "effect_fill_field", value] == "peak_direction",
  validation[metric == "effect_positive_fill", value] == "#EF8A62",
  validation[metric == "effect_negative_fill", value] == "#67A9CF",
  validation[metric == "figure4b_logx_rendered", value] == "TRUE",
  validation[metric == "figure4b_logx_axis_field", value] == "O2_pct",
  validation[metric == "figure4b_logx_transform", value] == "pseudo-log10",
  as.numeric(validation[metric == "figure4b_logx_sigma", value]) == 0.025,
  validation[metric == "figure4b_logx_zero_retained", value] == "TRUE",
  validation[
    metric == "figure4b_lite_endpoint_distribution_rendered", value
  ] == "FALSE",
  as.numeric(validation[metric == "figure4b_lite_output_width_in", value]) == 12,
  as.numeric(validation[metric == "figure4b_lite_output_height_in", value]) == 9,
  abs(as.numeric(validation[
    metric == "figure4b_lite_output_aspect_ratio", value
  ]) - 4 / 3) < 1e-12,
  validation[metric == "figure4b_lite_logx_rendered", value] == "TRUE",
  validation[
    metric == "figure4b_lite_logx_endpoint_distribution_rendered", value
  ] == "FALSE",
  as.numeric(validation[
    metric == "figure4b_lite_logx_output_width_in", value
  ]) == 12,
  as.numeric(validation[
    metric == "figure4b_lite_logx_output_height_in", value
  ]) == 9,
  abs(as.numeric(validation[
    metric == "figure4b_lite_logx_output_aspect_ratio", value
  ]) - 4 / 3) < 1e-12
)
cat("ranking_linear_logx_and_lite_validation_ok\n")
'

outputs=(
  "${PANEL_DIR}/parameter_continuous_ploidy_landscape.png"
  "${PANEL_DIR}/parameter_continuous_ploidy_landscape.pdf"
  "${PANEL_DIR}/parameter_continuous_ploidy_landscape.svg"
  "${PANEL_DIR}/Figure4B_lite.png"
  "${PANEL_DIR}/Figure4B_lite.pdf"
  "${PANEL_DIR}/Figure4B_lite.svg"
  "${PANEL_DIR}/parameter_continuous_ploidy_landscape_logx.png"
  "${PANEL_DIR}/parameter_continuous_ploidy_landscape_logx.pdf"
  "${PANEL_DIR}/parameter_continuous_ploidy_landscape_logx.svg"
  "${PANEL_DIR}/Figure4B_lite_logx.png"
  "${PANEL_DIR}/Figure4B_lite_logx.pdf"
  "${PANEL_DIR}/Figure4B_lite_logx.svg"
  "${DELIVERABLE_DIR}/Figure4B.png"
  "${DELIVERABLE_DIR}/Figure4B.pdf"
  "${DELIVERABLE_DIR}/Figure4B.svg"
  "${DELIVERABLE_DIR}/Figure4B_lite.png"
  "${DELIVERABLE_DIR}/Figure4B_lite.pdf"
  "${DELIVERABLE_DIR}/Figure4B_lite.svg"
  "${DELIVERABLE_DIR}/Figure4B_logx.png"
  "${DELIVERABLE_DIR}/Figure4B_logx.pdf"
  "${DELIVERABLE_DIR}/Figure4B_logx.svg"
  "${DELIVERABLE_DIR}/Figure4B_lite_logx.png"
  "${DELIVERABLE_DIR}/Figure4B_lite_logx.pdf"
  "${DELIVERABLE_DIR}/Figure4B_lite_logx.svg"
)
for path in "${outputs[@]}"; do
  [[ -s "${path}" ]] || {
    echo "Missing rendered output: ${path}" >&2
    exit 1
  }
done
printf 'sha256\tpath\n' > "${CHECKSUM_PATH}"
for path in "${outputs[@]}"; do
  printf '%s\t%s\n' "$(sha256sum "${path}" | awk '{print $1}')" "${path}" \
    >> "${CHECKSUM_PATH}"
done

CURRENT_STAGE="COMPLETE"
status COMPLETE "${CURRENT_STAGE}"
trap - ERR
echo "status_path=${STATUS_PATH}"
echo "checksum_path=${CHECKSUM_PATH}"
echo "Figure 4B linear and pseudo-log-x renders complete."
