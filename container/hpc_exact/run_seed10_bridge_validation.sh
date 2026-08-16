#!/usr/bin/env bash
# Run the canonical in-vitro seed10 fit inside the current SIF namespace while
# using the exact RED EasyBuild R stack. This is an intermediate validation;
# the resulting execution is not yet a standalone container runtime.

#SBATCH --job-name=o2ivt_s10_bridge
#SBATCH --qos=xxlarge
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=22
#SBATCH --mem=16G

set -euo pipefail

if [[ "$#" -ne 6 ]]; then
  echo "Usage: $0 PROJECT_ROOT RUN_ROOT SIF_PATH WORKER_DRIVER_TEMPLATE HEADLESS_PROFILE_TEMPLATE REFERENCE_SEED_DIR" >&2
  exit 64
fi

PROJECT_ROOT=$(realpath "$1")
RUN_ROOT=$(realpath -m "$2")
SIF_PATH=$(realpath "$3")
WORKER_DRIVER_TEMPLATE=$(realpath "$4")
HEADLESS_PROFILE_TEMPLATE=$(realpath "$5")
REFERENCE_SEED_DIR=$(realpath "$6")

HOST_R_LIBRARY=/home/4482173/R/x86_64-pc-linux-gnu-library/4.4
CONTAINER_R_LIBRARY=/opt/o2-host-r-library/4.4
HPC_RSCRIPT=/app/eb/software/R/4.4.2-gfbf-2024a/bin/Rscript
MODEL_CACHE_REL=oxygen/code/O2_supply_demand_MAP/model/.rcpp_cache_o2_supply_demand_map
MODEL_CACHE_HOST=${PROJECT_ROOT}/${MODEL_CACHE_REL}
ISOLATED_CACHE=${RUN_ROOT}/rcpp_cache
SEED_DIR=${RUN_ROOT}/seed10
LOG_DIR=${RUN_ROOT}/log
WORKER_DRIVER=${RUN_ROOT}/run_fit_invitro_worker_bindings.R
HEADLESS_PROFILE=${RUN_ROOT}/Rprofile_o2_headless.R

mkdir -p "$RUN_ROOT" "$SEED_DIR" "$LOG_DIR" "$ISOLATED_CACHE"
cp "$WORKER_DRIVER_TEMPLATE" "$WORKER_DRIVER"
cp "$HEADLESS_PROFILE_TEMPLATE" "$HEADLESS_PROFILE"

source /etc/profile.d/modules.sh
module use /app/eb/modules/all
module purge
ml R/4.4.2

HOST_MODULE_PATH=$PATH
HOST_MODULE_LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export LANG=C.UTF-8
export LC_ALL=C.UTF-8

bridge_exec() {
  env \
    APPTAINERENV_PATH="$HOST_MODULE_PATH" \
    APPTAINERENV_LD_LIBRARY_PATH="$HOST_MODULE_LD_LIBRARY_PATH" \
    APPTAINERENV_R_LIBS_USER="$CONTAINER_R_LIBRARY" \
    APPTAINERENV_O2SD_PROJECT_ROOT="$PROJECT_ROOT" \
    APPTAINERENV_OMP_NUM_THREADS=1 \
    APPTAINERENV_OPENBLAS_NUM_THREADS=1 \
    APPTAINERENV_MKL_NUM_THREADS=1 \
    APPTAINERENV_VECLIB_MAXIMUM_THREADS=1 \
    APPTAINERENV_LANG=C.UTF-8 \
    APPTAINERENV_LC_ALL=C.UTF-8 \
    APPTAINERENV_TZ=America/New_York \
    APPTAINERENV_MININGCLONEID_RCPP_REBUILD=TRUE \
    apptainer exec \
      --cleanenv \
      --containall \
      --pwd "$PROJECT_ROOT" \
      --bind /app/eb/software:/app/eb/software:ro \
      --bind "$HOST_R_LIBRARY:$CONTAINER_R_LIBRARY:ro" \
      --bind "$PROJECT_ROOT:$PROJECT_ROOT:rw" \
      --bind "$ISOLATED_CACHE:$MODEL_CACHE_HOST:rw" \
      "$SIF_PATH" \
      "$@"
}

echo "project_root=$PROJECT_ROOT"
echo "run_root=$RUN_ROOT"
echo "sif_path=$SIF_PATH"
sha256sum "$SIF_PATH"
echo "git_branch=$(git -C "$PROJECT_ROOT" branch --show-current)"
echo "git_head=$(git -C "$PROJECT_ROOT" rev-parse HEAD)"
echo "module_list_begin"
module -t list 2>&1
echo "module_list_end"

bridge_exec "$HPC_RSCRIPT" --vanilla -e '
args <- commandArgs(TRUE)
runtime_path <- args[[1L]]
packages_path <- args[[2L]]
ext <- extSoftVersion()
runtime <- data.frame(
  key = c(
    "r_version", "r_home", "lib_paths", "locale", "blas", "lapack",
    "rng_kind", "normal_kind", "sample_kind", "cc", "cxx17",
    "cxx17flags", "blas_libs"
  ),
  value = c(
    R.version.string,
    R.home(),
    paste(.libPaths(), collapse = ";"),
    paste(Sys.getlocale(), collapse = ";"),
    unname(ext["BLAS"]),
    paste(La_version(), collapse = "."),
    RNGkind(),
    system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CC"), stdout = TRUE),
    system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX17"), stdout = TRUE),
    system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX17FLAGS"), stdout = TRUE),
    system2(file.path(R.home("bin"), "R"), c("CMD", "config", "BLAS_LIBS"), stdout = TRUE)
  ),
  stringsAsFactors = FALSE
)
write.table(runtime, runtime_path, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), sub("runtime.tsv$", "sessionInfo.txt", runtime_path))
pkgs <- c("DEoptim", "dplyr", "ggplot2", "Matrix", "Rcpp", "RcppEigen", "tidyr")
rows <- lapply(pkgs, function(pkg) {
  path <- find.package(pkg)
  so <- list.files(file.path(path, "libs"), pattern = "\\.so$", full.names = TRUE, recursive = TRUE)
  data.frame(
    package = pkg,
    version = as.character(packageVersion(pkg)),
    built = as.character(packageDescription(pkg)$Built),
    library_path = path,
    native_files = paste(so, collapse = ";"),
    native_md5 = if (length(so)) paste(unname(tools::md5sum(so)), collapse = ";") else "",
    stringsAsFactors = FALSE
  )
})
write.table(do.call(rbind, rows), packages_path, sep = "\t", quote = FALSE, row.names = FALSE)
' "$RUN_ROOT/runtime.tsv" "$RUN_ROOT/runtime_packages.tsv"

bridge_exec "$HPC_RSCRIPT" --vanilla -e '
backend <- new.env(parent = globalenv())
original_command_args <- commandArgs
backend_path <- normalizePath(
  "oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_invitro_backend.R",
  mustWork = TRUE
)
assign(
  "commandArgs",
  function(trailingOnly = FALSE) if (isTRUE(trailingOnly)) character(0) else c("R", paste0("--file=", backend_path)),
  envir = .GlobalEnv
)
sys.source(backend_path, envir = backend, chdir = TRUE)
assign("commandArgs", original_command_args, envir = .GlobalEnv)
x <- backend$ivt_load_fit_objects_compat(
  backend$default_fit_objects_dir(must_exist = TRUE),
  backend$default_flow_density_path()
)
stopifnot(identical(x$death_enabled, FALSE), is.data.frame(x$death_data), nrow(x$death_data) == 0L)
write.table(
  data.frame(key = c("death_enabled", "death_rows"), value = c("FALSE", "0")),
  file = commandArgs(TRUE)[[1L]], sep = "\t", quote = FALSE, row.names = FALSE
)
' "$RUN_ROOT/death_preflight.tsv"

bridge_exec "$HPC_RSCRIPT" "$WORKER_DRIVER" \
  --config="${PROJECT_ROOT}/oxygen/config/O2_supply_demand.yaml" \
  --seed=10 \
  --itermax=500 \
  --itermax_max=500 \
  --de_reltol=1e-4 \
  --de_steptol=25 \
  --NP=80 \
  --n_cores=22 \
  --local_optim_maxit=200 \
  --dt=0.05 \
  --init_total_size=1e6 \
  --o2_upper_bound=21 \
  --out_dir="$SEED_DIR" \
  --parameter_table="${PROJECT_ROOT}/oxygen/data/O2_supply_demand/parameter_table_invitro_buffering.csv" \
  --fit_objects_dir="${PROJECT_ROOT}/oxygen/ploidyOxygen/data/fit_objects" \
  --flow_density_path="${PROJECT_ROOT}/oxygen/data/g0g1_ploidy_density_grid.csv" \
  --auto_viz=FALSE

APPTAINERENV_R_PROFILE_USER="$HEADLESS_PROFILE" \
  bridge_exec "$HPC_RSCRIPT" \
    "${PROJECT_ROOT}/oxygen/code/O2_supply_demand_MAP/runner/run_postfit_pipeline.R" \
    --fit_dir="$SEED_DIR" \
    --scope=invitro

test -s "$SEED_DIR/fit_summary.tsv"
test -s "$SEED_DIR/best_params.tsv"
test -s "$SEED_DIR/best_params_transformed.tsv"
test -s "$SEED_DIR/report/fit_report_seed10.html"

cmp "$REFERENCE_SEED_DIR/best_params.tsv" "$SEED_DIR/best_params.tsv"
cmp "$REFERENCE_SEED_DIR/best_params_transformed.tsv" "$SEED_DIR/best_params_transformed.tsv"

summary_value() {
  local summary_file=$1
  local metric=$2
  awk -F '\t' -v wanted="$metric" '$1 == wanted { print $2; exit }' "$summary_file"
}

{
  printf 'metric\treference\tbridge\texact_match\n'
  for metric in optimizer_deoptim_objective optimizer_local_objective objective_total; do
    reference_value=$(summary_value "$REFERENCE_SEED_DIR/fit_summary.tsv" "$metric")
    bridge_value=$(summary_value "$SEED_DIR/fit_summary.tsv" "$metric")
    exact_match=FALSE
    [[ "$reference_value" == "$bridge_value" ]] && exact_match=TRUE
    printf '%s\t%s\t%s\t%s\n' "$metric" "$reference_value" "$bridge_value" "$exact_match"
    [[ "$exact_match" == TRUE ]]
  done
} > "$RUN_ROOT/reference_comparison.tsv"

echo "bridge SIF seed10 fit, exact comparison, and HTML report completed"
