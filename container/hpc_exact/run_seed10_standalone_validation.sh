#!/usr/bin/env bash
# Validate that the HPC-exact SIF reproduces the canonical host seed10 result
# without binding the host EasyBuild tree or host R package library.

#SBATCH --job-name=o2ivt_s10_standalone
#SBATCH --qos=xxlarge
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=22
#SBATCH --mem=16G

set -euo pipefail

if [[ "$#" -ne 4 ]]; then
  echo "Usage: $0 PROJECT_ROOT RUN_ROOT SIF_PATH REFERENCE_SEED_DIR" >&2
  exit 64
fi

PROJECT_ROOT=$(realpath "$1")
RUN_ROOT=$(realpath -m "$2")
SIF_PATH=$(realpath "$3")
REFERENCE_SEED_DIR=$(realpath "$4")

CONTAINER_RSCRIPT=/usr/local/bin/o2-hpc-exact-rscript
WORKER_DRIVER=/opt/o2-hpc-exact/bin/run_fit_invitro_worker_bindings.R
HEADLESS_PROFILE=/opt/o2-hpc-exact/Rprofile_o2_headless.R
MODEL_CACHE_REL=oxygen/code/O2_supply_demand_MAP/model/.rcpp_cache_o2_supply_demand_map
MODEL_CACHE_HOST=${PROJECT_ROOT}/${MODEL_CACHE_REL}
ISOLATED_CACHE=${RUN_ROOT}/rcpp_cache
SEED_DIR=${RUN_ROOT}/seed10
LOG_DIR=${RUN_ROOT}/log

mkdir -p "$RUN_ROOT" "$SEED_DIR" "$LOG_DIR" "$ISOLATED_CACHE"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export LANG=C.UTF-8
export LC_ALL=C.UTF-8

standalone_exec() {
  env \
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

standalone_exec "$CONTAINER_RSCRIPT" --vanilla -e '
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

standalone_exec "$CONTAINER_RSCRIPT" "$WORKER_DRIVER" \
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
  standalone_exec "$CONTAINER_RSCRIPT" \
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
  printf 'metric\treference\tstandalone\texact_match\n'
  for metric in optimizer_deoptim_objective optimizer_local_objective objective_total; do
    reference_value=$(summary_value "$REFERENCE_SEED_DIR/fit_summary.tsv" "$metric")
    standalone_value=$(summary_value "$SEED_DIR/fit_summary.tsv" "$metric")
    exact_match=FALSE
    [[ "$reference_value" == "$standalone_value" ]] && exact_match=TRUE
    printf '%s\t%s\t%s\t%s\n' "$metric" "$reference_value" "$standalone_value" "$exact_match"
    [[ "$exact_match" == TRUE ]]
  done
} > "$RUN_ROOT/reference_comparison.tsv"

echo "standalone HPC-exact SIF seed10 fit, exact comparison, and HTML report completed"
