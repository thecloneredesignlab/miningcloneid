#!/usr/bin/env bash
# Build and verify a Julia structural-identifiability SIF from an immutable OCI digest.

#SBATCH --job-name=o2_structid_sif
#SBATCH --qos=xxlarge
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 zafiro/o2_supply_demand_map@sha256:DIGEST OUTPUT_SIF" >&2
  exit 64
fi

SOURCE_IMAGE=$1
OUTPUT_SIF=$2
if [[ ! "$SOURCE_IMAGE" =~ ^zafiro/o2_supply_demand_map@sha256:[a-f0-9]{64}$ ]]; then
  echo "An immutable zafiro/o2_supply_demand_map digest is required." >&2
  exit 64
fi
if [[ ! "$OUTPUT_SIF" = /* || "$OUTPUT_SIF" != *.sif ]]; then
  echo "OUTPUT_SIF must be an absolute .sif path." >&2
  exit 64
fi

OUTPUT_DIR=$(dirname "$OUTPUT_SIF")
OUTPUT_TMP="${OUTPUT_SIF}.new"
LOG_DIR="${OUTPUT_DIR}/build_logs"
LOG_PREFIX="${LOG_DIR}/$(basename "$OUTPUT_SIF" .sif)"
if [[ ! -d "$OUTPUT_DIR" ]]; then
  echo "Output directory does not exist: $OUTPUT_DIR" >&2
  exit 66
fi
if [[ -e "$OUTPUT_SIF" || -e "$OUTPUT_TMP" ]]; then
  echo "Refusing to overwrite an existing SIF or temporary candidate." >&2
  exit 73
fi

mkdir -p "$LOG_DIR"
LOCAL_TMP=$(mktemp -d "/tmp/o2-structid-sif-${SLURM_JOB_ID:-manual}.XXXXXX")
PUBLISHED=0
cleanup() {
  rm -rf "$LOCAL_TMP"
  if [[ "$PUBLISHED" -ne 1 && -e "$OUTPUT_TMP" ]]; then
    rm -f "$OUTPUT_TMP"
  fi
}
trap cleanup EXIT
export APPTAINER_TMPDIR="${LOCAL_TMP}/tmp"
export APPTAINER_CACHEDIR="${LOCAL_TMP}/cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

printf '%s\n' "$SOURCE_IMAGE" > "${LOG_PREFIX}.source.txt"
apptainer build "$OUTPUT_TMP" "docker://${SOURCE_IMAGE}"

{
  apptainer exec --cleanenv --containall "$OUTPUT_TMP" \
    /usr/bin/python3.9 -m pip check
  apptainer exec --cleanenv --containall "$OUTPUT_TMP" \
    /usr/bin/python3.9 \
    /opt/soft-coupling-environment/scripts/verify_python_environment.py \
    /opt/soft-coupling-environment/locks/requirements-efast-py39.lock.txt
  apptainer exec --cleanenv --containall "$OUTPUT_TMP" \
    /usr/bin/python3.9 \
    /opt/soft-coupling-environment/scripts/verify_efast_environment.py
  apptainer exec --cleanenv --containall "$OUTPUT_TMP" env \
    JULIA_DEPOT_PATH=/opt/julia-depot \
    JULIA_PROJECT=/opt/structural-identifiability \
    JULIA_PKG_OFFLINE=true \
    JULIA_NUM_THREADS=1 \
    /usr/local/julia/bin/julia --startup-file=no \
    /opt/soft-coupling-environment/scripts/verify_structural_identifiability.jl
  apptainer exec --cleanenv --containall "$OUTPUT_TMP" \
    /usr/local/bin/o2-hpc-exact-rscript --vanilla -e '
stopifnot(
  as.character(getRversion()) == "4.4.2",
  identical(R.home(), "/app/eb/software/R/4.4.2-gfbf-2024a/lib64/R"),
  identical(unname(extSoftVersion()["BLAS"]), "FlexiBLAS OPENBLAS"),
  identical(paste(La_version(), collapse = "."), "3.12.0"),
  identical(.libPaths()[[1L]], "/opt/o2-host-r-library/4.4")
)
pkgs <- c("DEoptim", "dplyr", "ggplot2", "Matrix", "Rcpp", "RcppEigen", "tidyr")
stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))
cat("HPC-exact R runtime verification: PASS\n")
'
} 2>&1 | tee "${LOG_PREFIX}.verification.log"

apptainer inspect "$OUTPUT_TMP" > "${LOG_PREFIX}.inspect.txt"
sha256sum "$OUTPUT_TMP" > "${LOG_PREFIX}.candidate.sha256"
mv "$OUTPUT_TMP" "$OUTPUT_SIF"
PUBLISHED=1
sha256sum "$OUTPUT_SIF" > "${LOG_PREFIX}.sha256"
echo "SIF_BUILD=PASS"
echo "output_sif=$OUTPUT_SIF"
