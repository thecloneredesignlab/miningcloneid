#!/usr/bin/env bash
# Resume an HPC-exact build after the sandbox and copied runtime are complete.
# APPTAINER_TMPDIR is placed on compute-node-local storage because the RED SMB
# share does not permit the lchown operations used by Apptainer's build packer.

#SBATCH --job-name=o2_sif_hpc_exact_resume
#SBATCH --qos=xxlarge
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "Usage: $0 BUILD_ROOT OUTPUT_SIF EXPECTED_PACKAGE_LIBRARY" >&2
  exit 64
fi

BUILD_ROOT=$(realpath "$1")
OUTPUT_SIF=$(realpath -m "$2")
EXPECTED_PACKAGE_LIBRARY=$(realpath "$3")

SANDBOX=${BUILD_ROOT}/sandbox
PROVENANCE=${BUILD_ROOT}/provenance
OUTPUT_TMP=${OUTPUT_SIF}.new
LOCAL_TMP=/tmp/o2_hpc_exact_build_${SLURM_JOB_ID}

if [[ ! -d "$SANDBOX" ]]; then
  echo "Completed sandbox is missing: $SANDBOX" >&2
  exit 66
fi
if [[ -e "$OUTPUT_TMP" || -e "$OUTPUT_SIF" || -e "$LOCAL_TMP" ]]; then
  echo "Refusing to overwrite an existing output or temporary path." >&2
  exit 73
fi

mkdir -p "$LOCAL_TMP"
trap 'rm -rf "$LOCAL_TMP"' EXIT
export APPTAINER_TMPDIR=$LOCAL_TMP

source_bytes=$(du -sb "$EXPECTED_PACKAGE_LIBRARY" | awk '{print $1}')
image_bytes=$(du -sb "${SANDBOX}/opt/o2-host-r-library/4.4" | awk '{print $1}')
printf 'source_package_library_bytes\t%s\n' "$source_bytes" > "${PROVENANCE}/package-library-size-check.tsv"
printf 'image_package_library_bytes\t%s\n' "$image_bytes" >> "${PROVENANCE}/package-library-size-check.tsv"

apptainer exec \
  --writable \
  --cleanenv \
  --containall \
  "$SANDBOX" \
  /usr/local/bin/o2-hpc-exact-rscript --vanilla -e '
stopifnot(
  identical(R.home(), "/app/eb/software/R/4.4.2-gfbf-2024a/lib64/R"),
  identical(unname(extSoftVersion()["BLAS"]), "FlexiBLAS OPENBLAS"),
  identical(paste(La_version(), collapse = "."), "3.12.0"),
  identical(.libPaths()[[1L]], "/opt/o2-host-r-library/4.4")
)
pkgs <- c("DEoptim", "dplyr", "ggplot2", "Matrix", "Rcpp", "RcppEigen", "tidyr")
stopifnot(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))
print(sessionInfo())
' > "${PROVENANCE}/sandbox-runtime-check.txt"

apptainer build "$OUTPUT_TMP" "$SANDBOX"
mv "$OUTPUT_TMP" "$OUTPUT_SIF"
sha256sum "$OUTPUT_SIF" | tee "${PROVENANCE}/output-sif.sha256"
apptainer inspect "$OUTPUT_SIF" > "${PROVENANCE}/output-sif.inspect.txt"

echo "standalone HPC-exact SIF resume build completed"
echo "output_sif=$OUTPUT_SIF"
