#!/usr/bin/env bash
# Build a standalone SIF by copying the exact EasyBuild prefixes and R package
# library used by the accepted RED R/4.4.2 seed10 run into an existing O2 SIF.

#SBATCH --job-name=o2_sif_hpc_exact
#SBATCH --qos=xxlarge
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

if [[ "$#" -ne 6 ]]; then
  echo "Usage: $0 PROJECT_ROOT BUILD_ROOT BASE_SIF OUTPUT_SIF GOLD_RUN_ROOT HOST_R_LIBRARY" >&2
  exit 64
fi

PROJECT_ROOT=$(realpath "$1")
BUILD_ROOT=$(realpath -m "$2")
BASE_SIF=$(realpath "$3")
OUTPUT_SIF=$(realpath -m "$4")
GOLD_RUN_ROOT=$(realpath "$5")
HOST_R_LIBRARY=$(realpath "$6")

MODULE_MANIFEST=${PROJECT_ROOT}/container/manifests/hpc-loaded-modules.tsv
ENVIRONMENT_SNAPSHOT=${PROJECT_ROOT}/container/hpc_snapshot/modules/selected-environment.tsv
WRAPPER_SOURCE=${PROJECT_ROOT}/container/hpc_exact/o2_hpc_exact_rscript
DRIVER_SOURCE=${GOLD_RUN_ROOT}/run_fit_invitro_worker_bindings_hostR442_20260814.R
PROFILE_SOURCE=${GOLD_RUN_ROOT}/Rprofile_o2_headless_hostR442_20260814.R
SANDBOX=${BUILD_ROOT}/sandbox
APPTAINER_TMPDIR=${BUILD_ROOT}/apptainer_tmp
PROVENANCE=${BUILD_ROOT}/provenance
OUTPUT_TMP=${OUTPUT_SIF}.new

for required in \
  "$MODULE_MANIFEST" \
  "$ENVIRONMENT_SNAPSHOT" \
  "$WRAPPER_SOURCE" \
  "$DRIVER_SOURCE" \
  "$PROFILE_SOURCE"; do
  if [[ ! -f "$required" ]]; then
    echo "Required build input is missing: $required" >&2
    exit 66
  fi
done

if [[ -e "$SANDBOX" || -e "$OUTPUT_TMP" || -e "$OUTPUT_SIF" ]]; then
  echo "Refusing to overwrite an existing build or image path." >&2
  echo "sandbox=$SANDBOX" >&2
  echo "output_tmp=$OUTPUT_TMP" >&2
  echo "output_sif=$OUTPUT_SIF" >&2
  exit 73
fi

mkdir -p "$BUILD_ROOT" "$APPTAINER_TMPDIR" "$PROVENANCE"
export APPTAINER_TMPDIR

source /etc/profile.d/modules.sh
module use /app/eb/modules/all
module purge
ml R/4.4.2

module -t list 2> "${PROVENANCE}/module-list.txt" || true
cp "$MODULE_MANIFEST" "${PROVENANCE}/hpc-loaded-modules.tsv"
cp "$ENVIRONMENT_SNAPSHOT" "${PROVENANCE}/selected-environment.tsv"
cp "${GOLD_RUN_ROOT}/runtime.tsv" "${PROVENANCE}/gold-runtime.tsv"
cp "${GOLD_RUN_ROOT}/runtime_packages.tsv" "${PROVENANCE}/gold-runtime-packages.tsv"
cp "${GOLD_RUN_ROOT}/sessionInfo.txt" "${PROVENANCE}/gold-sessionInfo.txt"
sha256sum "$BASE_SIF" > "${PROVENANCE}/base-sif.sha256"

apptainer build --sandbox "$SANDBOX" "$BASE_SIF"

mkdir -p "${SANDBOX}/app/eb/software"
while IFS=$'\t' read -r _order module_name; do
  [[ "$module_name" == "module" ]] && continue
  case "$module_name" in
    shared|DefaultModules|gfbf/*)
      continue
      ;;
  esac

  source_prefix=/app/eb/software/${module_name}
  destination_prefix=${SANDBOX}${source_prefix}
  if [[ ! -e "$source_prefix" ]]; then
    printf '%s\n' "$source_prefix" >> "${PROVENANCE}/missing-easybuild-prefixes.txt"
    continue
  fi

  resolved_prefix=$(realpath "$source_prefix")
  mkdir -p "$destination_prefix"
  rsync -a "${resolved_prefix}/" "${destination_prefix}/"
  printf '%s\t%s\n' "$source_prefix" "$resolved_prefix" >> "${PROVENANCE}/copied-easybuild-prefixes.tsv"
done < "$MODULE_MANIFEST"

if [[ -s "${PROVENANCE}/missing-easybuild-prefixes.txt" ]]; then
  echo "One or more EasyBuild prefixes were missing:" >&2
  cat "${PROVENANCE}/missing-easybuild-prefixes.txt" >&2
  exit 66
fi

mkdir -p "${SANDBOX}/opt/o2-host-r-library/4.4"
rsync -a \
  --exclude='00LOCK*' \
  "${HOST_R_LIBRARY}/" \
  "${SANDBOX}/opt/o2-host-r-library/4.4/"

mkdir -p \
  "${SANDBOX}/opt/o2-hpc-exact/bin" \
  "${SANDBOX}/opt/o2-hpc-exact/provenance" \
  "${SANDBOX}/usr/local/bin"
install -m 0755 "$WRAPPER_SOURCE" "${SANDBOX}/usr/local/bin/o2-hpc-exact-rscript"
install -m 0755 "$DRIVER_SOURCE" "${SANDBOX}/opt/o2-hpc-exact/bin/run_fit_invitro_worker_bindings.R"
install -m 0644 "$PROFILE_SOURCE" "${SANDBOX}/opt/o2-hpc-exact/Rprofile_o2_headless.R"
cp -a "${PROVENANCE}/." "${SANDBOX}/opt/o2-hpc-exact/provenance/"

if find "${SANDBOX}/app/eb/software" -xtype l -print -quit | grep -q .; then
  echo "Broken symlink found in copied EasyBuild closure:" >&2
  find "${SANDBOX}/app/eb/software" -xtype l -print >&2
  exit 65
fi

env \
  APPTAINERENV_OMP_NUM_THREADS=1 \
  APPTAINERENV_OPENBLAS_NUM_THREADS=1 \
  APPTAINERENV_MKL_NUM_THREADS=1 \
  apptainer exec \
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
'

apptainer build "$OUTPUT_TMP" "$SANDBOX"
mv "$OUTPUT_TMP" "$OUTPUT_SIF"
sha256sum "$OUTPUT_SIF" | tee "${PROVENANCE}/output-sif.sha256"
apptainer inspect "$OUTPUT_SIF" > "${PROVENANCE}/output-sif.inspect.txt"

echo "standalone HPC-exact SIF build completed"
echo "output_sif=$OUTPUT_SIF"
