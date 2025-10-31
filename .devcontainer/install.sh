#!/usr/bin/env bash
exec > >(tee .devcontainer/install.log) 2>&1
set -euo pipefail
echo "[install.sh] started $(date)"


# Resolve repo root robustly (works for any repo name)
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null || true)"
REPO_ROOT="${REPO_ROOT:-"$(cd "$SCRIPT_DIR/.." && pwd)"}"

if [ -n "${CODESPACES:-}" ]; then
  echo "[dotfiles] Bootstrapping Codespaces deps (incl. Java for rJava)…"

  # --- System deps ---
  sudo apt-get update
  sudo apt-get install -y \
    git-lfs \
    python3 python3-venv python3-dev \
    libnifti-dev \
    libmagick++-dev imagemagick gsfonts \
    libgdal-dev libproj-dev libgeos-dev libudunits2-dev libfftw3-dev gdal-bin\
    libmysqlclient-dev \
    libssl-dev libxml2-dev libcurl4-openssl-dev pkg-config \
    libpq-dev postgresql-client \
    default-jdk-headless \
    default-jdk \
    cmake pandoc libharfbuzz-dev libfribidi-dev

  # --- AWS CLI v2 ---
  # Installs the latest version of the AWS CLI from the official source
  echo "Installing AWS CLI v2..."
  curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
  unzip -q awscliv2.zip
  sudo ./aws/install
  rm -rf awscliv2.zip aws # Clean up installer files
  echo "AWS CLI installed successfully."
  
  # Make sure LFS is active (ok if already configured)
  git lfs install || true

  # --- Java env (best-effort detect) ---
  # This is the usual path on Ubuntu for JDK 17; adjust if `update-java-alternatives -l` differs.
  if [ -d "/usr/lib/jvm/java-21-openjdk-amd64" ]; then
    # Point to your JDK (21 on Ubuntu 24.04)
    export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
    echo "export JAVA_HOME=$JAVA_HOME" >> ~/.bashrc
    echo 'export PATH="$JAVA_HOME/bin:$PATH"' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH="$JAVA_HOME/lib/server:$LD_LIBRARY_PATH"' >> ~/.bashrc
  else
    # Fallback: pick first JDK from alternatives
    JAVA_HOME="$(update-java-alternatives -l | awk '{print $3; exit}')"
    export JAVA_HOME="${JAVA_HOME:-/usr/lib/jvm/default-java}"
    echo "export JAVA_HOME=${JAVA_HOME}" >> "${HOME}/.bashrc"
    echo 'export PATH="$JAVA_HOME/bin:$PATH"' >> "${HOME}/.bashrc"
  fi

  # Tell R where Java is and (re)configure JNI flags
  sudo R CMD javareconf

  # --- R deps (fast installer) ---
  R --quiet -e 'if (!requireNamespace("pak", quietly=TRUE)) install.packages("pak", repos="https://r-lib.github.io/p/pak/stable")'

  R --quiet -e 'pak::pkg_install(c(
    "glue","RNifti","data.table","magick","reticulate","raster","RMySQL","qualV",
    "gplots","gdata","RColorBrewer","gtools","flexclust",
    "Matrix","matlab","ape","rJava", "RPostgres", "yaml","dplyr","ggplot2","shiny",
    "devtools", "httr", "R.utils","plyr","survminer", "xlsx", "aws.s3"
  ))'

  R CMD INSTALL "$REPO_ROOT/code/cloneid_1.2.2.tar.gz"
  ##setup cloneid server access:
  Rscript -e 'cloneid::setupCLONEID(host=Sys.getenv()[["CLONEID_HOST"]],user=Sys.getenv()[["CLONEID_USER"]],password=Sys.getenv()[["CLONEID_PASSWORD"]])'


  # --- Clone supporting repositories ---
  echo "[dotfiles] Cloning dependent repositories..."

  # Ensure main workspace path exists
  WORKSPACE_PATH="$REPO_ROOT"
  # Clone cloneid-growthfit into /code
  mkdir -p "${WORKSPACE_PATH}/code"
  cd "${WORKSPACE_PATH}/code"
  if [ ! -d "cloneid-growthfit" ]; then
    echo "Cloning cloneid-growthfit..."
    git clone https://github.com/noemiandor/cloneid-growthfit.git
  else
    echo "cloneid-growthfit already present, skipping."
  fi

  # Clone Weakly-Supervised-Nuclei-Segmentation into /data/Feulgen_DNAPloidyAnalysis
  mkdir -p "${WORKSPACE_PATH}/data/Feulgen_DNAPloidyAnalysis"
  cd "${WORKSPACE_PATH}/data/Feulgen_DNAPloidyAnalysis"
  if [ ! -d "Weakly-Supervised-Nuclei-Segmentation" ]; then
    echo "Cloning Weakly-Supervised-Nuclei-Segmentation..."
    git clone https://github.com/CVIU-CSU/Weakly-Supervised-Nuclei-Segmentation.git
  else
    echo "Weakly-Supervised-Nuclei-Segmentation already present, skipping."
  fi

  # --- Test S3 Access ---
  # This relies on AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY being set as Codespaces secrets.
  if [ -z "${AWS_ACCESS_KEY_ID:-}" ] || [ -z "${AWS_SECRET_ACCESS_KEY:-}" ]; then
    echo "WARNING: AWS credentials not found in environment. Skipping S3 access test."
  else
    echo "Attempting to list S3 bucket to verify access..."
    if aws s3 ls "s3://cloneid4mysql8" >/dev/null; then
      echo "SUCCESS: Successfully connected to s3://cloneid4mysql8."
    else
      echo "ERROR: Failed to access s3://cloneid4mysql8. Check credentials and permissions."
      exit 1 # Exit with an error code
    fi
  fi
  
  ## pull the cell segmentation summary statistics
  aws s3 cp "s3://cloneid4mysql8/DetectionResults_InVitro.tar.gz" - | tar -xz -C data

  cd ~
  echo "[dotfiles] Done."
fi
