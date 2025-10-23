#!/usr/bin/env bash
set -euo pipefail

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
    default-jdk

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
    "Matrix","matlab","ape","rJava", "RPostgres"
  ))'

  R CMD INSTALL code/cloneid_1.2.2.tar.gz 

  # --- Clone supporting repositories ---
  echo "[dotfiles] Cloning dependent repositories..."

  # Ensure main workspace path exists
  WORKSPACE_PATH="/workspaces/IMO-workshop-2025"

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

  cd ~
  echo "[dotfiles] Done."
fi
