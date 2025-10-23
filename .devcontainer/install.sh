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
    libmagick++-dev imagemagick \
    libgdal-dev libproj-dev libgeos-dev libudunits2-dev libfftw3-dev \
    libmysqlclient-dev \
    libssl-dev libxml2-dev libcurl4-openssl-dev pkg-config \
    openjdk-17-jdk-headless

  # Make sure LFS is active (ok if already configured)
  git lfs install || true

  # --- Java env (best-effort detect) ---
  # This is the usual path on Ubuntu for JDK 17; adjust if `update-java-alternatives -l` differs.
  if [ -d "/usr/lib/jvm/java-17-openjdk-amd64" ]; then
    export JAVA_HOME="/usr/lib/jvm/java-17-openjdk-amd64"
  else
    # Fallback: pick first JDK from alternatives
    JAVA_HOME="$(update-java-alternatives -l | awk '{print $3; exit}')"
    export JAVA_HOME="${JAVA_HOME:-/usr/lib/jvm/default-java}"
  fi
  echo "export JAVA_HOME=${JAVA_HOME}" >> "${HOME}/.bashrc"
  echo 'export PATH="$JAVA_HOME/bin:$PATH"' >> "${HOME}/.bashrc"

  # Tell R where Java is and (re)configure JNI flags
  sudo R CMD javareconf

  # --- R deps (fast installer) ---
  R --quiet -e 'if (!requireNamespace("pak", quietly=TRUE)) install.packages("pak", repos="https://r-lib.github.io/p/pak/stable")'

  R --quiet -e 'pak::pkg_install(c(
    "glue","RNifti","data.table","magick","reticulate","raster","RMySQL","qualV",
    "gplots","gdata","RColorBrewer","gtools","flexclust",
    "Matrix","matlab","ape","rJava"
  ))'

  echo "[dotfiles] Done."
fi
