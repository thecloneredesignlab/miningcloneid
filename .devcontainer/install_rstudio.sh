#!/usr/bin/env bash
set -euo pipefail

# Install prerequisites
sudo apt-get update
sudo apt-get install -y --no-install-recommends gdebi-core wget

# Pick an RStudio Server build for Ubuntu 22.04 (Jammy).
RS_URL="https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2025.09.2-418-amd64.deb"

wget -O /tmp/rstudio-server.deb "$RS_URL"
sudo gdebi -n /tmp/rstudio-server.deb || sudo apt-get -f install -y


echo "$(whoami):rstudio" | sudo chpasswd

# Start RStudio Server
sudo rstudio-server stop || true
sudo rstudio-server start