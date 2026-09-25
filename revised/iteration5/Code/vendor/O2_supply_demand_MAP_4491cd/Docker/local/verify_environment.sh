#!/usr/bin/env bash

set -euo pipefail

LOCAL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=docker_runtime.sh
source "${LOCAL_ROOT}/docker_runtime.sh"

o2sd_docker_exec \
  Rscript /opt/soft-coupling-environment/scripts/verify_environment.R \
  /opt/soft-coupling-environment/locks/packages.lock.tsv target

o2sd_docker_exec \
  python3.9 /opt/soft-coupling-environment/scripts/verify_python_environment.py \
  /opt/soft-coupling-environment/locks/requirements-repository-all-target.lock.txt

o2sd_docker_exec bash -c \
  'printf "R="; Rscript --vanilla -e '\''cat(as.character(getRversion()), "\n")'\''; printf "Python="; python3.9 --version; printf "Git="; git --version; printf "aria2="; aria2c --version | head -n 1'
