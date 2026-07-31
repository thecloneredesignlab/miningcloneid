#!/usr/bin/env bash

#SBATCH --job-name=o2_env_verify
#SBATCH --qos=small
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -euo pipefail

if [[ "$#" -ne 1 ]]; then
  echo "Usage: verify_hpc_observed.sh AUDIT_ROOT" >&2
  exit 2
fi

AUDIT_ROOT="$1"
if [[ -f /etc/profile.d/modules.sh ]]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
fi
module use /app/eb/modules/all >/dev/null 2>&1 || true
if command -v ml >/dev/null 2>&1; then
  ml R/4.4
else
  module load R/4.4
fi

Rscript "${AUDIT_ROOT}/scripts/verify_environment.R" \
  "${AUDIT_ROOT}/locks/packages.lock.tsv" \
  hpc-observed \
  >"${AUDIT_ROOT}/hpc_snapshot/r/package-lock-verification.tsv"

