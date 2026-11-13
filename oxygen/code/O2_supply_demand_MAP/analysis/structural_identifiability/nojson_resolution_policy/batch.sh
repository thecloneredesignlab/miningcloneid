#!/usr/bin/env bash
set -euo pipefail
project="${1:?project}"
study="${2:?study}"
release="${3:?release}"
mode="${4:?mode}"
task="${5:?task}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ "${mode}" == controller ]]; then
  python3 "${release}/policy.py" controller "${project}" "${study}"
else
  bash "${release}/evidence.sh" "${project}" "${study}" "${release}" "${task}"
fi
