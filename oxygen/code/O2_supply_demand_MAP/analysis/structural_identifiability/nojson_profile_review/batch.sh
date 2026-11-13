#!/usr/bin/env bash
set -euo pipefail
project="${1:?project}"
study="${2:?study}"
mode="${3:?mode}"
task="${4:?task}"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ "${mode}" == controller || "${mode}" == validate || "${mode}" == validate_local ]]; then
  python3 "${study}/implementation/workflow.py" "${mode}" "${project}" "${study}"
else
  bash "${study}/implementation/run.sh" "${project}" "${study}" "${mode}" "${task}"
fi
