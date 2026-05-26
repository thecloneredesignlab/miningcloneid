#!/bin/bash
# Compatibility wrapper for the unified O2G submitter.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UNIFIED_SUBMIT="${UNIFIED_SUBMIT:-${SCRIPT_DIR}/submit_o2g_fit.sh}"

usage() {
  cat <<'EOF'
Usage:
  bash submit_invivo_invitro_joint_buffering.sh [submit_o2g_fit options]

Compatibility behavior:
  This wrapper now delegates to submit_o2g_fit.sh with:
    --fitting_mode=joint
    --joint_fitting_mode=${JOINT_FITTING_MODE:-JOINT}

Use this preferred entry directly for new runs:
  bash submit_o2g_fit.sh --help
EOF
}

for arg in "$@"; do
  case "${arg}" in
    --help|-h)
      usage
      echo
      bash "${UNIFIED_SUBMIT}" --help
      exit 0
      ;;
  esac
done

exec bash "${UNIFIED_SUBMIT}" \
  "--fitting_mode=joint" \
  "--joint_fitting_mode=${JOINT_FITTING_MODE:-JOINT}" \
  "$@"
