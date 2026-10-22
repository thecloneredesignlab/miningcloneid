#!/usr/bin/env bash
# Reject ambiguous legacy commands; explicit current options use current workers.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
alias_name="${1:?Compatibility filename}";shift
for arg in "$@"; do
  [[ "$arg" == --*=* ]] || { echo "$alias_name: legacy positional arguments/task indices are retired. Use run_Figure78.sh --figure=7|8 --stage=..." >&2; exit 2; }
done
[[ "$#" -gt 0 ]] || { echo "$alias_name: specify the current interface explicitly; no implicit execution" >&2; exit 2; }
case "$alias_name" in
  submit_*) exec bash "${script_dir}/submit_Figure78.sh" "$@" ;;
  *) exec bash "${script_dir}/run_Figure78.sh" "$@" ;;
esac
