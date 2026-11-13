#!/usr/bin/env bash
set -euo pipefail
study="${1:?study}"; mode="${2:?mode}"; task="${3:?task}"; release="${4:?release}"
exec bash "${study}/releases/${release}/run.sh" "${study}" "${mode}" "${task}" "${release}"
