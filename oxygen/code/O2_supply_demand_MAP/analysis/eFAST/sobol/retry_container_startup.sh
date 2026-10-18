#!/usr/bin/env bash
# Retry the RED identity lookup failure before a container starts its payload.
launch() {
  local attempt rc delay maximum base log
  maximum=${SOBOL_STARTUP_RETRIES:-6}
  base=${SOBOL_STARTUP_DELAY_SECONDS:-15}
  [[ "$maximum" =~ ^[1-9][0-9]*$ && "$base" =~ ^[0-9]+$ ]] || return 2
  for ((attempt=1; attempt<=maximum; attempt++)); do
    log="$scratch/launch.attempt${attempt}.log"
    set +e
    "$@" 2>&1 | tee "$log"
    rc=${PIPESTATUS[0]}
    set -e
    [[ "$rc" = 0 ]] && return 0
    # Both observed variants start with FATAL; ordinary analysis errors exit.
    if ! grep -Eq '^FATAL:.*user: unknown userid [0-9]+' "$log"; then
      return "$rc"
    fi
    echo "Container identity lookup failed before execution; attempt $attempt/$maximum" >&2
    if (( attempt < maximum )); then
      delay=$((base * (2 ** (attempt-1))))
      (( delay <= 120 )) || delay=120
      sleep "$delay"
    fi
  done
  return "$rc"
}
