#!/bin/bash
# Add wall-clock hours to running SLURM array tasks whose remaining time is low.

set -euo pipefail

DEFAULT_THRESHOLD_HOURS="6"
DEFAULT_DRY_RUN="TRUE"
DEFAULT_AUTO_CONFIRM="TRUE"

THRESHOLD_HOURS="${THRESHOLD_HOURS:-${DEFAULT_THRESHOLD_HOURS}}"
DRY_RUN="${DRY_RUN:-${DEFAULT_DRY_RUN}}"
AUTO_CONFIRM="${AUTO_CONFIRM:-${DEFAULT_AUTO_CONFIRM}}"

usage() {
  cat <<'USAGE'
Usage:
  add_hours_to_running_array_jobs.sh ARRAY_ID ADD_HOURS

Arguments:
  ARRAY_ID   Parent SLURM array job id, for example 17130014.
  ADD_HOURS  Number of hours passed to s-add-job-hours for each selected task.

Environment:
  THRESHOLD_HOURS  Only extend RUNNING tasks with TimeLeft below this many hours.
                   Default: 6
  DRY_RUN          TRUE/true/1/yes prints commands without running them.
                   Default: TRUE
  AUTO_CONFIRM     TRUE/true/1/yes sends "y" to each s-add-job-hours prompt.
                   Default: TRUE

Examples:
  ./add_hours_to_running_array_jobs.sh 17130014 12
  DRY_RUN=FALSE ./add_hours_to_running_array_jobs.sh 17130014 12
  AUTO_CONFIRM=FALSE DRY_RUN=FALSE ./add_hours_to_running_array_jobs.sh 17130014 12
  THRESHOLD_HOURS=4 DRY_RUN=FALSE ./add_hours_to_running_array_jobs.sh 17130014 12
USAGE
}

is_true() {
  case "${1}" in
    TRUE|true|1|yes|YES) return 0 ;;
    *) return 1 ;;
  esac
}

require_positive_integer() {
  local label="$1"
  local value="$2"

  if ! [[ "${value}" =~ ^[0-9]+$ ]]; then
    echo "${label} must be a positive integer, got: ${value}" >&2
    exit 1
  fi
  if (( 10#${value} <= 0 )); then
    echo "${label} must be > 0, got: ${value}" >&2
    exit 1
  fi
}

time_to_seconds() {
  local value="$1"
  local days="0"
  local rest="${value}"
  local first=""
  local second=""
  local third=""
  local hours="0"
  local minutes="0"
  local seconds="0"

  case "${value}" in
    ""|N/A|UNLIMITED|NOT_SET|INVALID)
      return 1
      ;;
  esac

  if [[ "${rest}" == *-* ]]; then
    days="${rest%%-*}"
    rest="${rest#*-}"
  fi

  IFS=: read -r first second third <<< "${rest}"
  if [[ -n "${third}" ]]; then
    hours="${first}"
    minutes="${second}"
    seconds="${third}"
  elif [[ -n "${second}" ]]; then
    minutes="${first}"
    seconds="${second}"
  else
    seconds="${first}"
  fi

  for part in "${days}" "${hours}" "${minutes}" "${seconds}"; do
    if ! [[ "${part}" =~ ^[0-9]+$ ]]; then
      return 1
    fi
  done

  echo $((10#${days} * 86400 + 10#${hours} * 3600 + 10#${minutes} * 60 + 10#${seconds}))
}

extract_numeric_job_id() {
  local array_task_id="$1"

  scontrol show job "${array_task_id}" |
    awk 'NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i ~ /^JobId=/) {
          sub(/^JobId=/, "", $i)
          print $i
          exit
        }
      }
    }'
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ "$#" -ne 2 ]]; then
  usage >&2
  exit 1
fi

ARRAY_ID="$1"
ADD_HOURS="$2"

require_positive_integer "ARRAY_ID" "${ARRAY_ID}"
require_positive_integer "ADD_HOURS" "${ADD_HOURS}"
require_positive_integer "THRESHOLD_HOURS" "${THRESHOLD_HOURS}"

for cmd in squeue scontrol awk; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Required command not found: ${cmd}" >&2
    exit 1
  fi
done

if ! is_true "${DRY_RUN}" && ! command -v s-add-job-hours >/dev/null 2>&1; then
  echo "Required command not found: s-add-job-hours" >&2
  exit 1
fi

threshold_seconds=$((10#${THRESHOLD_HOURS} * 3600))
matched_count=0
extended_count=0
skipped_count=0

echo "Scanning RUNNING tasks for array ${ARRAY_ID}"
echo "  threshold_hours: ${THRESHOLD_HOURS}"
echo "  add_hours: ${ADD_HOURS}"
echo "  dry_run: ${DRY_RUN}"
echo "  auto_confirm: ${AUTO_CONFIRM}"

while IFS='|' read -r array_task_id time_left; do
  if [[ -z "${array_task_id}" ]]; then
    continue
  fi

  matched_count=$((matched_count + 1))

  if ! left_seconds="$(time_to_seconds "${time_left}")"; then
    echo "SKIP: ${array_task_id} TimeLeft=${time_left} is not parseable"
    skipped_count=$((skipped_count + 1))
    continue
  fi

  if (( left_seconds >= threshold_seconds )); then
    echo "SKIP: ${array_task_id} TimeLeft=${time_left} >= ${THRESHOLD_HOURS}h"
    skipped_count=$((skipped_count + 1))
    continue
  fi

  numeric_job_id="$(extract_numeric_job_id "${array_task_id}")"
  if ! [[ "${numeric_job_id}" =~ ^[0-9]+$ ]]; then
    echo "SKIP: ${array_task_id} internal JobId='${numeric_job_id}' is not pure numeric"
    skipped_count=$((skipped_count + 1))
    continue
  fi

  echo "ADD: ${array_task_id} internal JobId=${numeric_job_id} TimeLeft=${time_left} +${ADD_HOURS}h"
  if is_true "${DRY_RUN}"; then
    if is_true "${AUTO_CONFIRM}"; then
      echo "  dry-run command: printf 'y\\n' | s-add-job-hours ${numeric_job_id} ${ADD_HOURS}"
    else
      echo "  dry-run command: s-add-job-hours ${numeric_job_id} ${ADD_HOURS}"
    fi
  else
    if is_true "${AUTO_CONFIRM}"; then
      printf 'y\n' | s-add-job-hours "${numeric_job_id}" "${ADD_HOURS}"
    else
      s-add-job-hours "${numeric_job_id}" "${ADD_HOURS}"
    fi
  fi
  extended_count=$((extended_count + 1))
done < <(squeue -h -r -j "${ARRAY_ID}" -t RUNNING -o "%i|%L")

echo "Done."
echo "  running_tasks_seen: ${matched_count}"
echo "  selected_for_add_hours: ${extended_count}"
echo "  skipped: ${skipped_count}"
