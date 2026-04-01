#!/usr/bin/env bash
#SBATCH --job-name=invitro_obj
#SBATCH --array=1-500
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --qos=normal
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
TASK_SCRIPT="${REPO_ROOT}/workflow/sample_invitro_objective_once.R"
SEED="${SEED:-${SLURM_ARRAY_TASK_ID:-1}}"
MAXIT="${MAXIT:-100}"

resolve_array_stamp() {
  if [[ -n "${ARRAY_STAMP:-}" ]]; then
    printf '%s\n' "$ARRAY_STAMP"
    return 0
  fi
  if [[ -n "${SLURM_ARRAY_JOB_ID:-}" ]]; then
    local submit_time
    submit_time=$(scontrol show job "${SLURM_ARRAY_JOB_ID}" | awk '{
      for (i = 1; i <= NF; i++) {
        if ($i ~ /^SubmitTime=/) {
          sub(/^SubmitTime=/, "", $i)
          print $i
          exit
        }
      }
    }')
    if [[ -n "$submit_time" && "$submit_time" != "Unknown" ]]; then
      printf '%s\n' "${submit_time//[-:T]/_}"
      return 0
    fi
  fi
  date +%Y%m%d_%H%M%S
}

ARRAY_STAMP=$(resolve_array_stamp)
OUT_TSV=${OUT_TSV:-${REPO_ROOT}/workflow/sampled_objective_draws/invitro_objective_array_${ARRAY_STAMP}.tsv}

mkdir -p "$(dirname "$OUT_TSV")"

"${REPO_ROOT}/scripts/agentRrunner.sh" "$TASK_SCRIPT" \
  --seed="${SEED}" \
  --maxit="${MAXIT}" \
  --out_tsv="${OUT_TSV}" \
  --quiet=TRUE
