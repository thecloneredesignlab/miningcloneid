#!/usr/bin/env bash

#SBATCH --job-name=o2_seed10_compare
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "Usage: $0 REFERENCE_SEED_DIR CANDIDATE_SEED_DIR OUTPUT_TSV" >&2
  exit 64
fi

REFERENCE_SEED_DIR=$(realpath "$1")
CANDIDATE_SEED_DIR=$(realpath "$2")
OUTPUT_TSV=$(realpath -m "$3")

for relative_path in \
  fit_summary.tsv \
  best_params.tsv \
  best_params_transformed.tsv; do
  test -s "${REFERENCE_SEED_DIR}/${relative_path}"
  test -s "${CANDIDATE_SEED_DIR}/${relative_path}"
done
test -s "${CANDIDATE_SEED_DIR}/report/fit_report_seed10.html"

cmp "${REFERENCE_SEED_DIR}/best_params.tsv" "${CANDIDATE_SEED_DIR}/best_params.tsv"
cmp \
  "${REFERENCE_SEED_DIR}/best_params_transformed.tsv" \
  "${CANDIDATE_SEED_DIR}/best_params_transformed.tsv"

summary_value() {
  local summary_file=$1
  local metric=$2
  awk -F '\t' -v wanted="$metric" '$1 == wanted { print $2; exit }' "$summary_file"
}

mkdir -p "$(dirname "$OUTPUT_TSV")"
{
  printf 'metric\treference\tcandidate\texact_match\n'
  for metric in optimizer_deoptim_objective optimizer_local_objective objective_total; do
    reference_value=$(summary_value "${REFERENCE_SEED_DIR}/fit_summary.tsv" "$metric")
    candidate_value=$(summary_value "${CANDIDATE_SEED_DIR}/fit_summary.tsv" "$metric")
    exact_match=FALSE
    [[ "$reference_value" == "$candidate_value" ]] && exact_match=TRUE
    printf '%s\t%s\t%s\t%s\n' "$metric" "$reference_value" "$candidate_value" "$exact_match"
    [[ "$exact_match" == TRUE ]]
  done
} > "$OUTPUT_TSV"

echo "seed10 numerical outputs and parameter tables match the reference exactly"
