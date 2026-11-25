#!/bin/bash
set -euo pipefail

source_dir="$(cd "$(dirname "$0")" && pwd)"
python="${CULTURE_PYTHON:-python3}"
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

usage() {
  printf '%s\n' 'Usage: bash run_cases.sh {all|prepare|certify|collect} OUTPUT_DIR [DESIGN_DIR]'
}

if [[ $# == 1 && "$1" == --help ]]; then
  usage
  exit 0
fi
if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage >&2
  exit 2
fi

stage="$1"
case "$stage" in
  all|prepare|certify|collect) ;;
  *) usage >&2; exit 2 ;;
esac
if [[ "$stage" == prepare && $# == 3 ]]; then
  usage >&2
  exit 2
fi

path_value() {
  "$python" -c 'import sys
sys.path.insert(0, sys.argv[1])
import paths
try:
    value = paths.output_directory(sys.argv[3]) if sys.argv[2] == "OUT" else getattr(paths, sys.argv[2])
except ValueError as error:
    sys.exit(str(error))
print(value)' "$source_dir" "$@"
}

out="$(path_value OUT "$2")"
design="${3:-$out/design}"
if [[ $# == 2 && "$stage" != all && "$stage" != prepare && ! -d "$design" ]]; then
  design="$(path_value FROZEN_DESIGN)"
fi
mkdir -p "$out"

prepare() {
  if [[ -f "$design/SHA256SUMS" ]]; then
    (cd "$design" && shasum -a 256 -c SHA256SUMS)
  else
    "$python" "$source_dir/freeze_design.py" --out "$design"
  fi
}

check_case_design() {
  local hashes receipt
  receipt="$out/case_design.sha256"
  hashes="$(cd "$design" && shasum -a 256 passage_graph.tsv event_mask.tsv case_register.tsv \
    parameter_roles_by_case.tsv contrast_register.tsv)"
  if [[ -f "$receipt" ]]; then
    if [[ "$hashes" != "$(cat "$receipt")" ]]; then
      printf 'The case design differs from the recorded run design in %s.\n' "$out" >&2
      return 1
    fi
  else
    printf '%s\n' "$hashes" > "$receipt"
  fi
}

run() {
  local output="$1"
  shift
  if [[ -s "$out/$output/summary.tsv" ]]; then
    return
  fi
  "$python" "$source_dir/joint_rank.py" --design "$design" --out "$out/$output" "$@" > "$out/$output.log" 2>&1
}

certify() {
  check_case_design
  if [[ ! -s "$out/archived_free_control/rank.tsv" ]]; then
    local prior finite_state
    prior="$(path_value STRUCTURAL)/20261118_culture_experimental_fixed_inoculum_01"
    finite_state="$(path_value FINITE_STATE)"
    PYTHONPATH="$finite_state${PYTHONPATH:+:$PYTHONPATH}" /usr/bin/time -p "$python" "$prior/certify_fixed.py" \
      --phase counts --visibility shared --bits 768 --out "$out/archived_free_control" \
      > "$out/archived_free_control.log" 2>&1
  fi
  "$python" "$source_dir/validate_joint.py" --design "$design" --out "$out" --controls --archive > "$out/controls.log" 2>&1
  run pilot_inherit_validated --case inherit --bits 384 --masks CDK --validate
  run pilot_free --case free --bits 384 --root SUM-159_NLS_2N_dp_seed --masks CDK --validate
  run free_384 --case free --bits 384
  run inherit_512 --case inherit --bits 512 --validate
  run free_512_same_root --case free --bits 512 --root SUM-159_NLS_2N_dp_seed --validate
  run free_512_second --case free --bits 512 --witness second --validate
  run inherit_384_second --case inherit --bits 384 --witness second
  run gaussian_384 --case gaussian --bits 384 --validate
  run tv_384 --case tv --bits 384 --masks CDK --validate
  run elapsed_inherit_384 --case elapsed_inherit --bits 384 --masks CDK --validate
  run elapsed_free_384 --case elapsed_free --bits 384 --masks CDK
  run free_seeds_384 --case free_seeds --bits 384 --masks CDK
}

collect() {
  check_case_design
  "$python" "$source_dir/collect_results.py" --design "$design" --out "$out" > "$out/collection.log" 2>&1
  Rscript "$source_dir/plot_joint.R" "$out" "$design" "$(path_value FIT)"
  "$python" "$source_dir/finalize_release.py" --design "$design" --out "$out"
}

case "$stage" in
  prepare) prepare ;;
  certify) certify ;;
  collect) collect ;;
  all)
    if [[ $# == 2 ]]; then prepare; fi
    certify
    collect
    ;;
esac
