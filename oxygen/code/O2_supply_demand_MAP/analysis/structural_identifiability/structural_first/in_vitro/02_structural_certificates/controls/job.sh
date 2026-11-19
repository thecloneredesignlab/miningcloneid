#!/bin/bash -l
set -eu
study=$1
task=$2
cd "$study"
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
export LD_PRELOAD="$study/runtime/io_guard.so"
python_bin=/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3
run_case() {
  mode=$1
  case_id=$2
  if [ "$mode" = preflight ]; then folder="$study/preflight/$case_id"; else folder="$study/controls/$case_id"; fi
  mkdir -p "$folder/tmp"
  export TMPDIR="$folder/tmp" CULTURE_GUARD_LOG="$folder/blocked_access.tsv"
  set +e
  timeout --kill-after=15s 300 /usr/bin/strace -ff -s 1024 -e trace=%file,%process -o "$folder/trace" \
    "$python_bin" -B -S "$study/releases/r01/worker.py" "$study" "$mode" "$case_id" > "$folder/console.log" 2>&1
  code=$?
  set -e
  printf 'application_exit_code\t%s\n' "$code" > "$folder/application_status.tsv"
  [ "$code" -eq 0 ] || return "$code"
  "$python_bin" -B -S "$study/releases/r01/policy_check.py" "$folder"
}
if [ "$task" = preflight ]; then
  for case_id in C02 C03 C04 C05; do run_case preflight "$case_id"; done
else
  case "$task" in C02|C03|C04|C05) run_case control "$task";; *) exit 2;; esac
fi
unset LD_PRELOAD
