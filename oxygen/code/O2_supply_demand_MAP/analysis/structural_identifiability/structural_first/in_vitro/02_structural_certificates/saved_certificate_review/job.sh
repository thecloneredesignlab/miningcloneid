#!/bin/bash -l
set -eu
study=$1
mode=$2
case "$mode" in
  review_C03) folder="$study/review/C03";;
  native_C04) folder="$study/controls/C04";;
  native_C05) folder="$study/controls/C05";;
  collect_C04) folder="$study/validation_retry/C04";;
  collect_C05) folder="$study/validation_retry/C05";;
  *) exit 2;;
esac
cd "$study"
mkdir -p "$folder/tmp"
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
export LD_PRELOAD="$study/runtime/io_guard.so"
export CULTURE_GUARD_LOG="$folder/blocked_access.tsv" TMPDIR="$folder/tmp"
set +e
timeout --kill-after=15s 300 /usr/bin/strace -ff -s 1024 -e trace=%file,%process -o "$folder/trace" \
 /app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S \
 "$study/releases/r01/worker.py" "$study" "$mode" > "$folder/console.log" 2>&1
code=$?
set -e
printf 'application_exit_code\t%s\n' "$code" > "$folder/application_status.tsv"
[ "$code" -eq 0 ] || exit "$code"
/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S "$study/releases/r01/policy_check.py" "$folder"
unset LD_PRELOAD
