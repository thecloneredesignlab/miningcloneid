#!/bin/bash -l
set -eu
study=$1
cd "$study"
mkdir -p verification/tmp
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
export LD_PRELOAD="$study/runtime/io_guard.so"
export CULTURE_GUARD_LOG="$study/verification/blocked_access.tsv" TMPDIR="$study/verification/tmp"
set +e
/usr/bin/strace -ff -s 1024 -e trace=%file,%process -o "$study/verification/trace" \
 /app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S \
 "$study/releases/r01/worker.py" "$study" > "$study/verification/console.log" 2>&1
code=$?
set -e
printf 'application_exit_code\t%s\n' "$code" > verification/application_status.tsv
/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S "$study/releases/r01/policy_check.py" "$study/verification"
unset LD_PRELOAD
exit "$code"
