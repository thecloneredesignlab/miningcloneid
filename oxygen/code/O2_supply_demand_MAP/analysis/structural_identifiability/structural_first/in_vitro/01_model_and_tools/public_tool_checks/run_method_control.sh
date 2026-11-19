#!/bin/bash -l
# Finite array of registered toy controls; no array throttle or project cases.
set -eu
study=$1
cases=(C02 C03 C04 C05 C01)
case_id=${cases[$SLURM_ARRAY_TASK_ID]}
cd "$study"
mkdir -p "controls/$case_id/tmp"
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
export CULTURE_GUARD_LOG="$study/controls/$case_id/blocked_access.tsv"
export CULTURE_DISABLE_DEBUGGER=1 PYTHONHASHSEED=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export TMPDIR="$study/controls/$case_id/tmp"
export LD_PRELOAD="$study/runtime/io_guard.so"
set +e
timeout --kill-after=15s 300 /usr/bin/strace -ff -s 1024 -e trace=%file,%process \
  -o "$study/controls/$case_id/trace" \
  /app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S \
  "$study/releases/r03/method_control.py" "$study" "$case_id" \
  > "$study/controls/$case_id/console.log" 2>&1
code=$?
set -e
printf 'application_exit_code\t%s\n' "$code" > "$study/controls/$case_id/application_status.tsv"
unset LD_PRELOAD
printf 'Registered method control ended: %s\n' "$case_id"
