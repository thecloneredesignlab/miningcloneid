#!/bin/bash -l
# Conditional fallback dependency gate, not an upstream structural algorithm.
set -eu
study=$1
cd "$study"
mkdir -p runtime/python/work runtime/python/tmp
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
module -t list > runtime/python/modules.txt 2>&1
printf '%s\n' "$PYTHONPATH" > runtime/python/module_pythonpath.txt
export CULTURE_GUARD_LOG="$study/runtime/python/blocked_access.tsv"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export TMPDIR="$study/runtime/python/tmp"
export LD_PRELOAD="$study/runtime/io_guard.so"
cd runtime/python/work
set +e
timeout --kill-after=15s 180 /usr/bin/strace -ff -s 1024 -e trace=%file,%process \
  -o "$study/runtime/python/trace" \
  /app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S \
  "$study/releases/r01/runtime_python.py" "$study" > "$study/runtime/python/console.log" 2>&1
code=$?
set -e
printf 'application_exit_code\t%s\n' "$code" > "$study/runtime/python/application_status.tsv"
unset LD_PRELOAD
printf 'Fallback dependency gate ended; inspect policy logs before any solver control.\n'
