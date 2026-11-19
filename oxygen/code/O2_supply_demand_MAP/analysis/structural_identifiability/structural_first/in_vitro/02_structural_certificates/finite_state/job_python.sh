#!/bin/bash -l
set -eu
study=$1
release=$2
task=$3
program=$4
out="$study/verification/$task"
mkdir -p "$out/tmp"
test "${SLURM_CPUS_PER_TASK}" = 1
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MAKEFLAGS=-j1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 SYMPY_GROUND_TYPES=python
export LD_PRELOAD="$study/runtime/io_guard.so" CULTURE_GUARD_LOG="$out/blocked_access.tsv" TMPDIR="$out/tmp"
/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S "$study/releases/$release/verify_files.py" "$study" "$release" > "$out/source_verification.txt"
set +e
/usr/bin/strace -ff -s 1024 -e trace=%file,%process -o "$out/trace" \
 /app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S \
 "$study/releases/$release/$program" "$study" "$task" > "$out/console.log" 2>&1
code=$?
set -e
printf 'application_exit_code\t%s\n' "$code" > "$out/application_status.tsv"
/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S "$study/releases/$release/audit_policy.py" "$out" Python
unset LD_PRELOAD
exit "$code"
