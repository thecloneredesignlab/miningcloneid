#!/bin/bash
# Runtime audit only; called within its own finite xxlarge allocation.
set -eu
study=$1
cd "$study"
mkdir -p runtime/matlab/preferences runtime/matlab/tmp runtime/matlab/work
gcc -shared -fPIC -O2 -Wall -Wextra releases/r01/io_guard.c -ldl -o runtime/io_guard.so
export CULTURE_GUARD_LOG="$study/runtime/matlab/blocked_access.tsv"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MATLAB_PREFDIR="$study/runtime/matlab/preferences"
export TMPDIR="$study/runtime/matlab/tmp"
export PYTHONDONTWRITEBYTECODE=1
export LD_PRELOAD="$study/runtime/io_guard.so"
# Verify that the policy guard denies a forbidden file before testing MATLAB.
/usr/bin/python3 -B -S -c 'import os; p=os.environ["TMPDIR"]+"/guard_negative.json"; failed=False
try: os.open(p,os.O_CREAT|os.O_WRONLY,0o600)
except PermissionError: failed=True
assert failed and not os.path.exists(p)
print("NATIVE_IO_GUARD_NEGATIVE_CONTROL_PASSED")' > runtime/guard_control.log 2>&1
mv runtime/matlab/blocked_access.tsv runtime/guard_expected_blocks.tsv
cd runtime/matlab/work
set +e
timeout --kill-after=15s 180 /usr/bin/strace -ff -s 1024 -e trace=%file,%process -o "$study/runtime/matlab/trace" \
  /app/other/MATLAB/R2022b/bin/matlab -nojvm -nodisplay -nodesktop -singleCompThread \
  -batch "run('$study/releases/r01/runtime_matlab.m')" > "$study/runtime/matlab/console.log" 2>&1
code=$?
set -e
printf 'application_exit_code\t%s\n' "$code" > "$study/runtime/matlab/application_status.tsv"
unset LD_PRELOAD
printf 'Runtime probe ended; inspect process/file traces before any application control.\n'
