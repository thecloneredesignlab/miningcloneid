#!/bin/bash -l
set -eu
study=$1
task=$2
out="$study/verification/$task"
mkdir -p "$out/tmp/home"
test "${SLURM_CPUS_PER_TASK}" = 1
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MAKEFLAGS=-j1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 SYMPY_GROUND_TYPES=python
py=/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3
"$py" -B -S "$study/release/verify_source.py" "$study" > "$out/source_before.txt"
if [ "$task" = P01 ]; then
  /usr/bin/gcc -shared -fPIC -O2 -Wall -Wextra "$study/release/io_guard.c" -ldl -o "$out/guard_b6.so" > "$out/guard_compile.txt" 2>&1
  /usr/bin/gcc --version > "$out/compiler_version.txt"
  chmod 555 "$out/guard_b6.so"
  guard_sha=$(sha256sum "$out/guard_b6.so")
  printf 'sha256\n%s\n' "${guard_sha%% *}" > "$out/guard_build.tsv"
fi
guard="$study/verification/P01/guard_b6.so"
test -f "$guard"
guard_sha=$(sha256sum "$guard")
recorded_guard=$(tail -n 1 "$study/verification/P01/guard_build.tsv")
test "${guard_sha%% *}" = "$recorded_guard"
set +e
if [ "$task" = P03 ]; then
  image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
  actual=$(sha256sum "$image")
  actual=${actual%% *}
  if [ "$actual" != dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372 ]; then exit 78; fi
  printf '%s\n' "$actual" > "$out/image_sha256.txt"
  apptainer exec --cleanenv --home "$out/tmp/home" --bind /share:/share \
   --bind "$study/inputs:$study/inputs:ro" --bind "$study/runtime:$study/runtime:ro" --bind "$study/release:$study/release:ro" \
   --bind /share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/20261116_culture_step04_stageB_01:/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/20261116_culture_step04_stageB_01:ro \
   --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 --env MAKEFLAGS=-j1 \
   --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null --env TMPDIR="$out/tmp" \
   "$image" "$study/runtime/trace_bundle/ld-linux-x86-64.so.2" --library-path "$study/runtime/trace_bundle" \
   "$study/runtime/trace_bundle/strace" -ff -s 1024 -e trace=%file,%process -o "$out/trace" \
   env LD_PRELOAD="$guard" CULTURE_GUARD_LOG="$out/blocked_access.tsv" CULTURE_ALLOWED_WRITE_ROOT="$out" \
   /usr/local/bin/o2-hpc-exact-rscript --vanilla "$study/release/passage.R" "$study" "$task" > "$out/console.log" 2>&1
  app=$?
  kind=R
else
  if [ "$task" = P01 ]; then program=preflight.py; else program=review.py; fi
  /usr/bin/strace -ff -s 1024 -e trace=%file,%process -o "$out/trace" \
   env LD_PRELOAD="$guard" CULTURE_GUARD_LOG="$out/blocked_access.tsv" CULTURE_ALLOWED_WRITE_ROOT="$out" TMPDIR="$out/tmp" \
   "$py" -B -S "$study/release/$program" "$study" "$task" > "$out/console.log" 2>&1
  app=$?
  kind=Python
fi
"$py" -B -S "$study/release/policy.py" "$out" "$kind" > "$out/policy_console.txt" 2>&1
policy=$?
"$py" -B -S "$study/release/verify_source.py" "$study" > "$out/source_after.txt" 2>&1
source_after=$?
set -e
printf 'application_exit_code\tpolicy_exit_code\tsource_after_exit_code\n%s\t%s\t%s\n' "$app" "$policy" "$source_after" > "$out/phase_exit_codes.tsv"
if [ "$app" != 0 ]; then exit "$app"; fi
if [ "$policy" != 0 ]; then exit "$policy"; fi
exit "$source_after"
