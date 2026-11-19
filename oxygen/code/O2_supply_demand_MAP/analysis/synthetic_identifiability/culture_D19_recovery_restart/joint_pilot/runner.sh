#!/bin/bash -l
set -eu
archive=$1
pilot=$2
task=$3
release=$4
out="$pilot/verification/$task"
test "${SLURM_CPUS_PER_TASK}" = 1
mkdir "$out"
mkdir -p "$out/tmp/home"
chmod 1777 "$out/tmp"
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MAKEFLAGS=-j1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
py=/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3
image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
verify() {
  "$py" -B -S "$release/admin.py" verify "$archive" "$pilot"
}
verify > "$out/source_before.txt"
actual=$(sha256sum "$image")
actual=${actual%% *}
test "$actual" = dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372
printf '%s\n' "$actual" > "$out/image_sha256.txt"
cd "$out"
set +e
apptainer exec --cleanenv --pwd "$out" --home "$out/tmp/home" --bind /share:/share \
  --bind "$archive:$archive:ro" --bind "$pilot/release:$pilot/release:ro" --bind "$release:$release:ro" --bind "$out/tmp:$out/tmp:rw" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 --env MAKEFLAGS=-j1 \
  --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null --env TMPDIR="$out/tmp" \
  "$image" "$archive/runtime/trace_bundle/ld-linux-x86-64.so.2" --library-path "$archive/runtime/trace_bundle" \
  "$archive/runtime/trace_bundle/strace" -ff -y -s 1024 -e trace=%file,%process,close,fchdir -o "$out/trace" \
  env LD_PRELOAD="$archive/runtime/guard.so" CULTURE_GUARD_LOG="$out/blocked_access.tsv" CULTURE_ALLOWED_WRITE_ROOT="$out" \
  /usr/local/bin/o2-hpc-exact-rscript --vanilla "$release/run.R" "$archive" "$pilot" "$task" "$release" > "$out/console.log" 2>&1
app=$?
"$py" -B -S "$archive/release/policy.py" "$out" R > "$out/policy_console.txt" 2>&1
policy=$?
"$py" -B -S "$release/admin.py" audit-reads "$archive" "$pilot" "$task" > "$out/read_audit_console.txt" 2>&1
reads=$?
verify > "$out/source_after.txt" 2>&1
source=$?
set -e
printf 'application_exit_code\tpolicy_exit_code\tread_audit_exit_code\tsource_after_exit_code\n%s\t%s\t%s\t%s\n' "$app" "$policy" "$reads" "$source" > "$out/phase_exit_codes.tsv"
"$py" -B -S "$archive/release/seal_phase.py" "$out" > "$out/collection_console.txt"
test "$app" = 0
test "$policy" = 0
test "$reads" = 0
test "$source" = 0
