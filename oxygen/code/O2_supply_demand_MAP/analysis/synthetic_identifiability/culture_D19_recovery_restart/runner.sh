#!/bin/bash -l
set -eu
study=$1
task=$2
out="$study/verification/$task"
mkdir -p "$out/tmp/home"
chmod 1777 "$out/tmp"
test "${SLURM_CPUS_PER_TASK}" = 1
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MAKEFLAGS=-j1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
py=/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3
"$py" -B -S "$study/release/verify_source.py" "$study" > "$out/source_before.txt"
image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
actual=$(sha256sum "$image");actual=${actual%% *}
test "$actual" = dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372
printf '%s\n' "$actual" > "$out/image_sha256.txt"
extra_guard=''
case "$task" in
 D_*|F_*|T02D|T03F)
  extra_guard=":$study/verification/T00/read_guard.so"
  export APPTAINERENV_REC_STUDY_ROOT="$study" APPTAINERENV_REC_READ_DENIAL_LOG="$out/read_denials.tsv"
  case "$task" in
   D_*) export APPTAINERENV_REC_MODE=D APPTAINERENV_REC_ORIGIN="${task:0:3}";;
   F_*) export APPTAINERENV_REC_MODE=F APPTAINERENV_REC_ORIGIN="${task:0:3}";;
   T02D) export APPTAINERENV_REC_MODE=D APPTAINERENV_REC_ORIGIN=D_A;;
   T03F) export APPTAINERENV_REC_MODE=F APPTAINERENV_REC_ORIGIN=F_E;;
  esac;;
esac
cd "$out"
set +e
apptainer exec --cleanenv --pwd "$out" --home "$out/tmp/home" --bind /share:/share \
 --bind "$out/tmp:$out/tmp:rw" --bind "$study/inputs:$study/inputs:ro" --bind "$study/runtime:$study/runtime:ro" --bind "$study/release:$study/release:ro" \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 --env MAKEFLAGS=-j1 \
 --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null --env TMPDIR="$out/tmp" \
 "$image" "$study/runtime/trace_bundle/ld-linux-x86-64.so.2" --library-path "$study/runtime/trace_bundle" \
 "$study/runtime/trace_bundle/strace" -ff -y -s 1024 -e trace=%file,%process,close,fchdir -o "$out/trace" \
 env LD_PRELOAD="$study/runtime/guard.so${extra_guard}" CULTURE_GUARD_LOG="$out/blocked_access.tsv" CULTURE_ALLOWED_WRITE_ROOT="$out" \
 /usr/local/bin/o2-hpc-exact-rscript --vanilla "$study/release/main.R" "$study" "$task" > "$out/console.log" 2>&1
app=$?
"$py" -B -S "$study/release/policy.py" "$out" R > "$out/policy_console.txt" 2>&1
policy=$?
"$py" -B -S "$study/release/verify_source.py" "$study" > "$out/source_after.txt" 2>&1
source_after=$?
set -e
printf 'application_exit_code\tpolicy_exit_code\tsource_after_exit_code\n%s\t%s\t%s\n' "$app" "$policy" "$source_after" > "$out/phase_exit_codes.tsv"
"$py" -B -S "$study/release/seal_phase.py" "$out" > "$out/collection_console.txt"
if [ "$app" != 0 ]; then exit "$app"; fi
if [ "$policy" != 0 ]; then exit "$policy"; fi
exit "$source_after"
