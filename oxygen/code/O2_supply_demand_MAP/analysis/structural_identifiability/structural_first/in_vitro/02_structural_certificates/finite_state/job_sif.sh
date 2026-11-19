#!/bin/bash -l
set -eu
study=$1
release=$2
phase=$3
program=${4:-semantic_probe.R}
project=/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability
image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
out="$study/verification/$phase"
mkdir -p "$out/tmp/home"
test "${SLURM_CPUS_PER_TASK}" = 1
actual=$(sha256sum "$image")
actual=${actual%% *}
test "$actual" = dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372
printf '%s\n' "$actual" > "$out/image_sha256.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MAKEFLAGS=-j1
set +e
if [ "$phase" = B0 ]
then
  set -- "$project" "$study"
else
  set -- "$project" "$study" "$phase"
fi
apptainer exec --cleanenv --home "$out/tmp/home" --bind /share:/share \
 --bind "$project/oxygen/results/structural_identifiability/20261109_seed144_wp0:$project/oxygen/results/structural_identifiability/20261109_seed144_wp0:ro" \
 --bind "$study/inputs:$study/inputs:ro" --bind "$study/releases:$study/releases:ro" \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 --env MAKEFLAGS=-j1 \
 --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null --env TMPDIR="$out/tmp" \
 "$image" "$study/runtime/trace_bundle/ld-linux-x86-64.so.2" \
 --library-path "$study/runtime/trace_bundle" "$study/runtime/trace_bundle/strace" \
 -ff -s 1024 -e trace=%file,%process -o "$out/trace" \
 env LD_PRELOAD="$study/runtime/io_guard.so" CULTURE_GUARD_LOG="$out/blocked_access.tsv" \
 /usr/local/bin/o2-hpc-exact-rscript --vanilla "$study/releases/$release/$program" "$@" > "$out/console.log" 2>&1
code=$?
set -e
printf 'application_exit_code\t%s\n' "$code" > "$out/application_status.tsv"
exit "$code"
