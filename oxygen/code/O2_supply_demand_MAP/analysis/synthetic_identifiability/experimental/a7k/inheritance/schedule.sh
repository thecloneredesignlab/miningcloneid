#!/bin/bash
set -euo pipefail
out=oxygen/results/synthetic_identifiability/20261119_A7K_inheritance_sensitivity_01
export OPENBLAS_NUM_THREADS=1
py=/Users/4482173/.pyenv/versions/3.13.0/bin/python3
remote_host=4482173@red.moffitt.org
remote_dir=$(cat "$out/hpc_directory.txt")
remote_code="$remote_dir/StructuralAndPracticalIdentifiability/$out"
mode="$1"; phase="$2"
shard_size=8
N=$("$py" "$out/recovery_tilts.py" prepare "$mode" "$phase")
if ((N == 0)); then
    printf '%s\t%s\tuntouched\n' "$mode" "$phase"
    exit
fi
last=$(( (N+shard_size-1)/shard_size-1 ))
staged_code=$(ssh -o BatchMode=yes "$remote_host" "mktemp -d '$remote_code/.update.XXXXXX'")
scp -q "$out/recovery_tilts.py" "$out/tilt_model.py" "$out/array.sh" "$remote_host:$staged_code/"
ssh -o BatchMode=yes "$remote_host" "mv '$staged_code'/* '$remote_code/' && rmdir '$staged_code'"
tar --no-xattrs -C "$out" -cf - "batch_${mode}_${phase}/manifest.tsv" "batch_${mode}_${phase}/frozen.npz" | \
    ssh -o BatchMode=yes "$remote_host" "tar -xf - -C '$remote_code'"
job=$(ssh -o BatchMode=yes "$remote_host" "cd '$remote_dir' && /cm/shared/apps/slurm/current/bin/sbatch --parsable --job-name=A7K-${mode#group_}-$phase --cpus-per-task=4 --mem=3G --time=00:30:00 --array=0-$last%28 --output='$remote_dir/slurm-%A_%a.out' '$remote_code/array.sh' '$remote_dir' '$mode' '$phase' '$shard_size'")
printf '%s\t%s\t%s\t%s\t%s\n' "$mode" "$phase" "$N" "$job" "$shard_size" | tee -a "$out/hpc_jobs.tsv"
