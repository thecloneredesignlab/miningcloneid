#!/bin/bash
set -euo pipefail
out=oxygen/results/synthetic_identifiability/20261119_experimental_death_A7K_comparison_01
export OPENBLAS_NUM_THREADS=1
py=/Users/4482173/.pyenv/versions/3.13.0/bin/python3
remote_host=4482173@red.moffitt.org
remote_dir=$(cat "$out/hpc_directory.txt")
remote_code="$remote_dir/StructuralAndPracticalIdentifiability/$out"
mode="$1"
phase="$2"
N=$("$py" "$out/array_run.py" freeze "$mode" "$phase")
width=4
last=$(( (N+width-1)/width-1 ))
scp -q "$out/array_run.py" "$remote_host:$remote_code/array_run.py"
tar --no-xattrs -C "$out" -cf - "arrays_${mode}_${phase}/frozen.npz" |   ssh -o BatchMode=yes "$remote_host" "tar -xf - -C '$remote_code'"
job=$(ssh -o BatchMode=yes "$remote_host" "cd '$remote_dir' && /cm/shared/apps/slurm/current/bin/sbatch --parsable --job-name=culture-${mode#stock_}-${phase} --cpus-per-task=2 --mem=2G --time=00:35:00 --array=0-${last}%40 --output='$remote_dir/slurm-%A_%a.out' '$remote_code/array.sh' '$remote_dir' '$mode' '$phase' '$width'")
printf '%s\t%s\t%s\t%s\n' "$mode" "$phase" "$N" "$job" | tee -a "$out/hpc_jobs.tsv"
