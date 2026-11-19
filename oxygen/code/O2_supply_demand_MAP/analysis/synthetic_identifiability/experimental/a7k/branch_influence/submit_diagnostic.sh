#!/bin/bash
set -euo pipefail
root=oxygen/results/synthetic_identifiability
d="$root/20261119_A7K_branch_influence_01"
prev="$root/20261119_A7K_inheritance_sensitivity_01"
host=4482173@red.moffitt.org
remote=$(cat "$prev/hpc_directory.txt")
folder="$remote/StructuralAndPracticalIdentifiability/$d"
phase="$1"
dependency=""
if [[ ${2:-} =~ ^[0-9]+$ ]]; then
    dependency="--dependency=afterany:$2"
fi
width=8
staged="$d/diag_$phase"
n=$(($(wc -l < "$staged/manifest.tsv")-1))
if [[ "$phase" == oracle ]]; then
    scp -q "$d/recovery_diagnostics.py" "$d/diagnostic_array.sh" "$host:$folder/"
fi
tar --no-xattrs -C "$d" -cf - "diag_$phase/manifest.tsv" "diag_$phase/frozen.npz" |
    ssh -o BatchMode=yes "$host" "tar -xf - -C '$folder'"
last=$(((n+width-1)/width-1))
job=$(ssh -o BatchMode=yes "$host" "/cm/shared/apps/slurm/current/bin/sbatch --parsable $dependency --job-name=A7K-diag-$phase --cpus-per-task=4 --mem=3G --time=00:40:00 --array=0-$last%20 --output='$remote/a7k-diag-%A_%a.out' '$folder/diagnostic_array.sh' '$remote' '$folder' '$phase' '$width'")
printf '%s\t%s\t%s\t%s\n' "$phase" "$n" "$job" "$width" | tee -a "$d/diagnostic_jobs.tsv"
