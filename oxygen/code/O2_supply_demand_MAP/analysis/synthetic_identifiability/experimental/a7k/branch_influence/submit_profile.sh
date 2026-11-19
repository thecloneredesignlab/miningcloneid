#!/bin/bash
set -euo pipefail
root=oxygen/results/synthetic_identifiability
d="$root/20261119_A7K_branch_influence_01"
prev="$root/20261119_A7K_inheritance_sensitivity_01"
prior="$root/20261119_experimental_death_A7K_comparison_01"
host=4482173@red.moffitt.org
remote=$(cat "$prev/hpc_directory.txt")
remote_repo="$remote/StructuralAndPracticalIdentifiability"
folder="$remote_repo/$d"
phase="$1"
list="profile_tasks_$phase.tsv"
shift
if [[ "$phase" == individual ]]; then
    ssh -o BatchMode=yes "$host" "mkdir -p '$folder'"
    scp -q "$prev"/best_strict.tsv "$prev"/best_group_050.tsv "$prev"/profile_points_group_050.tsv "$host:$remote_repo/$prev/"
fi
: > "$d/$list"
files=()
for m in strict group_050; do
    for k in "$@"; do
        printf '%s\t%s\n' "$m" "$k" >> "$d/$list"
        printf -v stem 'fit_ref_%s_%02d.tsv' "$m" "$k"
        files+=("$d/$stem")
    done
done
scp -q "$d/branch_influence.py" "$d/profile_array.sh" "$d/$list" "${files[@]}" "$host:$folder/"
n=$(wc -l < "$d/$list")
job=$(ssh -o BatchMode=yes "$host" "/cm/shared/apps/slurm/current/bin/sbatch --parsable --job-name=A7K-branch-$phase --cpus-per-task=1 --mem=2G --time=00:45:00 --array=0-$((n-1))%16 --output='$remote/a7k-branch-%A_%a.out' '$folder/profile_array.sh' '$remote' '$folder' '$list'")
printf '%s\t%s\t%s\n' "$phase" "$n" "$job" | tee -a "$d/profile_jobs.tsv"
