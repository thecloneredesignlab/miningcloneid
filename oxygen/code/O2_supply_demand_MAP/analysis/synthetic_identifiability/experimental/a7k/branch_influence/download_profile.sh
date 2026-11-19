#!/bin/bash
set -euo pipefail
root=oxygen/results/synthetic_identifiability
d="$root/20261119_A7K_branch_influence_01"
prev="$root/20261119_A7K_inheritance_sensitivity_01"
host=4482173@red.moffitt.org
remote=$(cat "$prev/hpc_directory.txt")
folder="$remote/StructuralAndPracticalIdentifiability/$d"
phase="$1"
read -r label n job <<< "$(awk -F '\t' -v phase="$phase" '$1==phase' "$d/profile_jobs.tsv" | tail -1)"
job=${job%%;*}
receipt="$d/slurm_profile_$phase.tsv"
ssh -o BatchMode=yes "$host" "/cm/shared/apps/slurm/current/bin/sacct -n -X -j '$job' --format JobID,State,ExitCode -P" > "$receipt"
awk -F '|' -v n="$n" '$2=="COMPLETED" && $3=="0:0" {count++} END {exit count!=n}' "$receipt" || exit 4
ssh -o BatchMode=yes "$host" "cd '$folder' && tar -cf - interval_*.tsv profile_ref_*.tsv profile_points_*.tsv" | tar -xf - -C "$d"
