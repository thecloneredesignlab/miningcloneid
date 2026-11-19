#!/bin/bash
set -euo pipefail
root=oxygen/results/synthetic_identifiability
d="$root/20261119_A7K_branch_influence_01"
prev="$root/20261119_A7K_inheritance_sensitivity_01"
host=4482173@red.moffitt.org
remote=$(cat "$prev/hpc_directory.txt")
folder="$remote/StructuralAndPracticalIdentifiability/$d"
phase="$1"
read -r label n job width <<< "$(awk -F '\t' -v phase="$phase" '$1==phase' "$d/diagnostic_jobs.tsv" | tail -1)"
job=${job%%;*}
expected=$(((n+width-1)/width))
receipt="$d/slurm_diagnostic_$phase.tsv"
ssh -o BatchMode=yes "$host" "/cm/shared/apps/slurm/current/bin/sacct -n -X -j '$job' --format JobID,State,ExitCode -P" > "$receipt"
awk -F '|' -v n="$expected" '$2=="COMPLETED" && $3=="0:0" {count++} END {exit count!=n}' "$receipt" || exit 4
ssh -o BatchMode=yes "$host" "cd '$folder' && tar -cf - diag_$phase/part_*.tsv" | tar -xf - -C "$d"
export OPENBLAS_NUM_THREADS=1
/Users/4482173/.pyenv/versions/3.13.0/bin/python3 "$d/recovery_diagnostics.py" collect "$phase"
