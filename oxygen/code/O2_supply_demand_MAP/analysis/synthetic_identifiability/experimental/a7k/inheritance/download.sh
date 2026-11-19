#!/bin/bash
set -euo pipefail
out=oxygen/results/synthetic_identifiability/20261119_A7K_inheritance_sensitivity_01
export OPENBLAS_NUM_THREADS=1
py=/Users/4482173/.pyenv/versions/3.13.0/bin/python3
remote_host=4482173@red.moffitt.org
remote_dir=$(cat "$out/hpc_directory.txt")
remote_code="$remote_dir/StructuralAndPracticalIdentifiability/$out"
mode="$1"; phase="$2"
entry=$(awk -F '\t' -v m="$mode" -v p="$phase" '$1==m && $2==p {print $3, $4, $5}' "$out/hpc_jobs.tsv" | tail -1)
read -r N job width <<< "$entry"
job=${job%%;*}
expected=$(( (N+width-1)/width ))
receipt="$out/slurm_${mode}_${phase}.tsv"
ssh -o BatchMode=yes "$remote_host" "/cm/shared/apps/slurm/current/bin/sacct -n -X -j '$job' --format=JobID,State,ExitCode -P" > "$receipt"
awk -F '|' -v n="$expected" '$2=="COMPLETED" && $3=="0:0" {count++} END {exit count != n}' "$receipt" || exit 4
data="batch_${mode}_${phase}"
ssh -o BatchMode=yes "$remote_host" "cd '$remote_code' && tar -cf - '$data'/part_*.tsv" | tar -C "$out" -xf -
"$py" "$out/recovery_tilts.py" collect "$mode" "$phase"
