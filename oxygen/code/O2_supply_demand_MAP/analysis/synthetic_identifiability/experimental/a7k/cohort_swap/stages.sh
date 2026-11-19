#!/bin/bash
set -euo pipefail
base=oxygen/results/synthetic_identifiability
d="$base/20261119_A7K_cohort_swap_01"
previous="$base/20261119_A7K_inheritance_sensitivity_01"
py=/Users/4482173/.pyenv/versions/3.13.0/bin/python3
export OPENBLAS_NUM_THREADS=1
host=4482173@red.moffitt.org
remote=$(cat "$previous/hpc_directory.txt")
dest="$remote/StructuralAndPracticalIdentifiability/$d"
mode="$1"
phase="$2"
width=8
batch="batch_${mode}_$phase"
manifest="$d/$batch/manifest.tsv"
if [[ "$mode:$phase" == group_050:pilot ]]; then
    ssh -o BatchMode=yes "$host" "mkdir -p '$dest'"
    scp -q "$d/cohort_swap.py" "$d/cohort_array.sh" "$host:$dest/"
fi
if [[ ! -f "$d/$batch/frozen.npz" ]]; then
    "$py" "$d/cohort_swap.py" prepare "$mode" "$phase"
fi
n=$(($(wc -l < "$manifest")-1))
tasks=$(((n+width-1)/width))
receipt="$d/slurm_${mode}_$phase.tsv"
tar --no-xattrs -C "$d" -cf - "$batch/manifest.tsv" "$batch/frozen.npz" |
    ssh -o BatchMode=yes "$host" "tar -xf - -C '$dest'"
job=$(ssh -o BatchMode=yes "$host" "/cm/shared/apps/slurm/current/bin/sbatch --parsable --job-name=A7K-cohort-${mode#group_}-$phase --cpus-per-task=4 --mem=3G --time=00:35:00 --array=0-$((tasks-1))%20 --output='$remote/a7k-cohort-%A_%a.out' '$dest/cohort_array.sh' '$remote' '$dest' '$mode' '$phase' '$width'")
printf '%s\t%s\t%s\t%s\n' "$mode" "$phase" "$n" "$job" | tee -a "$d/slurm_jobs.tsv"
job=${job%%;*}
for ((;;)); do
    if ssh -o BatchMode=yes -o ConnectTimeout=15 "$host" "/cm/shared/apps/slurm/current/bin/sacct -n -X -j '$job' --format JobID,State,ExitCode -P" 2>> "$d/poll.log" > "$receipt"; then
        if awk -F '|' -v n="$tasks" '$2=="COMPLETED" && $3=="0:0" {done++} END {exit done!=n}' "$receipt"; then
            break
        fi
        if awk -F '|' '$2 ~ /^(FAILED|TIMEOUT|CANCELLED|OUT_OF_MEMORY|NODE_FAIL)/ {bad=1} END {exit !bad}' "$receipt"; then
            cat "$receipt"
            exit 5
        fi
    fi
    sleep 22
done
ssh -o BatchMode=yes "$host" "cd '$dest' && tar -cf - '$batch'/part_*.tsv" | tar -xf - -C "$d"
"$py" "$d/cohort_swap.py" collect "$mode" "$phase"
"$py" "$d/cohort_swap.py" repair "$mode" "$phase"
if [[ "$phase" == expanded ]]; then
    "$py" "$d/crossfit_audit.py" "$mode"
    "$py" "$d/cohort_swap.py" summarize "$mode"
else
    "$py" - "$d" "$mode" <<'PY'
import csv,sys
from pathlib import Path
p=Path(sys.argv[1])
mode=sys.argv[2]
rows=list(csv.DictReader((p/f"fits_{mode}_pilot.tsv").open(),delimiter="\t"))
for g in ("2N","4N"):
    for j in ("1","2"):
        sub=[r for r in rows if r["branch"]==g and r["case_id"]==j]
        stable=sum(r["validated"]=="True" for r in sub)
        assert len(sub)==20 and stable>=18
        print(mode,g,j,"pilot",stable,sum(r["recovered"]=="True" for r in sub),flush=True)
PY
fi
