#!/bin/bash
set -euo pipefail
root=oxygen/results/synthetic_identifiability
d="$root/20261119_A7K_4N_branch_swap_01"
prior="$root/20261119_A7K_inheritance_sensitivity_01"
py=/Users/4482173/.pyenv/versions/3.13.0/bin/python3
export OPENBLAS_NUM_THREADS=1
host=4482173@red.moffitt.org
remote=$(cat "$prior/hpc_directory.txt")
dest="$remote/StructuralAndPracticalIdentifiability/$d"
variant="$1"; mode="$2"; phase="$3"
shift 3
width=8
batch="batch_${variant}_${mode}_$phase"
if [[ "$variant:$mode:$phase" == shared:group_050:pilot ]]; then
    ssh -o BatchMode=yes "$host" "mkdir -p '$dest'"
    scp -q "$d/branch_swap.py" "$d/branch_array.sh" "$host:$dest/"
fi
if [[ ! -f "$d/$batch/frozen.npz" ]]; then
    "$py" "$d/branch_swap.py" prepare "$variant" "$mode" "$phase" --background "$@"
fi
n=$(($(wc -l < "$d/$batch/manifest.tsv")-1))
tasks=$(((n+width-1)/width))
terminal="$d/slurm_${variant}_${mode}_$phase.tsv"
tar --no-xattrs -C "$d" -cf - "$batch/manifest.tsv" "$batch/frozen.npz" |
    ssh -o BatchMode=yes "$host" "tar -xf - -C '$dest'"
job=$(ssh -o BatchMode=yes "$host" "/cm/shared/apps/slurm/current/bin/sbatch --parsable --job-name=A7K-4N-$variant-${mode#group_}-$phase --cpus-per-task=4 --mem=3G --time=00:40:00 --array=0-$((tasks-1))%20 --output='$remote/a7k-4N-%A_%a.out' '$dest/branch_array.sh' '$remote' '$dest' '$variant' '$mode' '$phase' '$width'")
printf '%s\t%s\t%s\t%s\t%s\n' "$variant" "$mode" "$phase" "$n" "$job" | tee -a "$d/slurm_jobs.tsv"
job=${job%%;*}
for ((;;)); do
    if ssh -o BatchMode=yes -o ConnectTimeout=15 "$host" "/cm/shared/apps/slurm/current/bin/sacct -n -X -j '$job' --format JobID,State,ExitCode -P" 2>> "$d/poll.log" > "$terminal"; then
        if awk -F '|' -v n="$tasks" '$2=="COMPLETED" && $3=="0:0" {done++} END {exit done!=n}' "$terminal"; then
            break
        fi
        if awk -F '|' '$2 ~ /^(FAILED|TIMEOUT|CANCELLED|OUT_OF_MEMORY|NODE_FAIL)/ {bad=1} END {exit !bad}' "$terminal"; then
            cat "$terminal"
            exit 5
        fi
    fi
    sleep 24
done
ssh -o BatchMode=yes "$host" "cd '$dest' && tar -cf - '$batch'/part_*.tsv" | tar -xf - -C "$d"
"$py" "$d/branch_swap.py" collect "$variant" "$mode" "$phase"
"$py" "$d/branch_swap.py" repair "$variant" "$mode" "$phase"
if [[ "$phase" == expanded ]]; then
    "$py" "$d/branch_swap.py" summarize "$variant" "$mode"
else
    "$py" - "$d" "$variant" "$mode" <<'PY'
import csv,sys
from pathlib import Path
d,v,m=sys.argv[1:]
rows=list(csv.DictReader((Path(d)/f"fits_{v}_{m}_pilot.tsv").open(),delimiter="\t"))
parts=sorted({(r["background"],r["cell"],r["case_id"]) for r in rows})
for key in parts:
    subset=[r for r in rows if tuple(r[k] for k in ("background","cell","case_id"))==key]
    valid=sum(r["validated"]=="True" for r in subset)
    assert len(subset)==20 and valid>=18
    print(v,m,*key,"pilot",valid,sum(r["recovered"]=="True" for r in subset),flush=True)
PY
fi
