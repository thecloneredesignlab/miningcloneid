#!/bin/bash
set -euo pipefail
d=oxygen/results/synthetic_identifiability/20261119_A7K_branch_influence_01
host=4482173@red.moffitt.org
phase="$1"
read -r label n job width <<< "$(awk -F '\t' -v phase="$phase" '$1==phase' "$d/diagnostic_jobs.tsv" | tail -1)"
job=${job%%;*}
expected=$(((n+width-1)/width))
snapshot="$d/snapshot_$phase.tsv"
for ((;;)); do
    if ssh -o BatchMode=yes -o ConnectTimeout=15 "$host" "/cm/shared/apps/slurm/current/bin/sacct -n -X -j '$job' --format JobID,State,ExitCode -P" > "$snapshot"; then
        if awk -F '|' -v n="$expected" '$2=="COMPLETED" && $3=="0:0" {done++} END {exit done!=n}' "$snapshot"; then
            break
        fi
        if awk -F '|' '$2 ~ /^(FAILED|TIMEOUT|CANCELLED|OUT_OF_MEMORY|NODE_FAIL)/ {bad=1} END {exit !bad}' "$snapshot"; then
            awk -F '|' '$2 ~ /^(FAILED|TIMEOUT|CANCELLED|OUT_OF_MEMORY|NODE_FAIL)/ { print }' "$snapshot"
            exit 5
        fi
    fi
    sleep 25
done
bash "$d/download_diagnostic.sh" "$phase"
