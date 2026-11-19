#!/bin/bash
set -euo pipefail
d=oxygen/results/synthetic_identifiability/20261119_A7K_4N_branch_swap_01
p=/Users/4482173/.pyenv/versions/3.13.0/bin/python3
export OPENBLAS_NUM_THREADS=1
stage() {
    local variant="$1" mode="$2" phase="$3"
    shift 3
    local log="$d/${variant}_${mode}_$phase.log"
    if ! bash "$d/stage.sh" "$variant" "$mode" "$phase" "$@" > "$log" 2>&1; then
        tail -20 "$log"
        return 2
    fi
    tail -8 "$log"
}
for mode in group_050 group_025; do
    for phase in pilot expanded; do
        stage shared "$mode" "$phase" up down
    done
done
"$p" "$d/branch_swap.py" decide
for mode in group_050 group_025; do
    selected=$(awk -F '\t' -v m="$mode" '$1==m && $5=="True" {print $2}' "$d/independent_fit_gate.tsv")
    if [[ -n "$selected" ]]; then
        for phase in pilot expanded; do
            stage split4N "$mode" "$phase" $selected
        done
    fi
done
