#!/bin/bash -l
set -euo pipefail
module load Anaconda3/2024.02-1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export CULTURE_A7K_TV_LIB="$1/native.so"
cd "$2"
python branch_swap.py run "$3" "$4" "$5" --shard "$SLURM_ARRAY_TASK_ID" --width "$6" --cores "$SLURM_CPUS_PER_TASK"
