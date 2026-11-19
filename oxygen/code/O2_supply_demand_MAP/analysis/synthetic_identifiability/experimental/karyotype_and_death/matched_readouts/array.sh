#!/bin/bash -l
set -euo pipefail
module load Anaconda3/2024.02-1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export CULTURE_DEATH_A7K_LIB="$1/native.so"
code="$1/StructuralAndPracticalIdentifiability/oxygen/results/synthetic_identifiability/20261119_experimental_death_A7K_comparison_01"
python "$code/array_run.py" run "$2" "$3" --shard "$SLURM_ARRAY_TASK_ID" --width "$4" --cores "$SLURM_CPUS_PER_TASK"
