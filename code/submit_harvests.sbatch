#!/bin/bash
# ──────────────────────────────────────────────────────────────────────────
# submit_harvests.sbatch
# ──────────────────────────────────────────────────────────────────────────
# Slurm job array: one task = one harvest. Slurm will spread these across
# as many nodes as the partition allows, running multiple in parallel.
#
# This script requests 32 cores per task because N_CHAINS = 32 in
# mcmc_fit.py. On a 64-core node Slurm will pack 2 tasks per node.
#
# Submit with an array size equal to the number of harvests in harvests.txt:
#
#     python submit_harvests.py            # writes harvests.txt, prints cmd
#     sbatch --array=0-49%25 submit_harvests.sbatch    # e.g. 50 harvests, 25 at once
#
# The %N suffix caps simultaneous running tasks — useful if you want to
# stay polite on a shared cluster.
# ──────────────────────────────────────────────────────────────────────────
#SBATCH --job-name=harvest_mcmc
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=12G
#SBATCH --time=600:00:00
#SBATCH --output=harvest_zips/logs/slurm_%A_%a.out
#SBATCH --error=harvest_zips/logs/slurm_%A_%a.out
#SBATCH --partition=red
#SBATCH --qos=normal
# ────────────────────────────────────────────────────────────────────────

set -euo pipefail

mkdir -p harvest_zips/logs

# ── Activate the personal miniconda env ─────────────────────────────────
# Sourcing the user's own install directly is more robust than `module load`
# because it doesn't depend on the cluster's module naming.
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate myenv
# ────────────────────────────────────────────────────────────────────────

# Pick the harvest for this array task (1-indexed into harvests.txt)
HARVEST=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" harvests.txt)
if [ -z "${HARVEST:-}" ]; then
    echo "ERROR: no harvest at array index $SLURM_ARRAY_TASK_ID in harvests.txt" >&2
    exit 1
fi

echo "[$(date '+%F %T')] array=$SLURM_ARRAY_JOB_ID task=$SLURM_ARRAY_TASK_ID harvest=$HARVEST host=$(hostname) cores=$SLURM_CPUS_PER_TASK python=$(which python)"

python -u run_one_harvest.py "$HARVEST" --cores "$SLURM_CPUS_PER_TASK"
