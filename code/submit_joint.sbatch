#!/bin/bash
# ────────────────────────────────────────────────────────────────────────────
# submit_joint.sbatch
# ────────────────────────────────────────────────────────────────────────────
# Slurm SINGLE JOB (not an array) for the joint hierarchical fit.
#
# All mice are fitted simultaneously in one MCMC run, so only one job is
# needed regardless of how many harvests exist.
#
# Submit with:
#     sbatch submit_joint.sbatch
#
# Requires JOINT_FIT_MODE = True in beam_search_flip_rate_wgd.py.
# Do NOT use submit_harvests.sbatch / submit_harvests.py for joint mode.
# ────────────────────────────────────────────────────────────────────────────
#SBATCH --job-name=harvest_mcmc_joint
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-0:00:00
#SBATCH --output=harvest_zips/logs/slurm_joint_%j.out
#SBATCH --error=harvest_zips/logs/slurm_joint_%j.out
#SBATCH --partition=red
#SBATCH --qos=normal
# ────────────────────────────────────────────────────────────────────────────

set -euo pipefail

mkdir -p harvest_zips/logs

# ── Activate conda environment ────────────────────────────────────────────────
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate myenv
# ─────────────────────────────────────────────────────────────────────────────

echo "[$(date '+%F %T')] joint fit  host=$(hostname)  cores=$SLURM_CPUS_PER_TASK  python=$(which python)"

python -u run_joint.py --cores "$SLURM_CPUS_PER_TASK"
