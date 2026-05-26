"""
run_one_harvest.py
──────────────────
Runs a single harvest end-to-end (patch script -> run -> zip results).

This is the per-task worker for multi-node Slurm job arrays. Each array task
picks up one harvest name and processes it independently, so Slurm can spread
the work across as many nodes as the cluster has available. No MPI — the
parallelism is at the harvest level, not inside the simulation.

Usage:
    python run_one_harvest.py HARVEST_NAME [--cores N] [--no-skip]

Env vars consulted:
    SLURM_CPUS_PER_TASK   - overrides --cores if --cores not given
    SLURM_ARRAY_JOB_ID    - used to make tmp-script names unique per submission
    SLURM_ARRAY_TASK_ID   - used to make tmp-script names unique per task
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import zipfile


BEAM_SCRIPT = "beam_search_flip_rate_wgd.py"
OUTPUT_DIR = "harvest_zips"


def patch_script(script_path: str, harvest_name: str, tmp_path: str,
                 results_dir: str) -> None:
    """Copy beam-search script, patching SAMPLE_NAME and RESULTS_DIR in place."""
    with open(script_path) as f:
        src = f.read()

    src = re.sub(
        r'^SAMPLE_NAME\s*=\s*"[^"]*"',
        f'SAMPLE_NAME = "{harvest_name}"',
        src, count=1, flags=re.MULTILINE,
    )
    src = re.sub(
        r'^(\s*)RESULTS_DIR\s*=\s*"[^"]*"',
        rf'\1RESULTS_DIR = "{results_dir}"',
        src, count=1, flags=re.MULTILINE,
    )

    with open(tmp_path, "w") as f:
        f.write(src)


def zip_directory(folder: str, zip_path: str) -> None:
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for root, _dirs, files in os.walk(folder):
            for fname in files:
                full = os.path.join(root, fname)
                arcname = os.path.relpath(full, folder)
                zf.write(full, arcname)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("harvest", help="Harvest name (SAMPLE_NAME) to run")
    parser.add_argument("--cores", type=int, default=0,
                        help="Cores available to this harvest. "
                             "Default 0 = read from $SLURM_CPUS_PER_TASK or "
                             "fall back to 32.")
    parser.add_argument("--no-skip", dest="skip_completed",
                        action="store_false", default=True,
                        help="Re-run even if a zip already exists for this "
                             "harvest (default: skip completed).")
    args = parser.parse_args()

    harvest = args.harvest

    # Resolve core count: CLI flag > Slurm env > 32
    cores = args.cores
    if cores <= 0:
        cores = int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or 32

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    log_dir = os.path.join(OUTPUT_DIR, "logs")
    os.makedirs(log_dir, exist_ok=True)

    zip_name = os.path.join(OUTPUT_DIR, f"{harvest}.zip")
    if args.skip_completed and os.path.isfile(zip_name):
        print(f"⊘ {harvest} (already done, skipping)")
        return 0

    # Unique tmp-script name so concurrent array tasks on the same node,
    # or re-submissions, never overwrite each other.
    job_id = os.environ.get("SLURM_ARRAY_JOB_ID", str(os.getpid()))
    task_id = os.environ.get("SLURM_ARRAY_TASK_ID", "local")
    results_dir = f"results_{harvest}"
    tmp_script = f"_tmp_beam_{harvest}_{job_id}_{task_id}.py"

    if os.path.isdir(results_dir):
        shutil.rmtree(results_dir)

    patch_script(BEAM_SCRIPT, harvest, tmp_script, results_dir)

    env = os.environ.copy()
    # Cap thread-pool libraries so this task stays within its Slurm allocation
    env["OMP_NUM_THREADS"] = str(cores)
    env["MKL_NUM_THREADS"] = str(cores)
    env["OPENBLAS_NUM_THREADS"] = str(cores)
    env["NUMEXPR_NUM_THREADS"] = str(cores)
    # Read by mcmc_fit.py to size the inner multiprocessing.Pool
    env["HARVEST_MAX_WORKERS"] = str(cores)
    env["MPLBACKEND"] = "Agg"
    env["PYTORCH_ALLOC_CONF"] = "expandable_segments:True,max_split_size_mb:256"

    log_path = os.path.join(log_dir, f"{harvest}.log")
    print(f"▶ {harvest} (cores={cores}, log={log_path})", flush=True)

    try:
        with open(log_path, "w") as log_fh:
            result = subprocess.run(
                [sys.executable, "-u", tmp_script],
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                env=env,
            )
            rc = result.returncode
    except Exception as e:
        print(f"✗ {harvest} — subprocess error: {e}")
        return 1
    finally:
        try:
            os.remove(tmp_script)
        except OSError:
            pass

    if rc != 0:
        reason = "OOM killed" if rc == -9 else f"exit code {rc}"
        print(f"✗ {harvest} — {reason}")
        return rc if rc > 0 else 1

    if not os.path.isdir(results_dir):
        print(f"✗ {harvest} — no results dir produced")
        return 1

    zip_directory(results_dir, zip_name)
    shutil.rmtree(results_dir)
    print(f"✓ {harvest} -> {zip_name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
