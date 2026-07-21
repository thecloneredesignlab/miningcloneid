"""
run_all_harvests_parallel.py
────────────────────────────
Runs beam-search harvests in parallel, splitting CPU cores evenly across
concurrent harvest jobs.  Each subprocess is CPU-pinned so inner Pool /
ProcessPoolExecutor workers don't fight for the same cores.

Usage:
    python run_all_harvests_parallel.py              # auto-detect cores
    python run_all_harvests_parallel.py --workers 4   # run 4 harvests at once
"""

import argparse
import csv
import gc
import os
import re
import shutil
import subprocess
import sys
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

# ── Paths (edit if your layout differs) ──────────────────────────────────
CSV_PATH = "../data/InVivoData_Gemcitabine/harvest_ploidy_mapping.csv"
BEAM_SCRIPT = "beam_search_flip_rate_wgd.py"
OUTPUT_DIR = "harvest_zips"
RESULTS_DIR_PREFIX = "results"      # each harvest writes to results_<harvest>

# ── Tuning knobs ─────────────────────────────────────────────────────────
SKIP_COMPLETED = True                # skip harvests that already have a zip


def get_matching_harvests(csv_path: str) -> list[str]:
    harvests = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row["has_match"].strip() == "True":
                harvests.append(row["harvest"].strip())
    return harvests


def patch_script(script_path: str, harvest_name: str, tmp_path: str,
                 results_dir: str):
    """Copy beam-search script, patching SAMPLE_NAME and results dir."""
    with open(script_path) as f:
        src = f.read()

    # Patch sample name
    src = re.sub(
        r'^SAMPLE_NAME\s*=\s*"[^"]*"',
        f'SAMPLE_NAME = "{harvest_name}"',
        src, count=1, flags=re.MULTILINE,
    )
    # Patch results directory — beam_search now uses a single RESULTS_DIR
    # variable, so one substitution covers all output paths.
    src = re.sub(
        r'^(\s*)RESULTS_DIR\s*=\s*"[^"]*"',
        rf'\1RESULTS_DIR = "{results_dir}"',
        src, count=1, flags=re.MULTILINE,
    )

    with open(tmp_path, "w") as f:
        f.write(src)


def zip_directory(folder: str, zip_path: str):
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for root, _dirs, files in os.walk(folder):
            for fname in files:
                full = os.path.join(root, fname)
                arcname = os.path.relpath(full, folder)
                zf.write(full, arcname)


def run_one_harvest(harvest: str, idx: int, total: int,
                    cores_per_job: int) -> tuple[str, str]:
    """Run a single harvest in an isolated subprocess.

    Returns (harvest_name, "ok" | "fail:<reason>" | "skip").
    """
    zip_name = os.path.join(OUTPUT_DIR, f"{harvest}.zip")
    if SKIP_COMPLETED and os.path.isfile(zip_name):
        return harvest, "skip"

    # Each harvest gets its own results dir and temp script
    results_dir = f"results_{harvest}"
    tmp_script = f"_tmp_beam_{harvest}.py"

    if os.path.isdir(results_dir):
        shutil.rmtree(results_dir)

    patch_script(BEAM_SCRIPT, harvest, tmp_script, results_dir)

    # Per-harvest log file — keeps the main log clean
    log_dir = os.path.join(OUTPUT_DIR, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"{harvest}.log")

    # Limit inner parallelism: tell Python / OpenBLAS / MKL how many cores
    env = os.environ.copy()
    env["PYTORCH_ALLOC_CONF"] = (
        "expandable_segments:True,max_split_size_mb:256"
    )
    # These env vars cap internal thread/process counts so concurrent
    # harvests don't over-subscribe the node.
    env["OMP_NUM_THREADS"] = str(cores_per_job)
    env["MKL_NUM_THREADS"] = str(cores_per_job)
    env["OPENBLAS_NUM_THREADS"] = str(cores_per_job)
    env["NUMEXPR_NUM_THREADS"] = str(cores_per_job)
    # Used by mcmc_fit.py and beam_search — we patch those to read this.
    env["HARVEST_MAX_WORKERS"] = str(cores_per_job)
    env["MPLBACKEND"] = "Agg"  # headless matplotlib

    try:
        with open(log_path, "w") as log_fh:
            result = subprocess.run(
                [sys.executable, "-u", tmp_script],   # -u: unbuffered output
                stdout=log_fh,
                stderr=subprocess.STDOUT,             # merge stderr into log
                env=env,
            )
            rc = result.returncode
    except Exception as e:
        return harvest, f"fail:subprocess error {e}"

    # Cleanup temp script
    try:
        os.remove(tmp_script)
    except OSError:
        pass

    if rc != 0:
        reason = "OOM killed" if rc == -9 else f"exit code {rc}"
        return harvest, f"fail:{reason}"

    if os.path.isdir(results_dir):
        zip_directory(results_dir, zip_name)
        shutil.rmtree(results_dir)
        return harvest, "ok"
    else:
        return harvest, "fail:no results dir"


def main():
    # Force line-buffered stdout so nohup.out updates in real time
    sys.stdout.reconfigure(line_buffering=True)

    # ── Joint mode guard ───────────────────────────────────────────────────────
    # run_all_harvests.py spawns one beam-search process per harvest, which is
    # incompatible with JOINT_FIT_MODE = True.  In joint mode use run_joint.py.
    try:
        with open(BEAM_SCRIPT) as fh:
            _src = fh.read()
        _joint = "JOINT_FIT_MODE    = True" in _src or "JOINT_FIT_MODE = True" in _src
    except FileNotFoundError:
        _joint = False

    if _joint:
        print(
            "run_all_harvests.py: joint mode is ON — cannot run per-harvest jobs.\n"
            "Joint fitting must run exactly once for all mice together.\n"
            "Use:  python run_joint.py          (local)\n"
            "      sbatch submit_joint.sbatch   (Slurm)",
            file=sys.stderr,
        )
        return 1

    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=0,
                        help="Number of harvests to run in parallel. "
                             "0 = auto (total_cores / 32).")
    parser.add_argument("--skip", nargs="+", default=[], metavar="HARVEST",
                        help="Harvest names to skip (space-separated). "
                             "Useful for excluding harvests known to fail or hang.")
    parser.add_argument("--skip-file", type=str, default=None,
                        help="Path to a file with harvest names to skip, "
                             "one per line. Blank lines and lines starting "
                             "with '#' are ignored.")
    args = parser.parse_args()

    # Collect manually-skipped harvests from --skip and --skip-file
    manual_skip = set(args.skip)
    if args.skip_file:
        with open(args.skip_file) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    manual_skip.add(line)

    total_cores = os.cpu_count() or 1
    if args.workers > 0:
        n_parallel = args.workers
    else:
        # Each MCMC internally wants ~32 workers (N_CHAINS).
        # Divide available cores evenly.
        n_parallel = max(1, total_cores // 32)

    cores_per_job = max(1, total_cores // n_parallel)

    harvests = get_matching_harvests(CSV_PATH)
    n_total = len(harvests)

    # Filter out manually-skipped harvests
    if manual_skip:
        unknown = manual_skip - set(harvests)
        if unknown:
            print(f"Warning: --skip names not found in CSV: {sorted(unknown)}")
        harvests = [h for h in harvests if h not in manual_skip]
        n_skipped_manually = n_total - len(harvests)
        print(f"Manually skipping {n_skipped_manually} harvest(s)")

    print(f"Found {len(harvests)} harvests with matching ploidy files"
          + (f" (of {n_total} total)" if manual_skip else ""))
    print(f"CPU cores available : {total_cores}")
    print(f"Parallel harvests   : {n_parallel}")
    print(f"Cores per harvest   : {cores_per_job}")
    print()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    succeeded, failed, skipped = [], [], []

    with ProcessPoolExecutor(max_workers=n_parallel) as pool:
        futures = {
            pool.submit(run_one_harvest, h, i, len(harvests), cores_per_job): h
            for i, h in enumerate(harvests, 1)
        }
        done_count = 0
        for future in as_completed(futures):
            harvest = futures[future]
            done_count += 1
            ts = datetime.now().strftime("%H:%M:%S")
            try:
                _, status = future.result()
            except Exception as e:
                status = f"fail:exception {e}"
            if status == "ok":
                print(f"  [{ts}] [{done_count}/{len(harvests)}] ✓ {harvest}")
                succeeded.append(harvest)
            elif status == "skip":
                print(f"  [{ts}] [{done_count}/{len(harvests)}] ⊘ {harvest} (skipped)")
                skipped.append(harvest)
            else:
                print(f"  [{ts}] [{done_count}/{len(harvests)}] ✗ {harvest} — {status}")
                failed.append(harvest)

    print(f"\n{'='*72}")
    print("SUMMARY")
    print(f"{'='*72}")
    print(f"  Succeeded : {len(succeeded)}/{len(harvests)}")
    print(f"  Skipped   : {len(skipped)}/{len(harvests)}")
    print(f"  Failed    : {len(failed)}/{len(harvests)}")
    if failed:
        print(f"  Failed harvests: {failed}")
    print(f"\n  Zip files saved in : {os.path.abspath(OUTPUT_DIR)}/")
    print(f"  Per-harvest logs in: {os.path.abspath(os.path.join(OUTPUT_DIR, 'logs'))}/")


if __name__ == "__main__":
    main()