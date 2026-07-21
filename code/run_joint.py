"""
run_joint.py
────────────
Single-process runner for joint hierarchical fitting (JOINT_FIT_MODE = True).

All mice are fitted simultaneously in one MCMC run.  Results are written to
results_joint/<harvest_name>/ and then zipped into harvest_zips/<harvest_name>.zip,
so the output format matches what the original per-harvest pipeline produced.

Usage (local):
    python run_joint.py
    python run_joint.py --cores 32

Usage (Slurm, via submit_joint.sbatch):
    sbatch submit_joint.sbatch

Do NOT run this via run_one_harvest.py / run_all_harvests.py / submit_harvests.sbatch —
those scripts are for per-harvest (JOINT_FIT_MODE = False) runs only.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import zipfile
from time import time

BEAM_SCRIPT = "beam_search_flip_rate_wgd.py"
OUTPUT_DIR  = "harvest_zips"
RESULTS_DIR = "results_joint"   # must match the value hardcoded in beam_search


def _check_joint_mode() -> bool:
    """Return True iff JOINT_FIT_MODE = True is set in the beam-search script."""
    try:
        with open(BEAM_SCRIPT) as fh:
            src = fh.read()
        return "JOINT_FIT_MODE    = True" in src or "JOINT_FIT_MODE = True" in src
    except FileNotFoundError:
        return False


def zip_directory(folder: str, zip_path: str) -> None:
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for root, _dirs, files in os.walk(folder):
            for fname in files:
                full    = os.path.join(root, fname)
                arcname = os.path.relpath(full, folder)
                zf.write(full, arcname)


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(
        description="Run the joint hierarchical MCMC fit for all mice.",
    )
    parser.add_argument(
        "--cores", type=int, default=0,
        help="Number of CPU cores to use.  0 = auto (SLURM_CPUS_PER_TASK or cpu_count).",
    )
    parser.add_argument(
        "--no-skip", dest="skip_completed", action="store_false", default=True,
        help="Re-run even if all per-mouse zips already exist.",
    )
    args = parser.parse_args()

    # ── Guard: require JOINT_FIT_MODE = True ──────────────────────────────────
    if not _check_joint_mode():
        print(
            "ERROR: JOINT_FIT_MODE is not True in beam_search_flip_rate_wgd.py.\n"
            "Set JOINT_FIT_MODE = True before running this script.",
            file=sys.stderr,
        )
        return 1

    # ── Resolve core count ─────────────────────────────────────────────────────
    cores = args.cores
    if cores <= 0:
        cores = int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or (os.cpu_count() or 4)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    log_dir = os.path.join(OUTPUT_DIR, "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, "joint_fit.log")

    print(f"▶ Joint fit  (cores={cores}, results={RESULTS_DIR}, log={log_path})")
    t0 = time()

    # ── Run the beam-search script once, unpatched ─────────────────────────────
    env = os.environ.copy()
    env["HARVEST_MAX_WORKERS"]  = str(cores)
    env["OMP_NUM_THREADS"]      = str(cores)
    env["MKL_NUM_THREADS"]      = str(cores)
    env["OPENBLAS_NUM_THREADS"] = str(cores)
    env["NUMEXPR_NUM_THREADS"]  = str(cores)
    env["MPLBACKEND"]           = "Agg"
    env["PYTORCH_ALLOC_CONF"]   = "expandable_segments:True,max_split_size_mb:256"

    try:
        with open(log_path, "w") as log_fh:
            result = subprocess.run(
                [sys.executable, "-u", BEAM_SCRIPT],
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                env=env,
            )
            rc = result.returncode
    except Exception as exc:
        print(f"✗ Joint fit — subprocess error: {exc}", file=sys.stderr)
        return 1

    if rc != 0:
        reason = "OOM killed" if rc == -9 else f"exit code {rc}"
        print(f"✗ Joint fit failed ({reason}), see {log_path}", file=sys.stderr)
        return rc if rc > 0 else 1

    elapsed = time() - t0
    print(f"✓ Joint fit completed in {elapsed:.1f}s")

    # ── Zip each per-mouse results subdirectory ────────────────────────────────
    if not os.path.isdir(RESULTS_DIR):
        print(f"✗ Expected results directory '{RESULTS_DIR}' not found.", file=sys.stderr)
        return 1

    zipped = []
    skipped = []
    for name in sorted(os.listdir(RESULTS_DIR)):
        mouse_dir = os.path.join(RESULTS_DIR, name)
        if not (os.path.isdir(mouse_dir) and name.endswith("_harvest")):
            continue
        zip_path = os.path.join(OUTPUT_DIR, f"{name}.zip")
        if args.skip_completed and os.path.isfile(zip_path):
            skipped.append(name)
            print(f"  ⊘ {name} (zip already exists, skipping)")
            continue
        zip_directory(mouse_dir, zip_path)
        zipped.append(name)
        print(f"  ✓ {name} → {zip_path}")

    # Also package the global params JSON on its own
    global_json = os.path.join(RESULTS_DIR, "joint_global_params.json")
    if os.path.isfile(global_json):
        shutil.copy(global_json, os.path.join(OUTPUT_DIR, "joint_global_params.json"))
        print(f"  ✓ joint_global_params.json → {OUTPUT_DIR}/joint_global_params.json")

    print(f"\n{'='*60}")
    print(f"  Zipped : {len(zipped)}")
    print(f"  Skipped: {len(skipped)}")
    print(f"  Output : {os.path.abspath(OUTPUT_DIR)}/")
    print(f"  Log    : {os.path.abspath(log_path)}")
    print(f"{'='*60}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
