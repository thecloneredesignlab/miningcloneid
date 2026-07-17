"""
submit_harvests.py
──────────────────
Prepares a multi-node Slurm job array for the harvest pipeline:

  1. Reads the harvest_ploidy_mapping.csv
  2. Writes harvests.txt (one harvest name per line) — consumed by
     submit_harvests.sbatch via sed on $SLURM_ARRAY_TASK_ID
  3. Prints the exact sbatch command to submit the array
  4. Optionally runs sbatch itself with --submit

Usage:
    python submit_harvests.py
    python submit_harvests.py --max-concurrent 20         # cap simultaneous tasks
    python submit_harvests.py --submit                     # prepare AND submit
    python submit_harvests.py --submit --notify-email you@x   # + email when done
    python submit_harvests.py --skip HARVEST1 HARVEST2     # exclude some harvests

Design note — why an array, not MPI or a single huge allocation:
    Each harvest is independent and already CPU-parallel inside (32 chains
    via multiprocessing.Pool). The only axis left to scale is "how many
    harvests run at once". A Slurm array at 32 cores/task is the native
    way to express that; Slurm handles packing tasks onto nodes for us.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys


CSV_PATH = "../data/InVivoData_Gemcitabine/harvest_ploidy_mapping.csv"
HARVESTS_TXT = "harvests.txt"
SBATCH_SCRIPT = "submit_harvests.sbatch"


def get_matching_harvests(csv_path: str) -> list[str]:
    harvests = []
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["has_match"].strip() == "True":
                harvests.append(row["harvest"].strip())
    return harvests


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default=CSV_PATH,
                        help=f"Path to harvest mapping CSV (default: {CSV_PATH})")
    parser.add_argument("--out", default=HARVESTS_TXT,
                        help=f"Output task list path (default: {HARVESTS_TXT})")
    parser.add_argument("--sbatch-script", default=SBATCH_SCRIPT,
                        help=f"Sbatch script to submit (default: {SBATCH_SCRIPT})")
    parser.add_argument("--max-concurrent", type=int, default=0,
                        help="Cap on simultaneously-running array tasks. "
                             "0 = no cap (let Slurm schedule as widely as possible).")
    parser.add_argument("--skip", nargs="*", default=[], metavar="HARVEST",
                        help="Harvest names to exclude.")
    parser.add_argument("--skip-file", type=str, default=None,
                        help="File with harvest names to skip, one per line.")
    parser.add_argument("--skip-completed", action="store_true",
                        help="Exclude harvests whose zip already exists in harvest_zips/.")
    parser.add_argument("--submit", action="store_true",
                        help="Also run the sbatch command after writing the task list.")
    parser.add_argument("--notify-email", type=str, default=None, metavar="ADDR",
                        help="Email to notify when the whole array finishes. "
                             "Implies --submit. Chains a tiny dependent job "
                             "(afterany) that sends a single summary email "
                             "instead of one-per-task spam.")
    args = parser.parse_args()

    skip = set(args.skip)
    if args.skip_file:
        with open(args.skip_file) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    skip.add(line)

    harvests = get_matching_harvests(args.csv)
    n_total = len(harvests)

    if skip:
        unknown = skip - set(harvests)
        if unknown:
            print(f"Warning: --skip names not in CSV: {sorted(unknown)}")
        harvests = [h for h in harvests if h not in skip]

    if args.skip_completed:
        before = len(harvests)
        harvests = [
            h for h in harvests
            if not os.path.isfile(os.path.join("harvest_zips", f"{h}.zip"))
        ]
        print(f"Skipping {before - len(harvests)} already-completed harvest(s)")

    n = len(harvests)
    if n == 0:
        print("Nothing to do — no harvests remain after filtering.")
        return 0

    with open(args.out, "w") as f:
        for h in harvests:
            f.write(h + "\n")

    print(f"Matching harvests in CSV : {n_total}")
    print(f"Tasks to submit          : {n}")
    print(f"Task list written to     : {args.out}")

    # Build array spec: "0-(n-1)" optionally with "%K" concurrency cap
    array_spec = f"0-{n - 1}"
    if args.max_concurrent and args.max_concurrent < n:
        array_spec += f"%{args.max_concurrent}"

    cmd = ["sbatch", f"--array={array_spec}", args.sbatch_script]
    print("\nSubmit with:")
    print("  " + " ".join(cmd))

    # --notify-email implies --submit; we need to actually submit to know the jobid
    should_submit = args.submit or bool(args.notify_email)

    if should_submit:
        print("\nSubmitting array...")
        # Capture stdout so we can parse the job id Slurm prints.
        result = subprocess.run(cmd, capture_output=True, text=True)
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        if result.returncode != 0:
            return result.returncode

        # Slurm prints: "Submitted batch job <jobid>"
        jobid = None
        for line in result.stdout.splitlines():
            line = line.strip()
            if line.startswith("Submitted batch job"):
                jobid = line.split()[-1]
                break

        if args.notify_email:
            if not jobid:
                print("WARNING: could not parse array jobid from sbatch output; "
                      "not chaining the summary email job.", file=sys.stderr)
                return 0
            # One tiny job, depends on the whole array (afterany = runs once
            # the array is done regardless of per-task success/failure), sends
            # exactly one email when it completes.
            summary_cmd = [
                "sbatch",
                f"--dependency=afterany:{jobid}",
                "--job-name=harvest_mcmc_notify",
                "--time=00:02:00",
                "--mem=100M",
                "--output=harvest_zips/logs/slurm_notify_%j.out",
                "--mail-type=END",
                f"--mail-user={args.notify_email}",
                "--wrap",
                (f"echo 'Harvest array {jobid} complete.'; "
                 f"echo ''; "
                 f"echo 'Per-task status:'; "
                 f"sacct -j {jobid} --format=JobID,State,ExitCode,Elapsed -P || true; "
                 f"echo ''; "
                 f"echo 'Zips produced:'; "
                 f"ls -1 harvest_zips/*.zip 2>/dev/null | wc -l"),
            ]
            print("\nChaining summary-email job...")
            print("  " + " ".join(summary_cmd))
            notify_result = subprocess.run(summary_cmd)
            return notify_result.returncode

        return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())