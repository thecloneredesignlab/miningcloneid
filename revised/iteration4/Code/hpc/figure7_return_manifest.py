#!/usr/bin/env python3
"""Seal/verify a Figure7-only transfer; never includes manuscript text."""
import argparse
import csv
import hashlib
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--root", type=Path, required=True)
parser.add_argument("--run-id", required=True)
parser.add_argument("--verify", action="store_true")
args = parser.parse_args()
root = args.root.resolve()
run_id = args.run_id
if not run_id or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-" for c in run_id):
    raise SystemExit("Invalid run id")
audit = root / "audit" / "hpc_figure7_full_range" / run_id
manifest = audit / "return_sha256.tsv"
file_list = audit / "return_files.txt"

def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

if args.verify:
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        path = root / row["path"]
        if root not in path.resolve().parents or path.stat().st_size != int(row["size_bytes"]) or sha(path) != row["sha256"]:
            raise SystemExit(f"Transfer verification failed: {path}")
    print(f"Verified {len(rows)} files, {sum(int(row['size_bytes']) for row in rows)} bytes")
    raise SystemExit(0)

with (audit / "status.tsv").open() as handle:
    status = next(csv.DictReader(handle, delimiter="\t"))
if status["status"] != "COMPLETE" or status["stage"] != "COMPLETE":
    raise SystemExit("Only a complete calculation-and-render run may be exported")
data = root / "data" / "Figures" / "Figure7" / "fixed_pmisseg_v1"
pointer = data / "finite_time_full_q10_current.tsv"
with pointer.open() as handle:
    current = next(csv.DictReader(handle, delimiter="\t"))
if current["run_id"] != run_id:
    raise SystemExit("Current data pointer does not match requested run")
files = {pointer}
files.update(p for p in (data / "finite_time_full_q10_runs" / run_id).rglob("*") if p.is_file())
files.update(p for p in audit.rglob("*") if p.is_file() and p not in (manifest, file_list))
names = {
    8: "supp_fig7-8_steady_state_full_oxygen_range",
    9: "supp_fig7-9_inverse_response",
    10: "supp_fig7-10_invivo_continuous_full_range",
    11: "supp_fig7-11_invitro_continuous_full_range",
    12: "supp_fig7-12_invitro_passage_full_range",
}
for name in ["assembled_fig7", *names.values()]:
    for directory in (root / "Figures", root / "manuscript" / "Figures"):
        for suffix in (".pdf", ".png"):
            files.add(directory / (name + suffix))
for number, name in names.items():
    directory = root / "data" / "Figures" / f"Supp_Figure7_{number}"
    files.update(directory / (name + suffix) for suffix in (".pdf", ".png", "_render_validation.tsv"))
rows = []
for path in sorted(files):
    if not path.is_file():
        raise SystemExit(f"Missing return artifact: {path}")
    rows.append({"path": str(path.relative_to(root)), "size_bytes": path.stat().st_size, "sha256": sha(path)})
with manifest.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "sha256"], delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
with file_list.open("w") as handle:
    for row in rows:
        handle.write(row["path"] + "\n")
    handle.write(str(manifest.relative_to(root)) + "\n")
    handle.write(str(file_list.relative_to(root)) + "\n")
print(f"Sealed {len(rows)} files, {sum(row['size_bytes'] for row in rows)} bytes: {manifest}")
