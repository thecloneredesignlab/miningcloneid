#!/usr/bin/env python3
"""Verify and archive the compact deliverables, leaving raw evaluations on HPC."""
import argparse
import csv
import hashlib
import json
import tarfile
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def collect(root, archive):
    status = json.loads((root / "slurm/status.json").read_text())
    assert status["completed_tasks"] == 35000 and status["completed_seed_summaries"] == 500
    assert status["failed_tasks"] == 0 and status["missing_or_running_tasks"] == 0
    audit = json.loads((root / "slurm/numerical_recovery_audit.json").read_text())
    assert audit["status"] == "passed" and audit["completed_tasks"] == 35000
    with (root / "collection_manifest.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            path = root / row["path"]
            assert path.is_file() and path.stat().st_size == int(row["bytes"]), path
            assert sha256(path) == row["sha256"], path
    # Directory depth is deliberate: trajectory outputs, caches and logs stay on HPC.
    selected = set()
    for directory in (root, root / "summaries", root / "figures", root / "slurm", root / "pilot"):
        for path in directory.iterdir() if directory.exists() else []:
            if path.is_file() and path.suffix in (".tsv", ".gz", ".json", ".md", ".npz", ".pdf", ".png"):
                if directory == root / "slurm" and path.name.endswith(".tar.gz"):
                    continue
                if path.name != "artifact_collection_manifest.tsv":
                    selected.add(path)
    with (root / "artifact_collection_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("path", "sha256", "bytes"), delimiter="\t")
        writer.writeheader()
        for path in sorted(selected):
            assert path.stat().st_size < 95 * 1024 ** 2, path
            writer.writerow(dict(path=str(path.relative_to(root)), sha256=sha256(path), bytes=path.stat().st_size))
    selected.add(root / "artifact_collection_manifest.tsv")
    archive.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(str(archive), "w") as output:
        for path in sorted(selected):
            output.add(str(path), arcname=str(path.relative_to(root)), recursive=False)
    print(json.dumps(dict(status="passed", files=len(selected), bytes=archive.stat().st_size,
                          archive=str(archive), archive_sha256=sha256(archive)), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--archive", required=True)
    args = parser.parse_args()
    collect(Path(args.out_dir).resolve(), Path(args.archive).resolve())
