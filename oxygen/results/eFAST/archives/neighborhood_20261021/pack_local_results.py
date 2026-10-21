#!/usr/bin/env python3
"""Copy the sealed local inventory to a portable tar; never run a model."""
import argparse
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class Reader:
    def __init__(self, stream):
        self.stream = stream
        self.digest = hashlib.sha256()

    def read(self, size):
        block = self.stream.read(size)
        self.digest.update(block)
        return block


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True)
    parser.add_argument("--archive", required=True)
    args = parser.parse_args()
    repo, archive = Path(args.repo).resolve(), Path(args.archive).resolve()
    manifest = archive / "local_results.jsonl.gz"
    with gzip.open(str(manifest), "rt") as handle:
        rows = [json.loads(line) for line in handle]
    payload = archive / "payload"
    payload.mkdir(exist_ok=True)
    target = payload / "local_results.tar"
    if target.exists():
        raise RuntimeError("Portable result snapshot already exists")
    if shutil.disk_usage(payload).free < sum(r["bytes"] for r in rows) + 4 * 1024**3:
        raise RuntimeError("Insufficient space for independent local result copy")
    partial = target.with_name(target.name + ".partial")
    with tarfile.open(str(partial), "w") as handle:
        for row in rows:
            rel = Path(row["path"])
            if rel.is_absolute() or ".." in rel.parts:
                raise RuntimeError("Unsafe path in manifest")
            source = repo / rel
            info = handle.gettarinfo(str(source), arcname=str(rel))
            if row["kind"] == "symlink":
                if info.linkname != row["target"]:
                    raise RuntimeError("Symlink changed: " + str(rel))
                handle.addfile(info)
                continue
            if info.size != row["bytes"] or source.is_symlink():
                raise RuntimeError("Source changed: " + str(rel))
            with source.open("rb") as stream:
                reader = Reader(stream)
                handle.addfile(info, reader)
                if reader.digest.hexdigest() != row["sha256"]:
                    raise RuntimeError("Source SHA256 changed: " + str(rel))
    # Read the independent copy back, checking every archived file's content.
    expected = {r["path"]: r for r in rows}
    seen = set()
    with tarfile.open(str(partial), "r") as handle:
        for member in handle:
            row = expected[member.name]
            if member.issym():
                assert member.linkname == row["target"]
            else:
                digest = hashlib.sha256()
                stream = handle.extractfile(member)
                for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
                    digest.update(block)
                if digest.hexdigest() != row["sha256"]:
                    raise RuntimeError("Result snapshot checksum failed: " + member.name)
            seen.add(member.name)
    if seen != set(expected):
        raise RuntimeError("Result snapshot membership failed")
    os.replace(str(partial), str(target))
    receipt = dict(status="complete", file="payload/local_results.tar", bytes=target.stat().st_size,
                   sha256=sha(target), source_manifest_sha256=sha(manifest), files=len(rows),
                   all_member_sha256_checked=True, independent_copy=True,
                   completed_utc=datetime.datetime.utcnow().isoformat() + "Z")
    (archive / "local_payload.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
