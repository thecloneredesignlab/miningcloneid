#!/usr/bin/env python3
"""Seal existing neighborhood work without running models or changing results.

Python 3.6+ / standard library only. Large HPC payloads stay in place: the
manifest distinguishes a fresh SHA256 from a size/mtime inventory explicitly.
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile


RESULT = "oxygen/results/eFAST/"
SCOPES = [
    ("efast_neighborhood", RESULT + "neighborhood10pct"),
    ("sobol_neighborhood_pilot", RESULT + "sobol/neighborhood10pct"),
    ("sobol_neighborhood_N4096", RESULT + "sobol/full_N4096/neighborhood10pct"),
    ("sobol_neighborhood_refinement", RESULT + "sobol/neighborhood_refinement"),
    ("shared_sobol_pilot", RESULT + "sobol/resolution_pilot"),
    ("shared_sobol_inputs", RESULT + "sobol/inputs"),
    ("shared_sobol_runtime", RESULT + "sobol/slurm"),
    ("shared_sobol_N4096_inputs", RESULT + "sobol/full_N4096/inputs"),
    ("shared_sobol_N4096_runtime", RESULT + "sobol/full_N4096/slurm"),
    ("shared_sobol_N4096_recovery", RESULT + "sobol/full_N4096/recovery_20261018"),
    ("shared_sobol_N4096_diagnostics", RESULT + "sobol/full_N4096/diagnostics"),
    ("shared_sobol_refinement_control", RESULT + "sobol/refinement_20261018"),
]
CODE = ["oxygen/code/O2_supply_demand_MAP", "container"]
SKIP_NAMES = {".DS_Store", "NEIGHBORHOOD_ARCHIVED.json"}
SMALL = 256 * 1024


def now():
    return datetime.datetime.utcnow().isoformat() + "Z"


def sha(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save(path, obj):
    path = Path(path)
    tmp = path.with_name(path.name + ".partial")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
    os.replace(str(tmp), str(path))


def walk(root):
    if not root.exists():
        return
    if root.is_file() or root.is_symlink():
        yield root
        return
    for base, dirs, names in os.walk(str(root), followlinks=False):
        dirs[:] = sorted(d for d in dirs if d not in {"__pycache__", ".git"})
        # Preserve symlink identities, including links to directories.
        for d in list(dirs):
            p = Path(base) / d
            if p.is_symlink():
                yield p
                dirs.remove(d)
        for name in sorted(names):
            if name not in SKIP_NAMES:
                yield Path(base) / name


def record(path, relative, content):
    st = path.lstat()
    row = dict(path=str(relative), bytes=st.st_size, mtime_ns=st.st_mtime_ns)
    if path.is_symlink():
        row.update(kind="symlink", target=os.readlink(str(path)), integrity="link_target")
    elif path.is_file():
        row.update(kind="file", integrity="size_mtime_only")
        if content:
            row.update(sha256=sha(path), integrity="fresh_sha256")
            after = path.stat()
            if (st.st_size, st.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise RuntimeError("File changed while hashing: " + str(path))
    else:
        raise RuntimeError("Unsupported filesystem object: " + str(path))
    return row


def should_hash(path, location):
    if location == "local":
        return True
    st = path.lstat()
    # Hash all small files plus control metadata and compact final tables/plots.
    # Bulk trajectories, NPZ checkpoints and cache archives remain in place.
    parts = set(path.parts)
    return (st.st_size <= SMALL or path.suffix in {".json", ".md"}
            or "code_snapshot" in parts or "summaries" in parts
            or "figures" in parts)


def inventory(repo, dest, location):
    output = dest / (location + "_results.jsonl.gz")
    partial = output.with_name(output.name + ".partial")
    counts = collections.defaultdict(collections.Counter)
    absent = []
    total = 0
    scopes = list(SCOPES)
    for rel in [RESULT + "sobol", RESULT + "sobol/full_N4096"]:
        if (repo / rel).exists():
            scopes.extend(("shared_sobol_stage_control", str(p.relative_to(repo)))
                          for p in sorted((repo / rel).iterdir())
                          if p.is_file() and p.name not in SKIP_NAMES)
    with gzip.open(str(partial), "wt", encoding="utf-8", compresslevel=6) as handle:
        for scope, rel in scopes:
            root = repo / rel
            if not root.exists():
                absent.append(rel)
                continue
            print("Inventory " + scope, flush=True)
            for path in walk(root):
                row = record(path, path.relative_to(repo), should_hash(path, location))
                row["scope"] = scope
                handle.write(json.dumps(row, sort_keys=True) + "\n")
                counts[scope]["files"] += 1
                counts[scope]["bytes"] += row["bytes"]
                counts[scope][row["integrity"]] += 1
                if row["integrity"] == "fresh_sha256":
                    counts[scope]["hashed_bytes"] += row["bytes"]
                total += 1
                if total % 25000 == 0:
                    print("Inventoried %d files" % total, flush=True)
    os.replace(str(partial), str(output))
    return dict(file=output.name, sha256=sha(output), scopes=dict(counts),
                absent_roots=absent, policy="all content on local; on HPC all <=256KiB, JSON, Markdown, code snapshots, summaries and figures; others size/mtime only")


def snapshot(repo, dest, label, entries):
    """Make an independent copy of small code/input dependencies and read it back."""
    archive = dest / (label + ".tar.gz")
    manifest = dest / (label + "_manifest.json")
    partial = archive.with_name(archive.name + ".partial")
    rows = []
    with tarfile.open(str(partial), "w:gz", compresslevel=6) as handle:
        for path, member in sorted(entries, key=lambda x: x[1]):
            row = record(path, member, True)
            handle.add(str(path), arcname=member, recursive=False)
            after = path.lstat()
            if (row["bytes"], row["mtime_ns"]) != (after.st_size, after.st_mtime_ns):
                raise RuntimeError("Source changed while packaging: " + str(path))
            rows.append(row)
    os.replace(str(partial), str(archive))
    expected = {r["path"]: r for r in rows}
    with tarfile.open(str(archive), "r:gz") as handle:
        seen = set()
        for member in handle:
            row = expected[member.name]
            if member.issym():
                if member.linkname != row["target"]:
                    raise RuntimeError("Snapshot symlink mismatch")
            else:
                h = hashlib.sha256()
                stream = handle.extractfile(member)
                for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
                    h.update(block)
                if h.hexdigest() != row["sha256"]:
                    raise RuntimeError("Snapshot checksum mismatch: " + member.name)
            seen.add(member.name)
        if seen != set(expected):
            raise RuntimeError("Snapshot file list mismatch")
    save(manifest, dict(archive=archive.name, archive_sha256=sha(archive),
                        readback_verified=True, files=rows))
    return dict(file=archive.name, sha256=sha(archive), files=len(rows),
                bytes=archive.stat().st_size, manifest=manifest.name,
                manifest_sha256=sha(manifest), readback_verified=True)


def code_snapshot(repo, dest, label):
    names = subprocess.check_output(["git", "ls-files", "-z", "--"] + CODE,
                                    cwd=str(repo)).decode().split("\0")
    # Include the new archiver/policy even before their first Git commit.
    names = set(n for n in names if n)
    for n in ["oxygen/code/O2_supply_demand_MAP/analysis/eFAST/archive_neighborhood.py",
              "oxygen/code/O2_supply_demand_MAP/analysis/eFAST/AGENTS.md"]:
        if (repo / n).exists():
            names.add(n)
    entries = [(repo / n, n) for n in names if (repo / n).exists()]
    return snapshot(repo, dest, label, entries)


def fit_snapshot(repo, dest, fit):
    manifest = repo / RESULT / "neighborhood10pct/seed_manifest.tsv"
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 500:
        raise RuntimeError("Expected all 500 fitted endpoints")
    entries = []
    verified = []
    for row in rows:
        for name in ("best_params.tsv", "fit_config.rds", "parameter_table.csv", "fit_summary.tsv"):
            path = fit / row["fit_seed"] / name
            if not path.is_file():
                raise RuntimeError("Missing fitted input: " + str(path))
            key = {"best_params.tsv": "best_params_sha256", "fit_config.rds": "fit_config_sha256",
                   "parameter_table.csv": "parameter_table_sha256"}.get(name)
            if key and sha(path) != row[key]:
                raise RuntimeError("Fitted input differs from certified seed manifest: " + str(path))
            entries.append((path, "original_fit/" + row["fit_seed"] + "/" + name))
        verified.append(row["fit_seed"])
    result = snapshot(repo, dest, "hpc_fit_inputs", entries)
    result.update(original_root=str(fit), verified_endpoints=len(verified),
                  seed_manifest_sha256=sha(manifest))
    return result


def create(args):
    repo, dest = Path(args.repo).resolve(), Path(args.dest).resolve()
    dest.mkdir(parents=True, exist_ok=True)
    receipt = dest / (args.location + "_seal.json")
    if receipt.exists():
        raise RuntimeError("Seal already exists; verify it rather than overwriting")
    metadata = dict(status="building", started_utc=now(), repo=str(repo),
                    location=args.location, archive_mode="in_place_results_with_independent_code_and_fit_snapshots",
                    source_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(repo)).decode().strip(),
                    source_branch=subprocess.check_output(["git", "branch", "--show-current"], cwd=str(repo)).decode().strip(),
                    archive_tool_sha256=sha(Path(__file__)))
    metadata["code"] = code_snapshot(repo, dest, args.location + "_code")
    if args.fit_root:
        metadata["fit_inputs"] = fit_snapshot(repo, dest, Path(args.fit_root))
    metadata["results"] = inventory(repo, dest, args.location)
    metadata.update(status="sealed", completed_utc=now(),
                    bulk_payloads_copied=False, original_scientific_status_modified=False)
    save(receipt, metadata)
    print(json.dumps({k: v for k, v in metadata.items() if k != "results"}, indent=2), flush=True)


def verify(args):
    repo, path = Path(args.repo).resolve(), Path(args.manifest).resolve()
    counts = collections.Counter()
    failures = []
    with gzip.open(str(path), "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            p = repo / row["path"]
            # Manifests are data, never executable paths outside the checkout.
            if Path(row["path"]).is_absolute() or ".." in Path(row["path"]).parts:
                raise RuntimeError("Unsafe inventory path")
            try:
                st = p.lstat()
                if row["kind"] == "symlink":
                    if not p.is_symlink() or os.readlink(str(p)) != row["target"]:
                        raise RuntimeError("symlink target changed")
                    counts["link_target_checked"] += 1
                elif not p.is_file() or p.is_symlink():
                    raise RuntimeError("file type changed")
                elif st.st_size != row["bytes"]:
                    raise RuntimeError("size changed")
                elif "sha256" in row:
                    if sha(p) != row["sha256"]:
                        raise RuntimeError("SHA256 changed")
                    counts["sha256_checked"] += 1
                else:
                    if st.st_mtime_ns != row["mtime_ns"]:
                        raise RuntimeError("mtime changed; no baseline content hash")
                    counts["size_mtime_only_checked"] += 1
                counts["files_checked"] += 1
            except (OSError, RuntimeError) as error:
                failures.append(dict(path=row["path"], error=str(error)))
    result = dict(status="pass" if not failures else "fail", checked_utc=now(),
                  manifest_sha256=sha(path), counts=dict(counts), failures=failures,
                  note="Size/mtime checks do not certify content. Revalidate original trajectory/checkpoint identity receipts before reuse.")
    if args.receipt:
        save(Path(args.receipt), result)
    print(json.dumps(result, indent=2))
    return 0 if not failures else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command")
    p = sub.add_parser("create")
    p.add_argument("--repo", required=True)
    p.add_argument("--dest", required=True)
    p.add_argument("--location", choices=["local", "hpc"], required=True)
    p.add_argument("--fit-root")
    p = sub.add_parser("verify")
    p.add_argument("--repo", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--receipt")
    args = parser.parse_args()
    if args.command == "create":
        create(args)
        return 0
    if args.command == "verify":
        return verify(args)
    parser.error("Choose create or verify")


if __name__ == "__main__":
    raise SystemExit(main())
