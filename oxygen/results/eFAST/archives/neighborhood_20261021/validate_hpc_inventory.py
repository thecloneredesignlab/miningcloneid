#!/usr/bin/env python3
"""Validate the HPC inventory and recheck critical control/runtime files.

Bulk stat-only entries remain explicitly uncertified by content. No model calls.
"""
import collections
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time


def main():
    repo, archive, helper = map(Path, sys.argv[1:4])
    deadline = time.monotonic() + 3600
    while not (archive / "hpc_seal.json").exists():
        if time.monotonic() > deadline:
            raise RuntimeError("HPC inventory seal did not finish within one hour")
        time.sleep(10)
    spec = importlib.util.spec_from_file_location("archiver", str(helper))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    seal = json.loads((archive / "hpc_seal.json").read_text())
    manifest = archive / seal["results"]["file"]
    if module.sha(manifest) != seal["results"]["sha256"]:
        raise RuntimeError("Inventory manifest checksum changed")
    important = {"status.json", "plan.json", "submission_plan.json", "completion.json",
                 "input_manifest.json", "seed_manifest.tsv", "environment.tsv",
                 "parameter_bounds.tsv", "reference_best_params.tsv", "neighborhood_bounds.tsv",
                 "parameter_ranges.tsv", "selection.json", "targets.tsv.gz", "execution_policy.json",
                 "operator_input_audit.json", "rcpp_template.tar.gz", "compiled_review_template.tar.gz",
                 "runtime_cache_manifest.json", "rcpp_template_manifest.json", "preparation.json"}
    counts = collections.defaultdict(collections.Counter)
    seen = set()
    critical = archive / "hpc_critical_files.jsonl.gz"
    critical_count = 0
    with gzip.open(str(manifest), "rt") as source, gzip.open(str(critical), "wt") as output:
        for line in source:
            row = json.loads(line)
            rel = Path(row["path"])
            if rel.is_absolute() or ".." in rel.parts or str(rel) in seen:
                raise RuntimeError("Unsafe or duplicated inventory path")
            seen.add(str(rel))
            c = counts[row["scope"]]
            c["files"] += 1
            c["bytes"] += row["bytes"]
            c[row["integrity"]] += 1
            if "sha256" in row:
                if len(row["sha256"]) != 64 or any(x not in "0123456789abcdef" for x in row["sha256"]):
                    raise RuntimeError("Malformed content digest")
                c["hashed_bytes"] += row["bytes"]
            if rel.name in important or "code_snapshot" in rel.parts:
                new = module.record(repo / rel, rel, True)
                if (new["bytes"], new["mtime_ns"]) != (row["bytes"], row["mtime_ns"]):
                    raise RuntimeError("Critical source changed during archive: " + str(rel))
                if "sha256" in row and new["sha256"] != row["sha256"]:
                    raise RuntimeError("Critical content changed: " + str(rel))
                new["scope"] = row["scope"]
                output.write(json.dumps(new, sort_keys=True) + "\n")
                critical_count += 1
    if dict(counts) != seal["results"]["scopes"]:
        raise RuntimeError("Inventory counts do not match seal")
    result = dict(status="pass", checked_utc=module.now(), inventory_sha256=module.sha(manifest),
                  files=len(seen), bytes=sum(c["bytes"] for c in counts.values()),
                  fresh_sha256_files=sum(c["fresh_sha256"] for c in counts.values()),
                  size_mtime_only_files=sum(c["size_mtime_only"] for c in counts.values()),
                  duplicate_paths=0, scope_counts_match=True, critical_files=critical_count,
                  bulk_payload_content_rechecked=False,
                  note="Inventory digest/schema/counts checked; critical sources rehashed. Bulk stat-only records are not content-certified.")
    module.save(archive / "hpc_inventory_validation.json", result)
    class Args:
        pass
    args = Args()
    args.repo = str(repo)
    args.manifest = str(critical)
    args.receipt = str(archive / "hpc_critical_verification.json")
    if module.verify(args):
        raise RuntimeError("Critical source readback failed")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
