#!/usr/bin/env python3
"""Check collected archive receipts and write the final, nonrecursive index."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True)
    args = parser.parse_args()
    root = Path(args.archive).resolve()
    receipts = {name: json.loads((root / name).read_text()) for name in [
        "local_seal.json", "hpc_seal.json", "local_payload.json",
        "evidence/local_verification.json", "hpc_inventory_validation.json",
        "hpc_critical_verification.json", "hpc_runtime.json", "sobol_refinement_resume_state.json"]}
    for name in ["local_seal.json", "hpc_seal.json"]:
        if receipts[name]["status"] != "sealed":
            raise RuntimeError("Unsealed inventory: " + name)
        result = receipts[name]["results"]
        if sha(root / result["file"]) != result["sha256"]:
            raise RuntimeError("Collected inventory changed: " + name)
    for name in ["evidence/local_verification.json", "hpc_inventory_validation.json", "hpc_critical_verification.json"]:
        if receipts[name]["status"] != "pass":
            raise RuntimeError("Failed archive check: " + name)
    for path in root.glob("*_manifest.json"):
        m = json.loads(path.read_text())
        if not m["readback_verified"] or sha(root / m["archive"]) != m["archive_sha256"]:
            raise RuntimeError("Collected code/input snapshot changed: " + str(path))
    runtime = receipts["hpc_runtime.json"]
    if not runtime["expected_hash_matches"]:
        raise RuntimeError("SIF digest differs")
    payload = receipts["local_payload.json"]
    if (not payload["all_member_sha256_checked"] or not payload["independent_copy"]
            or (root / payload["file"]).stat().st_size != payload["bytes"]):
        raise RuntimeError("Incomplete local portable result copy")
    members = []
    for path in sorted(root.rglob("*")):
        rel = path.relative_to(root)
        if (not path.is_file() or "payload" in rel.parts or "__pycache__" in rel.parts
                or path.name in {".DS_Store", "archive.json"}):
            continue
        if path.name.endswith(".partial"):
            raise RuntimeError("Unfinished archive artifact: " + str(rel))
        members.append(dict(path=str(rel), bytes=path.stat().st_size, sha256=sha(path)))
    output = dict(format_version=1, state="archived", scope=["eFAST neighborhoods", "Sobol neighborhoods"],
                  archived_utc=datetime.datetime.utcnow().isoformat() + "Z",
                  resume_requires_explicit_user_request=True, global_analysis_paused=False,
                  original_results_preserved_in_place=True, bulk_hpc_payloads_copied=False,
                  git_tag="archive/neighborhood-20261021", local_source_head=receipts["local_seal.json"]["source_head"],
                  hpc_source_head=receipts["hpc_seal.json"]["source_head"],
                  hpc_inventory_summary=receipts["hpc_inventory_validation.json"],
                  local_inventory_summary=receipts["evidence/local_verification.json"]["counts"],
                  local_portable_backup=payload, sif=runtime,
                  stopped_jobs=["21140548", "21140549", "21140562", "21140563"],
                  sobol_refinement_state="incomplete_no_retry", sobol_failed_model_blocks=32,
                  files=members)
    (root / "archive.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({k: v for k, v in output.items() if k not in {"files", "sif", "local_portable_backup", "hpc_inventory_summary"}}, indent=2))


if __name__ == "__main__":
    main()
