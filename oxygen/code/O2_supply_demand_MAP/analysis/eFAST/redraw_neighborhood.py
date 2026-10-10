#!/usr/bin/env python3
"""Redraw completed neighborhood indices and record independent figure provenance."""
import argparse
import datetime
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import efast


def redraw(root):
    status = json.loads((root / "status.json").read_text())
    assert status["stage"] == "full_slurm" and status["status"] == "complete"
    plan = json.loads((root / "slurm/submission_plan.json").read_text())
    source_paths = (root / "convergence.tsv", root / "efast_o2_sensitivity_classification.tsv")
    source_hashes = {p.name: efast.sha256(p) for p in source_paths}
    efast.plot(SimpleNamespace(out_dir=str(root), figure4_dir=str(root), figure4_layout_dir=None,
        combined_only=True, classification_file=str(source_paths[1]), statistic="median",
        summary_label="median across 500 fit seeds of five-phase means"))
    assert source_hashes == {p.name: efast.sha256(p) for p in source_paths}
    commit = subprocess.check_output(["git", "-C", str(Path(__file__).parent), "rev-parse", "HEAD"],
                                     universal_newlines=True).strip()
    receipt = dict(recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        calculation_git_commit=plan["git_commit"], rendering_git_commit=commit,
        calculation_submission_plan_sha256=efast.sha256(root / "slurm/submission_plan.json"),
        rendering_efast_sha256=efast.sha256(Path(efast.__file__)),
        rendering_script_sha256=efast.sha256(Path(__file__)), source_sha256=source_hashes,
        S1_palette=["#FFFFFF", "#3F007D"], ST_palette=["#FFFFFF", "#E6550D"],
        statistical_results_changed=False)
    receipt_path = root / "figure_rendering_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    manifest_path = root / "collection_manifest.tsv"
    rows = {r["path"]: r for r in efast.read_table(manifest_path)}
    for path in (root / "figures/efast_four_panel.pdf", root / "figures/efast_four_panel.png",
                 root / "figure_redraw_manifest.tsv", receipt_path):
        name = str(path.relative_to(root))
        rows[name] = dict(path=name, sha256=efast.sha256(path), bytes=path.stat().st_size)
    efast.write_table(manifest_path, ["path", "sha256", "bytes"], rows.values())
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", required=True)
    redraw(Path(parser.parse_args().out_dir).resolve())
