#!/usr/bin/env python3
"""Stage versioned analysis sources and archived plan inputs in an isolated checkout."""

import argparse
import csv
from hashlib import sha256
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[5]
DOC_INDEX = Path("oxygen/results/analysis_documentation/index.tsv")
COUNT_PROVENANCE = Path(
    "oxygen/results/structural_identifiability/"
    "20261118_culture_experimental_counts_local_01/source_provenance.tsv"
)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def workspace_path(workspace, relative):
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Expected a repository-relative path: {relative}")
    absolute = workspace / relative
    if absolute.resolve() != absolute:
        raise ValueError(f"Workspace path crosses a symbolic link: {absolute}")
    return absolute


def copy_verified(source, target):
    if target.exists():
        if source.read_bytes() != target.read_bytes():
            raise ValueError(f"Source differs from the existing frozen copy: {target}")
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def relocate_count_provenance(workspace):
    manifest = workspace_path(workspace, COUNT_PROVENANCE)
    if not manifest.is_file():
        return 0
    rows = read_tsv(manifest)
    for row in rows:
        recorded = Path(row["path"])
        parts = recorded.parts
        if "soft_couping_org" in parts:
            target = workspace_path(workspace.parent, Path(*parts[parts.index("soft_couping_org"):]))
        elif "oxygen" in parts:
            target = workspace_path(workspace, Path(*parts[parts.index("oxygen"):]))
        else:
            raise ValueError(f"Unrecognized count provenance origin: {recorded}")
        if sha256(target.read_bytes()).hexdigest() != row["sha256"]:
            raise ValueError(f"Count provenance checksum mismatch: {target}")
        row["path"] = str(target)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("path", "sha256"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def prepare(workspace, rebase_count_provenance=False):
    workspace = workspace.resolve()
    if workspace == REPOSITORY or workspace.is_relative_to(REPOSITORY):
        raise ValueError("Provide a separate isolated checkout outside the source repository")
    if not (workspace / "oxygen").is_dir():
        raise ValueError(f"Missing oxygen directory in the isolated checkout: {workspace}")
    sources = read_tsv(HERE / "replay_layout.tsv")
    docs_path = workspace_path(workspace, DOC_INDEX)
    documents = read_tsv(docs_path) if docs_path.is_file() else []
    for row in sources:
        copy_verified(REPOSITORY / row["source"], workspace_path(workspace, Path(row["workspace_path"])))
    for row in documents:
        source = workspace_path(workspace, Path(row["archived_path"]))
        if sha256(source.read_bytes()).hexdigest() != row["sha256"]:
            raise ValueError(f"Archived plan checksum mismatch: {source}")
        copy_verified(source, workspace_path(workspace, Path(row["previous_path"])))
    rebased = relocate_count_provenance(workspace) if rebase_count_provenance else 0
    return len(sources), len(documents), rebased


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path, help="Separate checkout with the required archived inputs")
    parser.add_argument("--rebase-count-provenance", action="store_true", help="Verify and rebase C1 input paths within the isolated checkout")
    args = parser.parse_args()
    scripts, plans, paths = prepare(args.workspace, args.rebase_count_provenance)
    print(f"Staged {scripts} source files and {plans} archived documents. Rebased {paths} verified C1 input paths.")


if __name__ == "__main__":
    main()
