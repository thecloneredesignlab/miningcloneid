#!/usr/bin/env python3
"""Hash generated intermediates and verify published figure copies."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


SCRIPT = Path(__file__).resolve()


def find_workspace_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (
            (candidate / "manager.sh").is_file()
            and (candidate / "Code" / "Figures").is_dir()
        ):
            return candidate.resolve()
    raise RuntimeError(f"Could not locate figure workspace from: {start}")


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Hash generated intermediates and verify published figures."
    )
    parser.add_argument(
        "--only",
        default="",
        help="Optional comma-separated artifact IDs from figure_entrypoints.tsv.",
    )
    arguments = parser.parse_args()
    root = find_workspace_root(SCRIPT.parent)
    data_root = root / "data" / "Figures"
    figure_root = root / "Figures"
    published_root = root / "manuscript" / "Figures"
    entrypoint_path = root / "audit" / "manifests" / "figure_entrypoints.tsv"
    with entrypoint_path.open(encoding="utf-8", newline="") as handle:
        entrypoints = list(csv.DictReader(handle, delimiter="\t"))
    requested = {
        value.strip() for value in arguments.only.split(",") if value.strip()
    }
    if requested:
        available = {row["figure"] for row in entrypoints}
        unknown = sorted(requested - available)
        if unknown:
            raise RuntimeError("Unknown artifact ID(s): " + ", ".join(unknown))
        entrypoints = [row for row in entrypoints if row["figure"] in requested]

    required_names = sorted({
        Path(row["main_output"]).with_suffix(suffix).name
        for row in entrypoints
        for suffix in (".png", ".pdf")
    })
    failures: list[str] = []
    records: list[dict[str, object]] = []

    intermediate_roots = (
        [root / row["data_directory"] for row in entrypoints]
        if requested
        else [data_root]
    )
    intermediate_paths = sorted({
        path
        for intermediate_root in intermediate_roots
        for path in intermediate_root.rglob("*")
        if path.is_file()
        and path.name != ".DS_Store"
        and path.suffix.lower() != ".pyc"
    })
    for path in intermediate_paths:
        records.append({
            "category": "figure_intermediate",
            "relative_path": path.relative_to(root).as_posix(),
            "size_bytes": path.stat().st_size,
            "md5": md5(path),
        })

    for name in required_names:
        source = figure_root / name
        published = published_root / name
        if not source.is_file():
            failures.append(f"missing figure output: {source}")
            continue
        if not published.is_file():
            failures.append(f"missing published figure: {published}")
            continue
        source_md5 = md5(source)
        published_md5 = md5(published)
        for category, path, digest in (
            ("figure_output", source, source_md5),
            ("published_figure", published, published_md5),
        ):
            records.append({
                "category": category,
                "relative_path": path.relative_to(root).as_posix(),
                "size_bytes": path.stat().st_size,
                "md5": digest,
            })
        if source_md5 != published_md5:
            failures.append(f"published MD5 mismatch: {name}")

    scope_suffix = ""
    if requested:
        scope_suffix = "_" + "_".join(sorted(requested)).lower()
    manifest_path = (
        root / "audit" / "md5" / f"generated_artifact_md5{scope_suffix}.tsv"
    )
    write_tsv(
        manifest_path,
        ("category", "relative_path", "size_bytes", "md5"),
        records,
    )
    report_rows = [
        {
            "check": "figure_intermediate_files_hashed",
            "value": len(intermediate_paths),
            "status": "PASS",
        },
        {
            "check": "required_figure_pairs",
            "value": len(required_names) // 2,
            "status": "PASS" if not failures else "FAIL",
        },
        {
            "check": "published_md5_failures",
            "value": len(failures),
            "status": "PASS" if not failures else "FAIL",
        },
    ]
    write_tsv(
        root
        / "audit"
        / "reports"
        / f"artifact_md5_validation{scope_suffix}.tsv",
        ("check", "value", "status"),
        report_rows,
    )
    if failures:
        raise RuntimeError("\n".join(failures))
    print(
        f"Hashed {len(intermediate_paths)} intermediates and verified "
        f"{len(required_names)} published figure files."
    )


if __name__ == "__main__":
    main()
