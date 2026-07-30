#!/usr/bin/env python3
"""Validate runtime scientific inputs against the fixed MD5 baseline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path


SCRIPT = Path(__file__).resolve()
ROOT_ENVIRONMENT = {
    "invitro_result_root": "FIGURE_INVITRO_RESULT_ROOT",
    "invivo_result_root": "FIGURE_INVIVO_RESULT_ROOT",
    "joint_result_root": "FIGURE_JOINT_RESULT_ROOT",
    "gemcitabine_data_root": "FIGURE_GEMCITABINE_DATA_ROOT",
    "ltee_data_root": "FIGURE_LTEE_DATA_ROOT",
}
BASELINE_COLUMNS = (
    "input_root_id",
    "relative_path",
    "size_bytes",
    "expected_md5",
)


def find_workspace_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (
            (candidate / "manager.sh").is_file()
            and (candidate / "Code" / "Figures").is_dir()
        ):
            return candidate.resolve()
    raise RuntimeError(f"Could not locate figure workspace from: {start}")


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(
    path: Path,
    columns: tuple[str, ...],
    rows: list[dict[str, object]],
) -> None:
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


def resolve_input_roots() -> dict[str, Path]:
    roots: dict[str, Path] = {}
    for root_id, environment_name in ROOT_ENVIRONMENT.items():
        value = os.environ.get(environment_name, "").strip()
        if not value:
            raise RuntimeError(
                f"Missing required environment variable: {environment_name}"
            )
        path = Path(value).expanduser()
        if not path.is_dir():
            raise RuntimeError(
                f"Scientific input root does not exist: {environment_name}={path}"
            )
        roots[root_id] = path.resolve()
    return roots


def read_baseline(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise RuntimeError(f"Missing scientific-input MD5 baseline: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != BASELINE_COLUMNS:
            raise RuntimeError(
                "Invalid scientific-input MD5 baseline columns: "
                f"expected {BASELINE_COLUMNS}, observed {reader.fieldnames}"
            )
        rows = list(reader)
    if not rows:
        raise RuntimeError("Scientific-input MD5 baseline is empty.")

    keys: set[tuple[str, str]] = set()
    for line_number, row in enumerate(rows, start=2):
        root_id = row["input_root_id"].strip()
        relative_path = row["relative_path"].strip()
        expected_md5 = row["expected_md5"].strip().lower()
        if root_id not in ROOT_ENVIRONMENT:
            raise RuntimeError(
                f"Unknown input_root_id at baseline line {line_number}: {root_id}"
            )
        relative = Path(relative_path)
        if (
            not relative_path
            or relative.is_absolute()
            or ".." in relative.parts
        ):
            raise RuntimeError(
                f"Unsafe relative_path at baseline line {line_number}: "
                f"{relative_path}"
            )
        try:
            size_bytes = int(row["size_bytes"])
        except ValueError as error:
            raise RuntimeError(
                f"Invalid size_bytes at baseline line {line_number}: "
                f"{row['size_bytes']}"
            ) from error
        if size_bytes < 0:
            raise RuntimeError(
                f"Negative size_bytes at baseline line {line_number}."
            )
        if (
            len(expected_md5) != 32
            or any(character not in "0123456789abcdef" for character in expected_md5)
        ):
            raise RuntimeError(
                f"Invalid expected_md5 at baseline line {line_number}: "
                f"{row['expected_md5']}"
            )
        key = (root_id, relative_path)
        if key in keys:
            raise RuntimeError(
                f"Duplicate scientific-input baseline row: {root_id}/{relative_path}"
            )
        keys.add(key)
        row["input_root_id"] = root_id
        row["relative_path"] = relative_path
        row["size_bytes"] = str(size_bytes)
        row["expected_md5"] = expected_md5
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate scientific inputs against the Code/config MD5 baseline."
    )
    parser.add_argument(
        "--phase",
        choices=("preflight", "postflight", "manual"),
        default="manual",
        help="Label recorded in the validation summary.",
    )
    arguments = parser.parse_args()

    workspace_root = find_workspace_root(SCRIPT.parent)
    baseline_path = (
        workspace_root
        / "Code"
        / "config"
        / "manifests"
        / "expected_scientific_input_md5.tsv"
    )
    audit_table_path = (
        workspace_root / "audit" / "md5" / "scientific_input_md5.tsv"
    )
    report_path = (
        workspace_root / "audit" / "reports" / "input_md5_validation.tsv"
    )
    roots = resolve_input_roots()
    baseline = read_baseline(baseline_path)

    audit_rows: list[dict[str, object]] = []
    failures: list[str] = []
    for row in baseline:
        root_id = row["input_root_id"]
        relative_path = row["relative_path"]
        expected_size = int(row["size_bytes"])
        expected_md5 = row["expected_md5"]
        root = roots[root_id]
        path = root / relative_path

        exists = path.is_file()
        safely_contained = False
        if exists:
            resolved_path = path.resolve()
            try:
                resolved_path.relative_to(root)
                safely_contained = True
            except ValueError:
                safely_contained = False
        else:
            resolved_path = path

        actual_size = path.stat().st_size if exists and safely_contained else None
        actual_md5 = file_md5(path) if exists and safely_contained else ""
        size_matches = actual_size == expected_size
        md5_matches = actual_md5 == expected_md5
        passed = exists and safely_contained and size_matches and md5_matches

        if not passed:
            reasons: list[str] = []
            if not exists:
                reasons.append("missing")
            elif not safely_contained:
                reasons.append("resolves outside input root")
            else:
                if not size_matches:
                    reasons.append("size mismatch")
                if not md5_matches:
                    reasons.append("MD5 mismatch")
            failures.append(
                f"{root_id}/{relative_path}: {', '.join(reasons)}"
            )

        audit_rows.append(
            {
                "input_root_id": root_id,
                "relative_path": relative_path,
                "absolute_path": str(resolved_path),
                "expected_size_bytes": expected_size,
                "actual_size_bytes": "" if actual_size is None else actual_size,
                "expected_md5": expected_md5,
                "actual_md5": actual_md5,
                "exists": str(exists).upper(),
                "contained_in_root": str(safely_contained).upper(),
                "size_matches": str(size_matches).upper(),
                "md5_matches": str(md5_matches).upper(),
                "status": "PASS" if passed else "FAIL",
            }
        )

    audit_columns = (
        "input_root_id",
        "relative_path",
        "absolute_path",
        "expected_size_bytes",
        "actual_size_bytes",
        "expected_md5",
        "actual_md5",
        "exists",
        "contained_in_root",
        "size_matches",
        "md5_matches",
        "status",
    )
    write_tsv(audit_table_path, audit_columns, audit_rows)

    report_rows = [
        {"check": "validation_phase", "value": arguments.phase, "status": "INFO"},
        {
            "check": "baseline_file",
            "value": str(baseline_path),
            "status": "INFO",
        },
        {
            "check": "expected_input_files",
            "value": len(baseline),
            "status": "PASS",
        },
        {
            "check": "validated_input_files",
            "value": sum(row["status"] == "PASS" for row in audit_rows),
            "status": "PASS" if not failures else "FAIL",
        },
        {
            "check": "input_validation_failures",
            "value": len(failures),
            "status": "PASS" if not failures else "FAIL",
        },
    ]
    write_tsv(report_path, ("check", "value", "status"), report_rows)

    if failures:
        preview = "\n".join(failures[:50])
        remainder = len(failures) - 50
        if remainder > 0:
            preview += f"\n... and {remainder} additional failure(s)"
        raise RuntimeError(
            "Scientific-input MD5 validation failed:\n" + preview
        )
    print(
        f"Scientific-input MD5 validation passed ({arguments.phase}): "
        f"{len(baseline)} files"
    )


if __name__ == "__main__":
    main()
