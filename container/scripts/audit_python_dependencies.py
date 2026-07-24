#!/usr/bin/env python3
"""Inventory Python scripts, imports, local-module edges, and syntax floors."""

from __future__ import annotations

import ast
import csv
import os
import sys
from pathlib import Path


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def classify(path: str) -> str:
    if path.startswith("oxygen/code/O2_supply_demand_MAP/analysis/"):
        return "analysis"
    if path.startswith("oxygen/code/O2_supply_demand_MAP/"):
        return "o2_workflow"
    if path.startswith("oxygen/"):
        return "oxygen_auxiliary"
    return "repository_python_other"


def syntax_floor(source: str) -> str:
    explicit_floor = 7 if "from __future__ import annotations" in source else 6
    current_minor = sys.version_info.minor
    for minor in range(explicit_floor, current_minor + 1):
        try:
            ast.parse(source, feature_version=(3, minor))
            return f"3.{minor}"
        except (SyntaxError, ValueError):
            continue
    return f">{sys.version_info.major}.{current_minor}"


def imported_modules(tree: ast.AST) -> list[tuple[int, str, str]]:
    rows: list[tuple[int, str, str]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                rows.append((node.lineno, alias.name, "import"))
        elif isinstance(node, ast.ImportFrom):
            prefix = "." * node.level
            rows.append((node.lineno, prefix + (node.module or ""), "from"))
    return rows


def main() -> int:
    if len(sys.argv) < 3:
        raise SystemExit("Usage: audit_python_dependencies.py PROJECT_ROOT OUT_DIR")
    root = Path(sys.argv[1]).resolve()
    out_dir = Path(sys.argv[2])

    py_files: list[Path] = []
    pruned_names = {".git", "container", "results", "__pycache__", ".venv", "node_modules"}
    for current, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(name for name in dirnames if name not in pruned_names)
        current_path = Path(current)
        for filename in sorted(filenames):
            if filename.endswith(".py"):
                py_files.append(current_path / filename)
    py_files.sort()
    module_names = {path.stem for path in py_files}
    package_names = {path.parent.name for path in py_files if (path.parent / "__init__.py").exists()}
    local_names = module_names | package_names
    stdlib = set(getattr(sys, "stdlib_module_names", ()))

    script_rows: list[dict[str, object]] = []
    import_rows: list[dict[str, object]] = []
    for path in py_files:
        rel = path.relative_to(root).as_posix()
        source = path.read_text(encoding="utf-8")
        lines = source.splitlines()
        shebang = lines[0] if lines and lines[0].startswith("#!") else ""
        syntax_error = ""
        try:
            tree = ast.parse(source, filename=rel)
            imports = imported_modules(tree)
        except SyntaxError as exc:
            imports = []
            syntax_error = f"{exc.msg} at line {exc.lineno}"

        script_rows.append(
            {
                "file": rel,
                "scope": classify(rel),
                "shebang": shebang,
                "syntax_floor": syntax_floor(source),
                "parse_status": "ok" if not syntax_error else "error",
                "parse_error": syntax_error,
                "size_bytes": path.stat().st_size,
            }
        )
        for line, module, style in imports:
            top = module.lstrip(".").split(".", 1)[0]
            if module.startswith(".") or top in local_names:
                kind = "repository_local"
            elif top in stdlib:
                kind = "stdlib"
            else:
                kind = "third_party_or_unresolved"
            import_rows.append(
                {
                    "file": rel,
                    "line": line,
                    "scope": classify(rel),
                    "module": module,
                    "top_level_module": top,
                    "import_style": style,
                    "classification": kind,
                }
            )

    write_tsv(
        out_dir / "python-script-inventory.tsv",
        ["file", "scope", "shebang", "syntax_floor", "parse_status", "parse_error", "size_bytes"],
        script_rows,
    )
    (out_dir / "python-files.txt").write_text(
        "".join(f"{row['file']}\n" for row in script_rows),
        encoding="utf-8",
    )
    write_tsv(
        out_dir / "python-imports-by-file.tsv",
        ["file", "line", "scope", "module", "top_level_module", "import_style", "classification"],
        import_rows,
    )

    direct_rows: list[dict[str, object]] = []
    for scope_group, scopes in (
        ("python_o2_workflow", {"analysis", "o2_workflow", "oxygen_auxiliary"}),
        ("python_repository_all", {row["scope"] for row in script_rows}),
    ):
        modules = sorted(
            {
                str(row["top_level_module"])
                for row in import_rows
                if row["scope"] in scopes and row["classification"] == "third_party_or_unresolved"
            }
        )
        for module in modules:
            files = sorted(
                {
                    str(row["file"])
                    for row in import_rows
                    if row["scope"] in scopes and row["top_level_module"] == module
                }
            )
            direct_rows.append(
                {
                    "dependency_scope": scope_group,
                    "top_level_module": module,
                    "referencing_files": ";".join(files),
                }
            )
    write_tsv(
        out_dir / "python-direct-modules.tsv",
        ["dependency_scope", "top_level_module", "referencing_files"],
        direct_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
