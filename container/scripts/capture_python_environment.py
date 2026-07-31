#!/usr/bin/env python3
"""Capture a Python environment using syntax compatible with Python 3.6."""

import csv
import json
import os
import platform
import site
import subprocess
import sys


def write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def command_output(command):
    try:
        return subprocess.check_output(command, stderr=subprocess.STDOUT, universal_newlines=True)
    except Exception as exc:
        return "ERROR: {0}\n".format(exc)


def main():
    if len(sys.argv) < 3:
        raise SystemExit("Usage: capture_python_environment.py OUT_DIR DIRECT_MODULES_TSV")
    out_dir = sys.argv[1]
    direct_modules_path = sys.argv[2]
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    details = [
        ("sys.executable", sys.executable),
        ("sys.version", sys.version.replace("\n", " ")),
        ("implementation", platform.python_implementation()),
        ("implementation_version", platform.python_version()),
        ("platform", platform.platform()),
        ("machine", platform.machine()),
        ("prefix", sys.prefix),
        ("base_prefix", getattr(sys, "base_prefix", "")),
        ("real_prefix", getattr(sys, "real_prefix", "")),
        ("user_site", site.getusersitepackages()),
    ]
    write_tsv(
        os.path.join(out_dir, "python-version-details.tsv"),
        ["field", "value"],
        [{"field": key, "value": value} for key, value in details],
    )
    write_tsv(
        os.path.join(out_dir, "python-sys-path.tsv"),
        ["order", "path"],
        [{"order": i + 1, "path": value} for i, value in enumerate(sys.path)],
    )
    with open(os.path.join(out_dir, "python-version.txt"), "w", encoding="utf-8") as handle:
        handle.write(sys.version)
        handle.write("\n")
    with open(os.path.join(out_dir, "pip-freeze.txt"), "w", encoding="utf-8") as handle:
        handle.write(command_output([sys.executable, "-m", "pip", "freeze", "--all"]))
    with open(os.path.join(out_dir, "pip-list.json"), "w", encoding="utf-8") as handle:
        handle.write(command_output([sys.executable, "-m", "pip", "list", "--format=json"]))

    distributions = []
    import_map = {}
    try:
        import pkg_resources

        for dist in sorted(pkg_resources.working_set, key=lambda item: item.project_name.lower()):
            requirements = []
            try:
                requirements = [str(req) for req in dist.requires()]
            except Exception:
                pass
            top_levels = []
            try:
                if dist.has_metadata("top_level.txt"):
                    top_levels = [
                        line.strip()
                        for line in dist.get_metadata_lines("top_level.txt")
                        if line.strip()
                    ]
            except Exception:
                pass
            for module in top_levels:
                import_map.setdefault(module, []).append(dist.project_name)
            distributions.append(
                {
                    "distribution": dist.project_name,
                    "version": dist.version,
                    "location": dist.location,
                    "requires": ";".join(requirements),
                    "top_level_modules": ";".join(top_levels),
                }
            )
    except Exception as exc:
        distributions.append(
            {
                "distribution": "ERROR",
                "version": "",
                "location": "",
                "requires": str(exc),
                "top_level_modules": "",
            }
        )

    write_tsv(
        os.path.join(out_dir, "python-installed-distributions.tsv"),
        ["distribution", "version", "location", "requires", "top_level_modules"],
        distributions,
    )

    direct_rows = []
    with open(direct_modules_path, "r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            module = row["top_level_module"]
            mapped = import_map.get(module, [])
            import_ok = False
            import_error = ""
            try:
                __import__(module)
                import_ok = True
            except Exception as exc:
                import_error = str(exc)
            direct_rows.append(
                {
                    "dependency_scope": row["dependency_scope"],
                    "top_level_module": module,
                    "mapped_distributions": ";".join(sorted(set(mapped))),
                    "import_ok": str(import_ok),
                    "import_error": import_error,
                }
            )
    write_tsv(
        os.path.join(out_dir, "python-direct-module-resolution.tsv"),
        [
            "dependency_scope",
            "top_level_module",
            "mapped_distributions",
            "import_ok",
            "import_error",
        ],
        direct_rows,
    )
    with open(os.path.join(out_dir, "python-environment.json"), "w", encoding="utf-8") as handle:
        json.dump(
            {
                "executable": sys.executable,
                "version": sys.version,
                "path": sys.path,
                "prefix": sys.prefix,
                "base_prefix": getattr(sys, "base_prefix", ""),
            },
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

