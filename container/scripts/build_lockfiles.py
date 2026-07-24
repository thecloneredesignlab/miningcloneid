#!/usr/bin/env python3
"""Build compact R/Python lock artifacts from the captured HPC snapshot."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from collections import defaultdict, deque
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def clean(value: str | None) -> str:
    if value is None or value in {"", "NA", "nan", "None"}:
        return ""
    return value


def parse_runtime_edges(rows: list[dict[str, str]]) -> dict[str, set[str]]:
    edges: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["dependency_type"] in {"Depends", "Imports", "LinkingTo"}:
            edges[row["parent"]].add(row["dependency"])
    return edges


def dependency_closure(roots: set[str], edges: dict[str, set[str]]) -> set[str]:
    seen = set(roots)
    queue = deque(sorted(roots))
    while queue:
        parent = queue.popleft()
        for child in sorted(edges.get(parent, ())):
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return seen


def normalized_distribution_name(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def requirements_from_report(
    report_path: Path,
    output_path: Path,
    required_distributions: set[str] | None = None,
) -> None:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    lines = [
        "# Resolved for CPython 3.9, Linux x86_64, manylinux2014.",
        f"# Source report: {report_path.name}",
    ]
    installs = sorted(
        report.get("install", []),
        key=lambda item: item.get("metadata", {}).get("name", "").lower(),
    )
    installed_distributions = {
        normalized_distribution_name(item.get("metadata", {}).get("name", ""))
        for item in installs
    }
    missing = {
        normalized_distribution_name(name)
        for name in (required_distributions or set())
    } - installed_distributions
    if missing:
        raise RuntimeError(
            f"{report_path.name} is missing target-platform distributions: "
            + ", ".join(sorted(missing))
        )
    for item in installs:
        metadata = item.get("metadata", {})
        name = metadata.get("name", "")
        version = metadata.get("version", "")
        hashes = item.get("download_info", {}).get("archive_info", {}).get("hashes", {})
        sha256 = hashes.get("sha256", "")
        line = f"{name}=={version}"
        if sha256:
            line += f" --hash=sha256:{sha256}"
        lines.append(line)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def python_closure_rows(report_path: Path, scope: str) -> list[dict[str, object]]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    rows = []
    for item in report.get("install", []):
        metadata = item.get("metadata", {})
        archive = item.get("download_info", {}).get("archive_info", {})
        rows.append(
            {
                "dependency_scope": scope,
                "distribution": metadata.get("name", ""),
                "version": metadata.get("version", ""),
                "direct_request": str(bool(item.get("requested", False))).upper(),
                "requires_dist": ";".join(metadata.get("requires_dist") or []),
                "download_url": item.get("download_info", {}).get("url", ""),
                "sha256": archive.get("hashes", {}).get("sha256", ""),
            }
        )
    return sorted(rows, key=lambda row: str(row["distribution"]).lower())


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: build_lockfiles.py CONTAINER_ROOT")
    root = Path(sys.argv[1]).resolve()
    locks = root / "locks"
    locks.mkdir(parents=True, exist_ok=True)
    r_root = root / "hpc_snapshot" / "r"

    selected_rows = read_tsv(r_root / "r-installed-packages-selected.tsv")
    selected = {row["Package"]: row for row in selected_rows}
    direct_rows = read_tsv(root / "manifests" / "r-direct-packages.tsv")
    direct_scopes: dict[str, set[str]] = defaultdict(set)
    for row in direct_rows:
        direct_scopes[row["package"]].add(row["dependency_scope"])

    runtime_rows = read_tsv(r_root / "r-runtime-dependency-closure.tsv")
    analysis_rows = read_tsv(r_root / "r-analysis-runtime-dependency-closure.tsv")
    runtime_packages = {row["package"] for row in runtime_rows}
    analysis_packages = {row["package"] for row in analysis_rows}
    edge_rows = read_tsv(r_root / "r-installed-package-dependency-edges.tsv")
    runtime_edges = parse_runtime_edges(edge_rows)

    # magick is called by current report code but is absent on the captured HPC node.
    # Its current CRAN 2.9.1 runtime imports are added explicitly for the Docker target.
    magick_dependencies = {"Rcpp", "magrittr", "curl"}
    target_packages = set(runtime_packages)
    target_packages.add("magick")
    target_packages.update(dependency_closure(magick_dependencies, runtime_edges))

    package_rows: list[dict[str, object]] = []
    for package in sorted(target_packages, key=str.lower):
        observed = selected.get(package)
        is_magick = package == "magick"
        version = "2.9.1" if is_magick else clean(observed.get("Version") if observed else "")
        package_rows.append(
            {
                "package": package,
                "observed_hpc_version": clean(observed.get("Version") if observed else ""),
                "target_version": version,
                "status": "target_only_missing_on_hpc" if is_magick else "observed_on_hpc",
                "direct_scopes": ";".join(sorted(direct_scopes.get(package, ()))),
                "analysis_runtime": str(package in analysis_packages).upper(),
                "o2_runtime": str(package in runtime_packages).upper(),
                "source": "Repository" if is_magick else (
                    "R" if clean(observed.get("Priority") if observed else "") else
                    ("GitHub" if clean(observed.get("RemoteType") if observed else "").lower() == "github" else "Repository")
                ),
                "repository": "CRAN" if is_magick else clean(observed.get("Repository") if observed else ""),
                "remote_type": clean(observed.get("RemoteType") if observed else ""),
                "remote_repo": clean(observed.get("RemoteRepo") if observed else ""),
                "remote_ref": clean(observed.get("RemoteRef") if observed else ""),
                "remote_sha": clean(observed.get("RemoteSha") if observed else ""),
                "priority": clean(observed.get("Priority") if observed else ""),
                "built": clean(observed.get("Built") if observed else ""),
                "needs_compilation": "yes" if is_magick else clean(observed.get("NeedsCompilation") if observed else ""),
                "system_requirements": "ImageMagick++: Magick++-config (rpm: ImageMagick-c++-devel, deb: libmagick++-dev)" if is_magick else clean(observed.get("SystemRequirements") if observed else ""),
                "source_md5": "012ecba2c6bffc0d79f56658b3cb808c" if is_magick else "",
                "notes": "Referenced by current report code; not installed in captured HPC R libraries." if is_magick else "",
            }
        )

    package_fields = [
        "package", "observed_hpc_version", "target_version", "status", "direct_scopes",
        "analysis_runtime", "o2_runtime", "source", "repository", "remote_type",
        "remote_repo", "remote_ref", "remote_sha", "priority", "built",
        "needs_compilation", "system_requirements", "source_md5", "notes",
    ]
    write_tsv(locks / "packages.lock.tsv", package_fields, package_rows)

    renv_packages: dict[str, dict[str, object]] = {}
    for row in package_rows:
        package = str(row["package"])
        version = str(row["target_version"])
        if not version:
            continue
        record: dict[str, object] = {
            "Package": package,
            "Version": version,
            "Source": row["source"],
        }
        if row["source"] == "Repository":
            record["Repository"] = "CRAN"
        requirements = set(runtime_edges.get(package, ()))
        if package == "magick":
            requirements.update(magick_dependencies)
        requirements.intersection_update(target_packages)
        if requirements:
            record["Requirements"] = sorted(requirements)
        observed = selected.get(package, {})
        for source_key, target_key in (
            ("RemoteType", "RemoteType"),
            ("RemoteHost", "RemoteHost"),
            ("RemoteRepo", "RemoteRepo"),
            ("RemoteUsername", "RemoteUsername"),
            ("RemoteRef", "RemoteRef"),
            ("RemoteSha", "RemoteSha"),
        ):
            value = clean(observed.get(source_key))
            if value:
                record[target_key] = value
        renv_packages[package] = record

    renv_lock = {
        "R": {
            "Version": "4.4.2",
            "Repositories": [
                {"Name": "CRAN", "URL": "https://cloud.r-project.org"},
            ],
        },
        "Packages": dict(sorted(renv_packages.items(), key=lambda item: item[0].lower())),
    }
    (locks / "renv.lock").write_text(
        json.dumps(renv_lock, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )

    observed_python = read_tsv(root / "hpc_snapshot" / "python" / "python-installed-distributions.tsv")
    observed_lines = [
        "# Exact distributions observed in /usr/bin/python3 on the HPC compute node.",
    ]
    for row in sorted(observed_python, key=lambda item: item["distribution"].lower()):
        if row["distribution"] != "ERROR":
            observed_lines.append(f"{row['distribution']}=={row['version']}")
    (locks / "requirements-hpc-observed.txt").write_text(
        "\n".join(observed_lines) + "\n",
        encoding="utf-8",
    )

    requirements_from_report(
        locks / "python-o2-resolution-report.json",
        locks / "requirements-o2-target.lock.txt",
    )
    requirements_from_report(
        locks / "python-repository-resolution-report.json",
        locks / "requirements-repository-all-target.lock.txt",
        # pip's --python-version selects compatible wheels but evaluates
        # dependency environment markers against the resolver interpreter.
        # These Python <3.10 requirements therefore must be explicit inputs
        # when the report is resolved from a newer Python.
        required_distributions={"importlib-resources", "zipp"},
    )
    python_closure = python_closure_rows(
        locks / "python-o2-resolution-report.json", "python_o2_workflow"
    ) + python_closure_rows(
        locks / "python-repository-resolution-report.json", "python_repository_all"
    )
    write_tsv(
        root / "manifests" / "python-dependency-closure.tsv",
        [
            "dependency_scope", "distribution", "version", "direct_request",
            "requires_dist", "download_url", "sha256",
        ],
        python_closure,
    )

    analysis_unique: dict[str, dict[str, object]] = {}
    for row in analysis_rows:
        package = row["package"]
        depth = int(row["depth"])
        current = analysis_unique.get(package)
        if current is None or depth < int(current["minimum_depth"]):
            analysis_unique[package] = {
                "package": package,
                "version": row["version"],
                "minimum_depth": depth,
                "direct_analysis_reference": str(package in direct_scopes and "analysis_direct" in direct_scopes[package]).upper(),
                "installed": row["installed"],
                "library_path": row["library_path"],
            }
    write_tsv(
        root / "manifests" / "r-analysis-transitive-packages.tsv",
        [
            "package", "version", "minimum_depth", "direct_analysis_reference",
            "installed", "library_path",
        ],
        [analysis_unique[key] for key in sorted(analysis_unique, key=str.lower)],
    )

    gaps = [
        {
            "component": "R package",
            "name": "magick",
            "observed_hpc_state": "missing",
            "target_state": "pin CRAN 2.9.1",
            "impact": "Report image conversion paths cannot use magick on the captured HPC environment.",
        },
        {
            "component": "Python package set",
            "name": "O2 workflow",
            "observed_hpc_state": "openpyxl missing",
            "target_state": "openpyxl and et-xmlfile pinned in requirements-o2-target.lock.txt",
            "impact": "Warm-start XLSX output fails unless --no-xlsx is used.",
        },
        {
            "component": "Python package set",
            "name": "repository-wide scripts",
            "observed_hpc_state": "numpy scipy matplotlib PyYAML openpyxl missing",
            "target_state": "full resolved lock in requirements-repository-all-target.lock.txt",
            "impact": "Non-O2 repository Python analysis scripts are not runnable in default HPC python3.",
        },
        {
            "component": "External executable",
            "name": "ImageMagick magick",
            "observed_hpc_state": "missing",
            "target_state": "install ImageMagick in Docker",
            "impact": "oxygen/figures/assemble_iteration_panels.py is not runnable on captured HPC PATH.",
        },
        {
            "component": "External executable",
            "name": "Pandoc",
            "observed_hpc_state": "missing",
            "target_state": "install Pandoc in Docker",
            "impact": "R Markdown rendering requires a separately available Pandoc.",
        },
    ]
    write_tsv(
        root / "manifests" / "environment-gaps.tsv",
        ["component", "name", "observed_hpc_state", "target_state", "impact"],
        gaps,
    )

    module_text = (root / "hpc_snapshot" / "modules" / "module-list-after-R-4.4.txt").read_text(
        encoding="utf-8"
    )
    module_rows = []
    for line in module_text.splitlines():
        match = re.match(r"\s*(\d+)\)\s+(\S+)", line)
        if match:
            module_rows.append({"order": match.group(1), "module": match.group(2)})
    write_tsv(root / "manifests" / "hpc-loaded-modules.tsv", ["order", "module"], module_rows)

    owner_rows = read_tsv(
        root / "hpc_snapshot" / "shared_libraries" / "shared-library-rpm-owners.tsv"
    )
    rpm_owners = {row["library_path"]: row["rpm_owner"] for row in owner_rows}
    runtime_library_rows = []
    ldd_rows = read_tsv(
        root / "hpc_snapshot" / "shared_libraries" / "r-package-shared-libraries-ldd.tsv"
    )
    for row in ldd_rows:
        consumer_path = row["path"]
        consumer_package = ""
        for package in target_packages:
            if f"/{package}/" in consumer_path:
                consumer_package = package
                break
        if not consumer_package:
            continue
        line = row["ldd_output"].strip()
        dependency = ""
        resolved_path = ""
        status = "resolved"
        if " => " in line:
            dependency, remainder = line.split(" => ", 1)
            first = remainder.split(" ", 1)[0]
            if first == "not":
                status = "not_found"
            elif first.startswith("/"):
                resolved_path = first
        else:
            first = line.split(" ", 1)[0] if line else ""
            dependency = Path(first).name
            if first.startswith("/"):
                resolved_path = first
        runtime_library_rows.append(
            {
                "r_package": consumer_package,
                "consumer_path": consumer_path,
                "dependency": dependency.strip(),
                "resolved_path": resolved_path,
                "status": status,
                "rpm_owner": rpm_owners.get(resolved_path, ""),
            }
        )
    write_tsv(
        root / "manifests" / "system-runtime-libraries.tsv",
        ["r_package", "consumer_path", "dependency", "resolved_path", "status", "rpm_owner"],
        runtime_library_rows,
    )

    metadata = dict(
        (row["field"], row["value"])
        for row in read_tsv(root / "hpc_snapshot" / "python" / "python-version-details.tsv")
    )
    capture = dict(
        line.rstrip("\n").split("\t", 1)
        for line in (root / "hpc_snapshot" / "capture_metadata.tsv").read_text(encoding="utf-8").splitlines()
        if "\t" in line
    )
    os_release = {}
    for line in (root / "hpc_snapshot" / "system" / "os-release.txt").read_text(
        encoding="utf-8"
    ).splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            os_release[key] = value.strip('"')
    code_summary = {
        row["metric"]: row["value"]
        for row in read_tsv(root / "manifests" / "code-audit-summary.tsv")
    }
    environment_summary = [
        {"field": "capture_utc", "value": capture.get("capture_utc", "")},
        {"field": "slurm_job_id", "value": capture.get("slurm_job_id", "")},
        {"field": "compute_hostname", "value": capture.get("hostname", "")},
        {"field": "operating_system", "value": os_release.get("PRETTY_NAME", "")},
        {"field": "architecture", "value": metadata.get("machine", "")},
        {"field": "glibc", "value": (root / "hpc_snapshot" / "system" / "glibc.txt").read_text(encoding="utf-8").splitlines()[-1]},
        {"field": "requested_r_module", "value": capture.get("requested_r_module", "")},
        {"field": "resolved_r_module", "value": next((row["module"] for row in module_rows if row["module"].startswith("R/")), "")},
        {"field": "r_version", "value": (root / "hpc_snapshot" / "r" / "R-version.txt").read_text(encoding="utf-8").strip()},
        {"field": "python_executable", "value": metadata.get("sys.executable", "")},
        {"field": "python_version", "value": metadata.get("implementation_version", "")},
        {"field": "analysis_R_files", "value": code_summary.get("analysis_R_files", "")},
        {"field": "analysis_Python_files", "value": code_summary.get("analysis_Python_files", "")},
        {"field": "repository_Python_files", "value": str(len(read_tsv(root / "manifests" / "python-script-inventory.tsv")))},
        {"field": "analysis_runtime_R_packages", "value": str(sum(row["analysis_runtime"] == "TRUE" for row in package_rows))},
        {"field": "o2_runtime_R_packages", "value": str(sum(row["o2_runtime"] == "TRUE" for row in package_rows))},
        {"field": "loaded_modules", "value": str(len(module_rows))},
    ]
    write_tsv(root / "manifests" / "environment-summary.tsv", ["field", "value"], environment_summary)

    code_hash_rows = []
    for row in read_tsv(root / "manifests" / "code-scope.tsv"):
        path = root.parent / row["file"]
        if not path.is_file():
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        code_hash_rows.append(
            {
                "sha256": digest,
                "file": row["file"],
                "scope": row["scope"],
            }
        )
    write_tsv(root / "manifests" / "code-sha256.tsv", ["sha256", "file", "scope"], code_hash_rows)

    checksum_lines = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name == "SHA256SUMS" or "__pycache__" in path.parts:
            continue
        checksum_lines.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(root).as_posix()}")
    (root / "SHA256SUMS").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
