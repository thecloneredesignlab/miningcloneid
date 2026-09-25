#!/usr/bin/env python3
"""Build an auditable, shareable bundle from joint and separate O2 fits.

The script is a consume-only exporter. It validates recorded objectives and
source relationships, copies existing numerical/provenance outputs, and writes
file/checksum manifests. It never reads an RDS payload or recalculates a fit.
"""

import argparse
import csv
import hashlib
import math
import os
import re
import shutil
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[4]

DEFAULT_JOINT_ROOT = (
    PROJECT_ROOT
    / "oxygen/results/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540"
)
DEFAULT_INVIVO_ROOT = PROJECT_ROOT / "oxygen/results/fit_invivo_O2_buffering_500seed"
DEFAULT_INVITRO_ROOT = PROJECT_ROOT / "oxygen/results/fit_invitro_O2_buffering_500seed"
DEFAULT_POOLED_EMBEDDING_ROOT = (
    PROJECT_ROOT
    / "oxygen/results/analysis/best_fit_parameter_feature/04_combine_parameter_landscape"
    / "pooled_embedding_curve_class"
)
DEFAULT_FIXO2_EIGEN_EMBEDDING_ROOT = (
    PROJECT_ROOT
    / "oxygen/results/analysis/best_fit_parameter_feature/06_combine_FixO2_eigen_attractor"
    / "pooled_embedding_curve_class"
)
DEFAULT_FIXED_O2_CLASSIFICATION_TABLES = (
    PROJECT_ROOT
    / "oxygen/results/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification"
    / "monotonicity_classification/dense-grid_monotonicity_classification/tables"
)
DEFAULT_FIXED_O2_REGRESSION_CLASSIFICATION_TABLES = (
    PROJECT_ROOT
    / "oxygen/results/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification"
    / "monotonicity_classification/dense-grid_monotonicity_regression_classification/tables"
)
DEFAULT_FIXED_O2_ATTRACTOR_TABLES = (
    PROJECT_ROOT
    / "oxygen/results/analysis/best_fit_parameter_feature/01_fixed_o2"
    / "FixO2_invivo_500seed/attractors/tables"
)
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "oxygen/results/fitting_output_bundle_20260722"

POOLED_EMBEDDING_REQUIRED_FILES = (
    "pooled_embedding_curve_class_report.html",
    "tables/pooled_embedding_curve_class_manifest.tsv",
)

FIXO2_EIGEN_EMBEDDING_REQUIRED_FILES = (
    "fixo2_eigen_attractor_embedding_curve_class_report.html",
    "tables/pooled_embedding_curve_class_analysis_manifest.tsv",
    "tables/pooled_embedding_curve_class_manifest.tsv",
)

FIXED_O2_CLASSIFICATION_REQUIRED_FILES = (
    "fixed_o2_ploidy_monotonicity_regression_curves.tsv",
    "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv",
    "fixed_o2_ploidy_monotonicity_regression_run_arguments.tsv",
    "fixed_o2_ploidy_monotonicity_regression_validation.tsv",
)

FIXED_O2_ATTRACTOR_REQUIRED_FILES = (
    "fixed_o2_attractors_by_seed.tsv",
    "fixed_o2_attractor_mode_by_seed_o2.tsv",
    "analysis_run_arguments.tsv",
)

SUPPORTING_SOURCE_COLUMNS = (
    "artifact",
    "report",
    "direct_report_input",
    "source_root",
    "bundle_root",
    "file_count",
    "total_bytes",
)

JOINT_OVERVIEW_FILES = (
    "multi_warmup_best_seed_summary.tsv",
    "multi_warmup_manifest.tsv",
    ".multi_warmup_task_table_runs.tsv",
    "multi_warmup_integrated_seed_manifest.tsv",
    "multi_warmup_seed_plan_mode.tsv",
    "multi_warmup_task_table_jobs.tsv",
    "landscape_subcluster_selected_representatives.tsv",
    "landscape_subcluster_pair_settings.tsv",
    "landscape_subcluster_invitro_best_anchor.tsv",
    "multi_warmup_final_basin_assignments.tsv",
    "multi_warmup_integrated_objective_tradeoff.tsv",
    "multi-warm-up_results.html",
)

JOINT_REQUIRED_OVERVIEW_FILES = (
    "multi_warmup_best_seed_summary.tsv",
    "multi_warmup_manifest.tsv",
    "multi_warmup_integrated_seed_manifest.tsv",
)

JOINT_RUN_FILES = (
    "config.input.yaml",
    "config.resolved.yaml",
    "run_command.txt",
    "run_effective_args.tsv",
    "run_provenance.tsv",
    "fit_array_manifest.tsv",
    "parameter_table_input.csv",
)

JOINT_REQUIRED_RUN_FILES = (
    "config.input.yaml",
    "config.resolved.yaml",
    "run_effective_args.tsv",
    "run_provenance.tsv",
    "fit_array_manifest.tsv",
    "parameter_table_input.csv",
)

EXTRA_RESULT_FILES = (
    "seed_summary.tsv",
    "seed_objective_simple.tsv",
    "objective_components_long.tsv",
    "convergence_summary.tsv",
    "extra_results_report.html",
)

SEPARATE_RUN_FILES = (
    "config.input.yaml",
    "config.resolved.yaml",
    "fit_array_manifest.tsv",
    "parameter_table.csv",
    "parameter_table_input.csv",
    "extra_results_run.log",
)

SEPARATE_REQUIRED_FILES = (
    "fit_array_manifest.tsv",
    "extra_results/seed_summary.tsv",
)

SEED_REQUIRED_FILES = (
    "best_params.tsv",
    "fit_summary.tsv",
    "parameter_table_input.csv",
)

JOINT_SEED_REQUIRED_FILES = SEED_REQUIRED_FILES + (
    "best_params_transformed.tsv",
    "joint_components.tsv",
    "joint_soft_coupling.tsv",
    "run_effective_args.tsv",
    "run_provenance.tsv",
)

SEED_COPY_SUFFIXES = {
    ".csv",
    ".html",
    ".json",
    ".log",
    ".rda",
    ".rds",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}

JOINT_SUMMARY_COLUMNS = (
    "warmup_label",
    "invivo_seed",
    "invivo_seed_dir",
    "invitro_seed",
    "invitro_seed_dir",
    "joint_run_dir",
    "joint_soft_coupling_parameters_table",
    "best_joint_seed",
    "objective",
    "objective_invivo",
    "objective_invitro",
    "objective_soft_coupling",
)

SELECTED_RESULT_COLUMNS = (
    "record_type",
    "warmup_label",
    "invivo_seed",
    "invitro_seed",
    "selected_seed",
    "objective_metric",
    "objective",
    "objective_invivo",
    "objective_invitro",
    "objective_soft_coupling",
    "source_dir",
    "bundle_dir",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--joint-root", type=Path, default=DEFAULT_JOINT_ROOT)
    parser.add_argument("--invivo-root", type=Path, default=DEFAULT_INVIVO_ROOT)
    parser.add_argument("--invitro-root", type=Path, default=DEFAULT_INVITRO_ROOT)
    parser.add_argument(
        "--pooled-embedding-root", type=Path, default=DEFAULT_POOLED_EMBEDDING_ROOT
    )
    parser.add_argument(
        "--fixo2-eigen-embedding-root", type=Path, default=DEFAULT_FIXO2_EIGEN_EMBEDDING_ROOT
    )
    parser.add_argument(
        "--fixed-o2-classification-tables",
        type=Path,
        default=DEFAULT_FIXED_O2_CLASSIFICATION_TABLES,
    )
    parser.add_argument(
        "--fixed-o2-regression-classification-tables",
        type=Path,
        default=DEFAULT_FIXED_O2_REGRESSION_CLASSIFICATION_TABLES,
    )
    parser.add_argument(
        "--fixed-o2-attractor-tables", type=Path, default=DEFAULT_FIXED_O2_ATTRACTOR_TABLES
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing output directory/archive at the exact output path.",
    )
    parser.add_argument(
        "--no-archive",
        action="store_true",
        help="Build and validate the directory but do not create a tar.gz archive.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    require_file(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[Dict[str, object]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required file is missing: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Required file is empty: {path}")


def require_columns(rows: Sequence[Dict[str, str]], columns: Sequence[str], path: Path) -> None:
    if not rows:
        raise ValueError(f"Required table has no data rows: {path}")
    present = set(rows[0])
    missing = [column for column in columns if column not in present]
    if missing:
        raise ValueError(f"Missing columns in {path}: {', '.join(missing)}")


def as_float(value: object) -> float:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return math.nan


def normalized_seed(value: object) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError("Encountered an empty seed identifier")
    match = re.fullmatch(r"(?:seed)?(\d+)", text, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Invalid seed identifier: {text}")
    return f"seed{int(match.group(1))}"


def ensure_within(path: Path, root: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    root_resolved = root.expanduser().resolve()
    try:
        resolved.relative_to(root_resolved)
    except ValueError as exc:
        raise ValueError(f"{label} is outside the expected root: {resolved} not under {root_resolved}") from exc
    return resolved


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class BundleBuilder:
    def __init__(
        self,
        joint_root: Path,
        invivo_root: Path,
        invitro_root: Path,
        pooled_embedding_root: Path,
        fixo2_eigen_embedding_root: Path,
        fixed_o2_classification_tables: Path,
        fixed_o2_regression_classification_tables: Path,
        fixed_o2_attractor_tables: Path,
        output_dir: Path,
        force: bool,
        create_archive: bool,
    ) -> None:
        self.joint_root = joint_root.expanduser().resolve()
        self.invivo_root = invivo_root.expanduser().resolve()
        self.invitro_root = invitro_root.expanduser().resolve()
        self.pooled_embedding_root = pooled_embedding_root.expanduser().resolve()
        self.fixo2_eigen_embedding_root = fixo2_eigen_embedding_root.expanduser().resolve()
        self.fixed_o2_classification_tables = fixed_o2_classification_tables.expanduser().resolve()
        self.fixed_o2_regression_classification_tables = (
            fixed_o2_regression_classification_tables.expanduser().resolve()
        )
        self.fixed_o2_attractor_tables = fixed_o2_attractor_tables.expanduser().resolve()
        self.output_dir = output_dir.expanduser().resolve()
        self.archive_path = Path(f"{self.output_dir}.tar.gz")
        self.archive_checksum_path = Path(f"{self.archive_path}.sha256")
        self.force = force
        self.create_archive = create_archive
        self.stage_dir = self.output_dir.with_name(f".{self.output_dir.name}.tmp-{os.getpid()}")
        self.copied: List[Dict[str, object]] = []
        self.validation_rows: List[Dict[str, object]] = []
        self.selected_rows: List[Dict[str, object]] = []
        self.supporting_rows: List[Dict[str, object]] = []

    def build(self) -> None:
        self.validate_roots()
        self.prepare_destination()
        try:
            joint_rows = self.read_joint_summary()
            self.copy_joint_overview()
            self.copy_joint_pairs(joint_rows)
            self.copy_separate_fits(joint_rows)
            self.copy_supporting_analysis()
            self.write_generated_tables(joint_rows)
            self.write_readme(joint_rows)
            self.write_file_manifest()
            self.write_checksums()
            self.finalize_directory()
            if self.create_archive:
                self.write_archive()
        except Exception:
            if self.stage_dir.exists():
                shutil.rmtree(self.stage_dir)
            raise

    def validate_roots(self) -> None:
        for label, root in (
            ("joint root", self.joint_root),
            ("in-vivo root", self.invivo_root),
            ("in-vitro root", self.invitro_root),
            ("pooled embedding report root", self.pooled_embedding_root),
            ("FixO2 eigen embedding report root", self.fixo2_eigen_embedding_root),
            ("fixed-O2 classification tables", self.fixed_o2_classification_tables),
            (
                "fixed-O2 regression-classification tables",
                self.fixed_o2_regression_classification_tables,
            ),
            ("fixed-O2 attractor tables", self.fixed_o2_attractor_tables),
        ):
            if not root.is_dir():
                raise NotADirectoryError(f"{label} does not exist: {root}")
        for name in JOINT_REQUIRED_OVERVIEW_FILES:
            require_file(self.joint_root / name)
        for root in (self.invivo_root, self.invitro_root):
            for name in SEPARATE_REQUIRED_FILES:
                require_file(root / name)
        for root, names in (
            (self.pooled_embedding_root, POOLED_EMBEDDING_REQUIRED_FILES),
            (self.fixo2_eigen_embedding_root, FIXO2_EIGEN_EMBEDDING_REQUIRED_FILES),
            (self.fixed_o2_classification_tables, FIXED_O2_CLASSIFICATION_REQUIRED_FILES),
            (
                self.fixed_o2_regression_classification_tables,
                FIXED_O2_CLASSIFICATION_REQUIRED_FILES,
            ),
            (self.fixed_o2_attractor_tables, FIXED_O2_ATTRACTOR_REQUIRED_FILES),
        ):
            for name in names:
                require_file(root / name)

    def prepare_destination(self) -> None:
        self.output_dir.parent.mkdir(parents=True, exist_ok=True)
        existing = [
            path
            for path in (self.output_dir, self.archive_path, self.archive_checksum_path, self.stage_dir)
            if path.exists()
        ]
        if existing and not self.force:
            joined = "\n".join(str(path) for path in existing)
            raise FileExistsError(f"Output already exists; use --force to replace only these paths:\n{joined}")
        for path in existing:
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
        self.stage_dir.mkdir(parents=True)

    def read_joint_summary(self) -> List[Dict[str, str]]:
        summary_path = self.joint_root / "multi_warmup_best_seed_summary.tsv"
        rows = read_tsv(summary_path)
        require_columns(rows, JOINT_SUMMARY_COLUMNS, summary_path)
        labels = [row["warmup_label"].strip() for row in rows]
        if any(not label for label in labels) or len(set(labels)) != len(labels):
            raise ValueError("Joint warmup labels must be non-empty and unique")
        return rows

    def copy_one(self, source: Path, destination: Path, category: str, required: bool = False) -> bool:
        if not source.is_file():
            if required:
                raise FileNotFoundError(f"Required file is missing: {source}")
            return False
        if source.stat().st_size == 0 and required:
            raise ValueError(f"Required file is empty: {source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        self.copied.append(
            {
                "bundle_path": str(destination.relative_to(self.stage_dir)),
                "source_path": str(source),
                "category": category,
                "bytes": destination.stat().st_size,
                "sha256": sha256_file(destination),
            }
        )
        return True

    def copy_named_files(
        self,
        source_root: Path,
        destination_root: Path,
        names: Iterable[str],
        category: str,
        required_names: Iterable[str] = (),
    ) -> None:
        required = set(required_names)
        for name in names:
            self.copy_one(
                source_root / name,
                destination_root / name,
                category,
                required=name in required,
            )

    def copy_seed_tree(self, source_seed: Path, destination_seed: Path, category: str, joint: bool) -> None:
        required = JOINT_SEED_REQUIRED_FILES if joint else SEED_REQUIRED_FILES
        for name in required:
            require_file(source_seed / name)
        copied_count = 0
        for source in sorted(source_seed.rglob("*")):
            if not source.is_file() or source.suffix.lower() not in SEED_COPY_SUFFIXES:
                continue
            relative = source.relative_to(source_seed)
            self.copy_one(source, destination_seed / relative, category, required=False)
            copied_count += 1
        if copied_count == 0:
            raise ValueError(f"No seed-level files selected from {source_seed}")

    def copy_full_tree(
        self,
        source_root: Path,
        destination_root: Path,
        category: str,
    ) -> Tuple[int, int]:
        source_files = [
            path
            for path in sorted(source_root.rglob("*"))
            if path.is_file() and path.name != ".DS_Store"
        ]
        if not source_files:
            raise ValueError(f"No files found in supporting-analysis root: {source_root}")
        total_bytes = 0
        for source in source_files:
            relative = source.relative_to(source_root)
            self.copy_one(
                source,
                destination_root / relative,
                category,
                required=True,
            )
            total_bytes += source.stat().st_size

        destination_count = sum(1 for path in destination_root.rglob("*") if path.is_file())
        if destination_count != len(source_files):
            raise ValueError(
                "Supporting-analysis copy count mismatch for {}: source={}, destination={}".format(
                    source_root, len(source_files), destination_count
                )
            )
        self.validation_rows.append(
            {
                "check": "supporting_analysis_copy_count",
                "record": category,
                "status": "PASS",
                "expected": str(len(source_files)),
                "observed": str(destination_count),
                "details": str(source_root),
            }
        )
        return len(source_files), total_bytes

    def validate_report_manifest_references(
        self,
        report_root: Path,
        label: str,
        expected_tables: Sequence[str],
    ) -> None:
        manifests = sorted(report_root.rglob("*manifest*.tsv"))
        if not manifests:
            raise FileNotFoundError(f"No report manifest tables found under: {report_root}")
        manifest_text = "\n".join(
            path.read_text(encoding="utf-8", errors="replace") for path in manifests
        )
        missing = [name for name in expected_tables if name not in manifest_text]
        if missing:
            raise ValueError(
                "{} report manifests do not reference required fixed-O2 tables: {}".format(
                    label, ", ".join(missing)
                )
            )
        self.validation_rows.append(
            {
                "check": "report_fixed_o2_table_references",
                "record": label,
                "status": "PASS",
                "expected": str(len(expected_tables)),
                "observed": str(len(expected_tables)),
                "details": ", ".join(str(path) for path in manifests),
            }
        )

    def copy_supporting_analysis(self) -> None:
        self.validate_report_manifest_references(
            self.pooled_embedding_root,
            "pooled_embedding_curve_class",
            ("fixed_o2_ploidy_monotonicity_regression_by_seed.tsv",),
        )
        self.validate_report_manifest_references(
            self.fixo2_eigen_embedding_root,
            "fixo2_eigen_attractor_embedding_curve_class",
            (
                "fixed_o2_ploidy_monotonicity_regression_by_seed.tsv",
                "fixed_o2_ploidy_regression_curve_average_slope_by_seed.tsv",
            ),
        )
        specifications = (
            (
                "pooled_embedding_curve_class",
                "pooled_embedding_curve_class_report.html",
                "no",
                self.pooled_embedding_root,
                self.stage_dir / "supporting_analysis/pooled_embedding_curve_class",
                "supporting_pooled_embedding_report",
            ),
            (
                "fixo2_eigen_attractor_embedding_curve_class",
                "fixo2_eigen_attractor_embedding_curve_class_report.html",
                "no",
                self.fixo2_eigen_embedding_root,
                self.stage_dir
                / "supporting_analysis/fixo2_eigen_attractor_embedding_curve_class",
                "supporting_fixo2_eigen_report",
            ),
            (
                "dense_grid_monotonicity_classification_tables",
                "pooled_embedding_curve_class_report.html",
                "yes",
                self.fixed_o2_classification_tables,
                self.stage_dir
                / "supporting_analysis/fixed_o2_analysis_tables"
                / "dense_grid_monotonicity_classification",
                "supporting_fixed_o2_classification_tables",
            ),
            (
                "analytical_attractor_tables",
                "fixo2_eigen_attractor_embedding_curve_class_report.html",
                "no",
                self.fixed_o2_attractor_tables,
                self.stage_dir
                / "supporting_analysis/fixed_o2_analysis_tables/analytical_attractors",
                "supporting_fixed_o2_attractor_tables",
            ),
            (
                "dense_grid_monotonicity_regression_classification_tables",
                "fixo2_eigen_attractor_embedding_curve_class_report.html",
                "yes",
                self.fixed_o2_regression_classification_tables,
                self.stage_dir
                / "supporting_analysis/fixed_o2_analysis_tables"
                / "dense_grid_monotonicity_regression_classification",
                "supporting_fixed_o2_regression_classification_tables",
            ),
        )
        for artifact, report, direct_input, source_root, destination_root, category in specifications:
            file_count, total_bytes = self.copy_full_tree(
                source_root, destination_root, category
            )
            self.supporting_rows.append(
                {
                    "artifact": artifact,
                    "report": report,
                    "direct_report_input": direct_input,
                    "source_root": str(source_root),
                    "bundle_root": str(destination_root.relative_to(self.stage_dir)),
                    "file_count": file_count,
                    "total_bytes": total_bytes,
                }
            )

    def copy_joint_overview(self) -> None:
        destination = self.stage_dir / "joint_fit/overview"
        self.copy_named_files(
            self.joint_root,
            destination,
            JOINT_OVERVIEW_FILES,
            "joint_overview",
            JOINT_REQUIRED_OVERVIEW_FILES,
        )

    def copy_joint_pairs(self, rows: Sequence[Dict[str, str]]) -> None:
        for row in rows:
            label = row["warmup_label"].strip()
            pair_destination = self.stage_dir / "joint_fit/pairs" / label
            run_dir = ensure_within(Path(row["joint_run_dir"]), self.joint_root, "joint_run_dir")
            if not run_dir.is_dir():
                raise NotADirectoryError(f"Joint run directory does not exist: {run_dir}")
            best_seed = normalized_seed(row["best_joint_seed"])
            best_seed_dir = ensure_within(run_dir / best_seed, run_dir, "best joint seed directory")

            self.copy_named_files(
                run_dir,
                pair_destination / "run",
                JOINT_RUN_FILES,
                "joint_pair_run",
                JOINT_REQUIRED_RUN_FILES,
            )
            self.copy_named_files(
                run_dir / "extra_results",
                pair_destination / "run/extra_results",
                EXTRA_RESULT_FILES,
                "joint_pair_summary",
            )

            parameter_table = ensure_within(
                Path(row["joint_soft_coupling_parameters_table"]),
                self.joint_root,
                "joint parameter table",
            )
            self.copy_one(
                parameter_table,
                pair_destination / "run" / parameter_table.name,
                "joint_pair_run",
                required=True,
            )

            self.copy_seed_tree(
                best_seed_dir,
                pair_destination / "best_joint_seed" / best_seed,
                "joint_best_seed",
                joint=True,
            )
            self.validate_joint_objective(row, best_seed_dir)
            self.selected_rows.append(
                {
                    "record_type": "joint_pair_best",
                    "warmup_label": label,
                    "invivo_seed": normalized_seed(row["invivo_seed"]),
                    "invitro_seed": normalized_seed(row["invitro_seed"]),
                    "selected_seed": best_seed,
                    "objective_metric": "objective",
                    "objective": row["objective"],
                    "objective_invivo": row["objective_invivo"],
                    "objective_invitro": row["objective_invitro"],
                    "objective_soft_coupling": row["objective_soft_coupling"],
                    "source_dir": str(best_seed_dir),
                    "bundle_dir": str(
                        (Path("joint_fit/pairs") / label / "best_joint_seed" / best_seed)
                    ),
                }
            )

    def validate_joint_objective(self, row: Dict[str, str], best_seed_dir: Path) -> None:
        components = [
            as_float(row["objective_invivo"]),
            as_float(row["objective_invitro"]),
            as_float(row["objective_soft_coupling"]),
        ]
        objective = as_float(row["objective"])
        if not math.isfinite(objective) or not all(math.isfinite(value) for value in components):
            raise ValueError(f"Non-finite joint objective in summary row: {row['warmup_label']}")
        component_sum = sum(components)
        tolerance = 1e-8 * max(1.0, abs(objective))
        if abs(objective - component_sum) > tolerance:
            raise ValueError(
                f"Joint objective mismatch for {row['warmup_label']}: "
                f"objective={objective}, component_sum={component_sum}"
            )

        metrics = read_metric_table(best_seed_dir / "fit_summary.tsv")
        fit_objective = first_finite_metric(metrics, ("objective", "objective_total"))
        if not math.isfinite(fit_objective) or abs(objective - fit_objective) > tolerance:
            raise ValueError(
                f"Best-seed fit_summary objective mismatch for {row['warmup_label']}: "
                f"summary={objective}, fit_summary={fit_objective}"
            )
        self.validation_rows.append(
            {
                "check": "joint_objective_components",
                "record": row["warmup_label"],
                "status": "PASS",
                "expected": format(objective, ".17g"),
                "observed": format(component_sum, ".17g"),
                "details": str(best_seed_dir / "fit_summary.tsv"),
            }
        )

    def copy_separate_fits(self, joint_rows: Sequence[Dict[str, str]]) -> None:
        for row in joint_rows:
            for mode, root, seed_column, path_column in (
                ("invivo", self.invivo_root, "invivo_seed", "invivo_seed_dir"),
                ("invitro", self.invitro_root, "invitro_seed", "invitro_seed_dir"),
            ):
                recorded = ensure_within(Path(row[path_column]), root, path_column)
                expected = (root / normalized_seed(row[seed_column])).resolve()
                if recorded != expected:
                    raise ValueError(
                        "Recorded {} anchor does not match {} under the supplied root: {} != {}".format(
                            mode, row[seed_column], recorded, expected
                        )
                    )

        anchors = {
            "invivo": sorted({normalized_seed(row["invivo_seed"]) for row in joint_rows}, key=seed_number),
            "invitro": sorted({normalized_seed(row["invitro_seed"]) for row in joint_rows}, key=seed_number),
        }
        for mode, root in (("invivo", self.invivo_root), ("invitro", self.invitro_root)):
            destination = self.stage_dir / "separate_fits" / mode
            self.copy_named_files(
                root,
                destination / "run",
                SEPARATE_RUN_FILES,
                f"separate_{mode}_run",
                ("fit_array_manifest.tsv",),
            )
            self.copy_named_files(
                root / "extra_results",
                destination / "run/extra_results",
                EXTRA_RESULT_FILES,
                f"separate_{mode}_summary",
                ("seed_summary.tsv",),
            )

            best_seed, objective_metric, objective = select_separate_best(root, mode)
            selected = sorted(set(anchors[mode] + [best_seed]), key=seed_number)
            for seed in selected:
                source_seed = ensure_within(root / seed, root, f"{mode} seed directory")
                if not source_seed.is_dir():
                    raise NotADirectoryError(f"Separate-fit seed directory does not exist: {source_seed}")
                self.copy_seed_tree(
                    source_seed,
                    destination / "selected_seeds" / seed,
                    f"separate_{mode}_seed",
                    joint=False,
                )
            self.selected_rows.append(
                {
                    "record_type": f"separate_{mode}_global_best",
                    "warmup_label": "",
                    "invivo_seed": best_seed if mode == "invivo" else "",
                    "invitro_seed": best_seed if mode == "invitro" else "",
                    "selected_seed": best_seed,
                    "objective_metric": objective_metric,
                    "objective": format(objective, ".17g"),
                    "objective_invivo": "",
                    "objective_invitro": "",
                    "objective_soft_coupling": "",
                    "source_dir": str(root / best_seed),
                    "bundle_dir": str(Path("separate_fits") / mode / "selected_seeds" / best_seed),
                }
            )
            self.validation_rows.append(
                {
                    "check": "separate_fit_global_best",
                    "record": mode,
                    "status": "PASS",
                    "expected": best_seed,
                    "observed": best_seed,
                    "details": f"{objective_metric}={format(objective, '.17g')}",
                }
            )

    def write_generated_tables(self, joint_rows: Sequence[Dict[str, str]]) -> None:
        write_tsv(
            self.stage_dir / "selected_results.tsv",
            self.selected_rows,
            SELECTED_RESULT_COLUMNS,
        )
        write_tsv(
            self.stage_dir / "supporting_analysis_sources.tsv",
            self.supporting_rows,
            SUPPORTING_SOURCE_COLUMNS,
        )
        write_tsv(
            self.stage_dir / "VALIDATION.tsv",
            self.validation_rows,
            ("check", "record", "status", "expected", "observed", "details"),
        )

        anchor_rows = []
        for row in joint_rows:
            anchor_rows.append(
                {
                    "warmup_label": row["warmup_label"],
                    "invivo_seed": normalized_seed(row["invivo_seed"]),
                    "invivo_seed_dir": row["invivo_seed_dir"],
                    "invitro_seed": normalized_seed(row["invitro_seed"]),
                    "invitro_seed_dir": row["invitro_seed_dir"],
                    "joint_run_dir": row["joint_run_dir"],
                    "best_joint_seed": normalized_seed(row["best_joint_seed"]),
                }
            )
        write_tsv(
            self.stage_dir / "source_relationships.tsv",
            anchor_rows,
            (
                "warmup_label",
                "invivo_seed",
                "invivo_seed_dir",
                "invitro_seed",
                "invitro_seed_dir",
                "joint_run_dir",
                "best_joint_seed",
            ),
        )

    def write_readme(self, joint_rows: Sequence[Dict[str, str]]) -> None:
        generated = datetime.now(timezone.utc).isoformat()
        invivo_anchors = sorted(
            {normalized_seed(row["invivo_seed"]) for row in joint_rows}, key=seed_number
        )
        invitro_anchors = sorted(
            {normalized_seed(row["invitro_seed"]) for row in joint_rows}, key=seed_number
        )
        text = f"""# Fitting-output bundle

Generated: `{generated}`

This bundle is a consume-only export of recorded joint, in-vivo, and in-vitro
fitting outputs. No model was rerun and no fitted quantity was recalculated.

## Exact source roots

- Joint fit: `{self.joint_root}`
- In-vivo separate fit: `{self.invivo_root}`
- In-vitro separate fit: `{self.invitro_root}`
- Pooled embedding/curve-class report outputs: `{self.pooled_embedding_root}`
- Fixed-O2 eigen-attractor embedding/curve-class report outputs: `{self.fixo2_eigen_embedding_root}`
- Fixed-O2 dense-grid classification tables: `{self.fixed_o2_classification_tables}`
- Fixed-O2 dense-grid regression-classification tables: `{self.fixed_o2_regression_classification_tables}`
- Fixed-O2 analytical attractor tables: `{self.fixed_o2_attractor_tables}`

## Selection represented here

- Joint warm-up pairs: `{len(joint_rows)}`
- In-vivo anchors used by those pairs: `{', '.join(invivo_anchors)}`
- In-vitro anchors used by those pairs: `{', '.join(invitro_anchors)}`

`selected_results.tsv` records every pair-level best joint seed and the global
best separate-fit seeds. `source_relationships.tsv` preserves the exact source
paths connecting separate-fit anchors to joint runs. `VALIDATION.tsv` records
objective, selection, report-reference, and supporting-tree copy checks.
`supporting_analysis_sources.tsv` maps the five supporting-analysis directories
to their source roots and records exact copied file counts and byte totals.

For joint fits, `objective` is checked against
`objective_invivo + objective_invitro + objective_soft_coupling` and against the
selected seed's `fit_summary.tsv`. Separate in-vivo seeds are ranked by
`objective`; separate in-vitro seeds are ranked by `objective_total`.

## Bundle layout

- `joint_fit/overview/`: multi-warmup summaries and run manifests.
- `joint_fit/pairs/<warmup_label>/run/`: resolved/input configuration,
  effective arguments, provenance, array manifest, parameter inputs, and
  all-seed summary tables for that pair.
- `joint_fit/pairs/<warmup_label>/best_joint_seed/<seed>/`: existing numerical,
  configuration, RDS, log, and HTML outputs for the selected joint seed.
- `separate_fits/<mode>/run/`: separate-fit configuration, manifest, and
  all-seed summary tables.
- `separate_fits/<mode>/selected_seeds/<seed>/`: exact outputs for anchors used
  by joint fitting plus the global best separate-fit seed.
- `supporting_analysis/pooled_embedding_curve_class/`: complete output directory
  underlying `pooled_embedding_curve_class_report.html`.
- `supporting_analysis/fixo2_eigen_attractor_embedding_curve_class/`: complete
  output directory underlying
  `fixo2_eigen_attractor_embedding_curve_class_report.html`.
- `supporting_analysis/fixed_o2_analysis_tables/dense_grid_monotonicity_classification/`:
  direct fixed-O2 analysis-table inputs referenced by the pooled report.
- `supporting_analysis/fixed_o2_analysis_tables/dense_grid_monotonicity_regression_classification/`:
  direct fixed-O2 analysis-table inputs referenced by the FixO2 eigen report.
- `supporting_analysis/fixed_o2_analysis_tables/analytical_attractors/`:
  foundational fixed-O2 analytical attractor tables used by the eigen-attractor
  analysis.

Raster/PDF visualization files from selected fitting seeds and non-selected
joint seed directories are not duplicated. The two complete supporting report
directories do include their existing HTML, PDF, PNG, and tabular artifacts.
Existing fit/config RDS and numerical/provenance files for selected seeds are
included.

`FILE_MANIFEST.tsv` maps every copied file to its original absolute source and
records byte size and SHA-256. `SHA256SUMS.txt` verifies all bundle contents.
"""
        (self.stage_dir / "README.md").write_text(text, encoding="utf-8")

    def write_file_manifest(self) -> None:
        generated_paths = (
            "README.md",
            "selected_results.tsv",
            "source_relationships.tsv",
            "supporting_analysis_sources.tsv",
            "VALIDATION.tsv",
        )
        rows = list(self.copied)
        for relative in generated_paths:
            path = self.stage_dir / relative
            rows.append(
                {
                    "bundle_path": relative,
                    "source_path": "GENERATED_BY_BUNDLE_BUILDER",
                    "category": "bundle_metadata",
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
        rows.sort(key=lambda row: str(row["bundle_path"]))
        write_tsv(
            self.stage_dir / "FILE_MANIFEST.tsv",
            rows,
            ("bundle_path", "source_path", "category", "bytes", "sha256"),
        )

    def write_checksums(self) -> None:
        checksum_path = self.stage_dir / "SHA256SUMS.txt"
        paths = sorted(
            path for path in self.stage_dir.rglob("*") if path.is_file() and path != checksum_path
        )
        lines = [f"{sha256_file(path)}  {path.relative_to(self.stage_dir)}" for path in paths]
        checksum_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def finalize_directory(self) -> None:
        self.stage_dir.rename(self.output_dir)

    def write_archive(self) -> None:
        archive_tmp = Path(f"{self.archive_path}.tmp-{os.getpid()}")
        try:
            with tarfile.open(str(archive_tmp), "w:gz") as handle:
                handle.add(str(self.output_dir), arcname=self.output_dir.name)
            archive_tmp.rename(self.archive_path)
            digest = sha256_file(self.archive_path)
            self.archive_checksum_path.write_text(
                f"{digest}  {self.archive_path.name}\n", encoding="utf-8"
            )
        finally:
            if archive_tmp.exists():
                archive_tmp.unlink()


def read_metric_table(path: Path) -> Dict[str, str]:
    rows = read_tsv(path)
    if not rows:
        return {}
    if "metric" in rows[0] and "value" in rows[0]:
        return {
            row["metric"].strip(): row["value"].strip()
            for row in rows
            if row.get("metric", "").strip()
        }
    if len(rows) == 1:
        return {key: value for key, value in rows[0].items()}
    raise ValueError(f"Unsupported fit-summary layout: {path}")


def first_finite_metric(metrics: Dict[str, str], candidates: Sequence[str]) -> float:
    for candidate in candidates:
        value = as_float(metrics.get(candidate))
        if math.isfinite(value):
            return value
    return math.nan


def select_separate_best(root: Path, mode: str) -> Tuple[str, str, float]:
    summary_path = root / "extra_results/seed_summary.tsv"
    rows = read_tsv(summary_path)
    objective_candidates = ("objective", "optimizer_local_objective") if mode == "invivo" else (
        "objective_total",
        "objective",
        "optimizer_local_objective",
    )
    require_columns(rows, ("seed",), summary_path)
    ranked: List[Tuple[float, int, str, str]] = []
    for row in rows:
        seed = normalized_seed(row["seed"])
        for metric in objective_candidates:
            value = as_float(row.get(metric))
            if math.isfinite(value):
                ranked.append((value, seed_number(seed), seed, metric))
                break
    if not ranked:
        raise ValueError(f"No finite objective found in {summary_path}")
    objective, _, seed, metric = min(ranked)

    metrics = read_metric_table(root / seed / "fit_summary.tsv")
    seed_objective = first_finite_metric(metrics, (metric,) + objective_candidates)
    tolerance = 1e-8 * max(1.0, abs(objective))
    if not math.isfinite(seed_objective) or abs(seed_objective - objective) > tolerance:
        raise ValueError(
            f"Separate-fit best objective mismatch for {mode} {seed}: "
            f"seed_summary={objective}, fit_summary={seed_objective}"
        )
    return seed, metric, objective


def seed_number(seed: str) -> int:
    return int(normalized_seed(seed)[4:])


def main() -> None:
    args = parse_args()
    builder = BundleBuilder(
        joint_root=args.joint_root,
        invivo_root=args.invivo_root,
        invitro_root=args.invitro_root,
        pooled_embedding_root=args.pooled_embedding_root,
        fixo2_eigen_embedding_root=args.fixo2_eigen_embedding_root,
        fixed_o2_classification_tables=args.fixed_o2_classification_tables,
        fixed_o2_regression_classification_tables=(
            args.fixed_o2_regression_classification_tables
        ),
        fixed_o2_attractor_tables=args.fixed_o2_attractor_tables,
        output_dir=args.output_dir,
        force=args.force,
        create_archive=not args.no_archive,
    )
    builder.build()
    print(f"output_dir={builder.output_dir}")
    print(f"copied_files={len(builder.copied)}")
    print(f"selected_records={len(builder.selected_rows)}")
    print(f"supporting_sources={len(builder.supporting_rows)}")
    print(f"validation_checks={len(builder.validation_rows)}")
    if builder.create_archive:
        print(f"archive={builder.archive_path}")
        print(f"archive_sha256={sha256_file(builder.archive_path)}")


if __name__ == "__main__":
    main()
