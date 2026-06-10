#!/usr/bin/env python3
"""Rank separate-fit seeds by objective and collect the top report HTML files."""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]

DEFAULT_INVIVO_DIR = PROJECT_ROOT / "oxygen/results/fit_invivo_O2_buffering_200seed"
DEFAULT_INVITRO_DIR = PROJECT_ROOT / "oxygen/results/fit_invitro_O2_buffering_200seed"
DEFAULT_OUT_DIR = PROJECT_ROOT / "oxygen/results/warm_start_from_separate_200seed/best_result"

CSV_METADATA_COLUMNS = [
    "rank",
    "seed",
    "objective",
    "objective_metric",
    "optimizer_local_objective",
    "optimizer_deoptim_objective",
    "fit_mode",
]

LOG10_PARAMS = {
    "lam_max",
    "p_mis_base",
    "p_misseg",
    "k_o_mis",
    "buffer_beta",
    "buffer_n_exp",
    "p_wgd",
    "o2_S0",
    "kappa_O",
    "eta_o2",
    "rho_2N",
    "alpha_o2",
    "mu_hp",
    "O2_crit",
    "k_clear",
    "sigma_burden",
    "tau_O2",
    "c_vol_2N_eff_mm3",
    "sigma_growth",
    "sigma_kary",
    "init_sd_2N",
    "init_sd_4N",
    "o2_Nref",
}

LOGIT_PARAMS: set[str] = set()

IDENTITY_PARAMS = {
    "o2_min",
    "buffer_smax",
    "gamma_growth",
    "gamma_mu",
    "n_O",
    "beta_size",
    "init_mean_2N",
    "init_mean_4N",
    "alpha",
    "gamma",
    "ratio_4N_2N",
    "o2_S0_upper_bound",
}


def transform_label(parameter: str) -> str:
    if parameter in LOG10_PARAMS:
        return "log10"
    if parameter in LOGIT_PARAMS:
        return "logit"
    return "identity"


def transformed_column_name(parameter: str) -> str:
    label = transform_label(parameter)
    return f"{label}_{parameter}"


def transform_value(parameter: str, value: float) -> float:
    label = transform_label(parameter)
    if not math.isfinite(value):
        return math.nan
    if label == "identity":
        return value
    if label == "log10":
        return math.log10(value) if value > 0 else math.nan
    if label == "logit":
        return math.log(value / (1.0 - value)) if 0 < value < 1 else math.nan
    raise ValueError(f"Unsupported transform label: {label}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--invivo-dir", type=Path, default=DEFAULT_INVIVO_DIR)
    parser.add_argument("--invitro-dir", type=Path, default=DEFAULT_INVITRO_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--top-n", type=int, default=10)
    return parser.parse_args()


def read_metric_table(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {
            str(row["metric"]).strip(): str(row["value"]).strip()
            for row in reader
            if str(row.get("metric", "")).strip()
        }


def read_best_params(path: Path) -> dict[str, float]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {
            str(row["parameter"]).strip(): as_float(row.get("value"))
            for row in reader
            if str(row.get("parameter", "")).strip()
        }


def as_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def report_path(seed_dir: Path, seed: int) -> Path | None:
    preferred = seed_dir / "report" / f"fit_report_seed{seed}.html"
    if preferred.exists():
        return preferred
    report_dir = seed_dir / "report"
    if not report_dir.is_dir():
        return None
    candidates = sorted(report_dir.glob("*.html"))
    return candidates[0] if candidates else None


def collect_results(parent: Path, mode: str) -> list[dict[str, object]]:
    objective_key = "objective" if mode == "invivo" else "objective_total"
    fallback_keys = [objective_key, "optimizer_local_objective", "optimizer_deoptim_objective"]
    pattern = re.compile(r"^seed(\d+)$")
    rows: list[dict[str, object]] = []

    for seed_dir in sorted(parent.iterdir(), key=lambda p: p.name):
        if not seed_dir.is_dir():
            continue
        match = pattern.match(seed_dir.name)
        if not match:
            continue
        seed = int(match.group(1))
        summary_path = seed_dir / "fit_summary.tsv"
        params_path = seed_dir / "best_params.tsv"
        if not summary_path.exists() or not params_path.exists():
            continue

        metrics = read_metric_table(summary_path)
        params = read_best_params(params_path)
        objective = math.nan
        objective_metric = ""
        for key in fallback_keys:
            value = as_float(metrics.get(key))
            if math.isfinite(value):
                objective = value
                objective_metric = key
                break

        report = report_path(seed_dir, seed)
        rows.append(
            {
                "rank": "",
                "seed": seed,
                "objective": objective,
                "objective_metric": objective_metric,
                "optimizer_local_objective": as_float(metrics.get("optimizer_local_objective")),
                "optimizer_deoptim_objective": as_float(metrics.get("optimizer_deoptim_objective")),
                "fit_mode": metrics.get("fit_mode", ""),
                "seed_dir": str(seed_dir),
                "report_path": str(report) if report else "",
                "copied_report_path": "",
                "_params": params,
            }
        )

    rows.sort(
        key=lambda row: (
            float(row["objective"]) if math.isfinite(float(row["objective"])) else math.inf,
            int(row["seed"]),
        )
    )
    for idx, row in enumerate(rows, start=1):
        row["rank"] = idx
    return rows


def copy_top_reports(rows: list[dict[str, object]], target_dir: Path, top_n: int) -> None:
    target_dir.mkdir(parents=True, exist_ok=True)
    for row in rows[:top_n]:
        src_str = str(row["report_path"])
        if not src_str:
            continue
        src = Path(src_str)
        if not src.exists():
            continue
        objective = float(row["objective"])
        objective_label = f"{objective:.6g}" if math.isfinite(objective) else "NA"
        dest = target_dir / (
            f"rank{int(row['rank']):02d}_"
            f"seed{int(row['seed'])}_"
            f"objective{objective_label}_"
            f"{src.name}"
        )
        shutil.copy2(src, dest)
        row["copied_report_path"] = str(dest)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("")
        return
    param_names = []
    seen_params = set()
    for row in rows:
        for parameter in row.get("_params", {}):
            if parameter not in seen_params:
                seen_params.add(parameter)
                param_names.append(parameter)

    fieldnames = (
        CSV_METADATA_COLUMNS
        + param_names
        + [transformed_column_name(parameter) for parameter in param_names]
    )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = {key: row.get(key, "") for key in CSV_METADATA_COLUMNS}
            params = row.get("_params", {})
            for parameter in param_names:
                value = params.get(parameter, math.nan)
                out[parameter] = value
                out[transformed_column_name(parameter)] = transform_value(parameter, value)
            writer.writerow(out)


def main() -> None:
    args = parse_args()
    if args.top_n <= 0:
        raise SystemExit("--top-n must be > 0")

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    invivo_rows = collect_results(args.invivo_dir.resolve(), "invivo")
    invitro_rows = collect_results(args.invitro_dir.resolve(), "invitro")

    copy_top_reports(invivo_rows, out_dir / "invivo", args.top_n)
    copy_top_reports(invitro_rows, out_dir / "invitro", args.top_n)
    write_csv(out_dir / "invivo_objective_ranking.csv", invivo_rows)
    write_csv(out_dir / "invitro_objective_ranking.csv", invitro_rows)

    print(f"out_dir={out_dir}")
    print(f"invivo_results={len(invivo_rows)} top1_seed={invivo_rows[0]['seed']} top1_objective={invivo_rows[0]['objective']}")
    print(f"invitro_results={len(invitro_rows)} top1_seed={invitro_rows[0]['seed']} top1_objective={invitro_rows[0]['objective']}")
    print(f"invivo_reports_copied={sum(1 for _ in (out_dir / 'invivo').glob('*.html'))}")
    print(f"invitro_reports_copied={sum(1 for _ in (out_dir / 'invitro').glob('*.html'))}")


if __name__ == "__main__":
    main()
