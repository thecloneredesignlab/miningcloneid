#!/usr/bin/env python3
"""
Plot fitted 2N-vs-4N parameter differences from an invitro_fitting output folder.

Usage:
    python code/plot_ploidy_parameter_comparison.py code/invitro_fitting_outputs/20260513T164159
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_PARAMETERS = [
    "r",
    "K",
    "k_tr",
    "k_kill",
    "k_clear",
    "k_cyto",
    "beta_dose",
    "mu_base_death",
    "mu_confluence_death",
]


DISPLAY_NAMES = {
    "r": "r",
    "K": "K",
    "k_tr": "k_tr",
    "k_kill": "k_kill",
    "k_clear": "k_clear",
    "k_cyto": "k_cyto",
    "beta_dose": "beta_dose",
    "mu_base_death": "mu_base_death",
    "mu_confluence_death": "mu_confluence_death",
}


def choose_best_row(summary_df: pd.DataFrame) -> pd.Series:
    if "posterior_objective" in summary_df.columns:
        values = pd.to_numeric(summary_df["posterior_objective"], errors="coerce")
        if values.notna().any():
            return summary_df.loc[values.idxmin()]
    if len(summary_df) == 0:
        raise ValueError("joint_fit_summary.tsv has no rows")
    return summary_df.iloc[0]


def build_comparison_table(best_row: pd.Series, parameters: Iterable[str]) -> pd.DataFrame:
    rows = []
    for parameter in parameters:
        col_2n = f"2N_{parameter}"
        col_4n = f"4N_{parameter}"
        if col_2n not in best_row.index or col_4n not in best_row.index:
            continue
        value_2n = pd.to_numeric(pd.Series([best_row[col_2n]]), errors="coerce").iloc[0]
        value_4n = pd.to_numeric(pd.Series([best_row[col_4n]]), errors="coerce").iloc[0]
        if not np.isfinite(value_2n) or not np.isfinite(value_4n):
            continue
        if value_2n <= 0 or value_4n <= 0:
            continue
        ratio = float(value_4n / value_2n)
        rows.append(
            {
                "parameter": parameter,
                "display_name": DISPLAY_NAMES.get(parameter, parameter),
                "value_2N": float(value_2n),
                "value_4N": float(value_4n),
                "ratio_4N_over_2N": ratio,
                "log2_ratio_4N_over_2N": float(np.log2(ratio)),
            }
        )
    if not rows:
        raise ValueError("No positive paired 2N/4N parameter columns were found")
    return pd.DataFrame(rows)


def plot_fold_change(comparison_df: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(comparison_df))
    values = comparison_df["log2_ratio_4N_over_2N"].to_numpy(dtype=float)
    colors = np.where(values >= 0.0, "#4C78A8", "#F58518")

    ax.bar(x, values, color=colors, edgecolor="black", linewidth=0.6)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(comparison_df["display_name"], rotation=35, ha="right")
    ax.set_ylabel("log2(4N / 2N)")
    ax.set_title("Fitted Parameter Fold-Change: 4N vs 2N")
    ax.grid(axis="y", linestyle="--", alpha=0.35)

    for idx, value in enumerate(values):
        label_y = value + (0.06 if value >= 0 else -0.06)
        va = "bottom" if value >= 0 else "top"
        ax.text(idx, label_y, f"{2 ** value:.2g}x", ha="center", va=va, fontsize=8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_paired_values(comparison_df: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5.5))
    x_2n = np.zeros(len(comparison_df))
    x_4n = np.ones(len(comparison_df))

    for idx, row in comparison_df.iterrows():
        y_values = [row["value_2N"], row["value_4N"]]
        ax.plot([0, 1], y_values, color="#888888", linewidth=1.2, alpha=0.8)
        ax.scatter([0], [row["value_2N"]], color="#4C78A8", s=45, zorder=3)
        ax.scatter([1], [row["value_4N"]], color="#F58518", s=45, zorder=3)
        geometric_midpoint = np.sqrt(row["value_2N"] * row["value_4N"])
        ax.text(1.05, geometric_midpoint, row["display_name"], va="center", fontsize=9)

    ax.set_yscale("log")
    ax.set_xlim(-0.25, 1.85)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["2N", "4N"])
    ax.set_ylabel("Fitted parameter value (log scale)")
    ax.set_title("Paired Fitted Parameters by Ploidy")
    ax.grid(axis="y", which="both", linestyle="--", alpha=0.35)

    # Keep explicit arrays in the plot code so future readers can see that all
    # parameters share the same two ploidy positions.
    _ = x_2n, x_4n

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot 2N-vs-4N fitted parameter comparisons from joint_fit_summary.tsv."
    )
    parser.add_argument(
        "output_folder",
        type=Path,
        help="Folder containing joint_fit_summary.tsv.",
    )
    parser.add_argument(
        "--parameters",
        nargs="+",
        default=DEFAULT_PARAMETERS,
        help="Parameter names to plot, without the 2N_/4N_ prefixes.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_folder = args.output_folder
    summary_path = output_folder / "joint_fit_summary.tsv"
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing summary file: {summary_path}")

    summary_df = pd.read_csv(summary_path, sep="\t")
    best_row = choose_best_row(summary_df)
    comparison_df = build_comparison_table(best_row, args.parameters)

    table_path = output_folder / "ploidy_parameter_comparison.tsv"
    fold_change_path = output_folder / "ploidy_parameter_log2_fold_change.png"
    paired_path = output_folder / "ploidy_parameter_paired_values.png"

    comparison_df.to_csv(table_path, sep="\t", index=False)
    plot_fold_change(comparison_df, fold_change_path)
    plot_paired_values(comparison_df, paired_path)

    print(f"Wrote {table_path}")
    print(f"Wrote {fold_change_path}")
    print(f"Wrote {paired_path}")


if __name__ == "__main__":
    main()
