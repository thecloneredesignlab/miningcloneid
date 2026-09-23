#!/usr/bin/env python3
"""Reproducible SALib FAST design, analysis, and Figure 4 companion plots."""

import argparse
import csv
import gzip
import hashlib
from importlib.metadata import version
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
from SALib.analyze import fast
from SALib.sample import fast_sampler

salib_version = version("SALib")


ACTIVE = (
    "lam_max", "p_mis_base", "p_misseg", "k_o_mis", "buffer_smax",
    "buffer_beta", "buffer_n_exp", "p_wgd", "alpha_o2", "gamma_growth",
    "mu_hp", "gamma_mu", "O2_crit", "n_O",
)
STRUCTURAL = ("o2_S0", "kappa_O", "eta_o2", "k_clear")
GROUP = {
    "lam_max": "growth", "alpha_o2": "growth", "gamma_growth": "growth",
    "p_mis_base": "missegregation", "p_misseg": "missegregation",
    "k_o_mis": "missegregation", "p_wgd": "genome_doubling",
    "buffer_smax": "buffering", "buffer_beta": "buffering", "buffer_n_exp": "buffering",
    "mu_hp": "death", "gamma_mu": "death",
    "O2_crit": "shared_oxygen_stress", "n_O": "shared_oxygen_stress",
}
OUTPUTS = ("dominant_mean_ploidy", "dominant_growth_rate")


def read_table(path, delimiter="\t"):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_table(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def get_ranges(fit_root):
    path = Path(fit_root) / "parameter_table.csv"
    rows = read_table(path, ",")
    by_prototype = {row["param_prototype"]: row for row in rows if row["estimate"].upper() == "TRUE"}
    absent = set(ACTIVE) - set(by_prototype)
    if absent:
        raise ValueError("Missing estimated operator parameters: " + ", ".join(sorted(absent)))
    ranges = []
    for name in ACTIVE:
        row = by_prototype[name]
        encoded = row["param_name"]
        transform = "log10" if encoded == "log10_" + name else "identity"
        if transform == "identity" and encoded != name:
            raise ValueError("Unexpected parameter transform: " + encoded)
        lo, hi = float(row["lower_bound"]), float(row["upper_bound"])
        if not math.isfinite(lo) or not math.isfinite(hi) or lo >= hi:
            raise ValueError("Invalid bounds for " + name)
        ranges.append({
            "parameter": name, "transform": transform,
            "sampling_distribution": "log_uniform" if transform == "log10" else "uniform",
            "encoded_lower": lo, "encoded_upper": hi,
            "natural_lower": 10 ** lo if transform == "log10" else lo,
            "natural_upper": 10 ** hi if transform == "log10" else hi,
            "group": GROUP[name], "fit_table_name": encoded,
        })
    return path, ranges


def oxygen_grid(figure4_dir, full_grid, subset):
    path = Path(figure4_dir) / "fixed_o2_dominant_ploidy_201grid.tsv"
    with open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        values = sorted({float(row["O2_pct"]) for row in reader})
    if len(values) != 201 or not np.allclose(values, np.linspace(0, 5, 201), atol=1e-12):
        raise ValueError("Figure 4 oxygen grid differs from 0:0.025:5 (201 points)")
    if full_grid:
        return values, path
    wanted = [float(item) for item in subset.split(",")]
    if not wanted or any(min(abs(x - y) for y in values) > 1e-12 for x in wanted):
        raise ValueError("Pilot oxygen values must be on the Figure 4 grid")
    return sorted(set(wanted)), path


def prepare(args):
    if args.n <= 4 * args.m * args.m:
        raise ValueError("FAST requires N > 4*M^2; choose N > %d" % (4 * args.m * args.m))
    fit_table, ranges = get_ranges(args.fit_root)
    oxygen, grid_path = oxygen_grid(args.figure4_dir, args.full_grid, args.oxygen)
    output_root = Path(args.out_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    fields = list(ranges[0])
    write_table(output_root / "parameter_ranges.tsv", fields, ranges)
    correlation_path = Path(args.figure4_dir) / "continuous_ploidy_spearman_by_o2.tsv"
    correlation_copy = output_root / "figure4b_spearman_source.tsv"
    if not correlation_copy.exists():
        shutil.copyfile(correlation_path, correlation_copy)
    elif sha256(correlation_copy) != sha256(correlation_path):
        raise ValueError("Figure 4B correlation source changed between designs")
    group_path = Path(args.figure4_dir) / "parameter_function_groups.tsv"
    group_copy = output_root / "figure4_parameter_groups_source.tsv"
    if not group_copy.exists():
        shutil.copyfile(group_path, group_copy)
    elif sha256(group_copy) != sha256(group_path):
        raise ValueError("Figure 4 parameter order/group source changed between designs")
    problem = {
        "num_vars": len(ACTIVE), "names": list(ACTIVE),
        "bounds": [[r["encoded_lower"], r["encoded_upper"]] for r in ranges],
    }
    for replicate in range(1, args.replicates + 1):
        run_dir = output_root / "runs" / ("N%d_R%d" % (args.n, replicate))
        run_dir.mkdir(parents=True, exist_ok=True)
        seed = args.seed + replicate - 1
        encoded = fast_sampler.sample(problem, args.n, M=args.m, seed=seed)
        natural = encoded.copy()
        for j, spec in enumerate(ranges):
            if spec["transform"] == "log10":
                natural[:, j] = np.power(10.0, encoded[:, j])
        sample_rows = (
            dict(sample_id=i + 1, **{name: "%.17g" % natural[i, j]
                                     for j, name in enumerate(ACTIVE)})
            for i in range(len(natural))
        )
        sample_path = run_dir / "samples.tsv.gz"
        write_table(sample_path, ["sample_id"] + list(ACTIVE), sample_rows)
        metadata = {
            "salib_version": salib_version, "method": "eFAST", "M": args.m,
            "N": args.n, "replicate": replicate, "seed": seed,
            "n_parameters": len(ACTIVE), "n_samples": len(natural),
            "fit_root": str(Path(args.fit_root).resolve()),
            "figure4_dir": str(Path(args.figure4_dir).resolve()),
            "fit_parameter_table_sha256": sha256(fit_table),
            "figure4_grid_sha256": sha256(grid_path),
            "figure4b_spearman_sha256": sha256(correlation_copy),
            "figure4_parameter_groups_sha256": sha256(group_copy),
            "oxygen_pct": oxygen,
            "sampling": "independent uniform over fitted transformed bounds; log-uniform where log10",
            "samples_sha256": sha256(sample_path),
        }
        with open(run_dir / "metadata.json", "w") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
        print("Prepared %s: %d parameter vectors x %d oxygen points" %
              (run_dir, len(natural), len(oxygen)), flush=True)


def summarize(args):
    root = Path(args.out_dir)
    range_rows = read_table(root / "parameter_ranges.tsv")
    bounds = [[float(r["encoded_lower"]), float(r["encoded_upper"])] for r in range_rows]
    problem = {"num_vars": len(ACTIVE), "names": list(ACTIVE), "bounds": bounds}
    all_rows = []
    gap_qc_rows = []
    for metadata_path in sorted((root / "runs").glob("*/metadata.json")):
        run_dir = metadata_path.parent
        metadata = json.loads(metadata_path.read_text())
        outputs_path = run_dir / "outputs.tsv.gz"
        if not outputs_path.exists():
            raise FileNotFoundError(outputs_path)
        n_samples = metadata["n_samples"]
        grid = metadata["oxygen_pct"]
        data = {o2: {name: np.full(n_samples, np.nan) for name in OUTPUTS} for o2 in grid}
        seen = {o2: np.zeros(n_samples, dtype=bool) for o2 in grid}
        gap_stats = {band: {"n": 0, "lt_1e-6": 0, "lt_1e-4": 0,
                            "lt_1e-3": 0, "minimum": float("inf")}
                     for band in ("low_0_to_1pct", "middle_1_to_3pct", "high_3_to_5pct")}
        with gzip.open(outputs_path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                sample_id = int(row["sample_id"])
                o2 = float(row["O2_pct"])
                if o2 not in data or sample_id < 1 or sample_id > n_samples:
                    raise ValueError("Unexpected sample/O2 in " + str(outputs_path))
                index = sample_id - 1
                if seen[o2][index]:
                    raise ValueError("Duplicate sample/O2 in " + str(outputs_path))
                seen[o2][index] = True
                if row["status"] != "ok" or row["eigenvector_nonnegative"] != "TRUE":
                    raise ValueError("Nonvalid model evaluation at sample %d, O2 %s: %s" %
                                     (sample_id, o2, row["status"]))
                gap = float(row["spectral_gap"])
                if not math.isfinite(gap) or gap < 0:
                    raise ValueError("Invalid spectral gap at sample %d, O2 %s" % (sample_id, o2))
                band = "low_0_to_1pct" if o2 <= 1 else (
                    "high_3_to_5pct" if o2 >= 3 else "middle_1_to_3pct")
                stats = gap_stats[band]
                stats["n"] += 1
                stats["lt_1e-6"] += gap < 1e-6
                stats["lt_1e-4"] += gap < 1e-4
                stats["lt_1e-3"] += gap < 1e-3
                stats["minimum"] = min(stats["minimum"], gap)
                for name in OUTPUTS:
                    data[o2][name][index] = float(row[name])
        if any(not np.all(seen[o2]) for o2 in grid):
            raise ValueError("Incomplete FAST trajectories in " + str(outputs_path))
        for band, stats in gap_stats.items():
            if stats["n"]:
                gap_qc_rows.append({
                    "N": metadata["N"], "replicate": metadata["replicate"],
                    "oxygen_band": band, "n_evaluations": stats["n"],
                    "fraction_gap_lt_1e-6": "%.10g" % (stats["lt_1e-6"] / stats["n"]),
                    "fraction_gap_lt_1e-4": "%.10g" % (stats["lt_1e-4"] / stats["n"]),
                    "fraction_gap_lt_1e-3": "%.10g" % (stats["lt_1e-3"] / stats["n"]),
                    "minimum_gap": "%.10g" % stats["minimum"],
                })
        for o2 in grid:
            for name in OUTPUTS:
                y = data[o2][name]
                if not np.all(np.isfinite(y)) or np.var(y) <= 1e-16:
                    raise ValueError("Nonfinite or zero-variance output: %s at %g" % (name, o2))
                indices = fast.analyze(problem, y, M=metadata["M"], print_to_console=False,
                                       seed=metadata["seed"])
                for j, parameter in enumerate(ACTIVE):
                    all_rows.append({
                        "N": metadata["N"], "replicate": metadata["replicate"],
                        "seed": metadata["seed"], "O2_pct": "%.3f" % o2,
                        "output": name, "parameter": parameter, "group": GROUP[parameter],
                        "S1": "%.10g" % indices["S1"][j],
                        "ST": "%.10g" % indices["ST"][j],
                        "output_variance": "%.10g" % np.var(y),
                    })
        print("Analyzed " + str(run_dir), flush=True)
    fields = ["N", "replicate", "seed", "O2_pct", "output", "parameter", "group",
              "S1", "ST", "output_variance"]
    write_table(root / "indices.tsv", fields, all_rows)
    write_table(root / "spectral_gap_qc.tsv", list(gap_qc_rows[0]), gap_qc_rows)
    summarize_convergence(root, all_rows)
    summarize_band_replicates(root, all_rows)
    try:
        code_commit = subprocess.check_output(
            ["git", "-C", str(Path(__file__).resolve().parents[5]), "rev-parse", "HEAD"],
            text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        code_commit = "unavailable"
    manifest = [
        {"field": "summary_code_commit", "value": code_commit},
        {"field": "salib_version", "value": salib_version},
        {"field": "indices_sha256", "value": sha256(root / "indices.tsv")},
        {"field": "convergence_sha256", "value": sha256(root / "convergence.tsv")},
        {"field": "parameter_band_summary_sha256", "value": sha256(root / "parameter_band_summary.tsv")},
        {"field": "mechanism_band_summary_sha256", "value": sha256(root / "mechanism_band_summary.tsv")},
        {"field": "mechanism_band_replicates_sha256",
         "value": sha256(root / "mechanism_band_replicates.tsv")},
        {"field": "spectral_gap_qc_sha256", "value": sha256(root / "spectral_gap_qc.tsv")},
        {"field": "figure4_parameter_groups_sha256",
         "value": sha256(root / "figure4_parameter_groups_source.tsv")},
    ]
    write_table(root / "analysis_manifest.tsv", ["field", "value"], manifest)


def summarize_band_replicates(root, rows):
    highest_n = max(int(row["N"]) for row in rows)
    bands = {
        "low_0_to_1pct": lambda x: 0 <= x <= 1,
        "high_3_to_5pct": lambda x: 3 <= x <= 5,
    }
    selected = [row for row in rows if int(row["N"]) == highest_n]
    repeats = sorted(set(int(row["replicate"]) for row in selected))
    summary = []
    for output in OUTPUTS:
        for replicate in repeats:
            for band, predicate in bands.items():
                for group in sorted(set(GROUP.values())):
                    subset = [row for row in selected if row["output"] == output and
                              int(row["replicate"]) == replicate and row["group"] == group and
                              predicate(float(row["O2_pct"]))]
                    summary.append({
                        "N": highest_n, "output": output, "replicate": replicate,
                        "oxygen_band": band, "group": group,
                        "n_parameters": len(set(row["parameter"] for row in subset)),
                        "n_oxygen": len(set(row["O2_pct"] for row in subset)),
                        "S1_mean_per_parameter": "%.10g" % np.mean(
                            [float(row["S1"]) for row in subset]),
                        "ST_mean_per_parameter": "%.10g" % np.mean(
                            [float(row["ST"]) for row in subset]),
                    })
    write_table(root / "mechanism_band_replicates.tsv", list(summary[0]), summary)


def summarize_convergence(root, rows):
    by_key = {}
    for row in rows:
        key = (int(row["N"]), row["output"], row["parameter"], float(row["O2_pct"]))
        by_key.setdefault(key, []).append(row)
    conv = []
    for key, items in sorted(by_key.items()):
        n, output, parameter, o2 = key
        s1 = np.array([float(r["S1"]) for r in items])
        st = np.array([float(r["ST"]) for r in items])
        conv.append({"N": n, "output": output, "parameter": parameter,
                     "O2_pct": "%.3f" % o2, "replicates": len(items),
                     "S1_mean": "%.10g" % s1.mean(), "S1_range": "%.10g" % (s1.max() - s1.min()),
                     "ST_mean": "%.10g" % st.mean(), "ST_range": "%.10g" % (st.max() - st.min())})
    by_resolution = {(r["output"], r["parameter"], r["O2_pct"]): r
                     for r in conv if r["N"] == max(int(x["N"]) for x in conv)}
    for row in conv:
        high = by_resolution.get((row["output"], row["parameter"], row["O2_pct"]))
        row["S1_delta_vs_max_N"] = ("%.10g" % abs(float(row["S1_mean"]) - float(high["S1_mean"]))) if high else ""
        row["ST_delta_vs_max_N"] = ("%.10g" % abs(float(row["ST_mean"]) - float(high["ST_mean"]))) if high else ""
    write_table(root / "convergence.tsv", list(conv[0]), conv)
    summarize_mechanism_bands(root, conv)


def summarize_mechanism_bands(root, rows):
    highest_n = max(int(row["N"]) for row in rows)
    top = [row for row in rows if int(row["N"]) == highest_n]
    low_n = [row for row in rows if int(row["N"]) < highest_n]
    low_lookup = {(row["output"], row["parameter"], row["O2_pct"]): row for row in low_n}
    bands = {
        "low_0_to_1pct": lambda x: 0 <= x <= 1,
        "high_3_to_5pct": lambda x: 3 <= x <= 5,
    }
    param_rows = []
    for output in OUTPUTS:
        for band, predicate in bands.items():
            for parameter in ACTIVE:
                relevant = [row for row in top if row["output"] == output and
                            row["parameter"] == parameter and predicate(float(row["O2_pct"]))]
                if not relevant:
                    raise ValueError("Empty oxygen band for " + parameter)
                low_deltas = [abs(float(row["ST_mean"]) -
                                  float(low_lookup[(output, parameter, row["O2_pct"])]["ST_mean"]))
                              for row in relevant if (output, parameter, row["O2_pct"]) in low_lookup]
                param_rows.append({
                    "output": output, "oxygen_band": band, "parameter": parameter,
                    "group": GROUP[parameter], "N": highest_n, "n_oxygen": len(relevant),
                    "S1_mean": "%.10g" % np.mean([float(row["S1_mean"]) for row in relevant]),
                    "ST_mean": "%.10g" % np.mean([float(row["ST_mean"]) for row in relevant]),
                    "ST_replicate_range_p90": "%.10g" % np.percentile(
                        [float(row["ST_range"]) for row in relevant], 90),
                    "ST_resolution_delta_p90": "%.10g" % np.percentile(low_deltas, 90)
                    if low_deltas else "",
                })
    write_table(root / "parameter_band_summary.tsv", list(param_rows[0]), param_rows)
    group_rows = []
    for output in OUTPUTS:
        for band in bands:
            for group in sorted(set(GROUP.values())):
                subset = [row for row in param_rows if row["output"] == output and
                          row["oxygen_band"] == band and row["group"] == group]
                group_rows.append({
                    "output": output, "oxygen_band": band, "group": group,
                    "n_parameters": len(subset),
                    "S1_mean_per_parameter": "%.10g" % np.mean(
                        [float(row["S1_mean"]) for row in subset]),
                    "ST_mean_per_parameter": "%.10g" % np.mean(
                        [float(row["ST_mean"]) for row in subset]),
                    "ST_max_parameter": max(subset, key=lambda row: float(row["ST_mean"]))["parameter"],
                    "ST_max_parameter_mean": "%.10g" % max(float(row["ST_mean"]) for row in subset),
                })
    write_table(root / "mechanism_band_summary.tsv", list(group_rows[0]), group_rows)


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    root = Path(args.out_dir)
    rows = read_table(root / "convergence.tsv")
    max_n = max(int(row["N"]) for row in rows)
    rows = [row for row in rows if int(row["N"]) == max_n]
    o2_values = sorted(set(float(row["O2_pct"]) for row in rows))
    if len(o2_values) != 201:
        raise ValueError("Publication heatmaps require the full Figure 4 oxygen grid")
    fig_dir = root / "figures"
    fig_dir.mkdir(exist_ok=True)
    group_rows = sorted(read_table(root / "figure4_parameter_groups_source.tsv"),
                        key=lambda row: int(row["parameter_order"]))
    parameters = [row["parameter"] for row in group_rows]
    if len(parameters) != len(ACTIVE) + len(STRUCTURAL) or set(parameters) != set(ACTIVE + STRUCTURAL):
        raise ValueError("Figure 4 parameter order does not match fixed-O2 heatmap rows")
    group_boundaries = [i for i in range(1, len(parameters))
                        if group_rows[i]["parameter_group"] != group_rows[i - 1]["parameter_group"]]
    for output in OUTPUTS:
        for index in ("S1", "ST"):
            lookup = {(row["parameter"], float(row["O2_pct"])): float(row[index + "_mean"])
                      for row in rows if row["output"] == output}
            grid = np.full((len(parameters), len(o2_values)), np.nan)
            for i, parameter in enumerate(parameters):
                for j, o2 in enumerate(o2_values):
                    if parameter in ACTIVE:
                        grid[i, j] = lookup[(parameter, o2)]
            fig, ax = plt.subplots(figsize=(12, 6.5), constrained_layout=True)
            cmap = plt.colormaps["viridis"].copy()
            cmap.set_bad("#dddddd")
            im = ax.imshow(grid, aspect="auto", origin="upper", extent=[0, 5, len(parameters), 0],
                           vmin=0, vmax=1, cmap=cmap, interpolation="nearest")
            ax.set_yticks(np.arange(len(parameters)) + .5, parameters)
            ax.set_xticks(np.arange(0, 5.1, .5))
            ax.set_xlabel("Fixed oxygen (%)")
            ax.set_title("%s | %s | eFAST N=%d, mean across replicates" %
                         ("Dominant mean ploidy" if output == OUTPUTS[0] else "Asymptotic net live growth (day$^{-1}$)",
                          "First order S1" if index == "S1" else "Total effect ST", max_n))
            for boundary in group_boundaries:
                ax.axhline(boundary, color="white", linewidth=1)
            for row_index, parameter in enumerate(parameters):
                if parameter in STRUCTURAL:
                    ax.text(2.5, row_index + .5, "N/A for fixed oxygen", ha="center",
                            va="center", fontsize=8, color="#555555")
            fig.colorbar(im, ax=ax, label=index + " variance fraction")
            stem = "%s_%s" % (output, index)
            fig.savefig(fig_dir / (stem + ".png"), dpi=250)
            fig.savefig(fig_dir / (stem + ".pdf"))
            plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(19, 12), sharex=True, sharey=True,
                             constrained_layout=True)
    for i, output in enumerate(OUTPUTS):
        for j, index in enumerate(("S1", "ST")):
            ax = axes[i, j]
            lookup = {(row["parameter"], float(row["O2_pct"])): float(row[index + "_mean"])
                      for row in rows if row["output"] == output}
            grid = np.full((len(parameters), len(o2_values)), np.nan)
            for pi, parameter in enumerate(parameters):
                if parameter in ACTIVE:
                    for oi, o2 in enumerate(o2_values):
                        grid[pi, oi] = lookup[(parameter, o2)]
            cmap = plt.colormaps["viridis"].copy()
            cmap.set_bad("#dddddd")
            im = ax.imshow(grid, aspect="auto", origin="upper",
                           extent=[0, 5, len(parameters), 0], vmin=0, vmax=1,
                           cmap=cmap, interpolation="nearest")
            for boundary in group_boundaries:
                ax.axhline(boundary, color="white", linewidth=1)
            ax.set_yticks(np.arange(len(parameters)) + .5, parameters)
            ax.tick_params(axis="y", labelleft=(j == 0))
            ax.set_xticks(np.arange(0, 5.1, .5))
            ax.set_title(("Ploidy" if i == 0 else "Net live growth") + " | " + index)
            if i == 1:
                ax.set_xlabel("Fixed oxygen (%)")
            for row_index, parameter in enumerate(parameters):
                if parameter in STRUCTURAL:
                    ax.text(2.5, row_index + .5, "N/A for fixed oxygen", ha="center",
                            va="center", fontsize=7, color="#555555")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="eFAST variance fraction", shrink=.86)
    fig.suptitle("Independent in-vivo eFAST | N=%d | mean across phase replicates" % max_n)
    fig.savefig(fig_dir / "efast_four_panel.png", dpi=250)
    fig.savefig(fig_dir / "efast_four_panel.pdf")
    plt.close(fig)

    # Existing Figure 4B correlations retain the sign that FAST indices omit.
    correlation = read_table(root / "figure4b_spearman_source.tsv")
    for name in ("parameter", "O2_pct", "spearman_rho"):
        if name not in correlation[0]:
            raise ValueError("Missing Figure 4 correlation column " + name)
    lookup = {(row["parameter"], float(row["O2_pct"])): float(row["spearman_rho"])
              for row in correlation}
    grid = np.array([[lookup[(p, o2)] for o2 in o2_values] for p in parameters])
    fig, ax = plt.subplots(figsize=(12, 6.5), constrained_layout=True)
    im = ax.imshow(grid, aspect="auto", origin="upper", extent=[0, 5, len(parameters), 0],
                   norm=TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1), cmap="coolwarm", interpolation="nearest")
    ax.set_yticks(np.arange(len(parameters)) + .5, parameters)
    ax.set_xticks(np.arange(0, 5.1, .5))
    ax.set_xlabel("Fixed oxygen (%)")
    ax.set_title("Figure 4B fitted-solution Spearman correlation (directional context)")
    fig.colorbar(im, ax=ax, label="Spearman rho")
    fig.savefig(fig_dir / "figure4b_spearman_direction.png", dpi=250)
    fig.savefig(fig_dir / "figure4b_spearman_direction.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--fit-root", required=True)
    p.add_argument("--figure4-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--replicates", type=int, default=2)
    p.add_argument("--seed", type=int, default=20260923)
    p.add_argument("--m", type=int, default=4)
    p.add_argument("--oxygen", default="0,0.5,2.5,5")
    p.add_argument("--full-grid", action="store_true")
    a = sub.add_parser("summarize")
    a.add_argument("--out-dir", required=True)
    f = sub.add_parser("plot")
    f.add_argument("--out-dir", required=True)
    f.add_argument("--figure4-dir", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args)
    elif args.command == "summarize":
        summarize(args)
    else:
        plot(args)


if __name__ == "__main__":
    main()
