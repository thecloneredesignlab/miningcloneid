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


def salib_version():
    """Return the SALib version only for commands that require SALib."""
    return version("SALib")


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
    from SALib.sample import fast_sampler

    if args.n <= 4 * args.m * args.m:
        raise ValueError("FAST requires N > 4*M^2; choose N > %d" % (4 * args.m * args.m))
    fit_table, ranges = get_ranges(args.fit_root)
    ranges_file = getattr(args, "ranges_file", None)
    if ranges_file:
        ranges = read_table(ranges_file)
        for row in ranges:
            for field in ("encoded_lower", "encoded_upper", "natural_lower", "natural_upper"):
                row[field] = float(row[field])
        if [row["parameter"] for row in ranges] != list(ACTIVE):
            raise ValueError("Neighborhood bounds are not in FAST parameter order")
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
        metadata_path = run_dir / "metadata.json"
        sample_path = run_dir / "samples.tsv.gz"
        identity = {
            "salib_version": salib_version(), "method": "eFAST", "M": args.m,
            "N": args.n, "replicate": replicate, "seed": seed,
            "n_parameters": len(ACTIVE), "n_samples": args.n * len(ACTIVE),
            "fit_parameter_table_sha256": sha256(fit_table),
            "figure4_grid_sha256": sha256(grid_path), "oxygen_pct": oxygen,
        }
        fit_seed = getattr(args, "fit_seed", "seed25")
        if ranges_file:
            identity.update({
                "fit_seed": fit_seed, "scope": "neighborhood10pct",
                "bounds_sha256": sha256(ranges_file),
                "best_params_sha256": sha256(Path(args.fit_root) / fit_seed / "best_params.tsv"),
                "fit_config_sha256": sha256(Path(args.fit_root) / fit_seed / "fit_config.rds"),
            })
        if metadata_path.exists():
            previous = json.loads(metadata_path.read_text())
            if any(previous.get(k) != v for k, v in identity.items()):
                raise ValueError("Existing design metadata differs: " + str(run_dir))
            if not sample_path.exists() or sha256(sample_path) != previous["samples_sha256"]:
                raise ValueError("Existing design sample hash differs: " + str(run_dir))
            print("Retained verified design " + str(run_dir), flush=True)
            continue
        if (run_dir / "outputs.tsv.gz").exists() or sample_path.exists():
            raise ValueError("Orphan design files; inspect before resuming: " + str(run_dir))
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
        write_table(sample_path, ["sample_id"] + list(ACTIVE), sample_rows)
        metadata = {
            "salib_version": salib_version(), "method": "eFAST", "M": args.m,
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
        metadata.update(identity)
        if ranges_file:
            metadata["sampling"] = "independent natural-span +/-10% clipped bounds; inherited transforms"
        with open(metadata_path, "w") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
        print("Prepared %s: %d parameter vectors x %d oxygen points" %
              (run_dir, len(natural), len(oxygen)), flush=True)


def summarize(args):
    from SALib.analyze import fast

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
        {"field": "salib_version", "value": salib_version()},
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
    write_global_interpretation(root)


def write_global_interpretation(root):
    """Regenerate the narrative from completed phase summaries, never stale R1/R2 values."""
    conv = read_table(root / "convergence.tsv")
    top_n = max(int(r["N"]) for r in conv)
    top = [r for r in conv if int(r["N"]) == top_n]
    repeats = sorted({int(r["replicates"]) for r in top})
    groups = read_table(root / "mechanism_band_summary.tsv")
    lines = ["# Independent in-vivo fixed-oxygen eFAST", "",
        "The global analysis varies 14 parameters independently over the original in-vivo fit bounds, using the documented identity/log10 transforms. The 201 oxygen points span 0–5% by 0.025%.", "",
        "Highest resolution: N=%d; phase repetitions per cell: %s. S1/ST have no sign. The existing Figure 4B Spearman panel supplies direction." % (top_n, repeats), "",
        "## Mechanism comparison", ""]
    for output in OUTPUTS:
        lookup = {(r["oxygen_band"], r["group"]): r for r in groups if r["output"] == output}
        for index in ("S1", "ST"):
            field = index + "_mean_per_parameter"
            low = {g: float(lookup["low_0_to_1pct", g][field]) for g in set(GROUP.values())}
            high = {g: float(lookup["high_3_to_5pct", g][field]) for g in set(GROUP.values())}
            conditions = (low["death"] > max(v for g, v in low.items() if g != "death"),
                          high["missegregation"] > low["missegregation"], high["buffering"] > low["buffering"])
            lines.append("- %s, %s: low-O2 death group mean %.4f; low-O2 largest group %s (%.4f). Missegregation %.4f -> %.4f and buffering %.4f -> %.4f from low to high O2. All three proposed conditions: %s." %
                (output, index, low["death"], max(low, key=low.get), max(low.values()),
                 low["missegregation"], high["missegregation"], low["buffering"], high["buffering"],
                 "met in these point estimates" if all(conditions) else "not met in these point estimates"))
    lines += ["", "## Numerical diagnostics", ""]
    lower_ns = sorted({int(r["N"]) for r in conv if int(r["N"]) < top_n})
    for output in OUTPUTS:
        for index in ("S1", "ST"):
            repeat = np.percentile([float(r[index + "_range"]) for r in top if r["output"] == output], 90)
            delta = np.percentile([float(r[index + "_delta_vs_max_N"]) for r in conv
                    if lower_ns and int(r["N"]) == lower_ns[-1] and r["output"] == output], 90) if lower_ns else float("nan")
            lines.append("- %s, %s: p90 phase range %.4f; p90 change from N=%s to N=%d %.4f." %
                         (output, index, repeat, lower_ns[-1] if lower_ns else "NA", top_n, delta))
    lines += ["", "Five repetitions measure phase variability; they do not establish resolution convergence. Fine ST rankings require caution wherever phase ranges or resolution changes remain large. Near-degenerate leading modes are separately recorded in spectral_gap_qc.tsv.", "",
        "Group values are means of parameter indices, not additive group variance fractions. Total effects overlap through interactions. These indices are conditional on the chosen ranges and independent input distributions, not posterior uncertainty. The separate neighborhood10pct analysis examines robustness around each fitted endpoint.", ""]
    (root / "interpretation.md").write_text("\n".join(lines))


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


def benjamini_hochberg(p_values):
    """Benjamini-Hochberg adjusted p values in the original row order."""
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1 or len(values) == 0 or np.any(~np.isfinite(values)):
        raise ValueError("BH adjustment requires finite p values")
    order = np.argsort(values)
    ranked = values[order]
    adjusted_ranked = np.minimum.accumulate(
        (ranked * len(values) / np.arange(1, len(values) + 1))[::-1]
    )[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def normalized_trapezoid_weights(o2_values, lower, upper):
    """Normalized trapezoid weights for one closed oxygen window."""
    grid = np.asarray(o2_values, dtype=float)
    selected = np.flatnonzero((grid >= lower - 1e-12) & (grid <= upper + 1e-12))
    window_grid = grid[selected]
    if (len(window_grid) < 2 or abs(window_grid[0] - lower) > 1e-12 or
            abs(window_grid[-1] - upper) > 1e-12):
        raise ValueError("The oxygen grid does not span an exact window boundary")
    delta = np.diff(window_grid)
    local = np.zeros(len(window_grid))
    local[0] = delta[0] / 2
    local[-1] = delta[-1] / 2
    if len(window_grid) > 2:
        local[1:-1] = (delta[:-1] + delta[1:]) / 2
    weights = np.zeros(len(grid))
    weights[selected] = local / (upper - lower)
    if abs(weights.sum() - 1) > 1e-12:
        raise ValueError("Normalized oxygen-window weights do not sum to one")
    return weights


def derive_efast_o2_classification(root, max_n, o2_values, parameter_order,
                                   bootstrap_reps=5000, bootstrap_seed=5826,
                                   minimum_peak=0.3):
    """Classify eFAST rows with the Figure 4B oxygen-window decision rule.

    The classification curve uses dominant mean ploidy only. Complete
    phase-repeat curves are the bootstrap unit. S1 and ST are classified
    independently; the growth heatmap follows the corresponding ploidy order.
    """
    raw_rows = [row for row in read_table(Path(root) / "indices.tsv")
                if int(row["N"]) == max_n]
    replicates = sorted({int(row["replicate"]) for row in raw_rows})
    if len(replicates) < 2:
        raise ValueError("eFAST oxygen classification requires at least two phase repeats")
    if bootstrap_reps < 1000:
        raise ValueError("At least 1000 eFAST bootstrap replicates are required")
    low_weights = normalized_trapezoid_weights(o2_values, 0, 1)
    high_weights = normalized_trapezoid_weights(o2_values, 3, 5)
    o2_index = {value: index for index, value in enumerate(o2_values)}
    parameters = list(parameter_order)
    expected = len(replicates) * len(OUTPUTS) * len(parameters) * len(o2_values)
    if len(raw_rows) != expected:
        raise ValueError("Incomplete high-resolution eFAST replicate table")

    raw_lookup = {}
    for row in raw_rows:
        key = (int(row["replicate"]), row["output"], row["parameter"],
               float(row["O2_pct"]))
        if key in raw_lookup:
            raise ValueError("Duplicate high-resolution eFAST replicate row")
        raw_lookup[key] = row

    rng = np.random.default_rng(bootstrap_seed)
    bootstrap_index = rng.integers(
        0, len(replicates), size=(bootstrap_reps, len(replicates)))
    group_levels = ("High O2", "Low O2", "O2-independent")
    all_rows = []
    for index in ("S1", "ST"):
        provisional = []
        for parameter in parameters:
            phase_curves = []
            for replicate in replicates:
                curve = np.empty(len(o2_values))
                for o2 in o2_values:
                    curve[o2_index[o2]] = float(raw_lookup[
                        (replicate, "dominant_mean_ploidy", parameter, o2)][index])
                phase_curves.append(curve)
            phase_curves = np.asarray(phase_curves)
            observed_curve = phase_curves.mean(axis=0)
            low_score = float(observed_curve @ low_weights)
            high_score = float(observed_curve @ high_weights)
            delta = low_score - high_score
            peak = float(observed_curve.max())
            bootstrap_curves = phase_curves[bootstrap_index].mean(axis=1)
            bootstrap_delta = (bootstrap_curves @ low_weights -
                               bootstrap_curves @ high_weights)
            lower_tail = (np.count_nonzero(bootstrap_delta <= 0) + 1) / (
                bootstrap_reps + 1)
            upper_tail = (np.count_nonzero(bootstrap_delta >= 0) + 1) / (
                bootstrap_reps + 1)
            provisional.append({
                "index": index,
                "parameter": parameter,
                "low_o2_score": low_score,
                "high_o2_score": high_score,
                "low_minus_high": delta,
                "global_peak": peak,
                "bootstrap_sign_p_value": min(1.0, 2 * min(lower_tail, upper_tail)),
            })
        adjusted = benjamini_hochberg(
            [row["bootstrap_sign_p_value"] for row in provisional])
        for row, q_value in zip(provisional, adjusted):
            peak_passes = row["global_peak"] > minimum_peak
            if not peak_passes:
                group = "O2-independent"
                rule = "Global peak <= 0.3"
            elif q_value >= 0.05:
                group = "O2-independent"
                rule = "No Low-High difference at BH q < 0.05"
            elif row["low_minus_high"] < 0:
                group = "High O2"
                rule = "High [3,5] exceeds Low [0,1]"
            elif row["low_minus_high"] > 0:
                group = "Low O2"
                rule = "Low [0,1] exceeds High [3,5]"
            else:
                group = "O2-independent"
                rule = "Zero Low-High contrast"
            row.update({
                "o2_sensitivity_group": group,
                "o2_sensitivity_group_order": group_levels.index(group) + 1,
                "bh_adjusted_p_value": float(q_value),
                "minimum_global_peak": minimum_peak,
                "global_peak_passes": str(peak_passes).upper(),
                "decision_rule": rule,
                "bootstrap_reps": bootstrap_reps,
                "bootstrap_seed": bootstrap_seed,
                "bootstrap_unit": "complete phase-repeat curve",
                "n_phase_repeats": len(replicates),
                "classification_output": "dominant_mean_ploidy",
                "parameter_order": parameter_order[row["parameter"]],
            })
        provisional.sort(key=lambda row: (
            row["o2_sensitivity_group_order"], -row["global_peak"],
            row["parameter_order"]))
        group_rank = {group: 0 for group in group_levels}
        for display_order, row in enumerate(provisional, 1):
            group_rank[row["o2_sensitivity_group"]] += 1
            row["within_group_rank"] = group_rank[row["o2_sensitivity_group"]]
            row["display_order"] = display_order
            all_rows.append(row)

    fields = [
        "index", "parameter", "o2_sensitivity_group",
        "o2_sensitivity_group_order", "display_order", "within_group_rank",
        "low_o2_score", "high_o2_score", "low_minus_high", "global_peak",
        "minimum_global_peak", "global_peak_passes", "bootstrap_sign_p_value",
        "bh_adjusted_p_value", "decision_rule", "bootstrap_reps",
        "bootstrap_seed", "bootstrap_unit", "n_phase_repeats",
        "classification_output", "parameter_order",
    ]
    formatted = []
    for row in all_rows:
        formatted.append({
            key: ("%.10g" % value if isinstance(value, float) else value)
            for key, value in row.items()
        })
    write_table(Path(root) / "efast_o2_sensitivity_classification.tsv", fields, formatted)
    return all_rows


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
    from matplotlib.patches import Patch

    root = Path(args.out_dir)
    rows = read_table(root / "convergence.tsv")
    max_n = max(int(row["N"]) for row in rows)
    rows = [row for row in rows if int(row["N"]) == max_n]
    o2_values = sorted(set(float(row["O2_pct"]) for row in rows))
    if len(o2_values) != 201:
        raise ValueError("Publication heatmaps require the full Figure 4 oxygen grid")
    fig_dir = root / "figures"
    fig_dir.mkdir(exist_ok=True)

    # The row classifications and order are defined by the iteration5 Figure 4B
    # window analysis. Keep local copies beside the eFAST numerical results so
    # the figure can be reproduced without recomputing that independent analysis.
    source_names = {
        "continuous_ploidy_o2_window_classification.tsv":
            "figure4b_o2_classification_source.tsv",
        "continuous_ploidy_parameter_ranking.tsv":
            "figure4b_parameter_ranking_source.tsv",
        "parameter_function_groups.tsv": "figure4_parameter_groups_source.tsv",
        "parameter_function_group_palette.tsv":
            "figure4_parameter_group_palette_source.tsv",
    }
    if args.figure4_layout_dir:
        layout_dir = Path(args.figure4_layout_dir)
        for source_name, copy_name in source_names.items():
            source_path = layout_dir / source_name
            if not source_path.exists():
                raise FileNotFoundError(source_path)
            copy_path = root / copy_name
            if not copy_path.exists() or sha256(copy_path) != sha256(source_path):
                shutil.copyfile(source_path, copy_path)
    for copy_name in source_names.values():
        if not (root / copy_name).exists():
            raise FileNotFoundError(
                "%s; pass --figure4-layout-dir to copy the iteration5 Figure 4B sources" %
                (root / copy_name))

    figure4_classification = [
        row for row in read_table(root / "figure4b_o2_classification_source.tsv")
        if row["parameter"] in ACTIVE
    ]
    if (len(figure4_classification) != len(ACTIVE) or
            {row["parameter"] for row in figure4_classification} != set(ACTIVE)):
        raise ValueError("Figure 4B classification does not contain the 14 eFAST parameters exactly once")
    figure4_group_by_parameter = {
        row["parameter"]: row["o2_association_group"]
        for row in figure4_classification
    }
    o2_group_levels = ("High O2", "Low O2", "O2-independent")
    if any(group not in o2_group_levels for group in figure4_group_by_parameter.values()):
        raise ValueError("Unexpected Figure 4B O2 association group")

    function_rows = read_table(root / "figure4_parameter_groups_source.tsv")
    function_by_parameter = {row["parameter"]: row["parameter_group"]
                             for row in function_rows}
    configured_parameter_order = {
        row["parameter"]: int(row["parameter_order"])
        for row in function_rows if row["parameter"] in ACTIVE
    }
    if (not set(ACTIVE).issubset(function_by_parameter) or
            set(configured_parameter_order) != set(ACTIVE)):
        raise ValueError("Missing parameter-function annotation")

    classification_file = getattr(args, "classification_file", None)
    efast_classification = (read_table(classification_file) if classification_file else
        derive_efast_o2_classification(root, max_n, o2_values, configured_parameter_order))
    efast_rows_by_index = {
        index: sorted(
            [row for row in efast_classification if row["index"] == index],
            key=lambda row: int(row["display_order"]),
        )
        for index in ("S1", "ST")
    }
    panel_parameters = {
        index: [row["parameter"] for row in efast_rows_by_index[index]]
        for index in ("S1", "ST")
    }
    if any(set(panel_parameters[index]) != set(ACTIVE) for index in ("S1", "ST")):
        raise ValueError("eFAST oxygen classification does not contain all sampled parameters")
    efast_group_by_index = {
        index: {row["parameter"]: row["o2_sensitivity_group"]
                for row in efast_rows_by_index[index]}
        for index in ("S1", "ST")
    }

    palette_rows = sorted(
        read_table(root / "figure4_parameter_group_palette_source.tsv"),
        key=lambda row: int(row["group_order"]),
    )
    function_palette = {row["parameter_group"]: row["color"] for row in palette_rows}
    function_labels = {row["parameter_group"]: row["display_label"] for row in palette_rows}
    if any(group not in function_palette for group in function_by_parameter.values()):
        raise ValueError("Missing parameter-function color")
    o2_palette = {
        "High O2": "#B2182B",
        "Low O2": "#2166AC",
        "O2-independent": "#8A8A8A",
    }
    o2_labels = {
        "High O2": "High O2",
        "Low O2": "Low O2",
        "O2-independent": "O2-independent",
    }
    def grid_for(output, index, parameters):
        lookup = {(row["parameter"], float(row["O2_pct"])):
                  float(row[index + "_" + getattr(args, "statistic", "mean")])
                  for row in rows if row["output"] == output}
        expected = len(parameters) * len(o2_values)
        if len(lookup) != expected:
            raise ValueError("Incomplete %s %s grid: %d of %d cells" %
                             (output, index, len(lookup), expected))
        return np.array([[lookup[(parameter, o2)] for o2 in o2_values]
                         for parameter in parameters])

    grids = {(index, output): grid_for(output, index, panel_parameters[index])
             for index in ("S1", "ST") for output in OUTPUTS}
    output_max = {(index, output): float(np.nanmax(grids[index, output]))
                  for index in ("S1", "ST") for output in OUTPUTS}
    if any(not math.isfinite(value) or value <= 0 for value in output_max.values()):
        raise ValueError("Invalid eFAST color scale maximum")
    oxygen_array = np.asarray(o2_values)
    oxygen_edges = np.empty(len(oxygen_array) + 1)
    oxygen_edges[0] = 0
    oxygen_edges[1:-1] = (oxygen_array[:-1] + oxygen_array[1:]) / 2
    oxygen_edges[-1] = oxygen_array[-1] + (oxygen_array[-1] - oxygen_array[-2]) / 2
    oxygen_ticks = (0, 0.025, 0.1, 0.5, 1, 2, 5)
    oxygen_tick_labels = ("0", ".025", ".1", ".5", "1", "2", "5")
    oxygen_window_boundaries = (1.0, 3.0)

    # Two main panels: A contains S1 for both outputs and B contains ST for
    # both outputs. Ploidy and growth use separate data-driven color scales.
    index_styles = {
        "S1": {
            "letter": "A", "title": "First-order effects (S1)",
            "cmap": LinearSegmentedColormap.from_list(
                "white_to_deep_purple", ["#FFFFFF", "#3F007D"]),
            "color": "#3F007D",
        },
        "ST": {
            "letter": "B", "title": "Total effects (ST)",
            "cmap": LinearSegmentedColormap.from_list(
                "white_to_orange_st", ["#FFFFFF", "#E6550D"]),
            "color": "#E6550D",
        },
    }
    output_titles = {
        "dominant_mean_ploidy": "Dominant mean ploidy",
        "dominant_growth_rate": "Asymptotic net live-cell growth rate",
    }

    fig = plt.figure(figsize=(18, 13.5))
    outer = fig.add_gridspec(2, 1, left=.19, right=.965, top=.95, bottom=.15,
                             hspace=.30)
    panel_title_y = {"S1": .968, "ST": .514}
    for panel_index, index in enumerate(("S1", "ST")):
        parameters = panel_parameters[index]
        function_groups = [function_by_parameter[parameter] for parameter in parameters]
        figure4_groups = [figure4_group_by_parameter[parameter] for parameter in parameters]
        efast_groups = [efast_group_by_index[index][parameter] for parameter in parameters]
        group_boundaries = [i for i in range(1, len(parameters))
                            if efast_groups[i] != efast_groups[i - 1]]
        inner = outer[panel_index].subgridspec(
            2, 6, width_ratios=(.14, .14, .14, 3.5, .25, 3.5),
            height_ratios=(1, .055), wspace=.055, hspace=.28)
        function_ax = fig.add_subplot(inner[0, 0])
        figure4_ax = fig.add_subplot(inner[0, 1], sharey=function_ax)
        efast_ax = fig.add_subplot(inner[0, 2], sharey=function_ax)
        heat_axes = [fig.add_subplot(inner[0, column], sharey=function_ax)
                     for column in (3, 5)]
        color_axes = [fig.add_subplot(inner[1, column]) for column in (3, 5)]

        function_codes = np.array([
            list(function_palette).index(group) for group in function_groups
        ])[:, None]
        function_cmap = ListedColormap(list(function_palette.values()))
        function_ax.imshow(function_codes, aspect="auto", origin="upper",
                           extent=[0, 1, len(parameters), 0], cmap=function_cmap,
                           vmin=-.5, vmax=len(function_palette) - .5,
                           interpolation="nearest")
        function_ax.set_yticks(np.arange(len(parameters)) + .5, parameters)
        function_ax.tick_params(axis="y", labelsize=9.5, length=0, pad=7)
        function_ax.set_xticks([])
        function_ax.set_title("1", fontsize=8.5, pad=5)

        o2_cmap = ListedColormap([o2_palette[group] for group in o2_group_levels])
        for annotation_ax, groups, title in (
                (figure4_ax, figure4_groups, "2"),
                (efast_ax, efast_groups, "3")):
            codes = np.array([o2_group_levels.index(group) for group in groups])[:, None]
            annotation_ax.imshow(
                codes, aspect="auto", origin="upper",
                extent=[0, 1, len(parameters), 0], cmap=o2_cmap,
                vmin=-.5, vmax=len(o2_group_levels) - .5,
                interpolation="nearest")
            annotation_ax.tick_params(axis="y", labelleft=False, left=False)
            annotation_ax.set_xticks([])
            annotation_ax.set_title(title, fontsize=8.5, pad=5)

        for output, ax, color_ax in zip(OUTPUTS, heat_axes, color_axes):
            image = ax.pcolormesh(
                oxygen_edges, np.arange(len(parameters) + 1), grids[index, output],
                shading="flat", vmin=0,
                vmax=output_max[index, output], cmap=index_styles[index]["cmap"],
                rasterized=True,
            )
            ax.set_xscale("symlog", base=10, linthresh=.025, linscale=1)
            ax.set_xlim(0, oxygen_edges[-1])
            ax.set_ylim(len(parameters), 0)
            ax.set_yticks(np.arange(len(parameters)) + .5)
            ax.tick_params(axis="y", labelleft=False, left=False)
            ax.set_xticks(oxygen_ticks, oxygen_tick_labels)
            ax.tick_params(axis="x", labelsize=8.5)
            ax.set_xlabel("Fixed oxygen (%)", labelpad=4)
            ax.set_title(output_titles[output], fontsize=11, pad=8)
            for boundary_o2 in oxygen_window_boundaries:
                ax.axvline(
                    boundary_o2, color="#000000", linewidth=.9,
                    linestyle=(0, (4, 3)), zorder=3,
                )
            for boundary in group_boundaries:
                ax.axhline(boundary, color="#FFFFFF", linewidth=1.1)
            colorbar = fig.colorbar(image, cax=color_ax, orientation="horizontal")
            colorbar_ticks = np.linspace(0, output_max[index, output], 3)
            colorbar.set_ticks(colorbar_ticks)
            colorbar.set_ticklabels(["%.4f" % value for value in colorbar_ticks])
            colorbar.ax.tick_params(labelsize=8, length=2, pad=2)
            colorbar.outline.set_linewidth(.6)
        for annotation_ax in (function_ax, figure4_ax, efast_ax):
            for boundary in group_boundaries:
                annotation_ax.axhline(boundary, color="#FFFFFF", linewidth=1.1)
            for spine in annotation_ax.spines.values():
                spine.set_color("#444444")
                spine.set_linewidth(.6)

        fig.text(.032, panel_title_y[index], index_styles[index]["letter"],
                 fontsize=18, fontweight="bold", va="top")
        fig.text(.058, panel_title_y[index], index_styles[index]["title"],
                 fontsize=13, fontweight="bold", va="top")

    function_handles = [
        Patch(facecolor=function_palette[row["parameter_group"]], edgecolor="none",
              label=function_labels[row["parameter_group"]])
        for row in palette_rows
    ]
    o2_handles = [
        Patch(facecolor=o2_palette[group], edgecolor="none", label=o2_labels[group])
        for group in o2_group_levels
    ]
    process_legend = fig.legend(
        handles=function_handles, loc="lower left", bbox_to_anchor=(.19, .055),
        ncol=5, frameon=False, title="1  Process", fontsize=8.5,
        title_fontsize=9.5, handlelength=1.2, columnspacing=1.4,
    )
    fig.add_artist(process_legend)
    fig.legend(
        handles=o2_handles, loc="lower left", bbox_to_anchor=(.19, .015),
        ncol=3, frameon=False, title="2  O2 correlation     3  Sensitivity",
        fontsize=8.5,
        title_fontsize=9.5, handlelength=1.2, columnspacing=1.8,
    )

    combined_png = fig_dir / "efast_four_panel.png"
    combined_pdf = fig_dir / "efast_four_panel.pdf"
    fig.savefig(combined_png, dpi=250, facecolor="white")
    fig.savefig(combined_pdf, facecolor="white")
    plt.close(fig)

    figure_manifest = [
        {"field": "figure", "value": "efast_four_panel"},
        {"field": "max_N", "value": max_n},
        {"field": "replicate_summary", "value": getattr(args, "summary_label", None) or
         "mean across %s phase replicates" % rows[0]["replicates"]},
        {"field": "n_parameters", "value": len(ACTIVE)},
        {"field": "S1_parameter_order", "value": ",".join(panel_parameters["S1"])},
        {"field": "ST_parameter_order", "value": ",".join(panel_parameters["ST"])},
        {"field": "efast_o2_group_windows", "value": "Low [0,1]; High [3,5]"},
        {"field": "efast_o2_group_classification_output",
         "value": "dominant_mean_ploidy"},
        {"field": "efast_o2_group_bootstrap_unit", "value": efast_classification[0]["bootstrap_unit"]},
        {"field": "efast_o2_group_bootstrap_reps", "value": "5000"},
        {"field": "efast_o2_group_bootstrap_seed", "value": "5826"},
        {"field": "efast_o2_group_bh_scope", "value": "14 parameters separately for S1 and ST"},
        {"field": "efast_o2_group_peak_gate", "value": "global peak strictly greater than 0.3"},
        {"field": "oxygen_axis_scale", "value": "symlog base 10; linear threshold 0.025 percent"},
        {"field": "oxygen_window_boundary_lines",
         "value": "black dashed lines at 1 and 3 percent oxygen"},
        {"field": "colorbar_position", "value": "horizontal below each heatmap"},
        {"field": "S1_vmin", "value": "0"},
        {"field": "S1_ploidy_vmax",
         "value": "%.10g" % output_max["S1", "dominant_mean_ploidy"]},
        {"field": "S1_growth_vmax",
         "value": "%.10g" % output_max["S1", "dominant_growth_rate"]},
        {"field": "S1_color", "value": index_styles["S1"]["color"]},
        {"field": "ST_vmin", "value": "0"},
        {"field": "ST_ploidy_vmax",
         "value": "%.10g" % output_max["ST", "dominant_mean_ploidy"]},
        {"field": "ST_growth_vmax",
         "value": "%.10g" % output_max["ST", "dominant_growth_rate"]},
        {"field": "ST_color", "value": index_styles["ST"]["color"]},
        {"field": "convergence_sha256", "value": sha256(root / "convergence.tsv")},
        {"field": "classification_sha256",
         "value": sha256(root / "figure4b_o2_classification_source.tsv")},
        {"field": "ranking_sha256",
         "value": sha256(root / "figure4b_parameter_ranking_source.tsv")},
        {"field": "parameter_groups_sha256",
         "value": sha256(root / "figure4_parameter_groups_source.tsv")},
        {"field": "parameter_group_palette_sha256",
         "value": sha256(root / "figure4_parameter_group_palette_source.tsv")},
        {"field": "efast_o2_sensitivity_classification_sha256",
         "value": sha256(root / "efast_o2_sensitivity_classification.tsv")},
        {"field": "figure_png_sha256", "value": sha256(combined_png)},
        {"field": "figure_pdf_sha256", "value": sha256(combined_pdf)},
    ]
    write_table(root / "figure_redraw_manifest.tsv", ["field", "value"], figure_manifest)

    if args.combined_only:
        return

    for index in ("S1", "ST"):
        parameters = panel_parameters[index]
        efast_groups = [efast_group_by_index[index][parameter] for parameter in parameters]
        group_boundaries = [i for i in range(1, len(parameters))
                            if efast_groups[i] != efast_groups[i - 1]]
        for output in OUTPUTS:
            fig, ax = plt.subplots(figsize=(12, 6.5), constrained_layout=True)
            im = ax.pcolormesh(
                oxygen_edges, np.arange(len(parameters) + 1), grids[index, output],
                shading="flat", vmin=0, vmax=output_max[index, output],
                cmap=index_styles[index]["cmap"], rasterized=True)
            ax.set_xscale("symlog", base=10, linthresh=.025, linscale=1)
            ax.set_xlim(0, oxygen_edges[-1])
            ax.set_ylim(len(parameters), 0)
            ax.set_yticks(np.arange(len(parameters)) + .5, parameters)
            ax.set_xticks(oxygen_ticks, oxygen_tick_labels)
            ax.set_xlabel("Fixed oxygen (%)")
            ax.set_title("%s | %s | eFAST N=%d, mean across replicates" %
                         (output_titles[output], index_styles[index]["title"], max_n))
            for boundary in group_boundaries:
                ax.axhline(boundary, color="white", linewidth=1)
            fig.colorbar(im, ax=ax, orientation="horizontal", pad=.14,
                         label=index + " variance fraction")
            stem = "%s_%s" % (output, index)
            fig.savefig(fig_dir / (stem + ".png"), dpi=250)
            fig.savefig(fig_dir / (stem + ".pdf"))
            plt.close(fig)

    # Existing Figure 4B correlations retain the sign that FAST indices omit.
    correlation = read_table(root / "figure4b_spearman_source.tsv")
    for name in ("parameter", "O2_pct", "spearman_rho"):
        if name not in correlation[0]:
            raise ValueError("Missing Figure 4 correlation column " + name)
    lookup = {(row["parameter"], float(row["O2_pct"])): float(row["spearman_rho"])
              for row in correlation}
    parameters = [row["parameter"] for row in sorted(
        figure4_classification, key=lambda row: int(row["display_order"]))]
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
    p.add_argument("--replicates", type=int, default=5)
    p.add_argument("--ranges-file", help="Audited per-fit-seed neighborhood bounds")
    p.add_argument("--fit-seed", default="seed25")
    p.add_argument("--seed", type=int, default=20260923)
    p.add_argument("--m", type=int, default=4)
    p.add_argument("--oxygen", default="0,0.5,2.5,5")
    p.add_argument("--full-grid", action="store_true")
    a = sub.add_parser("summarize")
    a.add_argument("--out-dir", required=True)
    f = sub.add_parser("plot")
    f.add_argument("--out-dir", required=True)
    f.add_argument("--figure4-dir", required=True)
    f.add_argument("--figure4-layout-dir",
                   help="iteration5 Figure 4 directory containing classification/order sources")
    f.add_argument("--combined-only", action="store_true",
                   help="redraw only efast_four_panel.pdf/png")
    f.add_argument("--classification-file", help="Precomputed neighborhood oxygen classification")
    f.add_argument("--statistic", choices=("mean", "median"), default="mean")
    f.add_argument("--summary-label")
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args)
    elif args.command == "summarize":
        summarize(args)
    else:
        plot(args)


if __name__ == "__main__":
    main()
