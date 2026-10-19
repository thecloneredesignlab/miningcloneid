#!/usr/bin/env python3
"""Local, table-only Sobol redraw. No model evaluation or index estimation."""

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import itertools
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import NullLocator
import numpy as np


SCOPES = ("global", "neighborhood10pct")
OUTPUTS = ("dominant_mean_ploidy", "dominant_growth_rate")
FOCAL = ("p_mis_base", "mu_hp", "gamma_mu")
POSITIVE, NEGATIVE = "#54278f", "#2166ac"
OUT_LABEL = {OUTPUTS[0]: "Dominant mean ploidy", OUTPUTS[1]: "Net live-cell growth rate"}
SCOPE_LABEL = {"global": "Global full bounds", "neighborhood10pct": "500 neighborhoods"}
BANDS = (("low_0_1", 0., 1.), ("middle_1_3", 1., 3.), ("high_3_5", 3., 5.))
OXYGEN_LINTHRESH = .025
OXYGEN_TICKS = (0., .025, .1, .5, 1., 3., 5.)
OXYGEN_NOTE = ("O₂: log₁₀ above 0.025%; linear from 0 to 0.025% to retain zero. "
               "Low: 0–1%; intermediate: >1–<3%; high: 3–5%.")


def oxygen_axis(ax):
    """Keep real zero; place thresholds/labels in the actual oxygen coordinate."""
    ax.set_xscale("symlog", base=10, linthresh=OXYGEN_LINTHRESH, linscale=1)
    ax.set_xlim(0, 5)
    ax.set_xticks(OXYGEN_TICKS, ("0", "0.025", "0.1", "0.5", "1", "3", "5"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Fixed oxygen (%)", fontsize=8)
    transform = ax.xaxis.get_transform()
    for lo, hi, label in ((0., 1., "Low"), (1., 3., "Mid"), (3., 5., "High")):
        middle = float(transform.inverted().transform(np.mean(transform.transform(np.array([lo, hi])))))
        ax.text(middle, 1.012, label, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=7, color="#333333", clip_on=False)
    for x in (1., 3.):
        ax.axvline(x, color="black", lw=.6, ls="--", alpha=.8)


def oxygen_mesh(ax, data, oxygen, **kwargs):
    """Use real cell edges, not equal-width columns relabeled as logarithmic."""
    edges = np.r_[0., (oxygen[:-1] + oxygen[1:]) / 2., 5.]
    assert len(edges) == len(oxygen) + 1 and np.all(np.diff(edges) > 0)
    row_edges = np.arange(data.shape[0] + 1) - .5
    image = ax.pcolormesh(edges, row_edges, np.ma.masked_invalid(data), shading="flat",
                          rasterized=True, antialiased=False, **kwargs)
    ax.set_ylim(data.shape[0] - .5, -.5)
    oxygen_axis(ax)
    return image


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_table(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0])
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig, stem, rendered):
    stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".png"):
        path = stem.with_suffix(suffix)
        temporary = path.with_name(path.name + ".tmp")
        fig.savefig(temporary, format=suffix[1:], dpi=170)
        temporary.replace(path)
        rendered.append(path)
    plt.close(fig)


def load_sources(root):
    sources = {s: root / s / "summaries/indices_scope_summary.tsv.gz" for s in SCOPES}
    hashes = {s: sha(p) for s, p in sources.items()}
    rows = {s: read_table(p) for s, p in sources.items()}
    mapping = read_table(root / "model_structure/parameter_equation_map.tsv")
    active = [r["parameter"] for r in mapping if r["present_fixed_O2"] == "TRUE"]
    groups = {r["parameter"]: r["primary_group"] for r in mapping}
    assert len(active) == 14 and len(set(active)) == 14
    canonical = {frozenset(p): p for p in itertools.combinations(active, 2)}
    focal_pairs, seen = [], set()
    for focus in FOCAL:
        for other in active:
            if other == focus:
                continue
            key = frozenset((focus, other))
            if key in seen:
                continue
            seen.add(key)
            pi, pj = canonical[key]
            focal_pairs.append(dict(row=len(focal_pairs) + 1, focal_parameter=focus,
                                    partner=other, parameter_i=pi, parameter_j=pj,
                                    display_label=focus + " × " + other))
    assert len(focal_pairs) == 36
    oxygen = np.array(sorted({float(r["O2_pct"]) for r in rows["global"]}))
    assert len(oxygen) == 201 and np.allclose(oxygen, np.linspace(0, 5, 201), rtol=0, atol=1e-12)
    lookups = {}
    for scope, data in rows.items():
        lookup = {(r["index"], r["parameter_i"], r["parameter_j"],
                   r["output"], float(r["O2_pct"])): r for r in data}
        assert len(lookup) == len(data) == 47838
        assert {int(r["N"]) for r in data} == {4096}
        assert {int(r["n_endpoints"]) for r in data} == ({1} if scope == "global" else {500})
        assert all(int(r["n_defined"]) == int(r["n_endpoints"]) for r in data)
        lookups[scope] = lookup
    return sources, hashes, rows, active, groups, focal_pairs, oxygen, lookups


def metric_rows(index, active, focal_pairs):
    if index == "S2":
        return [(r["parameter_i"], r["parameter_j"]) for r in focal_pairs]
    return [(p, "") for p in active]


def matrix(lookup, keys, index, output, oxygen, field="value_mean"):
    return np.array([[float(lookup[index, pi, pj, output, float(o)][field])
                      for o in oxygen] for pi, pj in keys])


def palette(vmax, has_negative):
    # Symmetric signed scale: the same absolute estimate has the same saturation.
    # This avoids giving tiny negative sampling estimates half the whole scale.
    if has_negative:
        cmap = LinearSegmentedColormap.from_list("signed_purple", [NEGATIVE, "white", POSITIVE], N=513)
        norm = Normalize(-vmax, vmax)
    else:
        cmap = LinearSegmentedColormap.from_list("positive_purple", ["white", POSITIVE], N=257)
        norm = Normalize(0, vmax)
    cmap.set_bad("#d9d9d9")
    return cmap, norm


def draw(ax, data, oxygen, labels, index, limit=None, signed=None, show_labels=True,
         pair_breaks=(), diagnostic=False):
    finite = data[np.isfinite(data)]
    if diagnostic:
        cmap = LinearSegmentedColormap.from_list("diagnostic", ["#f2f2f2", "#333333"])
        norm = Normalize(0, 1)
    else:
        vmax = max(float(np.max(np.abs(finite))), 1e-9) if limit is None else limit
        cmap, norm = palette(vmax, bool(np.any(finite < 0)) if signed is None else signed)
    im = oxygen_mesh(ax, data, oxygen, cmap=cmap, norm=norm)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels if show_labels else [], fontsize=7.8)
    ax.tick_params(axis="both", length=2, width=.5, labelsize=7.5)
    ax.tick_params(axis="y", labelsize=7.8)
    for label in ax.get_yticklabels():
        if any(label.get_text().split(" [")[0] == p for p in FOCAL):
            label.set_fontweight("bold")
    for y in pair_breaks:
        ax.axhline(y - .5, color="#555555", lw=.65)
    for spine in ax.spines.values():
        spine.set_linewidth(.5)
    cb = ax.figure.colorbar(im, ax=ax, orientation="horizontal", pad=.045, fraction=.04, aspect=35)
    ticks = ([0, .5, 1] if diagnostic else
             ([-norm.vmax, 0., norm.vmax] if norm.vmin < 0 else [0, norm.vmax / 2, norm.vmax]))
    cb.set_ticks(ticks)
    cb.set_ticklabels([format(x, ".3g") for x in ticks])
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("Diagnostic pass fraction" if diagnostic else index + " estimate", fontsize=7.5)
    return dict(vmin=float(norm.vmin), vmax=float(norm.vmax), negative_blue=bool(norm.vmin < 0))


def draw_grid(root, review, lookups, active, groups, focal_pairs, oxygen, scopes, outputs,
              name, title, rendered, scales, diagnostic=False):
    columns = list(zip(scopes, outputs))
    cols = len(columns)
    fig, axes = plt.subplots(3, cols, figsize=(14.5 if cols == 2 else 25, 17.7),
                             gridspec_kw={"height_ratios": [14, 14, 36]}, layout="constrained",
                             squeeze=False)
    fig.set_constrained_layout_pads(w_pad=.05, h_pad=.07, hspace=.07, wspace=.06)
    breaks = [i for i in range(1, 36)
              if focal_pairs[i]["focal_parameter"] != focal_pairs[i - 1]["focal_parameter"]]
    for ri, index in enumerate(("S1", "ST", "S2")):
        keys = metric_rows(index, active, focal_pairs)
        labels = ([r["display_label"] for r in focal_pairs] if index == "S2" else
                  [p + " [" + groups[p] + "]" for p in active])
        for ci, (scope, output) in enumerate(columns):
            ax = axes[ri, ci]
            data = matrix(lookups[scope], keys, index, output, oxygen,
                          "pass_fraction_all_endpoints" if diagnostic else "value_mean")
            limit, signed = scales[index, output]
            draw(ax, data, oxygen, labels, index, limit, signed, ci == 0,
                 breaks if index == "S2" else (), diagnostic)
            panel = chr(ord("A") + ri * cols + ci)
            index_label = "S2: 36 focal pairs" if index == "S2" else index
            ax.set_title(panel + "  " + SCOPE_LABEL[scope] + " | " + index_label + "\n" + OUT_LABEL[output],
                         loc="left", fontsize=9.3, pad=23)
    suffix = ("Global: pass/fail per cell; neighborhoods: fraction of all 500 seeds passing." if diagnostic else
              "Same metric/output uses identical color limits across scopes; outputs and metrics have separate scales.")
    fig.suptitle(title + " | N4096, five repeats, 201 oxygen points\n" + suffix, fontsize=11)
    fig.supxlabel("Global = five-repeat mean; neighborhoods = median of 500 five-repeat means. "
                  "No convergence exclusions.\n" +
                  ("Diagnostic thresholds: repeat range ≤ 0.05 and N2048–N4096 mean change ≤ 0.02." if diagnostic else
                   "Purple = positive estimate; blue = negative sampling estimate (not antagonism); gray = undefined. "
                   "Dashed lines: 1% and 3% O₂.") + "\n" + OXYGEN_NOTE, fontsize=8)
    save_figure(fig, review / name, rendered)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    review = root / "figure_comparison_20261019"
    review.mkdir(parents=True, exist_ok=True)
    sources, hashes, rows, active, groups, pairs, oxygen, lookups = load_sources(root)
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.titlesize": 10, "axes.labelsize": 8, "savefig.facecolor": "white"})
    rendered, scale_rows, original_files = [], [], []
    # Keep collection-time files/receipts intact as provenance. This is a local figure overlay.
    for scope in SCOPES:
        for index in ("s1", "st"):
            for suffix in (".pdf", ".png"):
                path = root / scope / "figures" / (index + "_parameter_by_oxygen" + suffix)
                backup = review / "original_figures" / scope / path.name
                backup.parent.mkdir(parents=True, exist_ok=True)
                if not backup.exists():
                    shutil.copy2(path, backup)
                original_files.append(dict(path=str(path.relative_to(root)),
                                           original_sha256=sha(backup), backup=str(backup.relative_to(root))))
    write_table(review / "focal_pair_order.tsv", pairs)
    selected_rows = []
    for scope in SCOPES:
        for row in rows[scope]:
            if row["index"] != "S2" or any(row["parameter_i"] == p or row["parameter_j"] == p for p in FOCAL):
                selected_rows.append(row)
    write_table(review / "plotted_indices.tsv.gz", selected_rows)
    write_table(review / "focal_s2_source.tsv.gz", [r for r in selected_rows if r["index"] == "S2"])
    diagnostics, band_rows = [], []
    for scope, data in rows.items():
        for output in OUTPUTS:
            for index in ("S1", "ST", "S2"):
                subset = [r for r in data if r["output"] == output and r["index"] == index]
                total = sum(int(r["n_endpoints"]) for r in subset)
                passed = sum(int(r["n_passed"]) for r in subset)
                diagnostics.append(dict(scope=scope, output=output, index=index, total_cells=total,
                                        passed_cells=passed, pass_fraction=passed / total,
                                        negative_scope_summary_cells=sum(float(r["value_mean"]) < 0 for r in subset),
                                        scope_summary_cells=len(subset)))
                for pi, pj in ([(p, "") for p in active] if index != "S2" else list(itertools.combinations(active, 2))):
                    parameter_rows = [r for r in subset if r["parameter_i"] == pi and r["parameter_j"] == pj]
                    for band, lo, hi in BANDS:
                        keep = [r for r in parameter_rows if
                                (lo <= float(r["O2_pct"]) <= hi if band != "middle_1_3" else
                                 lo < float(r["O2_pct"]) < hi)]
                        band_rows.append(dict(scope=scope, output=output, index=index, parameter_i=pi,
                                              parameter_j=pj, band=band, n_oxygen=len(keep),
                                              mean_of_scope_summaries=float(np.mean([float(r["value_mean"]) for r in keep])),
                                              diagnostic_pass_fraction=sum(int(r["n_passed"]) for r in keep) /
                                              sum(int(r["n_endpoints"]) for r in keep)))
    write_table(review / "diagnostic_summary.tsv", diagnostics)
    write_table(review / "oxygen_band_comparison.tsv", band_rows)
    focal_points = [r for r in selected_rows if r["index"] in ("S1", "ST") and r["parameter_i"] in FOCAL
                    and float(r["O2_pct"]) in (0., .5, 2.5, 5.)]
    write_table(review / "focal_parameter_oxygen_points.tsv", focal_points)
    shared_scales = {}
    for index in ("S1", "ST", "S2"):
        keys = metric_rows(index, active, pairs)
        for output in OUTPUTS:
            both = np.concatenate([matrix(lookups[s], keys, index, output, oxygen).ravel() for s in SCOPES])
            shared_scales[index, output] = max(float(np.nanmax(np.abs(both))), 1e-9), bool(np.nanmin(both) < 0)
    for scope in SCOPES:
        breaks = [i for i in range(1, 36) if pairs[i]["focal_parameter"] != pairs[i - 1]["focal_parameter"]]
        for index in ("S1", "ST", "S2"):
            keys = metric_rows(index, active, pairs)
            labels = ([r["display_label"] for r in pairs] if index == "S2" else
                      [p + " [" + groups[p] + "]" for p in active])
            fig, axes = plt.subplots(1, 2, figsize=(15, 11.7 if index == "S2" else 5.8), layout="constrained")
            for k, (ax, output) in enumerate(zip(axes, OUTPUTS)):
                data = matrix(lookups[scope], keys, index, output, oxygen)
                limits = draw(ax, data, oxygen, labels, index, show_labels=k == 0,
                              pair_breaks=breaks if index == "S2" else ())
                ax.set_title(OUT_LABEL[output], fontsize=11, pad=23)
                scale_rows.append(dict(figure=scope + "/" + index, output=output, **limits))
            fig.suptitle(SCOPE_LABEL[scope] + " | " + ("36 focal parameter pairs" if index == "S2" else index)
                         + " | N4096, five repeats", fontsize=12)
            fig.supxlabel("Global: five-repeat mean; neighborhoods: median of 500 five-repeat means. All estimates retained.\n"
                          "Positive estimates: white to purple. Negative estimates: blue, not antagonism. Diagnostics in companion figures/tables.\n"
                          + OXYGEN_NOTE, fontsize=8)
            stem = (root / scope / "figures" / (index.lower() + "_parameter_by_oxygen") if index != "S2" else
                    root / scope / "figures/s2_focal_pairs_by_oxygen")
            save_figure(fig, stem, rendered)
        # All individual scope groups also use shared limits, so PDF-to-PDF comparison is valid.
        draw_grid(root, review, lookups, active, groups, pairs, oxygen,
                  [scope, scope], list(OUTPUTS), scope + "_s1_st_focal_s2",
                  SCOPE_LABEL[scope] + ": first, total and focal second-order sensitivity", rendered, shared_scales)
    draw_grid(root, review, lookups, active, groups, pairs, oxygen,
              ["global", "neighborhood10pct", "global", "neighborhood10pct"],
              [OUTPUTS[0], OUTPUTS[0], OUTPUTS[1], OUTPUTS[1]], "global_neighborhood_s1_st_focal_s2",
              "Global and neighborhood Sobol sensitivity", rendered, shared_scales)
    for output, short in zip(OUTPUTS, ("ploidy", "growth")):
        draw_grid(root, review, lookups, active, groups, pairs, oxygen, list(SCOPES), [output, output],
                  "global_neighborhood_" + short + "_s1_st_focal_s2",
                  "Global and neighborhood comparison: " + OUT_LABEL[output], rendered, shared_scales)
    draw_grid(root, review, lookups, active, groups, pairs, oxygen,
              ["global", "neighborhood10pct", "global", "neighborhood10pct"],
              [OUTPUTS[0], OUTPUTS[0], OUTPUTS[1], OUTPUTS[1]], "global_neighborhood_diagnostic_fraction",
              "Diagnostic coverage for the displayed sensitivity cells", rendered, shared_scales, diagnostic=True)
    for (index, output), (limit, signed) in shared_scales.items():
        scale_rows.append(dict(figure="all_assembled_groups/" + index, output=output,
                               vmin=-limit if signed else 0., vmax=limit, negative_blue=signed))
    write_table(review / "color_scales.tsv", scale_rows)
    after = {scope: sha(path) for scope, path in sources.items()}
    assert after == hashes, "Numerical source tables changed"
    assert len([r for r in selected_rows if r["index"] == "S2"]) == 2 * 2 * 201 * 36
    # The purple positive mapping is exactly the same with/without a negative half.
    signed_cmap, signed_norm = palette(1., True)
    positive_cmap, positive_norm = palette(1., False)
    assert np.allclose([signed_cmap(signed_norm(v)) for v in (0., .5, 1.)],
                       [positive_cmap(positive_norm(v)) for v in (0., .5, 1.)])
    receipt = dict(status="complete", created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                   method="local_table_only_redraw", N=4096, repeats=5, oxygen_points=201,
                   n_neighborhoods=500, focal_parameters=list(FOCAL), unique_focal_pairs=36,
                   input_sha256=hashes, numerical_sources_unchanged=True,
                   model_evaluations=0, sensitivity_estimations=0,
                   palette=dict(positive=POSITIVE, negative=NEGATIVE, zero="white"),
                   signed_scale="symmetric absolute range; negative estimates retained",
                   oxygen_axis=dict(scale="symlog", base=10, linthresh=OXYGEN_LINTHRESH,
                                    ticks=list(OXYGEN_TICKS), boundaries=[1., 3.], zero_retained=True),
                   shared_scales="same metric and output shared across global/neighborhood assembled groups",
                   original_figures=original_files, script_sha256=sha(__file__),
                   runtime=dict(python=__import__("sys").version, numpy=np.__version__, matplotlib=matplotlib.__version__),
                   rendered=[dict(path=str(p.relative_to(root)), sha256=sha(p), bytes=p.stat().st_size) for p in rendered])
    (review / "redraw_manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(dict(status="complete", n_rendered_files=len(rendered), unique_pairs=36,
                          sources_unchanged=True, model_evaluations=0, output=str(review)), indent=2))


if __name__ == "__main__":
    main()
