#!/usr/bin/env python3
"""Local annotated display overlay; reuses N4096 indices without model evaluation."""
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import csv
import datetime as dt
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import shutil

import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm, to_rgb
from matplotlib.patches import Patch

import redraw_full_n4096_local as h

AUDIT = "row_annotations_20261019"
GROUPS = ("High O2", "Low O2", "No resolved O2 preference")
COLORS = dict(zip(GROUPS, ("#B2182B", "#2166AC", "#8C8C8C")))
VERSION = "paired-whole-curve-bootstrap-v1"


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def save_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2) + "\n")
    tmp.replace(path)


def identifier(index, pi, pj):
    return index + ":" + pi + ":" + pj


def definitions(active):
    return [(ix, p, "") for ix in ("S1", "ST") for p in active] + [
        ("S2", a, b) for a, b in itertools.combinations(active, 2)]


def weights(oxygen, lo, hi):
    keep = np.flatnonzero((oxygen >= lo) & (oxygen <= hi))
    w = np.zeros(len(oxygen))
    delta = np.diff(oxygen[keep])
    w[keep[:-1]] += delta / 2 / (hi - lo)
    w[keep[1:]] += delta / 2 / (hi - lo)
    assert np.isclose(w.sum(), 1.)
    return w


def target_specs(root):
    audit = json.loads((root / "axis_redraw_20261019/completion.json").read_text())
    return [r for r in audit["outputs"] if "/full/" not in r["path"]
            and "standalone_N4096" not in r["path"]]


def prepare(root, reps, seed):
    out = root / AUDIT
    pp = out / "plan.json"
    if pp.exists():
        plan = json.loads(pp.read_text())
        assert (plan["bootstrap_reps"], plan["bootstrap_seed"]) == (reps, seed)
        return plan
    out.mkdir(parents=True, exist_ok=True)
    old = json.loads((root / "axis_redraw_20261019/plan.json").read_text())
    targets = target_specs(root)
    assert len(targets) == 30
    for r in targets:
        p = root / r["path"]
        assert h.sha(p) == r["sha256"], p
        backup = out / "original_figures" / r["path"]
        backup.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, backup)
        r["backup"] = str(backup.relative_to(root))
    annroot = root.parents[1]
    annotations = [annroot / n for n in ("figure4_parameter_groups_source.tsv",
                    "figure4_parameter_group_palette_source.tsv", "figure4b_o2_classification_source.tsv")]
    source_paths = [root / s / "summaries/indices_scope_summary.tsv.gz" for s in h.SCOPES]
    source_paths += [root / "global/standalone_N4096/indices_repeats.tsv.gz"]
    source_paths += [root / "neighborhood10pct/full" / ("seed" + str(i)) /
                     "indices_summary.tsv.gz" for i in range(1, 501)]
    source_paths += annotations
    plan = dict(created_utc=now(), bootstrap_reps=reps, bootstrap_seed=seed,
                sources={str(p): h.sha(p) for p in source_paths}, targets=targets,
                preserved_numerical_files=old["numerical_input_sha256"],
                preserved_images=[r for r in json.loads((root / "axis_redraw_20261019/completion.json").read_text())["outputs"]
                                  if r["path"] not in {t["path"] for t in targets}],
                preserved_responses=[r for r in h.read_table(root / "collection_manifest.tsv")
                                     if "/figures/responses/" in r["path"] and r["path"].endswith((".pdf", ".png"))],
                model_evaluations=0, sensitivity_estimations=0,
                annotation_sources=[str(p) for p in annotations],
                comparison_order_anchor="neighborhood10pct dominant_mean_ploidy",
                peak_threshold=None, classification_output=h.OUTPUTS[0])
    save_json(pp, plan)
    for p in annotations:
        shutil.copy2(p, out / p.name)
    for name in ("README.md", "redraw_manifest.json", "pdf_validation.json", "interpretation_zh.md"):
        src = root / "figure_comparison_20261019" / name
        dest = out / "documents_before" / name
        dest.parent.mkdir(exist_ok=True)
        shutil.copy2(src, dest)
    return plan


def read_seed(payload):
    path, defs = payload
    lookup = {tuple(row): i for i, row in enumerate(defs)}
    arr = np.full((len(defs), 201), np.nan)
    count = 0
    with gzip.open(path, "rt", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["N"] != "4096" or r["output"] != h.OUTPUTS[0]:
                continue
            arr[lookup[r["index"], r["parameter_i"], r["parameter_j"]],
                int(round(float(r["O2_pct"]) / .025))] = float(r["value_mean"])
            count += 1
    assert count == len(defs) * 201 and np.isfinite(arr).all(), path
    return arr


def curves(root, defs, lookups, oxygen, workers, plan):
    fingerprint = hashlib.sha256(json.dumps(plan["sources"], sort_keys=True).encode()).hexdigest()
    cache = Path("/private/tmp") / ("n4096_annotation_curves_" + fingerprint[:16] + ".npz")
    if cache.exists():
        with np.load(cache) as loaded:
            return {s: loaded[s] for s in h.SCOPES}
    paths = [root / "neighborhood10pct/full" / ("seed" + str(i)) / "indices_summary.tsv.gz"
             for i in range(1, 501)]
    result = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, a in enumerate(pool.map(read_seed, [(p, defs) for p in paths]), 1):
            result.append(a)
            if i % 50 == 0:
                print("READ seeds %d/500" % i, flush=True)
    neighborhood = np.stack(result)
    repeats = np.full((5, len(defs), len(oxygen)), np.nan)
    indices = {v: i for i, v in enumerate(defs)}
    count = 0
    with gzip.open(root / "global/standalone_N4096/indices_repeats.tsv.gz", "rt", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["N"] == "4096" and r["output"] == h.OUTPUTS[0]:
                repeats[int(r["repeat"]) - 1, indices[r["index"], r["parameter_i"], r["parameter_j"]],
                        int(round(float(r["O2_pct"]) / .025))] = float(r["value"])
                count += 1
    assert count == repeats.size and np.isfinite(repeats).all()
    for s, arr, reducer in (("global", repeats, np.mean), ("neighborhood10pct", neighborhood, np.median)):
        expected = np.stack([h.matrix(lookups[s], [(pi, pj)], ix, h.OUTPUTS[0], oxygen)[0]
                             for ix, pi, pj in defs])
        assert np.allclose(reducer(arr, axis=0), expected, rtol=1e-10, atol=1e-12), s
    arrays = {"global": repeats, "neighborhood10pct": neighborhood}
    np.savez(cache, **arrays)
    return arrays


def bootstrap_one(scope, curve, draws, low_w, high_w):
    reducer = np.mean if scope == "global" else np.median
    observed = reducer(curve, axis=0)
    low, high = float(observed @ low_w), float(observed @ high_w)
    used = np.flatnonzero(low_w + high_w)
    contrast_w = (low_w - high_w)[used]
    selected = np.ascontiguousarray(curve[:, used].T)
    delta = np.empty(len(draws))
    for start in range(0, len(draws), 16):
        sampled = selected[:, draws[start:start + 16]]
        summarized = (np.mean(sampled, axis=-1) if scope == "global" else
                      np.median(sampled, axis=-1, overwrite_input=True))
        delta[start:start + len(summarized.T)] = contrast_w @ summarized
    tail = min(np.count_nonzero(delta <= 0), np.count_nonzero(delta >= 0))
    return dict(low_score=low, high_score=high, low_minus_high=low - high,
                global_peak=float(observed.max()), bootstrap_sign_p=min(1., 2 * (tail + 1) / (len(delta) + 1)),
                bootstrap_delta_q025=float(np.quantile(delta, .025)),
                bootstrap_delta_q975=float(np.quantile(delta, .975)))


def bh(values):
    p = np.asarray(values)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.empty_like(p)
    q[order] = np.minimum(1., np.minimum.accumulate(ranked[::-1])[::-1])
    return q


def classify(root, plan, active, oxygen, lookups, workers):
    out = root / AUDIT
    defs = definitions(active)
    source_id = hashlib.sha256(json.dumps(plan["sources"], sort_keys=True).encode()).hexdigest()
    token = dict(version=VERSION, source_id=source_id, reps=plan["bootstrap_reps"], seed=plan["bootstrap_seed"])
    complete = out / "classification.json"
    if complete.exists():
        saved = json.loads(complete.read_text())
        if saved["config"] == token:
            return saved["rows"]
    arrays = curves(root, defs, lookups, oxygen, workers, plan)
    low_w, high_w = weights(oxygen, 0, 1), weights(oxygen, 3, 5)
    rng = np.random.default_rng(plan["bootstrap_seed"])
    samples = {s: rng.integers(0, len(arrays[s]), size=(plan["bootstrap_reps"], len(arrays[s])))
               for s in h.SCOPES}
    results = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {}
        for s in h.SCOPES:
            for i, (ix, pi, pj) in enumerate(defs):
                key = s + ":" + identifier(ix, pi, pj)
                checkpoint = out / "bootstrap_rows" / (hashlib.sha256(key.encode()).hexdigest()[:20] + ".json")
                if checkpoint.exists():
                    record = json.loads(checkpoint.read_text())
                    if record["config"] == token:
                        results[key] = record["scores"]
                        continue
                future = pool.submit(bootstrap_one, s, arrays[s][:, i], samples[s], low_w, high_w)
                futures[future] = key, checkpoint
        for n, future in enumerate(as_completed(futures), 1):
            key, checkpoint = futures[future]
            results[key] = future.result()
            save_json(checkpoint, dict(config=token, scores=results[key]))
            if n % 5 == 0 or n == len(futures):
                save_json(out / "progress.json", dict(phase="classification", completed=len(results), total=238, updated_utc=now()))
                print("CLASSIFY %d/238" % len(results), flush=True)
    rows = []
    for s in h.SCOPES:
        for ix in ("S1", "ST", "S2"):
            selected = [d for d in defs if d[0] == ix]
            qvals = bh([results[s + ":" + identifier(*d)]["bootstrap_sign_p"] for d in selected])
            for (_, pi, pj), q in zip(selected, qvals):
                scores = results[s + ":" + identifier(ix, pi, pj)]
                # Negative-only window estimates are sampling artifacts, not positive sensitivity.
                positive_window = max(scores["low_score"], scores["high_score"]) > 0
                group = (GROUPS[0] if scores["low_minus_high"] < 0 else GROUPS[1]) if (
                    positive_window and q < .05 and scores["low_minus_high"] != 0) else GROUPS[2]
                rows.append(dict(scope=s, index=ix, parameter_i=pi, parameter_j=pj,
                                 sensitivity_group=group, bh_q=float(q), **scores,
                                 positive_window=positive_window, peak_threshold="none",
                                 classification_output=h.OUTPUTS[0], n_units=len(arrays[s]),
                                 bootstrap_unit="whole scramble-repeat curve" if s == "global" else "whole fitted-seed curve (five-repeat mean)",
                                 bootstrap_reps=plan["bootstrap_reps"], bootstrap_seed=plan["bootstrap_seed"],
                                 bh_family_size=len(selected),
                                 low_diagnostic_fraction=float(h.matrix(lookups[s], [(pi, pj)], ix, h.OUTPUTS[0], oxygen,
                                                                       "pass_fraction_all_endpoints")[0] @ low_w),
                                 high_diagnostic_fraction=float(h.matrix(lookups[s], [(pi, pj)], ix, h.OUTPUTS[0], oxygen,
                                                                        "pass_fraction_all_endpoints")[0] @ high_w)))
    h.write_table(out / "classification.tsv", rows)
    save_json(complete, dict(config=token, rows=rows))
    return rows


def ordering(keys, scope, index, lookup):
    def key(pair):
        r = lookup[scope, index, *pair]
        group = GROUPS.index(r["sensitivity_group"])
        score = (-r["high_score"] if group == 0 else r["low_score"] if group == 1 else -r["global_peak"])
        return group, score, keys.index(pair)
    return sorted(keys, key=key)


def independent_bootstrap_check():
    """Compare the optimized bootstrap against full-curve, draw-by-draw integration."""
    rng = np.random.default_rng(937)
    oxygen = np.linspace(0, 5, 201)
    curve = rng.normal(size=(7, 201))
    draws = rng.integers(0, 7, size=(23, 7))
    low_w, high_w = weights(oxygen, 0, 1), weights(oxygen, 3, 5)
    def contrast(a):
        low = np.trapezoid(a[oxygen <= 1], oxygen[oxygen <= 1])
        high = np.trapezoid(a[oxygen >= 3], oxygen[oxygen >= 3]) / 2
        return float(low), float(high)
    for scope, reducer in (("global", np.mean), ("neighborhood10pct", np.median)):
        got = bootstrap_one(scope, curve, draws, low_w, high_w)
        low, high = contrast(reducer(curve, axis=0))
        deltas = np.array([np.subtract(*contrast(reducer(curve[d], axis=0))) for d in draws])
        assert np.allclose([got["low_score"], got["high_score"]], [low, high])
        assert np.allclose([got["bootstrap_delta_q025"], got["bootstrap_delta_q975"]], np.quantile(deltas, [.025, .975]))
        expected_p = min(1., 2 * (min(np.count_nonzero(deltas <= 0), np.count_nonzero(deltas >= 0)) + 1) / (len(draws) + 1))
        assert np.isclose(got["bootstrap_sign_p"], expected_p)
    keys = [(p, "") for p in ("A", "B", "C", "D", "E")]
    rows = {("x", "ST", p, ""): dict(sensitivity_group=g, high_score=hscore, low_score=lscore, global_peak=9.)
            for p, g, hscore, lscore in [("A", GROUPS[0], .2, .1), ("B", GROUPS[1], .1, .5),
                                       ("C", GROUPS[2], .3, .3), ("D", GROUPS[0], .8, .1), ("E", GROUPS[1], .1, .2)]}
    assert [p for p, _ in ordering(keys, "x", "ST", rows)] == ["D", "A", "E", "B", "C"]
    return dict(status="passed", full_curve_oracle_scopes=2, draws_each=23, physical_window_integration=True,
                pointwise_median_before_integration=True, requested_sort_directions=True)


def legend(fig, palettes):
    process = [Patch(facecolor=r["color"], label=r["display_label"]) for r in palettes]
    height = fig.get_figheight()
    fig.legend(handles=process, loc="lower center", bbox_to_anchor=(.5, .58 / height), ncol=5, frameon=False,
               title="1  Process", fontsize=8, title_fontsize=9)
    fig.legend(handles=[Patch(facecolor=COLORS[g], label=g.replace(" O2", " O₂")) for g in GROUPS],
               loc="lower center", bbox_to_anchor=(.5, .06 / height), ncol=3, frameon=False,
               title="2  Figure 4B O₂ correlation     3  Sobol O₂ sensitivity (ploidy-based)", fontsize=8, title_fontsize=9)


def annotation_colors(keys, members_list, scope, index, classes, process, correlation, pcolors):
    colors = []
    for (pi, pj), members in zip(keys, members_list):
        colors.append([to_rgb(pcolors[process[p]]) for p in members] +
                      [to_rgb(COLORS[correlation[p]]) for p in members] +
                      [to_rgb(COLORS[classes[scope, index, pi, pj]["sensitivity_group"]])] * 2)
    return np.asarray(colors)


def draw_panel(fig, slot, data, oxygen, keys, labels, members_list, scope, index, classes, annotation,
               title, limit, signed, diagnostic=False, legacy=False, show_labels=True):
    slot_height = slot.get_position(fig).height * fig.get_figheight()
    gs = slot.subgridspec(3, 2, width_ratios=[.72, 6],
                         height_ratios=[slot_height - .59, .47, .12], hspace=0, wspace=.07)
    aa, ax, cbax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[2, 1])
    rgb = annotation_colors(keys, members_list, scope, index, classes, *annotation)
    aa.imshow(rgb, extent=(0, 3, len(keys) - .5, -.5), aspect="auto", interpolation="nearest")
    aa.set_xticks([.5, 1.5, 2.5], ["1", "2", "3"])
    aa.xaxis.tick_top()
    aa.set_yticks(np.arange(len(keys)), labels if show_labels else [])
    aa.tick_params(axis="both", length=0, labelsize=8)
    for x in (1, 2):
        aa.axvline(x, color="white", lw=1)
    for spine in aa.spines.values():
        spine.set_color("#666666")
        spine.set_linewidth(.5)
    if diagnostic:
        cmap = (h.plt.get_cmap("viridis") if legacy else
                h.LinearSegmentedColormap.from_list("diagnostic", ["#f2f2f2", "#333333"]))
        norm = Normalize(0, 1)
    elif legacy:
        lo, hi = min(0., float(np.nanmin(data))), max(.001, float(np.nanmax(data)))
        cmap, _ = h.palette(max(abs(lo), hi), lo < 0)
        norm = TwoSlopeNorm(0, lo, hi) if lo < 0 else Normalize(0, hi)
    else:
        cmap, norm = h.palette(limit, signed)
    mesh = h.oxygen_mesh(ax, data, oxygen, cmap=cmap, norm=norm)
    assert np.array_equal(np.asarray(mesh.get_array()), data, equal_nan=True)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=8, length=2)
    ax.set_title(title, loc="left", fontsize=10, pad=23)
    for spine in ax.spines.values():
        spine.set_linewidth(.5)
    cb = fig.colorbar(mesh, cax=cbax, orientation="horizontal")
    ticks = ([0, .5, 1] if diagnostic else [norm.vmin, 0, norm.vmax] if norm.vmin < 0 else [0, norm.vmax / 2, norm.vmax])
    cb.set_ticks(ticks, labels=[format(x, ".3g") for x in ticks])
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("Diagnostic pass fraction" if diagnostic else index + " estimate", fontsize=8)
    return dict(vmin=float(norm.vmin), vmax=float(norm.vmax),
                annotation_axes=aa, heatmap_axes=ax, colorbar_axes=cbax, rgb=rgb)


def render(root, classes_rows, active, pairs, oxygen, lookups, only=None):
    audit, review = root / AUDIT, root / "figure_comparison_20261019"
    classes = {(r["scope"], r["index"], r["parameter_i"], r["parameter_j"]): r for r in classes_rows}
    process = {r["parameter"]: r["parameter_group"] for r in h.read_table(audit / "figure4_parameter_groups_source.tsv")}
    palettes = h.read_table(audit / "figure4_parameter_group_palette_source.tsv")
    pcolors = {r["parameter_group"]: r["color"] for r in palettes}
    correlation = {r["parameter"]: r["o2_association_group"].replace("O2-independent", GROUPS[2])
                   for r in h.read_table(audit / "figure4b_o2_classification_source.tsv")}
    assert all(p in process and p in correlation for p in active)
    annotation = process, correlation, pcolors
    focal_members = {(r["parameter_i"], r["parameter_j"]): (r["focal_parameter"], r["partner"]) for r in pairs}
    keysets = {"S1": [(p, "") for p in active], "ST": [(p, "") for p in active],
               "S2": [(r["parameter_i"], r["parameter_j"]) for r in pairs],
               "S2_all": list(itertools.combinations(active, 2))}
    scales = {}
    for ix in ("S1", "ST", "S2"):
        for output in h.OUTPUTS:
            data = np.concatenate([h.matrix(lookups[s], keysets[ix], ix, output, oxygen).ravel() for s in h.SCOPES])
            scales[ix, output] = max(float(np.max(np.abs(data))), 1e-9), bool(np.min(data) < 0)
    specs = []
    for scope in h.SCOPES:
        for ix in ("S1", "ST", "S2", "S2_all"):
            name = (ix.lower() + "_parameter_by_oxygen" if ix in ("S1", "ST") else
                    "s2_focal_pairs_by_oxygen" if ix == "S2" else "s2_pair_by_oxygen")
            specs.append((root / scope / "figures" / name, [(ix, scope, o) for o in h.OUTPUTS], 1, False, False))
        if scope == "neighborhood10pct":
            specs.append((root / scope / "figures/s2_diagnostic_fraction",
                          [("S2_all", scope, o) for o in h.OUTPUTS], 1, False, True))
        specs.append((review / (scope + "_s1_st_focal_s2"),
                      [(ix, scope, o) for ix in ("S1", "ST", "S2") for o in h.OUTPUTS], 3, False, False))
    for name, columns, diag in (
        ("global_neighborhood_s1_st_focal_s2", [(s, o) for o in h.OUTPUTS for s in h.SCOPES], False),
        ("global_neighborhood_ploidy_s1_st_focal_s2", [(s, h.OUTPUTS[0]) for s in h.SCOPES], False),
        ("global_neighborhood_growth_s1_st_focal_s2", [(s, h.OUTPUTS[1]) for s in h.SCOPES], False),
        ("global_neighborhood_diagnostic_fraction", [(s, o) for o in h.OUTPUTS for s in h.SCOPES], True)):
        specs.append((review / name, [(ix, s, o) for ix in ("S1", "ST", "S2") for s, o in columns], 3, True, diag))
    all_records, order_rows = [], []
    for scope in h.SCOPES:
        for kind, keys in keysets.items():
            ix = kind.split("_")[0]
            for i, (pi, pj) in enumerate(ordering(keys, scope, ix, classes), 1):
                rr = classes[scope, ix, pi, pj]
                order_rows.append(dict(scope=scope, kind=kind, display_order=i, parameter_i=pi, parameter_j=pj,
                                       sensitivity_group=rr["sensitivity_group"], high_score=rr["high_score"], low_score=rr["low_score"]))
    h.write_table(audit / "row_order.tsv", order_rows)
    h.plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42,
                          "savefig.facecolor": "white"})
    for stem, panels, nrows, comparison, diagnostic in specs:
        if only and stem.name not in only:
            continue
        cols = len(panels) // nrows
        longest = max(len(keysets[k]) for k, _, _ in panels)
        height = 25 if longest == 91 else 19.5 if nrows == 3 else 13.5 if longest == 36 else 8
        fig = h.plt.figure(figsize=(19 if cols == 2 else 33, height))
        grid = fig.add_gridspec(nrows, cols, left=.13 if cols == 2 else .08, right=.985,
                                top=.92 if nrows == 1 else .935, bottom=2.35 / height,
                                hspace=.48 if nrows == 3 else .1, wspace=.26 if cols == 2 else .17,
                                height_ratios=[14, 14, 36] if nrows == 3 else None)
        records, axes = [], []
        for pidx, (kind, scope, output) in enumerate(panels):
            ix = kind.split("_")[0]
            anchor = "neighborhood10pct" if comparison else scope
            keys = ordering(keysets[kind], anchor, ix, classes)
            members_list = [focal_members[pi, pj] if kind == "S2" else (pi, pj or pi) for pi, pj in keys]
            labels = [a + (" × " + b if pj else "") for (pi, pj), (a, b) in zip(keys, members_list)]
            data = h.matrix(lookups[scope], keys, ix, output, oxygen,
                            "pass_fraction_all_endpoints" if diagnostic else "value_mean")
            legacy = kind == "S2_all"
            limit, signed = scales[ix, output] if nrows == 3 else (max(float(np.max(np.abs(data))), 1e-9), bool(np.min(data) < 0))
            panel_id = chr(ord("A") + pidx)
            detail = "S2: 91 pairs" if kind == "S2_all" else "S2: 36 focal pairs" if ix == "S2" else ix
            res = draw_panel(fig, grid[pidx // cols, pidx % cols], data, oxygen, keys, labels, members_list, scope, ix,
                             classes, annotation, panel_id + "  " + h.SCOPE_LABEL[scope] + " | " + detail + "\n" + h.OUT_LABEL[output],
                             limit, signed, diagnostic, legacy, show_labels=(pidx % cols == 0))
            axes.append(res)
            anchor_groups = [classes[anchor, ix, *k]["sensitivity_group"] for k in keys]
            for y in range(1, len(keys)):
                if anchor_groups[y] != anchor_groups[y - 1]:
                    for axis in (res["annotation_axes"], res["heatmap_axes"]):
                        axis.axhline(y - .5, color="#555555", lw=.65)
            records.append(dict(panel=panel_id, scope=scope, output=output, index=ix, kind=kind,
                                ordering_scope=anchor, ordering_output=h.OUTPUTS[0],
                                keys=[list(k) for k in keys], labels=labels,
                                displayed_members=[list(m) for m in members_list],
                                sensitivity_groups=[classes[scope, ix, *k]["sensitivity_group"] for k in keys],
                                process_groups=[[process[p] for p in members] for members in members_list],
                                correlation_groups=[[correlation[p] for p in members] for members in members_list],
                                annotation_rgb=res["rgb"].tolist(), vmin=res["vmin"], vmax=res["vmax"],
                                displayed_data_sha256=hashlib.sha256(np.ascontiguousarray(data).tobytes()).hexdigest()))
        fig.suptitle("Sobol sensitivity | N4096, five repeats, 201 oxygen points\n" +
                     ("Cross-scope row order: neighborhood ploidy; each panel retains its own sensitivity annotation." if comparison else
                      "Row order and sensitivity annotations: ploidy-based; growth follows the same rows."), fontsize=12, y=.988)
        caption = fig.text(.5, 1.19 / height,
                 "High group: high-O₂ score descending; Low group: low-O₂ score ascending; unresolved last. "
                 "No 0.3 sensitivity gate.\n"
                 "S2 annotations 1/2: left/right half = first/second named parameter. "
                 "Gray = no resolved O₂ preference, not no output effect.\n"
                 "Cells: global five-repeat mean; neighborhood median of 500 five-repeat means. No diagnostic exclusions.\n"
                 + h.OXYGEN_NOTE, ha="center", fontsize=8)
        legend(fig, palettes)
        fig.canvas.draw()
        for record, res in zip(records, axes):
            aa, ax = res["annotation_axes"], res["heatmap_axes"]
            yy = np.arange(len(record["keys"]))
            left = aa.transData.transform(np.c_[np.ones(len(yy)), yy])[:, 1]
            right = ax.transData.transform(np.c_[np.ones(len(yy)), yy])[:, 1]
            record["row_alignment_error_points"] = float(np.max(np.abs(left - right)) * 72 / fig.dpi)
            assert record["row_alignment_error_points"] < 1e-8
            renderer = fig.canvas.get_renderer()
            xlabel_gap = ax.xaxis.label.get_window_extent(renderer).y0 - res["colorbar_axes"].get_window_extent(renderer).y1
            footer_gap = res["colorbar_axes"].xaxis.label.get_window_extent(renderer).y0 - caption.get_window_extent(renderer).y1
            record["xlabel_colorbar_gap_points"] = float(xlabel_gap * 72 / fig.dpi)
            record["colorbar_footer_gap_points"] = float(footer_gap * 72 / fig.dpi)
            assert record["xlabel_colorbar_gap_points"] > 2, record
            assert record["colorbar_footer_gap_points"] > 2, record
        legend_top, legend_bottom = [lg.get_window_extent(fig.canvas.get_renderer()) for lg in fig.legends]
        assert caption.get_window_extent(fig.canvas.get_renderer()).y0 > legend_top.y1 + 2 * fig.dpi / 72
        assert legend_top.y0 > legend_bottom.y1 + 2 * fig.dpi / 72
        # Manual layout is stable for PDF and PNG; use the common atomic exporter.
        h.save_figure(fig, stem, [])
        rec = dict(stem=str(stem.relative_to(root)), panels=records,
                   files=[dict(path=str(stem.with_suffix(s).relative_to(root)), sha256=h.sha(stem.with_suffix(s))) for s in (".pdf", ".png")])
        save_json(audit / "figure_receipts" / (stem.name + ("_" + panels[0][1] if nrows == 1 else "") + ".json"), rec)
        all_records.append(rec)
        print("RENDER " + str(stem.relative_to(root)), flush=True)
    return all_records


def validate(root, plan, active, pairs, oxygen, lookups, classes_rows):
    import pdfplumber
    audit = root / AUDIT
    independent_check = independent_bootstrap_check()
    records = [json.loads(p.read_text()) for p in (audit / "figure_receipts").glob("*.json")]
    assert len(records) == 15, len(records)
    assert {f["path"] for r in records for f in r["files"]} == {t["path"] for t in plan["targets"]}
    cl = {(r["scope"], r["index"], r["parameter_i"], r["parameter_j"]): r for r in classes_rows}
    process = {r["parameter"]: r["parameter_group"] for r in h.read_table(audit / "figure4_parameter_groups_source.tsv")}
    pcolors = {r["parameter_group"]: r["color"] for r in h.read_table(audit / "figure4_parameter_group_palette_source.tsv")}
    correlation = {r["parameter"]: r["o2_association_group"].replace("O2-independent", GROUPS[2])
                   for r in h.read_table(audit / "figure4b_o2_classification_source.tsv")}
    for scope in h.SCOPES:
        for ix in ("S1", "ST", "S2"):
            subset = [r for r in classes_rows if r["scope"] == scope and r["index"] == ix]
            assert len(subset) == (91 if ix == "S2" else 14)
            assert np.allclose(bh([r["bootstrap_sign_p"] for r in subset]), [r["bh_q"] for r in subset])
    pdf_checks = []
    for r in records:
        for f in r["files"]:
            assert h.sha(root / f["path"]) == f["sha256"]
        for panel in r["panels"]:
            keys = [tuple(k) for k in panel["keys"]]
            assert len(keys) == len(set(keys)) == (91 if panel["kind"] == "S2_all" else 36 if panel["index"] == "S2" else 14)
            assert keys == ordering(keys, panel["ordering_scope"], panel["index"], cl)
            data = h.matrix(lookups[panel["scope"]], keys, panel["index"], panel["output"], oxygen,
                            "pass_fraction_all_endpoints" if "diagnostic" in r["stem"] else "value_mean")
            assert hashlib.sha256(np.ascontiguousarray(data).tobytes()).hexdigest() == panel["displayed_data_sha256"]
            assert panel["sensitivity_groups"] == [cl[panel["scope"], panel["index"], *k]["sensitivity_group"] for k in keys]
            members = panel["displayed_members"]
            expected_rgb = annotation_colors(keys, members, panel["scope"], panel["index"], cl,
                                             process, correlation, pcolors)
            assert np.array_equal(expected_rgb, panel["annotation_rgb"])
            assert panel["row_alignment_error_points"] < 1e-8
            assert panel["xlabel_colorbar_gap_points"] > 2
            assert panel["colorbar_footer_gap_points"] > 2
            for (pi, pj), names, label in zip(keys, members, panel["labels"]):
                assert (set(names) == {pi, pj} if pj else names == [pi, pi])
                assert label == names[0] + (" × " + names[1] if pj else "")
            if "global_neighborhood" in r["stem"]:
                assert panel["ordering_scope"] == "neighborhood10pct"
        for i in range(0, len(r["panels"]), 2 if len(r["panels"]) <= 6 else 4):
            matching = [p for p in r["panels"] if p["kind"] == r["panels"][i]["kind"]]
            assert all(p["keys"] == matching[0]["keys"] for p in matching)
        pdf = root / (r["stem"] + ".pdf")
        with pdfplumber.open(pdf) as doc:
            assert len(doc.pages) == 1
            page = doc.pages[0]
            text = page.extract_text() or ""
            flat = " ".join(text.split())
            for term in ("Process", "correlation", "sensitivity", "0.025", "Mid", "High", "Low"):
                assert term in flat, (pdf, term)
            for p in r["panels"]:
                for label in p["labels"]:
                    assert label in flat, (pdf, label)
            outside = [c for c in page.chars if c["text"].strip() and
                       (c["x0"] < -.5 or c["x1"] > page.width + .5 or c["top"] < -.5 or c["bottom"] > page.height + .5)]
            assert not outside, (pdf, outside[:2])
            pdf_checks.append(dict(path=str(pdf.relative_to(root)), pages=1, outside_page=0))
    for path, sha in plan["sources"].items():
        assert h.sha(path) == sha, path
    for path, sha in plan["preserved_numerical_files"].items():
        assert h.sha(root / path) == sha, path
    for r in plan["preserved_images"]:
        assert h.sha(root / r["path"]) == r["sha256"], r["path"]
    assert len(plan["preserved_responses"]) == 108
    for r in plan["preserved_responses"]:
        assert h.sha(root / r["path"]) == r["sha256"], r["path"]
    for r in plan["targets"]:
        assert h.sha(root / r["backup"]) == r["sha256"]
    save_json(audit / "validation.json", dict(status="passed", checked_utc=now(), pdfs=pdf_checks,
              independent_bootstrap_and_sort_check=independent_check,
              classification_rows=len(classes_rows), bootstrap_reps=plan["bootstrap_reps"],
              displayed_values_verified=True, pair_completeness_verified=True, annotation_and_order_verified=True,
              unchanged_numerical_files=len(plan["preserved_numerical_files"]),
              unchanged_response_images=len(plan["preserved_responses"]),
              unchanged_other_oxygen_images=len(plan["preserved_images"]), original_backups_verified=True))
    completion = dict(status="complete", completed_utc=now(), figures=15, files=30,
                      model_evaluations=0, sensitivity_estimations=0,
                      script_sha256=h.sha(__file__), helper_sha256=h.sha(h.__file__),
                      figures_and_annotations=records, validation="validation.json")
    save_json(audit / "completion.json", completion)
    print("VALIDATED 15 figures; 30 files; 2107 numerical files unchanged", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--phase", choices=("all", "classify", "render", "validate"), default="all")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--bootstrap-reps", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=5826)
    parser.add_argument("--only", action="append")
    args = parser.parse_args()
    root = args.root.resolve()
    plan = prepare(root, args.bootstrap_reps, args.seed)
    _, _, _, active, _, pairs, oxygen, lookups = h.load_sources(root)
    rows = classify(root, plan, active, oxygen, lookups, args.workers)
    if args.phase in ("all", "render"):
        render(root, rows, active, pairs, oxygen, lookups, set(args.only) if args.only else None)
    if args.phase in ("all", "validate") and not args.only:
        validate(root, plan, active, pairs, oxygen, lookups, rows)


if __name__ == "__main__":
    main()
