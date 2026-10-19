#!/usr/bin/env python3
"""Resumable local redraw of every current N4096 oxygen-axis figure."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
import numpy as np

import redraw_full_n4096_local as h


AUDIT = "axis_redraw_20261019"
SCRIPT = Path(__file__).resolve()
HELPER = SCRIPT.with_name("redraw_full_n4096_local.py")


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def save_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(data, indent=2) + "\n")
    temp.replace(path)


def settings():
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42,
                         "ps.fonttype": 42, "savefig.facecolor": "white"})


def receipt_path(root, stem):
    relative = str(stem.relative_to(root))
    name = hashlib.sha256(relative.encode()).hexdigest()[:24]
    return root / AUDIT / "receipts" / (name + ".json")


def renderer():
    return dict(batch_sha256=h.sha(SCRIPT), helper_sha256=h.sha(HELPER),
                scale="symlog", base=10, linthresh=.025, linscale=1,
                oxygen_ticks=list(h.OXYGEN_TICKS), boundaries=[1., 3.],
                pdf_layout_frozen_before_export=True,
                band_label_offset_points=3,
                zero_retained=True, cell_edges="physical oxygen midpoints, clipped to 0 and 5",
                positive_palette=h.POSITIVE, model_evaluations=0, sensitivity_estimations=0)


def prepare(root, refresh=False):
    directory = root / AUDIT
    directory.mkdir(parents=True, exist_ok=True)
    pp = directory / "plan.json"
    if pp.exists():
        plan = json.loads(pp.read_text())
        if refresh:
            if plan["renderer"] != renderer() and (directory / "completion.json").exists():
                old_id = plan["renderer"]["helper_sha256"][:12]
                history = directory / "iterations" / ("completed_" + old_id)
                history.mkdir(parents=True, exist_ok=True)
                shutil.copy2(pp, history / "plan.json")
                (directory / "completion.json").replace(history / "completion.json")
            plan["renderer"] = renderer()
            save_json(pp, plan)
        elif plan["renderer"] != renderer():
            raise ValueError("Renderer changed; explicitly refresh the plan before resuming")
        return plan
    existing = sorted(p for p in root.rglob("*.pdf") if AUDIT not in p.parts
                      and "original_figures" not in p.parts and "responses" not in p.parts)
    assert len(existing) == 1521, ("Unexpected current oxygen-figure inventory", len(existing))
    backups = []
    for pdf in existing:
        for suffix in (".pdf", ".png"):
            source = pdf.with_suffix(suffix)
            target = directory / "original_figures" / source.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                try:
                    os.link(source, target)
                except OSError:
                    shutil.copy2(source, target)
            assert h.sha(source) == h.sha(target)
            backups.append(dict(path=str(source.relative_to(root)), original_sha256=h.sha(source),
                                backup=str(target.relative_to(root))))
    sources = [root / s / "summaries/indices_scope_summary.tsv.gz" for s in h.SCOPES]
    datasets = []
    paths = [(root / "global/full/global", "global_case", "Global full bounds; case global"),
             (root / "global/standalone_N4096", "global_standalone", "Global full bounds; standalone N4096")]
    paths += [(root / "neighborhood10pct/full" / ("seed" + str(i)), "seed" + str(i),
               "Neighborhood seed" + str(i) + "; mean of five repeats") for i in range(1, 501)]
    for path, case, title in paths:
        source = path / "indices_summary.tsv.gz"
        assert source.exists()
        sources.append(source)
        datasets.append(dict(case=case, title=title, source=str(source.relative_to(root)),
                             source_sha256=h.sha(source), directory=str((path / "figures").relative_to(root))))
    # Verify all numerical TSV/NPZ inputs, including derived comparison sources.
    numerical = sorted(p for p in root.rglob("*") if p.is_file() and AUDIT not in p.parts
                       and "original_figures" not in p.parts
                       and (p.name.endswith(".tsv") or p.name.endswith(".tsv.gz") or p.suffix == ".npz"))
    numeric_hashes = {str(p.relative_to(root)): h.sha(p) for p in numerical}
    plan = dict(created_utc=now(), status="prepared", renderer=renderer(),
                expected_figures=1521, expected_files=3042, original_files=backups,
                numerical_input_sha256=numeric_hashes, datasets=datasets,
                source_tables={str(p.relative_to(root)): h.sha(p) for p in sources},
                free_gib=shutil.disk_usage(root).free / 1024 ** 3,
                preserved_response_figures=108)
    save_json(pp, plan)
    for name in ("README.md", "interpretation_zh.md", "redraw_manifest.json", "pdf_validation.json"):
        source = root / "figure_comparison_20261019" / name
        if source.exists():
            dest = directory / "documents_before" / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            if not dest.exists():
                shutil.copy2(source, dest)
    print("PREPARED 1521 figures; 3042 old files backed up; %d numerical files hashed" % len(numerical), flush=True)
    return plan


def valid_receipt(root, stem, config):
    path = receipt_path(root, stem)
    if not path.exists():
        return False
    record = json.loads(path.read_text())
    if record.get("renderer") != config:
        return False
    return all((root / r["path"]).exists() and h.sha(root / r["path"]) == r["sha256"]
               for r in record["files"])


def perform(root, stem, config, source, callback):
    if valid_receipt(root, stem, config):
        return 0
    callback()
    files = [stem.with_suffix(s) for s in (".pdf", ".png")]
    assert all(p.is_file() and p.stat().st_size > 1000 for p in files)
    save_json(receipt_path(root, stem), dict(status="complete", stem=str(stem.relative_to(root)),
               renderer=config, source=source, completed_utc=now(),
               files=[dict(path=str(p.relative_to(root)), sha256=h.sha(p), bytes=p.stat().st_size) for p in files]))
    return 1


def legacy_figure(lookup, active, groups, oxygen, stem, index, title, diagnostic=False):
    keys = [(p, "") for p in active] if index != "S2" else list(h.itertools.combinations(active, 2))
    labels = [pi + (" × " + pj if pj else " [" + groups[pi] + "]") for pi, pj in keys]
    fig, axes = plt.subplots(1, 2, figsize=(16, 23 if index == "S2" else 6.5), layout="constrained")
    for k, (ax, output) in enumerate(zip(axes, h.OUTPUTS)):
        data = h.matrix(lookup, keys, index, output, oxygen,
                        "pass_fraction_all_endpoints" if diagnostic else "value_mean")
        finite = data[np.isfinite(data)]
        if diagnostic:
            cmap, norm = plt.get_cmap("viridis"), Normalize(0, 1)
        else:
            lo = min(0., float(finite.min())) if len(finite) else 0.
            hi = float(finite.max()) if len(finite) and finite.max() > 0 else .001
            cmap = LinearSegmentedColormap.from_list("old_limits_purple",
                [h.NEGATIVE, "white", h.POSITIVE] if lo < 0 else ["white", h.POSITIVE])
            cmap.set_bad("#d9d9d9")
            norm = TwoSlopeNorm(vmin=lo, vcenter=0, vmax=hi) if lo < 0 else Normalize(0, hi)
        h.oxygen_mesh(ax, data, oxygen, cmap=cmap, norm=norm)
        ax.set_yticks(np.arange(len(labels)), labels if k == 0 else [], fontsize=7)
        ax.tick_params(axis="x", labelsize=8, length=2)
        ax.tick_params(axis="y", length=2)
        ax.set_title(h.OUT_LABEL[output], fontsize=11, pad=23)
        for spine in ax.spines.values():
            spine.set_linewidth(.5)
        cb = fig.colorbar(ax.collections[0], ax=ax, orientation="horizontal", pad=.025, fraction=.025)
        ticks = [0., .5, 1.] if diagnostic else ([float(norm.vmin), 0., float(norm.vmax)] if norm.vmin < 0
                                                else [0., float(norm.vmax) / 2, float(norm.vmax)])
        cb.set_ticks(ticks, labels=[format(v, ".3g") for v in ticks])
        cb.ax.tick_params(labelsize=8)
        cb.set_label("Fraction of 500 seeds passing diagnostics" if diagnostic else
                     index + " (raw estimates; blue = negative sampling estimates)", fontsize=8)
    fig.suptitle(title + "; N4096\nAll estimates retained; gray = undefined; diagnostics in source tables", fontsize=11)
    fig.supxlabel(h.OXYGEN_NOTE, fontsize=8)
    h.save_figure(fig, stem, [])


def main_figures(root, plan):
    settings()
    source_paths, hashes, rows, active, groups, pairs, oxygen, lookups = h.load_sources(root)
    review = root / "figure_comparison_20261019"
    config = plan["renderer"]
    source = {str(p.relative_to(root)): hashes[s] for s, p in source_paths.items()}
    scales = {}
    for index in ("S1", "ST", "S2"):
        keys = h.metric_rows(index, active, pairs)
        for output in h.OUTPUTS:
            values = np.concatenate([h.matrix(lookups[s], keys, index, output, oxygen).ravel() for s in h.SCOPES])
            scales[index, output] = max(float(np.nanmax(np.abs(values))), 1e-9), bool(np.nanmin(values) < 0)
    created = 0
    for scope in h.SCOPES:
        for index in ("S1", "ST", "S2"):
            name = index.lower() + "_parameter_by_oxygen" if index != "S2" else "s2_focal_pairs_by_oxygen"
            stem = root / scope / "figures" / name
            def plot(scope=scope, index=index, stem=stem):
                keys = h.metric_rows(index, active, pairs)
                labels = ([r["display_label"] for r in pairs] if index == "S2" else
                          [p + " [" + groups[p] + "]" for p in active])
                breaks = [i for i in range(1, 36) if pairs[i]["focal_parameter"] != pairs[i - 1]["focal_parameter"]]
                fig, axes = plt.subplots(1, 2, figsize=(15, 11.7 if index == "S2" else 6.2), layout="constrained")
                for k, (ax, output) in enumerate(zip(axes, h.OUTPUTS)):
                    h.draw(ax, h.matrix(lookups[scope], keys, index, output, oxygen), oxygen, labels, index,
                           show_labels=k == 0, pair_breaks=breaks if index == "S2" else ())
                    ax.set_title(h.OUT_LABEL[output], fontsize=11, pad=23)
                fig.suptitle(h.SCOPE_LABEL[scope] + " | " + ("36 focal parameter pairs" if index == "S2" else index)
                             + " | N4096, five repeats", fontsize=12)
                fig.supxlabel("Global: five-repeat mean; neighborhoods: median of 500 five-repeat means. All estimates retained.\n"
                              "Positive estimates: white to purple. Negative estimates: blue, not antagonism.\n" + h.OXYGEN_NOTE, fontsize=8)
                h.save_figure(fig, stem, [])
            created += perform(root, stem, config, source, plot)
        stem = root / scope / "figures/s2_pair_by_oxygen"
        created += perform(root, stem, config, source, lambda s=scope, t=stem:
                           legacy_figure(lookups[s], active, groups, oxygen, t, "S2",
                                         h.SCOPE_LABEL[s] + ("; five-repeat mean" if s == "global" else "; median of 500 five-repeat means")))
        if scope == "neighborhood10pct":
            stem = root / scope / "figures/s2_diagnostic_fraction"
            created += perform(root, stem, config, source, lambda t=stem:
                               legacy_figure(lookups["neighborhood10pct"], active, groups, oxygen, t, "S2",
                                             "Second-order diagnostic pass fraction", diagnostic=True))
        name = scope + "_s1_st_focal_s2"
        created += perform(root, review / name, config, source, lambda s=scope, n=name:
                           h.draw_grid(root, review, lookups, active, groups, pairs, oxygen,
                                       [s, s], list(h.OUTPUTS), n, h.SCOPE_LABEL[s] + ": first, total and focal second-order sensitivity", [], scales))
    groups_to_draw = [
        ("global_neighborhood_s1_st_focal_s2", ["global", "neighborhood10pct", "global", "neighborhood10pct"],
         [h.OUTPUTS[0], h.OUTPUTS[0], h.OUTPUTS[1], h.OUTPUTS[1]], "Global and neighborhood Sobol sensitivity", False),
        ("global_neighborhood_ploidy_s1_st_focal_s2", list(h.SCOPES), [h.OUTPUTS[0]] * 2,
         "Global and neighborhood comparison: dominant mean ploidy", False),
        ("global_neighborhood_growth_s1_st_focal_s2", list(h.SCOPES), [h.OUTPUTS[1]] * 2,
         "Global and neighborhood comparison: net live-cell growth rate", False),
        ("global_neighborhood_diagnostic_fraction", ["global", "neighborhood10pct", "global", "neighborhood10pct"],
         [h.OUTPUTS[0], h.OUTPUTS[0], h.OUTPUTS[1], h.OUTPUTS[1]], "Diagnostic coverage for displayed sensitivity cells", True)]
    for name, scopes, outputs, title, diag in groups_to_draw:
        created += perform(root, review / name, config, source,
                           lambda n=name, ss=scopes, oo=outputs, tt=title, dd=diag:
                           h.draw_grid(root, review, lookups, active, groups, pairs, oxygen, ss, oo, n, tt, [], scales, dd))
    print("MAIN complete: %d newly rendered figures" % created, flush=True)


def case_worker(root_string, dataset, config, active, groups):
    root = Path(root_string)
    settings()
    assert config == renderer(), "Renderer changed during run"
    stems = [root / dataset["directory"] / (index.lower() + "_parameter_by_oxygen" if index != "S2" else "s2_pair_by_oxygen")
             for index in ("S1", "ST", "S2")]
    source = root / dataset["source"]
    assert h.sha(source) == dataset["source_sha256"], "Numerical input changed"
    if all(valid_receipt(root, stem, config) for stem in stems):
        return dict(case=dataset["case"], created=0)
    rows = h.read_table(source)
    assert len(rows) == 47838 and {int(r["N"]) for r in rows} == {4096}
    lookup = {(r["index"], r["parameter_i"], r["parameter_j"], r["output"], float(r["O2_pct"])): r for r in rows}
    assert len(lookup) == len(rows)
    oxygen = np.array(sorted({float(r["O2_pct"]) for r in rows}))
    assert np.allclose(oxygen, np.linspace(0, 5, 201), rtol=0, atol=1e-12)
    created = 0
    for index, stem in zip(("S1", "ST", "S2"), stems):
        created += perform(root, stem, config, {dataset["source"]: dataset["source_sha256"]},
                           lambda ix=index, st=stem: legacy_figure(lookup, active, groups, oxygen, st, ix, dataset["title"]))
    return dict(case=dataset["case"], created=created)


def cases(root, plan, workers, case_ids=None):
    mapping = h.read_table(root / "model_structure/parameter_equation_map.tsv")
    active = [r["parameter"] for r in mapping if r["present_fixed_O2"] == "TRUE"]
    groups = {r["parameter"]: r["primary_group"] for r in mapping}
    datasets = [d for d in plan["datasets"] if case_ids is None or d["case"] in case_ids]
    completed, failures = [], []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(case_worker, str(root), d, plan["renderer"], active, groups): d["case"] for d in datasets}
        for future in as_completed(futures):
            case = futures[future]
            try:
                completed.append(future.result())
            except Exception as exc:
                failures.append(dict(case=case, error=repr(exc)))
            save_json(root / AUDIT / "progress.json", dict(updated_utc=now(), completed_datasets=len(completed),
                      requested_datasets=len(datasets), failures=failures,
                      newly_rendered=sum(r["created"] for r in completed)))
            if len(completed) % 5 == 0 or failures:
                print("DATASETS %d/%d; new figures %d; failures %d" %
                      (len(completed), len(datasets), sum(r["created"] for r in completed), len(failures)), flush=True)
    if failures:
        raise RuntimeError("Figure datasets failed: " + json.dumps(failures))


def finish(root, plan):
    records = [json.loads(p.read_text()) for p in (root / AUDIT / "receipts").glob("*.json")]
    assert len(records) == plan["expected_figures"], (len(records), plan["expected_figures"])
    outputs = {r["path"]: r for item in records for r in item["files"]}
    assert set(outputs) == {r["path"] for r in plan["original_files"]}
    for path, item in outputs.items():
        assert h.sha(root / path) == item["sha256"]
    for relative, expected in plan["numerical_input_sha256"].items():
        assert h.sha(root / relative) == expected, "Numerical file changed: " + relative
    for item in plan["original_files"]:
        assert h.sha(root / item["backup"]) == item["original_sha256"]
    surface_records = [r for r in h.read_table(root / "collection_manifest.tsv")
                       if "/figures/responses/" in r["path"] and r["path"].endswith((".pdf", ".png"))]
    assert len(surface_records) == 108
    assert all(h.sha(root / r["path"]) == r["sha256"] for r in surface_records)
    result = dict(status="complete", completed_utc=now(), figures=len(records), files=len(outputs),
                  renderer=plan["renderer"], all_numerical_files_unchanged=True,
                  n_numerical_files_verified=len(plan["numerical_input_sha256"]),
                  original_figure_backups_verified=True, model_evaluations=0, sensitivity_estimations=0,
                  response_surfaces_unchanged=True, n_response_files_verified=len(surface_records),
                  outputs=list(outputs.values()))
    save_json(root / AUDIT / "completion.json", result)
    h.write_table(root / AUDIT / "figure_manifest.tsv", list(outputs.values()))
    print(json.dumps({k: v for k, v in result.items() if k not in ("outputs", "renderer")}, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--phase", choices=("prepare", "main", "cases", "finish", "all"), default="all")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--case", action="append")
    parser.add_argument("--refresh-renderer", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    plan = prepare(root, args.refresh_renderer)
    if args.phase in ("main", "all"):
        main_figures(root, plan)
    if args.phase in ("cases", "all"):
        cases(root, plan, args.workers, set(args.case) if args.case else None)
    if args.phase in ("finish", "all"):
        finish(root, plan)


if __name__ == "__main__":
    main()
