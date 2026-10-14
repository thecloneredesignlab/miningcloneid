#!/usr/bin/env python3
"""Six-endpoint N4097 pilot using the unchanged, audited eFAST evaluator."""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
import efast
import neighborhood as nh
import neighborhood_slurm as ns
import resolution_pilot as rp

HERE = Path(__file__).resolve().parent
SEEDS = ("seed25", "seed101", "seed138", "seed174", "seed66", "seed311")
N = 4097
BASELINE_NS = (513, 1025, 2049)
PAIRS = ((513, 1025), (1025, 2049), (2049, 4097))
REASONS = ("best fitting endpoint; unresolved key ST", "broad severe ST instability",
           "nonmonotonic ST and moving intermediate-oxygen peaks", "low-oxygen stress case",
           "repeat stability with remaining cross-N drift", "relatively stable comparison endpoint")
SOURCES = (*rp.SOURCES, "resolution_4097.py", "test_resolution_4097.py",
           "submit_resolution_4097.sh", "run_resolution_4097.sbatch")


def task_rows():
    tasks, summaries = [], []
    for seed in SEEDS:
        number = int(seed[4:])
        sid = len(summaries) + 1
        summaries.append(dict(summary_id=sid, fit_seed=seed, N=N))
        for repeat in range(1, 6):
            for p, parameter in enumerate(efast.ACTIVE):
                tasks.append(dict(task_id=len(tasks)+1, summary_id=sid, fit_seed=seed,
                    fit_seed_number=number, N=N, M=4, replicate=repeat, parameter_index=p,
                    focal_parameter=parameter, phase_seed=202610080+number*10+repeat-1,
                    parent_sample_id_start=p*N+1, parent_sample_id_end=(p+1)*N))
    return tasks, summaries


def checked_indices(data, seed, n):
    if (list(data["parameters"]) != list(efast.ACTIVE) or
        list(data["outputs"]) != list(efast.OUTPUTS) or
        list(data["index_names"]) != list(nh.INDICES) or
        not np.allclose(data["oxygen"], np.arange(201)*.025, rtol=0, atol=1e-12)):
        raise ValueError("Baseline axes/grid differ: %s N%d" % (seed, n))
    if str(data["fit_seed"].item()) != seed or int(data["N"]) != n:
        raise ValueError("Baseline endpoint/resolution differs")
    values = data["indices"].copy()
    if values.shape != (5, 2, 14, 201, 2):
        raise ValueError("Require five phases and complete axes")
    return values


def prepare(args):
    root = Path(args.out_dir).resolve()
    baseline = Path(args.baseline_dir).resolve()
    previous = baseline / "pilot/N1025_N2049"
    if root != baseline / "pilot/N4097":
        raise ValueError("Use the approved separate pilot/N4097 directory")
    root.mkdir(parents=True, exist_ok=True)
    if (root/"slurm/submission_plan.json").exists():
        context(root, verify_sources=True)
        print("Retained verified N4097 plan", flush=True)
        return
    old = rp.context(previous)
    completion = json.loads((previous/"completion.json").read_text())
    if completion["status"] != "complete" or completion["n_tasks"] != 1680:
        raise ValueError("Require the complete preceding pilot")
    if str(baseline) != old["baseline_dir"]:
        raise ValueError("Previous pilot belongs to different neighborhood inputs")
    inputs = json.loads((baseline/"input_manifest.json").read_text())
    for name, field in (("seed_manifest.tsv", "seed_manifest_sha256"),
                        ("neighborhood_bounds.tsv", "neighborhood_bounds_sha256"),
                        ("operator_input_audit.json", "operator_audit_sha256")):
        if efast.sha256(baseline/name) != inputs[field]:
            raise ValueError("Baseline audit changed: " + name)
    if efast.sha256(Path(old["figure4_dir"])/"fixed_o2_dominant_ploidy_201grid.tsv") != inputs["reference_grid_sha256"]:
        raise ValueError("Figure 4 oxygen grid changed")
    if efast.sha256(Path(old["sif"])) != ns.SIF_SHA256:
        raise ValueError("Required SIF hash differs")
    if shutil.disk_usage(root).free < 100*1024**3:
        raise ValueError("Require 100 GiB free for raw pilot outputs")
    endpoints = {r["fit_seed"]: r for r in efast.read_table(previous/"seed_selection.tsv")}
    selected = []
    for i, seed in enumerate(SEEDS):
        row = dict(endpoints[seed], pilot_order=i+1, selection_reason=REASONS[i])
        for field, path in (("best_params", Path(old["fit_root"])/seed/"best_params.tsv"),
                            ("fit_config", Path(old["fit_root"])/seed/"fit_config.rds"),
                            ("parameter_table", Path(old["fit_root"])/seed/"parameter_table.csv"),
                            ("bounds", baseline/"bounds"/(seed+".tsv"))):
            if efast.sha256(path) != row[field+"_sha256"]:
                raise ValueError("Fitted endpoint inputs changed: " + str(path))
        selected.append(row)
    sources = {str(previous/name): efast.sha256(previous/name) for name in
               ("slurm/submission_plan.json", "completion.json", "seed_selection.tsv", "baseline_N513.npz")}
    if sources[str(previous/"baseline_N513.npz")] != old["baseline_selected_indices_sha256"]:
        raise ValueError("Previous N513 baseline changed")
    values = {}
    with np.load(previous/"baseline_N513.npz") as data:
        if (list(data["parameters"]) != list(efast.ACTIVE) or
            list(data["outputs"]) != list(efast.OUTPUTS) or
            list(data["index_names"]) != list(nh.INDICES) or
            not np.allclose(data["oxygen"], np.arange(201)*.025, rtol=0, atol=1e-12)):
            raise ValueError("N513 baseline axes/grid differ")
        for seed, item in zip(data["fit_seeds"], data["indices"]):
            if str(seed) in SEEDS:
                if item.shape != (5, 2, 14, 201, 2):
                    raise ValueError("Incomplete N513 phase axes")
                values[(str(seed), 513)] = item.copy()
    for seed in SEEDS:
        for n in BASELINE_NS[1:]:
            prefix = previous/"summaries"/("%s_N%d" % (seed, n))
            receipt = json.loads(prefix.with_suffix(".receipt.json").read_text())
            files = ((prefix.with_suffix(".npz"), "indices_sha256"),
                     (prefix.with_name(prefix.name+"_repeat_diagnostics.tsv"), "diagnostics_sha256"))
            if receipt["status"] != "complete" or receipt["n_trajectories"] != 70:
                raise ValueError("Incomplete baseline summary")
            for path, field in files:
                digest = efast.sha256(path)
                if digest != receipt[field]:
                    raise ValueError("Baseline summary changed: " + str(path))
                sources[str(path)] = digest
            sources[str(prefix.with_suffix(".receipt.json"))] = efast.sha256(prefix.with_suffix(".receipt.json"))
            with np.load(prefix.with_suffix(".npz")) as data:
                values[(seed, n)] = checked_indices(data, seed, n)
    rp.save_npz(root/"baseline_indices.npz", indices=np.stack([
        np.stack([values[(s, n)] for n in BASELINE_NS]) for s in SEEDS]),
        fit_seeds=SEEDS, resolutions=BASELINE_NS, oxygen=np.arange(201)*.025,
        parameters=efast.ACTIVE, outputs=efast.OUTPUTS, index_names=nh.INDICES)
    efast.write_table(root/"seed_selection.tsv", list(selected[0]), selected)
    tasks, summaries = task_rows()
    directory = root/"slurm"
    directory.mkdir(exist_ok=True)
    for name, rows in (("task_manifest.tsv", tasks), ("summary_manifest.tsv", summaries)):
        efast.write_table(directory/name, list(rows[0]), rows)
    plan = copy.deepcopy(old)
    plan.update(prepared_utc=ns.now(), git_commit=subprocess.check_output(
        ["git", "-C", str(HERE), "rev-parse", "HEAD"], text=True).strip(),
        result_root=str(root), previous_pilot_dir=str(previous), resolutions=[N],
        comparison_resolutions=[*BASELINE_NS, N], primary_comparison=[2049, N],
        fit_seeds=list(SEEDS), n_seeds=6, n_tasks=len(tasks), n_summaries=len(summaries),
        total_operator_evaluations=sum(t["N"]*201 for t in tasks),
        comparison_source_sha256=sources,
        comparison_indices_sha256=efast.sha256(root/"baseline_indices.npz"),
        seed_selection_sha256=efast.sha256(root/"seed_selection.tsv"),
        task_manifest_sha256=efast.sha256(directory/"task_manifest.tsv"),
        summary_manifest_sha256=efast.sha256(directory/"summary_manifest.tsv"),
        code_sha256={name: efast.sha256(HERE/name) for name in SOURCES},
        scope="six preselected diagnostic endpoints; all 14 parameters; no population estimate",
        convergence_policy="diagnostic_only", qos="xxlarge", time="12:00:00",
        specified_compute_node=None, array_concurrency_limit=None)
    nh.atomic_json(directory/"submission_plan.json", plan)
    print(json.dumps({k: plan[k] for k in ("fit_seeds", "resolutions", "n_tasks", "n_summaries",
        "total_operator_evaluations", "git_commit", "qos", "time")}, indent=2), flush=True)


def context(root, verify_sources=False):
    plan = rp.context(root)
    if (plan["fit_seeds"] != list(SEEDS) or plan["resolutions"] != [N] or
        plan["n_tasks"] != 420 or plan["n_summaries"] != 6):
        raise ValueError("N4097 plan differs from approved scope")
    if efast.sha256(root/"baseline_indices.npz") != plan["comparison_indices_sha256"]:
        raise ValueError("Selected comparison arrays changed")
    if verify_sources:
        for name, digest in plan["comparison_source_sha256"].items():
            if efast.sha256(Path(name)) != digest:
                raise ValueError("Previous pilot source changed: " + name)
    return plan


def resolution_rows(seed, low, high, a, b, oxygen):
    rows = rp.resolution_diagnostics(seed, low, high, a, b, oxygen)
    for row in rows:
        row["primary_comparison"] = (low, high) == (2049, N)
    return rows


def stability_rows(bands):
    rows = rp.stability_rows(bands)
    lookup = {(r["fit_seed"], r["N"], r["index"], r["output"], r["parameter"]): r for r in bands}
    for seed in SEEDS:
        for index in nh.INDICES:
            for output in efast.OUTPUTS:
                for parameter in efast.ACTIVE:
                    a = lookup[(seed, 2049, index, output, parameter)]
                    b = lookup[(seed, N, index, output, parameter)]
                    defined = np.isfinite([a["high_minus_low"], b["high_minus_low"]]).all()
                    rows.append(dict(fit_seed=seed, N_low=2049, N_high=N, index=index, output=output,
                        parameter=parameter, key_parameter=parameter in rp.KEY_PARAMETERS,
                        high_minus_low_at_lower_N=a["high_minus_low"], high_minus_low_at_higher_N=b["high_minus_low"],
                        contrast_direction_unchanged=bool(np.sign(a["high_minus_low"]) == np.sign(b["high_minus_low"])) if defined else None,
                        low_band_rank_delta=b["low_band_rank"]-a["low_band_rank"] if a["low_band_rank"] is not None and b["low_band_rank"] is not None else None,
                        high_band_rank_delta=b["high_band_rank"]-a["high_band_rank"] if a["high_band_rank"] is not None and b["high_band_rank"] is not None else None,
                        noemi_support_at_lower_N=a["noemi_all_three_supported"], noemi_support_at_higher_N=b["noemi_all_three_supported"],
                        noemi_support_unchanged=a["noemi_all_three_supported"] == b["noemi_all_three_supported"] if a["noemi_all_three_supported"] is not None and b["noemi_all_three_supported"] is not None else None))
    return rows


def plot_comparison(root, values, oxygen):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {513: "#999999", 1025: "#7551A8", 2049: "#E6550D", N: "#087E8B"}
    dest = root/"figures"
    dest.mkdir(exist_ok=True)
    for label, k in (("ploidy_ST", 0), ("growth_ST", 1)):
        fig, axes = plt.subplots(6, 5, figsize=(17, 13), sharex=True, sharey=True)
        for i, seed in enumerate(SEEDS):
            for j, parameter in enumerate(rp.KEY_PARAMETERS):
                ax = axes[i, j]
                for n in (*BASELINE_NS, N):
                    y = values[(seed, n)][:, 1, efast.ACTIVE.index(parameter), :, k]
                    mean, _, _, count = rp.phase_stats(y)
                    valid = count == 5
                    ax.plot(oxygen, np.where(valid, mean, np.nan), color=colors[n], lw=1, label="N=%d" % n)
                    ax.fill_between(oxygen, np.min(np.where(np.isfinite(y), y, np.inf), axis=0),
                        np.max(np.where(np.isfinite(y), y, -np.inf), axis=0), where=valid,
                        color=colors[n], alpha=.10, lw=0)
                ax.axvspan(0, 1, color="#CFE2F3", alpha=.18)
                ax.axvspan(3, 5, color="#FCE5CD", alpha=.18)
                ax.set(xlim=(0, 5), ylim=(0, 1.02))
                if i == 0: ax.set_title(parameter)
                if j == 0: ax.set_ylabel(seed+"\nST")
                if i == 5: ax.set_xlabel("Oxygen (%)")
        axes[0, 0].legend(fontsize=7)
        fig.suptitle(label.replace("_", " ")+": five-phase means and full phase ranges\nSix diagnostic endpoints; shaded windows: 0-1% and 3-5% O2", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, .95))
        for ext in ("pdf", "png"): fig.savefig(dest/(label+"_resolution_comparison."+ext), dpi=150)
        plt.close(fig)
    # These cells show error metrics, not sensitivity indices.
    focus = {}
    for seed in SEEDS:
        for n in (2049, N):
            for row in rp.repeat_diagnostics(seed, n, values[(seed, n)], oxygen):
                if row["index"] == "ST" and row["output"] == efast.OUTPUTS[0] and row["oxygen_window"] == rp.WINDOWS[0][0]:
                    focus[(seed, n, row["parameter"])] = row
        for row in resolution_rows(seed, 2049, N, values[(seed, 2049)], values[(seed, N)], oxygen):
            if row["index"] == "ST" and row["output"] == efast.OUTPUTS[0] and row["oxygen_window"] == rp.WINDOWS[0][0]:
                focus[(seed, "delta", row["parameter"])] = row
    from matplotlib.colors import ListedColormap
    fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=True)
    for ax, key, title in zip(axes, (2049, N, "delta"), ("A  N2049 phase range", "B  N4097 phase range", "C  N2049 to N4097 mean change")):
        cells = np.zeros((6, 5), dtype=int)
        for i, seed in enumerate(SEEDS):
            for j, p in enumerate(rp.KEY_PARAMETERS):
                r = focus[(seed, key, p)]
                cells[i, j] = {"passed_by_threshold": 0, "localized_review": 1, "not_passed": 2, "undefined_or_incomplete": 3}[r["diagnostic_status"]]
        ax.imshow(cells, cmap=ListedColormap(["#eeeeee", "#d9c6ec", "#f6b37b", "#777777"]), vmin=0, vmax=3, aspect="auto")
        for i, seed in enumerate(SEEDS):
            for j, p in enumerate(rp.KEY_PARAMETERS):
                r = focus[(seed, key, p)]
                v = r["absolute_delta_p90" if key == "delta" else "phase_range_p90"]
                ax.text(j, i, "NA" if v is None else "%.3f%s" % (v, "*" if cells[i, j] == 1 else ""), ha="center", va="center")
        ax.set(xticks=range(5), yticks=range(6), xticklabels=rp.KEY_PARAMETERS, yticklabels=SEEDS, title=title)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.suptitle("Ploidy ST numerical stability across the full 0-5% O2 grid")
    fig.text(.05, .02, "Numbers: p90 across 201 O2 points. Grey: p90 and maximum pass; purple *: local maximum fails; orange: p90 fails; dark grey: incomplete.\nPhase-range threshold 0.10; mean-change threshold 0.05 (absolute index units). Five-phase ranges are not confidence intervals.", fontsize=9)
    fig.tight_layout(rect=(0, .12, 1, .94))
    for ext in ("pdf", "png"): fig.savefig(dest/("ploidy_ST_convergence_overview."+ext), dpi=150)
    plt.close(fig)


def finalize(root):
    plan = context(root, verify_sources=True)
    with ns.locked(root/"slurm/finalize.lock"):
        status = rp.status(root)
        if status["completed_tasks"] != 420 or status["completed_seed_N_summaries"] != 6:
            raise ValueError("Require all 420 trajectories and six validated summaries")
        values = {}
        with np.load(root/"baseline_indices.npz") as data:
            oxygen = data["oxygen"].copy()
            for i, seed in enumerate(SEEDS):
                for j, n in enumerate(BASELINE_NS): values[(seed, n)] = data["indices"][i, j].copy()
        for seed in SEEDS:
            prefix = root/"summaries"/(seed+"_N4097")
            receipt = json.loads(prefix.with_suffix(".receipt.json").read_text())
            if receipt["status"] != "complete" or receipt["n_trajectories"] != 70:
                raise ValueError("Incomplete N4097 summary")
            for path, field in ((prefix.with_suffix(".npz"), "indices_sha256"),
                               (prefix.with_name(prefix.name+"_repeat_diagnostics.tsv"), "diagnostics_sha256")):
                if efast.sha256(path) != receipt[field]: raise ValueError("Summary hash differs")
            with np.load(prefix.with_suffix(".npz")) as data: values[(seed, N)] = checked_indices(data, seed, N)
        repeats, resolutions, bands, curves = [], [], [], []
        for seed in SEEDS:
            for n in (*BASELINE_NS, N):
                v = values[(seed, n)]
                repeats.extend(rp.repeat_diagnostics(seed, n, v, oxygen))
                bands.extend(rp.comparison_rows(seed, n, v, oxygen))
                mean, sd, span, count = rp.phase_stats(v)
                for t, index in enumerate(nh.INDICES):
                    for p, param in enumerate(efast.ACTIVE):
                        for j, o in enumerate(oxygen):
                            for k, output in enumerate(efast.OUTPUTS):
                                curves.append(dict(fit_seed=seed, N=n, index=index, parameter=param, O2_pct=float(o), output=output,
                                    phase_mean=float(mean[t, p, j, k]), phase_sd=float(sd[t, p, j, k]),
                                    phase_range=float(span[t, p, j, k]), n_defined_phases=int(count[t, p, j, k])))
            for low, high in PAIRS:
                resolutions.extend(resolution_rows(seed, low, high, values[(seed, low)], values[(seed, high)], oxygen))
        for name, rows in (("repeat_diagnostics.tsv", repeats), ("resolution_diagnostics.tsv", resolutions),
                           ("band_rank_mechanism_comparison.tsv", bands), ("band_rank_mechanism_stability.tsv", stability_rows(bands)),
                           ("oxygen_curve_summary.tsv.gz", curves)):
            efast.write_table(root/name, list(rows[0]), rows)
        for name, rows in (("key_ploidy_ST_repeat_diagnostics.tsv", repeats), ("key_ploidy_ST_resolution_diagnostics.tsv", resolutions)):
            keep = [r for r in rows if r["key_parameter"] and r["index"] == "ST" and r["output"] == efast.OUTPUTS[0]]
            efast.write_table(root/name, list(keep[0]), keep)
        plot_comparison(root, values, oxygen)
        nh.atomic_json(root/"completion.json", dict(status="complete", completed_utc=ns.now(), n_tasks=420,
            n_seed_N_summaries=6, diagnostic_only=True, primary_comparison=[2049, N],
            interpretation="Practical repeat/resolution diagnostics; selected endpoints are not a population sample.",
            artifact_sha256={str(p.relative_to(root)): efast.sha256(p) for p in [*root.glob("*.tsv*"), *root.glob("figures/*")]}))
        print("Completed N4097 repeat/cross-resolution tables and ST figures", flush=True)


def smoke(root):
    plan = context(root, verify_sources=True)
    for pattern in ("test_resolution_pilot.py", "test_resolution_4097.py"):
        subprocess.run(["python3", "-m", "unittest", "discover", "-s", str(HERE), "-p", pattern], check=True)
    seed, p = "seed138", efast.ACTIVE.index("p_mis_base")
    target = rp.ensure_design(root, seed, N, plan, smoke=True)
    parent = target/"runs/N4097_R1"
    part = ns.split_trajectory(parent, p)
    begin = time.monotonic()
    result = rp.evaluate_part(root, part, plan)
    if result["indices"].shape != (2, 3, 2) or not np.isfinite(result["indices"]).all():
        raise ValueError("N4097 smoke indices are incomplete/nonfinite")
    subprocess.run(["Rscript", "--vanilla", str(HERE/"resolution_pilot_smoke.R"), str(root), str(part)], check=True)
    proof = json.loads((root/"smoke/perron_validation.json").read_text())
    if proof["status"] != "passed": raise ValueError("Perron smoke failed")
    nh.atomic_json(root/"smoke/smoke_receipt.json", dict(status="passed", completed_utc=ns.now(),
        submission_plan_sha256=efast.sha256(root/"slurm/submission_plan.json"),
        unit_tests=["test_resolution_pilot.py", "test_resolution_4097.py"],
        trajectory_checks=[dict(seed=seed, N=N, focal_parameter=efast.ACTIVE[p], oxygen_points=[0, 2.5, 5],
            n_rows=N*3, elapsed_seconds=time.monotonic()-begin, indices_shape=list(result["indices"].shape),
            output_sha256=efast.sha256(part/"outputs.tsv.gz"))], perron_validation=proof,
        host=os.uname().nodename, slurm_job_id=os.environ.get("SLURM_JOB_ID")))
    print("Passed N4097 SIF smoke, analytic FFT tests, and 50/100-digit matrix check", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "worker", "seed-summary", "finalize", "smoke", "status"))
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--baseline-dir")
    parser.add_argument("--task-id", type=int)
    parser.add_argument("--summary-id", type=int)
    args = parser.parse_args()
    root = Path(args.out_dir).resolve()
    if args.command == "prepare": prepare(args)
    elif args.command == "smoke": smoke(root)
    elif args.command == "finalize": finalize(root)
    else:
        context(root)
        if args.command == "worker": rp.worker(root, args.task_id)
        elif args.command == "seed-summary": rp.seed_summary(root, args.summary_id)
        else: rp.status(root)


if __name__ == "__main__": main()
