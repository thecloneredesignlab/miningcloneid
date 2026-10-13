#!/usr/bin/env python3
"""Audited 12-endpoint, five-phase N=1025/2049 eFAST pilot.

The N=513 ensemble remains read-only. One worker evaluates a complete focal
trajectory. Summaries use small trajectory index caches, never concatenated
multi-million-row phase outputs. Thresholds are practical diagnostics, not a
proof of convergence or a gate for retaining an endpoint.
"""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from types import SimpleNamespace

import numpy as np
import efast
import neighborhood as nh
import neighborhood_slurm as ns

HERE = Path(__file__).resolve().parent
RESOLUTIONS = (1025, 2049)
KEY_PARAMETERS = ("buffer_beta", "buffer_smax", "p_mis_base", "mu_hp", "O2_crit")
ORIGINAL_SEEDS = ("seed25", "seed464", "seed165", "seed2", "seed214", "seed138",
                  "seed66", "seed395", "seed101", "seed311")
REPEAT_THRESHOLD = .10
RESOLUTION_THRESHOLD = .05
WINDOWS = (("overall_0_to_5pct", 0., 5.), ("low_0_to_1pct", 0., 1.),
           ("high_3_to_5pct", 3., 5.))
SOURCES = ("resolution_pilot.py", "run_resolution_pilot.sbatch",
           "submit_resolution_pilot.sh", "test_resolution_pilot.py",
           "efast.py", "neighborhood.py", "neighborhood_slurm.py",
           "evaluate_fixed_o2.R", "repair_eigen_points_safe.R",
           "perron_high_precision_safe.cpp", "resolution_pilot_smoke.R",
           "../../simulation/o2/fixed_o2/run_fixed_o2_simulation.R",
           "../../model/model_O2_supply_demand_MAP.R",
           "../../model/model_O2_supply_demand_MAP.cpp")


def save_npz(path, **arrays):
    temporary = path.with_name(path.stem + ".incomplete.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def phase_stats(values):
    """Only complete five-phase cells support a repeat/convergence claim."""
    count = np.isfinite(values).sum(axis=0)
    complete = count == values.shape[0]
    # Avoid all-NaN warnings and never turn an undefined index into a zero.
    mean = np.divide(np.nansum(values, axis=0), count,
                     out=np.full(count.shape, np.nan), where=count > 0)
    low = np.min(np.where(np.isfinite(values), values, np.inf), axis=0)
    high = np.max(np.where(np.isfinite(values), values, -np.inf), axis=0)
    span = np.where(complete, high - low, np.nan)
    sd = np.sqrt(np.divide(np.nansum((values - mean) ** 2, axis=0), count,
                          out=np.full(count.shape, np.nan), where=count > 0))
    return mean, np.where(complete, sd, np.nan), span, count


def band_scores(curves, oxygen, lo, hi):
    weights = efast.normalized_trapezoid_weights(oxygen, lo, hi)
    selected = weights > 0
    valid = np.isfinite(curves[..., selected]).all(axis=-1)
    return np.where(valid, np.sum(curves[..., selected] * weights[selected], axis=-1), np.nan)


def metric(curve, oxygen, mask):
    valid = mask & np.isfinite(curve)
    if not valid.any():
        return dict(p90=None, maximum=None, maximum_O2_pct=None, n_valid_oxygen=0)
    ids = np.flatnonzero(valid)
    maximum_at = ids[np.argmax(curve[valid])]
    return dict(p90=float(np.quantile(curve[valid], .9)), maximum=float(curve[maximum_at]),
                maximum_O2_pct=float(oxygen[maximum_at]), n_valid_oxygen=int(valid.sum()))


def flag(stats, total, threshold):
    if stats["n_valid_oxygen"] != total:
        return "undefined_or_incomplete"
    if stats["p90"] > threshold:
        return "not_passed"
    if stats["maximum"] > threshold:
        return "localized_review"
    return "passed_by_threshold"


def repeat_diagnostics(seed, n, values, oxygen):
    # values: phase, index, parameter, oxygen, output
    _, sd, span, count = phase_stats(values)
    rows = []
    for t, index in enumerate(nh.INDICES):
        for p, parameter in enumerate(efast.ACTIVE):
            for k, output in enumerate(efast.OUTPUTS):
                curves = values[:, t, p, :, k]
                for window, lo, hi in WINDOWS:
                    mask = (oxygen >= lo) & (oxygen <= hi)
                    stats = metric(span[t, p, :, k], oxygen, mask)
                    bands = band_scores(curves, oxygen, lo, hi)
                    band_valid = np.isfinite(bands).all()
                    peaks = [float(np.max(c[mask])) if np.isfinite(c[mask]).all() else np.nan for c in curves]
                    rows.append(dict(fit_seed=seed, N=n, output=output, index=index,
                        parameter=parameter, key_parameter=parameter in KEY_PARAMETERS, oxygen_window=window,
                        phase_range_p90=stats["p90"], phase_range_max=stats["maximum"],
                        phase_range_max_O2_pct=stats["maximum_O2_pct"],
                        n_valid_oxygen=stats["n_valid_oxygen"], n_expected_oxygen=int(mask.sum()),
                        n_cells_with_five_defined_phases=int((count[t,p,mask,k] == 5).sum()),
                        phase_sd_p90=metric(sd[t,p,:,k], oxygen, mask)["p90"],
                        phase_band_mean=float(np.mean(bands)) if band_valid else None,
                        phase_band_range=float(np.ptp(bands)) if band_valid else None,
                        phase_band_sd=float(np.std(bands)) if band_valid else None,
                        phase_peak_range=float(np.ptp(peaks)) if np.isfinite(peaks).all() else None,
                        threshold_absolute_index_units=REPEAT_THRESHOLD,
                        p90_pass=stats["n_valid_oxygen"] == int(mask.sum()) and stats["p90"] <= REPEAT_THRESHOLD,
                        diagnostic_status=flag(stats, int(mask.sum()), REPEAT_THRESHOLD)))
    return rows


def resolution_diagnostics(seed, low_n, high_n, low_values, high_values, oxygen):
    a, _, _, ca = phase_stats(low_values)
    b, _, _, cb = phase_stats(high_values)
    difference = np.where((ca == 5) & (cb == 5), np.abs(b-a), np.nan)
    rows = []
    for t, index in enumerate(nh.INDICES):
        for p, parameter in enumerate(efast.ACTIVE):
            for k, output in enumerate(efast.OUTPUTS):
                for window, lo, hi in WINDOWS:
                    mask = (oxygen >= lo) & (oxygen <= hi)
                    stats = metric(difference[t,p,:,k], oxygen, mask)
                    x, y = a[t,p,:,k], b[t,p,:,k]
                    all_valid = np.isfinite(difference[t,p,mask,k]).all()
                    rows.append(dict(fit_seed=seed, N_low=low_n, N_high=high_n,
                        primary_comparison=(low_n,high_n) == RESOLUTIONS,
                        output=output, index=index, parameter=parameter,
                        key_parameter=parameter in KEY_PARAMETERS, oxygen_window=window,
                        absolute_delta_p90=stats["p90"], absolute_delta_max=stats["maximum"],
                        absolute_delta_max_O2_pct=stats["maximum_O2_pct"],
                        n_valid_oxygen=stats["n_valid_oxygen"], n_expected_oxygen=int(mask.sum()),
                        band_delta=float(band_scores(y,oxygen,lo,hi)-band_scores(x,oxygen,lo,hi)) if all_valid else None,
                        peak_delta=float(np.max(y[mask])-np.max(x[mask])) if all_valid else None,
                        peak_O2_low=float(oxygen[mask][np.argmax(x[mask])]) if all_valid else None,
                        peak_O2_high=float(oxygen[mask][np.argmax(y[mask])]) if all_valid else None,
                        threshold_absolute_index_units=RESOLUTION_THRESHOLD,
                        p90_pass=stats["n_valid_oxygen"] == int(mask.sum()) and stats["p90"] <= RESOLUTION_THRESHOLD,
                        diagnostic_status=flag(stats,int(mask.sum()),RESOLUTION_THRESHOLD)))
    return rows


def task_rows(seeds):
    tasks, summaries = [], []
    for seed in seeds:
        number = int(seed.removeprefix("seed"))
        for n in RESOLUTIONS:
            summary_id = len(summaries) + 1
            summaries.append(dict(summary_id=summary_id, fit_seed=seed, N=n))
            for repeat in range(1, 6):
                for p, name in enumerate(efast.ACTIVE):
                    tasks.append(dict(task_id=len(tasks)+1, summary_id=summary_id,
                        fit_seed=seed, fit_seed_number=number, N=n, M=4, replicate=repeat,
                        parameter_index=p, focal_parameter=name,
                        phase_seed=202610080+number*10+repeat-1,
                        parent_sample_id_start=p*n+1, parent_sample_id_end=(p+1)*n))
    return tasks, summaries


def baseline_selection(baseline):
    manifest = efast.read_table(baseline / "seed_manifest.tsv")
    original = [r["fit_seed"] for r in sorted((r for r in manifest if r["pilot_order"]),
                                             key=lambda r:int(r["pilot_order"]))]
    if tuple(original) != ORIGINAL_SEEDS:
        raise ValueError("Original representative selection differs from the approved ten endpoints")
    collection = {r["path"]:r["sha256"] for r in efast.read_table(baseline / "artifact_collection_manifest.tsv")}
    scores, saved, source_hashes = [], {}, {}
    for path in sorted((baseline / "summaries").glob("phase_indices_*.npz")):
        digest = efast.sha256(path)
        relative = str(path.relative_to(baseline))
        if collection.get(relative) != digest:
            raise ValueError("Baseline source differs from collected audit: " + relative)
        source_hashes[str(path)] = digest
        with np.load(path) as data:
            if list(data["parameters"]) != list(efast.ACTIVE) or list(data["index_names"]) != list(nh.INDICES) or list(data["outputs"]) != list(efast.OUTPUTS):
                raise ValueError("Baseline axis labels differ")
            oxygen = data["oxygen"].copy()
            if not np.allclose(oxygen,np.arange(201)*.025,rtol=0,atol=1e-12):
                raise ValueError("Baseline oxygen grid differs")
            for seed, values in zip(data["fit_seeds"], data["indices"]):
                if values.shape != (5,2,14,201,2):
                    raise ValueError("Baseline must contain five phases")
                _, _, span, _ = phase_stats(values)
                full = span[1,:,:,0]
                valid = np.isfinite(full)
                overall = float(np.quantile(full[valid],.9)) if valid.any() else -1.
                key_scores = [float(np.nanquantile(full[efast.ACTIVE.index(p)],.9)) for p in KEY_PARAMETERS]
                which = int(np.argmax(key_scores))
                scores.append(dict(fit_seed=str(seed), ploidy_ST_overall_phase_range_p90=overall,
                                   max_key_parameter_phase_range_p90=key_scores[which],
                                   worst_key_parameter=KEY_PARAMETERS[which]))
                saved[str(seed)] = values.copy()
    if len(saved) != 500 or set(saved) != {r["fit_seed"] for r in manifest}:
        raise ValueError("Require all 500 audited baseline endpoints")
    eligible = [r for r in scores if r["fit_seed"] not in original]
    extra1 = max(eligible, key=lambda r:(r["ploidy_ST_overall_phase_range_p90"],-int(r["fit_seed"][4:])))
    extra2 = max((r for r in eligible if r != extra1),
                 key=lambda r:(r["max_key_parameter_phase_range_p90"],-int(r["fit_seed"][4:])))
    selected = original+[extra1["fit_seed"],extra2["fit_seed"]]
    if selected[-2:] != ["seed174","seed494"]:
        raise ValueError("Baseline diagnostics differ from approved two stress-test endpoints")
    selected_rows = []
    lookup = {r["fit_seed"]:r for r in scores}
    endpoint = {r["fit_seed"]:r for r in manifest}
    for i, seed in enumerate(selected):
        reason = "original normalized parameter-space/fit/boundary selection" if i<10 else (
            "largest overall ploidy ST repeat-range p90 outside original ten" if i==10 else
            "largest key-parameter ploidy ST repeat-range p90 among remaining endpoints")
        selected_rows.append({**endpoint[seed], **lookup[seed], "pilot_order":i+1,
                              "selection_reason":reason})
    return selected_rows, scores, np.stack([saved[s] for s in selected]), oxygen, source_hashes


def prepare(args):
    root, baseline = Path(args.out_dir).resolve(), Path(args.baseline_dir).resolve()
    if root == baseline or root != baseline / "pilot/N1025_N2049":
        raise ValueError("Use the approved separate pilot/N1025_N2049 directory")
    root.mkdir(parents=True, exist_ok=True)
    if (root / "slurm/submission_plan.json").exists():
        context(root)
        print("Retained verified existing pilot plan",flush=True)
        return
    inputs = json.loads((baseline / "input_manifest.json").read_text())
    for name, field in (("seed_manifest.tsv","seed_manifest_sha256"),
                        ("neighborhood_bounds.tsv","neighborhood_bounds_sha256"),
                        ("operator_input_audit.json","operator_audit_sha256")):
        if efast.sha256(baseline / name) != inputs[field]:
            raise ValueError("Baseline input audit changed: " + name)
    if efast.sha256(Path(args.figure4_dir)/"fixed_o2_dominant_ploidy_201grid.tsv") != inputs["reference_grid_sha256"]:
        raise ValueError("Figure 4 grid changed")
    selected, scores, values, oxygen, source_hashes = baseline_selection(baseline)
    for endpoint in selected:
        seed = endpoint["fit_seed"]
        for field, path in (("best_params",Path(args.fit_root)/seed/"best_params.tsv"),
                            ("fit_config",Path(args.fit_root)/seed/"fit_config.rds"),
                            ("parameter_table",Path(args.fit_root)/seed/"parameter_table.csv"),
                            ("bounds",baseline/"bounds"/(seed+".tsv"))):
            if efast.sha256(path) != endpoint[field+"_sha256"]:
                raise ValueError("Endpoint input changed: "+str(path))
    if shutil.disk_usage(root).free < 100*1024**3:
        raise ValueError("Require 100 GiB free for pilot raw outputs")
    sif = Path(args.sif).resolve()
    if efast.sha256(sif) != ns.SIF_SHA256:
        raise ValueError("Required SALib 1.5.2 SIF hash differs")
    archives = {}
    template = json.loads((baseline/"slurm/rcpp_template_manifest.json").read_text())
    if template["sif_sha256"] != ns.SIF_SHA256 or template["model_cpp_sha256"] != efast.sha256(HERE/"../../model/model_O2_supply_demand_MAP.cpp"):
        raise ValueError("Model cache source/image changed")
    for name, expected in (("rcpp_template.tar.gz",template["archive_sha256"]),
                           ("perron_safe_template.tar.gz","427b6dd7b24e5fde65612c4523adf185bbd0e6c7fb0643d323ab7e803160c2f5")):
        path = baseline/"slurm"/name
        if efast.sha256(path) != expected:
            raise ValueError("Private cache archive hash differs: "+name)
        archives[name] = dict(path=str(path),sha256=expected,bytes=path.stat().st_size,
                              mtime_ns=path.stat().st_mtime_ns)
    if efast.sha256(HERE/"perron_high_precision_safe.cpp") != "806c444f83b06bcc167c78404ca380b6aae495216bab8e71bdf7640ff29347f9":
        raise ValueError("Safe Perron source changed")
    directory = root/"slurm"
    directory.mkdir(exist_ok=True)
    for name, rows in (("seed_selection.tsv",selected),("selection_diagnostics_all500.tsv",scores)):
        efast.write_table(root/name,list(rows[0]),rows)
    save_npz(root/"baseline_N513.npz",indices=values,oxygen=oxygen,
             fit_seeds=[r["fit_seed"] for r in selected],parameters=efast.ACTIVE,
             outputs=efast.OUTPUTS,index_names=nh.INDICES)
    tasks, summaries = task_rows([r["fit_seed"] for r in selected])
    efast.write_table(directory/"task_manifest.tsv",list(tasks[0]),tasks)
    efast.write_table(directory/"summary_manifest.tsv",list(summaries[0]),summaries)
    plan = dict(prepared_utc=ns.now(),git_commit=subprocess.check_output(
        ["git","-C",str(HERE),"rev-parse","HEAD"],text=True).strip(),
        baseline_dir=str(baseline),result_root=str(root),fit_root=str(Path(args.fit_root).resolve()),
        figure4_dir=str(Path(args.figure4_dir).resolve()),sif=str(sif),sif_sha256=ns.SIF_SHA256,
        sif_bytes=sif.stat().st_size,sif_mtime_ns=sif.stat().st_mtime_ns,
        archives=archives,baseline_source_sha256=source_hashes,
        baseline_input_manifest_sha256=efast.sha256(baseline/"input_manifest.json"),
        baseline_selected_indices_sha256=efast.sha256(root/"baseline_N513.npz"),
        seed_selection_sha256=efast.sha256(root/"seed_selection.tsv"),
        task_manifest_sha256=efast.sha256(directory/"task_manifest.tsv"),
        summary_manifest_sha256=efast.sha256(directory/"summary_manifest.tsv"),
        code_sha256={name:efast.sha256(HERE/name) for name in SOURCES},
        resolutions=RESOLUTIONS,M=4,n_phases=5,n_parameters=14,n_oxygen=201,n_seeds=12,
        n_tasks=len(tasks),n_summaries=len(summaries),
        total_operator_evaluations=sum(r["N"]*201 for r in tasks),
        key_parameters=KEY_PARAMETERS,repeat_range_p90_threshold=REPEAT_THRESHOLD,
        resolution_delta_p90_threshold=RESOLUTION_THRESHOLD,
        convergence_policy="diagnostic_only",qos="xxlarge",time="12:00:00",
        specified_compute_node=None,array_concurrency_limit=None,
        granularity="seed x N x phase x one complete focal FAST trajectory; all 201 oxygen points",
        tasks_cpus=1,tasks_mem="4G",summary_cpus=1,summary_mem="8G",final_cpus=4,final_mem="32G",
        scope="selected pilot endpoints; no replacement of original 500-seed results")
    nh.atomic_json(directory/"submission_plan.json",plan)
    print(json.dumps(plan,indent=2),flush=True)


def context(root):
    plan = json.loads((root/"slurm/submission_plan.json").read_text())
    for name, digest in plan["code_sha256"].items():
        if efast.sha256(HERE/name) != digest:
            raise ValueError("Submitted pilot code changed: "+name)
    for name, field in (("seed_selection.tsv","seed_selection_sha256"),
                        ("slurm/task_manifest.tsv","task_manifest_sha256"),
                        ("slurm/summary_manifest.tsv","summary_manifest_sha256")):
        if efast.sha256(root/name) != plan[field]:
            raise ValueError("Submitted manifest changed: "+name)
    sif = Path(plan["sif"])
    if sif.stat().st_size != plan["sif_bytes"] or sif.stat().st_mtime_ns != plan["sif_mtime_ns"]:
        raise ValueError("Validated SIF changed")
    for archive in plan["archives"].values():
        path = Path(archive["path"])
        if path.stat().st_size != archive["bytes"] or path.stat().st_mtime_ns != archive["mtime_ns"]:
            raise ValueError("Validated cache archive changed")
    if efast.sha256(Path(plan["baseline_dir"])/"input_manifest.json") != plan["baseline_input_manifest_sha256"]:
        raise ValueError("Baseline inputs changed")
    return plan


def ensure_design(root, seed, n, plan, smoke=False):
    target = root/("smoke/runs" if smoke else "runs")/seed/("N%d"%n)
    with ns.locked(target/"prepare.lock"):
        endpoint = next(r for r in efast.read_table(root/"seed_selection.tsv") if r["fit_seed"] == seed)
        bounds = Path(plan["baseline_dir"])/"bounds"/(seed+".tsv")
        for field,path in (("best_params",Path(plan["fit_root"])/seed/"best_params.tsv"),
                           ("fit_config",Path(plan["fit_root"])/seed/"fit_config.rds"),
                           ("parameter_table",Path(plan["fit_root"])/seed/"parameter_table.csv"),
                           ("bounds",bounds)):
            if efast.sha256(path) != endpoint[field+"_sha256"]:
                raise ValueError("Endpoint inputs changed: "+str(path))
        receipt = target/"prepared.json"
        signature = dict(bounds_sha256=endpoint["bounds_sha256"],N=n,smoke=smoke,
                         submission_plan_sha256=efast.sha256(root/"slurm/submission_plan.json"))
        if receipt.exists():
            if json.loads(receipt.read_text()) != signature:
                raise ValueError("Prepared design signature changed")
        else:
            efast.prepare(SimpleNamespace(fit_root=plan["fit_root"],figure4_dir=plan["figure4_dir"],
                out_dir=str(target),n=n,m=4,replicates=1 if smoke else 5,
                seed=202610080+int(seed[4:])*10,full_grid=not smoke,oxygen="0,2.5,5",
                fit_seed=seed,ranges_file=str(bounds)))
            nh.atomic_json(receipt,signature)
    return target


def read_rows(path):
    with gzip.open(path,"rt",newline="") as handle:
        yield from csv.DictReader(handle,delimiter="\t")


def verify_parent(parent, task, plan):
    meta=json.loads((parent/"metadata.json").read_text())
    endpoint=next(r for r in efast.read_table(Path(plan["result_root"])/"seed_selection.tsv")
                  if r["fit_seed"]==task["fit_seed"])
    expected=dict(fit_seed=task["fit_seed"],N=int(task["N"]),M=4,replicate=int(task["replicate"]),
        seed=int(task["phase_seed"]),n_samples=int(task["N"])*14,n_parameters=14,
        salib_version="1.5.2",bounds_sha256=endpoint["bounds_sha256"],
        best_params_sha256=endpoint["best_params_sha256"],fit_config_sha256=endpoint["fit_config_sha256"])
    if any(meta.get(k)!=v for k,v in expected.items()) or not np.allclose(meta["oxygen_pct"],np.arange(201)*.025,rtol=0,atol=1e-12):
        raise ValueError("Parent design identity/grid differs from submitted task")


def fallback(root, part):
    """Recover only negative-vector points; all other failures stop the task."""
    from SALib.analyze.fast import compute_orders
    metadata = json.loads((part/"metadata.json").read_text())
    original = part/"outputs.tsv.gz.part01.gz"
    grid = np.asarray(metadata["oxygen_pct"])
    lookup = {round(float(o),10):j for j,o in enumerate(grid)}
    data = np.full((metadata["N"],len(grid),2),np.nan)
    seen = np.zeros(data.shape[:2],dtype=bool)
    bad = []
    for row in read_rows(original):
        i, j = int(row["sample_id"])-1, lookup.get(round(float(row["O2_pct"]),10))
        if j is None or not 0<=i<metadata["N"] or seen[i,j] or row["status"] != "ok":
            raise ValueError("Canonical failure is incomplete/duplicate/status-invalid; no numerical fallback")
        values = [float(row[k]) for k in efast.OUTPUTS]
        if not np.isfinite(values).all() or not np.isfinite(float(row["spectral_gap"])) or float(row["spectral_gap"])<0:
            raise ValueError("Canonical failure has nonfinite outputs or invalid spectrum")
        data[i,j] = values
        seen[i,j] = True
        if row["eigenvector_nonnegative"] != "TRUE":
            bad.append(row)
    if not seen.all() or not bad:
        raise ValueError("No complete negative-vector trajectory to recover")
    recovery = part/"recovery"
    recovery.mkdir(exist_ok=True)
    nh.atomic_json(recovery/"failed_points.json",bad)
    subprocess.run(["Rscript","--vanilla",str(HERE/"repair_eigen_points_safe.R"),str(root),str(part)],check=True)
    corrections = json.loads((recovery/"corrections.json").read_text())
    if corrections["status"] != "passed" or corrections["n_points"] != len(bad):
        raise ValueError("Incomplete high-precision proof")
    replacements = {(int(r["sample_id"]),round(float(r["O2_pct"]),10)):r for r in corrections["reports"]}
    if set(replacements) != {(int(r["sample_id"]),round(float(r["O2_pct"]),10)) for r in bad} or len(replacements)!=len(bad):
        raise ValueError("Recovered points differ from failed points")
    def corrected_rows():
        for row in read_rows(original):
            r = replacements.get((int(row["sample_id"]),round(float(row["O2_pct"]),10)))
            if r is not None:
                if r["original"] != row:
                    raise ValueError("Recovery source row differs")
                for k in efast.OUTPUTS:
                    row[k] = repr(float(r[k]))
                row["eigenvector_nonnegative"] = "TRUE"
            yield row
    temporary = part/"outputs.recovery.incomplete.tsv.gz"
    efast.write_table(temporary,ns.OUTPUT_FIELDS,corrected_rows())
    temporary.replace(part/"outputs.tsv.gz")
    before = np.full((2,len(grid),2),np.nan)
    for j in range(len(grid)):
        for k in range(2):
            y=data[:,j,k]
            if np.var(y)>1e-20*max(1.,float(np.mean(y))**2):
                before[:,j,k]=compute_orders(y,metadata["N"],metadata["M"],(metadata["N"]-1)//8)
    after = ns.analyze_trajectory(part)["indices"]
    impact = []
    for t,index in enumerate(nh.INDICES):
        for k,output in enumerate(efast.OUTPUTS):
            valid=np.isfinite(before[t,:,k])&np.isfinite(after[t,:,k])
            impact.append(dict(index=index,output=output,
                max_absolute_index_difference=float(np.max(np.abs(before[t,valid,k]-after[t,valid,k]))) if valid.any() else None,
                n_compared_oxygen_points=int(valid.sum()),
                undefined_pattern_changed=bool(np.any(np.isfinite(before[t,:,k])!=np.isfinite(after[t,:,k])))))
    artifacts={name:efast.sha256(recovery/name) for name in ("failed_points.json","corrections.json")}
    artifacts["../outputs.tsv.gz.part01.gz"]=efast.sha256(original)
    for r in corrections["reports"]:
        artifacts[r["matrix_file"]]=r["matrix_sha256"]
    proof=dict(status="passed",recovered_utc=ns.now(),n_recomputed_points=len(bad),
        n_unchanged_good_rows=int(seen.sum())-len(bad),index_impact=impact,
        method="Noda positive inverse iteration; 50/100-digit agreement, residual and leading-root bound",
        outputs_sha256=efast.sha256(part/"outputs.tsv.gz"),artifact_sha256=artifacts,
        code_sha256={name:efast.sha256(HERE/name) for name in ("resolution_pilot.py","repair_eigen_points_safe.R","perron_high_precision_safe.cpp")},
        no_matrix_regularization=True,no_sampling_change=True)
    nh.atomic_json(recovery/"numerical_recovery.json",proof)
    nh.atomic_json(part/"outputs.tsv.gz.receipt.json",dict(fit_seed=metadata["fit_seed"],
        n_rows=int(seen.sum()),metadata_sha256=efast.sha256(part/"metadata.json"),
        samples_sha256=metadata["samples_sha256"],outputs_sha256=proof["outputs_sha256"],
        evaluator_sha256=efast.sha256(HERE/"evaluate_fixed_o2.R"),
        recovery_manifest_sha256=efast.sha256(recovery/"numerical_recovery.json"),
        evaluation_method="canonical fixed-O2 operator with certified targeted fallback",workers=1,
        completed_at=ns.now()))


def evaluate_part(root, part, plan):
    if not nh.verify_completion(part):
        begin=time.monotonic()
        try:
            subprocess.run(["Rscript","--vanilla",str(HERE/"evaluate_fixed_o2.R"),
                "--samples="+str(part/"samples.tsv.gz"),"--metadata="+str(part/"metadata.json"),
                "--fit_root="+plan["fit_root"],"--figure4_dir="+plan["figure4_dir"],
                "--out="+str(part/"outputs.tsv.gz"),"--workers=1","--validate=TRUE"],check=True)
        except subprocess.CalledProcessError:
            if not (part/"outputs.tsv.gz.part01.gz").exists():
                raise
            fallback(root,part)
            receipt=json.loads((part/"outputs.tsv.gz.receipt.json").read_text())
            receipt["elapsed_seconds"]=time.monotonic()-begin
            nh.atomic_json(part/"outputs.tsv.gz.receipt.json",receipt)
    if not nh.verify_completion(part):
        raise ValueError("Missing complete trajectory receipt")
    result=ns.analyze_trajectory(part)
    save_npz(part/"indices.npz",**result)
    nh.atomic_json(part/"indices.receipt.json",dict(
        outputs_sha256=efast.sha256(part/"outputs.tsv.gz"),metadata_sha256=efast.sha256(part/"metadata.json"),
        cache_sha256=efast.sha256(part/"indices.npz"),analyzer_sha256=efast.sha256(HERE/"resolution_pilot.py"),
        salib_version=efast.salib_version()))
    return result


def worker(root, task_id):
    plan=context(root)
    rows=efast.read_table(root/"slurm/task_manifest.tsv")
    if not 1<=task_id<=len(rows):
        raise ValueError("Task ID outside pilot manifest")
    task=rows[task_id-1]
    marker=root/"slurm/task_receipts"/("task_%05d.json"%task_id)
    try:
        target=ensure_design(root,task["fit_seed"],int(task["N"]),plan)
        parent=target/"runs"/("N%s_R%s"%(task["N"],task["replicate"]))
        with ns.locked(parent/("trajectory_%s.lock"%task["parameter_index"])):
            verify_parent(parent,task,plan)
            part=ns.split_trajectory(parent,int(task["parameter_index"]))
            evaluate_part(root,part,plan)
            nh.atomic_json(marker,dict(task,status="complete",completed_utc=ns.now(),
                trajectory_receipt_sha256=efast.sha256(part/"outputs.tsv.gz.receipt.json"),
                indices_sha256=efast.sha256(part/"indices.npz"),
                slurm_job_id=os.environ.get("SLURM_JOB_ID"),host=os.uname().nodename))
    except Exception as error:
        nh.atomic_json(marker,dict(task,status="failed",error=str(error),recorded_utc=ns.now(),
            slurm_job_id=os.environ.get("SLURM_JOB_ID"),host=os.uname().nodename))
        raise


def load_part(part):
    if not nh.verify_completion(part):
        raise ValueError("Missing complete trajectory: "+str(part))
    receipt=json.loads((part/"indices.receipt.json").read_text())
    meta=json.loads((part/"metadata.json").read_text())
    if efast.sha256(part/"samples.tsv.gz")!=meta["samples_sha256"]:
        raise ValueError("Trajectory samples changed")
    for name,field in (("outputs.tsv.gz","outputs_sha256"),("metadata.json","metadata_sha256"),("indices.npz","cache_sha256")):
        if efast.sha256(part/name)!=receipt[field]:
            raise ValueError("Trajectory cache receipt changed: "+str(part))
    if receipt["analyzer_sha256"]!=efast.sha256(HERE/"resolution_pilot.py"):
        raise ValueError("Trajectory analyzer changed")
    with np.load(part/"indices.npz") as data:
        return {k:data[k] for k in data.files}


def seed_summary(root, summary_id):
    plan=context(root)
    summaries=efast.read_table(root/"slurm/summary_manifest.tsv")
    if not 1<=summary_id<=len(summaries):
        raise ValueError("Summary ID outside pilot manifest")
    summary=summaries[summary_id-1]
    seed,n=summary["fit_seed"],int(summary["N"])
    target=ensure_design(root,seed,n,plan)
    with ns.locked(target/"summary.lock"):
        tasks=[r for r in efast.read_table(root/"slurm/task_manifest.tsv") if int(r["summary_id"])==summary_id]
        result=np.full((5,2,14,201,2),np.nan)
        variances=np.full((5,14,201,2),np.nan)
        sources=[]
        recovery_rows=[]
        verified_parents=set()
        for task in tasks:
            marker=root/"slurm/task_receipts"/("task_%05d.json"%int(task["task_id"]))
            if not marker.exists() or json.loads(marker.read_text()).get("status")!="complete":
                raise ValueError("Incomplete task: "+task["task_id"])
            parent=target/"runs"/("N%d_R%s"%(n,task["replicate"]))
            if parent not in verified_parents:
                verify_parent(parent,task,plan)
                verified_parents.add(parent)
            p,r=int(task["parameter_index"]),int(task["replicate"])-1
            part=parent/"trajectories"/("P%02d_%s"%(p+1,efast.ACTIVE[p]))
            part_metadata=json.loads((part/"metadata.json").read_text())
            if part_metadata["parent_metadata_sha256"]!=efast.sha256(parent/"metadata.json") or part_metadata["trajectory_parameter_index"]!=p:
                raise ValueError("Trajectory parent/order changed")
            data=load_part(part)
            task_receipt=json.loads(marker.read_text())
            if task_receipt["indices_sha256"]!=efast.sha256(part/"indices.npz") or task_receipt["trajectory_receipt_sha256"]!=efast.sha256(part/"outputs.tsv.gz.receipt.json"):
                raise ValueError("Task marker differs from current trajectory receipt/cache")
            if data["indices"].shape!=(2,201,2) or not np.allclose(data["oxygen"],np.arange(201)*.025,rtol=0,atol=1e-12):
                raise ValueError("Trajectory index shape/grid differs")
            result[r,:,p]=data["indices"]
            variances[r,p]=data["variance"]
            sources.append(dict(task_id=task["task_id"],indices_sha256=efast.sha256(part/"indices.npz"),
                                receipt_sha256=efast.sha256(part/"outputs.tsv.gz.receipt.json")))
            proof=part/"recovery/numerical_recovery.json"
            if proof.exists():
                recovery_rows.append(dict(task_id=task["task_id"],**json.loads(proof.read_text())))
        oxygen=np.arange(201)*.025
        destination=root/"summaries"
        destination.mkdir(exist_ok=True)
        path=destination/("%s_N%d.npz"%(seed,n))
        save_npz(path,indices=result,variance=variances,oxygen=oxygen,parameters=efast.ACTIVE,
                 outputs=efast.OUTPUTS,index_names=nh.INDICES,fit_seed=seed,N=n)
        rows=repeat_diagnostics(seed,n,result,oxygen)
        table=destination/("%s_N%d_repeat_diagnostics.tsv"%(seed,n))
        efast.write_table(table,list(rows[0]),rows)
        nh.atomic_json(destination/("%s_N%d.receipt.json"%(seed,n)),dict(status="complete",
            fit_seed=seed,N=n,n_trajectories=len(tasks),indices_sha256=efast.sha256(path),
            diagnostics_sha256=efast.sha256(table),completed_utc=ns.now(),sources=sources,
            numerical_recovery=recovery_rows))
        print("Completed %s N%d summary from 70 validated trajectory caches"%(seed,n),flush=True)


def comparison_rows(seed,n,values,oxygen):
    means,_,_,counts=phase_stats(values)
    means=np.where(counts==5,means,np.nan)
    rows=[]
    for t,index in enumerate(nh.INDICES):
        for k,output in enumerate(efast.OUTPUTS):
            bands={band:band_scores(means[t,:,:,k],oxygen,lo,hi) for band,lo,hi in WINDOWS[1:]}
            ranks={band:np.argsort(np.argsort(-np.where(np.isfinite(v),v,-np.inf),kind="stable"),kind="stable")+1 for band,v in bands.items()}
            groups={band:{group:float(np.mean([v[p] for p,name in enumerate(efast.ACTIVE) if efast.GROUP[name]==group]))
                          for group in set(efast.GROUP.values())} for band,v in bands.items()}
            low,high=groups["low_0_to_1pct"],groups["high_3_to_5pct"]
            defined=all(np.isfinite(list(g.values())).all() for g in (low,high))
            death_largest=defined and low["death"]>=max(low.values())
            mis_increase=defined and high["missegregation"]>low["missegregation"]
            buffer_increase=defined and high["buffering"]>low["buffering"]
            for p,parameter in enumerate(efast.ACTIVE):
                l,h=bands["low_0_to_1pct"][p],bands["high_3_to_5pct"][p]
                rows.append(dict(fit_seed=seed,N=n,index=index,output=output,parameter=parameter,
                    low_band_mean=float(l),high_band_mean=float(h),high_minus_low=float(h-l),
                    low_band_rank=int(ranks["low_0_to_1pct"][p]) if np.isfinite(l) else None,
                    high_band_rank=int(ranks["high_3_to_5pct"][p]) if np.isfinite(h) else None,
                    n_oxygen_with_five_phases=int((counts[t,p,:,k]==5).sum()),
                    low_death_largest_group=death_largest if defined else None,
                    high_missegregation_increase=mis_increase if defined else None,
                    high_buffering_increase=buffer_increase if defined else None,
                    noemi_all_three_supported=(death_largest and mis_increase and buffer_increase) if defined else None,
                    aggregation="mean sensitivity across parameters within a process; all parameters equally weighted"))
    return rows


def status(root):
    plan=context(root)
    failed,complete,missing=[],[],[]
    for task_id in range(1,plan["n_tasks"]+1):
        marker=root/"slurm/task_receipts"/("task_%05d.json"%task_id)
        if not marker.exists():missing.append(task_id)
        elif json.loads(marker.read_text()).get("status")=="complete":complete.append(task_id)
        else:failed.append(task_id)
    summaries=list((root/"summaries").glob("seed*_N*.receipt.json"))
    report=dict(checked_utc=ns.now(),total_tasks=plan["n_tasks"],completed_tasks=len(complete),
        failed_tasks=len(failed),missing_or_running_tasks=len(missing),failed_task_ids=failed,
        missing_task_ids=missing,completed_seed_N_summaries=len(summaries),convergence_policy="diagnostic_only")
    nh.atomic_json(root/"slurm/status.json",report)
    print(json.dumps({k:v for k,v in report.items() if not k.endswith("_ids")},indent=2),flush=True)
    return report


def stability_rows(bands):
    lookup={(r["fit_seed"],r["N"],r["index"],r["output"],r["parameter"]):r for r in bands}
    rows=[]
    for seed in dict.fromkeys(r["fit_seed"] for r in bands):
        for low,high in ((513,1025),(1025,2049)):
            for index in nh.INDICES:
                for output in efast.OUTPUTS:
                    for parameter in efast.ACTIVE:
                        a=lookup[(seed,low,index,output,parameter)]
                        b=lookup[(seed,high,index,output,parameter)]
                        defined=np.isfinite([a["high_minus_low"],b["high_minus_low"]]).all()
                        rows.append(dict(fit_seed=seed,N_low=low,N_high=high,index=index,output=output,
                            parameter=parameter,key_parameter=parameter in KEY_PARAMETERS,
                            high_minus_low_at_lower_N=a["high_minus_low"],
                            high_minus_low_at_higher_N=b["high_minus_low"],
                            contrast_direction_unchanged=bool(np.sign(a["high_minus_low"])==np.sign(b["high_minus_low"])) if defined else None,
                            low_band_rank_delta=b["low_band_rank"]-a["low_band_rank"] if a["low_band_rank"] is not None and b["low_band_rank"] is not None else None,
                            high_band_rank_delta=b["high_band_rank"]-a["high_band_rank"] if a["high_band_rank"] is not None and b["high_band_rank"] is not None else None,
                            noemi_support_at_lower_N=a["noemi_all_three_supported"],
                            noemi_support_at_higher_N=b["noemi_all_three_supported"],
                            noemi_support_unchanged=a["noemi_all_three_supported"]==b["noemi_all_three_supported"] if a["noemi_all_three_supported"] is not None and b["noemi_all_three_supported"] is not None else None))
    return rows


def plot_comparison(root, all_values, oxygen, selected):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors={513:"#999999",1025:"#7551A8",2049:"#E6550D"}
    figs=root/"figures"
    figs.mkdir(exist_ok=True)
    for output,k in (("ploidy_ST",0),("growth_ST",1)):
        fig,axes=plt.subplots(12,5,figsize=(16,24),sharex=True,sharey=True)
        for i,endpoint in enumerate(selected):
            seed=endpoint["fit_seed"]
            for j,param in enumerate(KEY_PARAMETERS):
                p=efast.ACTIVE.index(param)
                for n in (513,*RESOLUTIONS):
                    values=all_values[(seed,n)][:,1,p,:,k]
                    mean,_,_,count=phase_stats(values)
                    defined=count==5
                    minimum=np.min(np.where(np.isfinite(values),values,np.inf),axis=0)
                    maximum=np.max(np.where(np.isfinite(values),values,-np.inf),axis=0)
                    axes[i,j].plot(oxygen,np.where(defined,mean,np.nan),color=colors[n],lw=1.1,label="N=%d"%n)
                    axes[i,j].fill_between(oxygen,minimum,maximum,where=defined,color=colors[n],alpha=.12,lw=0)
                axes[i,j].axvspan(0,1,color="#CFE2F3",alpha=.18)
                axes[i,j].axvspan(3,5,color="#FCE5CD",alpha=.18)
                if i==0:axes[i,j].set_title(param)
                if j==0:axes[i,j].set_ylabel(seed+"\nST")
                if i==11:axes[i,j].set_xlabel("Oxygen (%)")
                axes[i,j].set_xlim(0,5)
                axes[i,j].set_ylim(0,1.02)
        axes[0,0].legend(fontsize=7,loc="upper right")
        fig.suptitle(output.replace("_"," ")+": five-phase means and full phase ranges\nSelected pilot endpoints; shaded oxygen windows: 0-1% and 3-5%",fontsize=13)
        fig.tight_layout(rect=(0,0,1,.975))
        for ext in ("pdf","png"):fig.savefig(figs/(output+"_resolution_comparison."+ext),dpi=150)
        plt.close(fig)


def finalize(root):
    plan=context(root)
    with ns.locked(root/"slurm/finalize.lock"):
        report=status(root)
        if report["completed_tasks"]!=1680 or report["completed_seed_N_summaries"]!=24:
            raise ValueError("Pilot aggregation requires all 1,680 tasks and 24 summaries; see status.json")
        if efast.sha256(root/"baseline_N513.npz")!=plan["baseline_selected_indices_sha256"]:
            raise ValueError("Selected N513 baseline changed")
        selected=efast.read_table(root/"seed_selection.tsv")
        all_values={}
        with np.load(root/"baseline_N513.npz") as baseline:
            oxygen=baseline["oxygen"]
            for seed,values in zip(baseline["fit_seeds"],baseline["indices"]):all_values[(str(seed),513)]=values.copy()
        for endpoint in selected:
            seed=endpoint["fit_seed"]
            for n in RESOLUTIONS:
                prefix=root/"summaries"/("%s_N%d"%(seed,n))
                receipt=json.loads(prefix.with_suffix(".receipt.json").read_text())
                if receipt["status"]!="complete" or efast.sha256(prefix.with_suffix(".npz"))!=receipt["indices_sha256"]:
                    raise ValueError("Seed summary cache changed")
                if efast.sha256(prefix.with_name(prefix.name+"_repeat_diagnostics.tsv"))!=receipt["diagnostics_sha256"]:
                    raise ValueError("Seed diagnostic table changed")
                with np.load(prefix.with_suffix(".npz")) as data:all_values[(seed,n)]=data["indices"].copy()
        repeats,resolutions,bands=[],[],[]
        curve_rows=[]
        for endpoint in selected:
            seed=endpoint["fit_seed"]
            for n in (513,*RESOLUTIONS):
                values=all_values[(seed,n)]
                repeats.extend(repeat_diagnostics(seed,n,values,oxygen))
                bands.extend(comparison_rows(seed,n,values,oxygen))
                mean,sd,span,count=phase_stats(values)
                for t,index in enumerate(nh.INDICES):
                    for p,param in enumerate(efast.ACTIVE):
                        for j,o in enumerate(oxygen):
                            for k,output in enumerate(efast.OUTPUTS):
                                curve_rows.append(dict(fit_seed=seed,N=n,index=index,parameter=param,O2_pct=float(o),output=output,
                                    phase_mean=float(mean[t,p,j,k]),phase_sd=float(sd[t,p,j,k]),
                                    phase_range=float(span[t,p,j,k]),n_defined_phases=int(count[t,p,j,k])))
            for low,high in ((513,1025),(1025,2049)):
                resolutions.extend(resolution_diagnostics(seed,low,high,all_values[(seed,low)],all_values[(seed,high)],oxygen))
        for name,rows in (("repeat_diagnostics.tsv",repeats),("resolution_diagnostics.tsv",resolutions),
                          ("band_rank_mechanism_comparison.tsv",bands),
                          ("band_rank_mechanism_stability.tsv",stability_rows(bands)),
                          ("oxygen_curve_summary.tsv.gz",curve_rows)):
            efast.write_table(root/name,list(rows[0]),rows)
        for name,rows in (("key_ploidy_ST_repeat_diagnostics.tsv",repeats),("key_ploidy_ST_resolution_diagnostics.tsv",resolutions)):
            keep=[r for r in rows if r["key_parameter"] and r["index"]=="ST" and r["output"]=="dominant_mean_ploidy"]
            efast.write_table(root/name,list(keep[0]),keep)
        plot_comparison(root,all_values,oxygen,selected)
        nh.atomic_json(root/"completion.json",dict(status="complete",completed_utc=ns.now(),
            n_tasks=1680,n_seed_N_summaries=24,diagnostic_only=True,
            interpretation="Threshold passes are practical checks; localized failures and undefined cells require review. Selected pilot is not a population estimate.",
            artifact_sha256={str(p.relative_to(root)):efast.sha256(p) for p in root.glob("*.tsv*")}))
        print("Completed cross-resolution diagnostic tables and comparative ST figures",flush=True)


def smoke(root):
    plan=context(root)
    subprocess.run(["python3","-m","unittest","discover","-s",str(HERE),"-p","test_resolution_pilot.py"],check=True)
    results=[]
    for n,seed,p in ((1025,"seed174",10),(2049,"seed494",1)):
        target=ensure_design(root,seed,n,plan,smoke=True)
        parent=target/"runs"/("N%d_R1"%n)
        part=ns.split_trajectory(parent,p)
        begin=time.monotonic()
        result=evaluate_part(root,part,plan)
        results.append(dict(seed=seed,N=n,focal_parameter=efast.ACTIVE[p],
            oxygen_points=[0,2.5,5],n_rows=n*3,elapsed_seconds=time.monotonic()-begin,
            indices_shape=list(result["indices"].shape),nonfinite_indices=int((~np.isfinite(result["indices"])).sum()),
            output_sha256=efast.sha256(part/"outputs.tsv.gz")))
    # Test the existing 50/100-digit solver and targeted recovery on one actual
    # canonical matrix, independent of whether smoke trajectories happen to fail.
    part=ns.split_trajectory(parent,p)
    subprocess.run(["Rscript","--vanilla",str(HERE/"resolution_pilot_smoke.R"),str(root),str(part)],check=True)
    proof=json.loads((root/"smoke/perron_validation.json").read_text())
    if proof["status"]!="passed":raise ValueError("Smoke Perron validation failed")
    nh.atomic_json(root/"smoke/smoke_receipt.json",dict(status="passed",completed_utc=ns.now(),
        submission_plan_sha256=efast.sha256(root/"slurm/submission_plan.json"),
        unit_tests="test_resolution_pilot.py",trajectory_checks=results,perron_validation=proof,
        host=os.uname().nodename,slurm_job_id=os.environ.get("SLURM_JOB_ID")))
    print("Passed SIF unit tests, full-N three-oxygen trajectories, and actual-matrix 50/100-digit validation",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",choices=("prepare","worker","seed-summary","status","finalize","smoke"))
    parser.add_argument("--out-dir",required=True)
    parser.add_argument("--baseline-dir")
    parser.add_argument("--fit-root")
    parser.add_argument("--figure4-dir")
    parser.add_argument("--sif")
    parser.add_argument("--task-id",type=int)
    parser.add_argument("--summary-id",type=int)
    args=parser.parse_args()
    root=Path(args.out_dir).resolve()
    if args.command=="prepare":prepare(args)
    elif args.command=="worker":worker(root,args.task_id)
    elif args.command=="seed-summary":seed_summary(root,args.summary_id)
    elif args.command=="status":status(root)
    elif args.command=="smoke":smoke(root)
    else:finalize(root)


if __name__=="__main__":main()
