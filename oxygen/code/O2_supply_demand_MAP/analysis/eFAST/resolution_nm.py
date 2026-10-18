#!/usr/bin/env python3
"""Four-endpoint eFAST pilot with independently labelled N and M comparisons.

The canonical model, sample transformations, FFT estimator, and certified
Perron solver remain imported from the frozen preceding pilot. This adapter
carries M through preparation, cache identity, recovery diagnostics and reports.
"""
import argparse
import copy
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np
import efast
import neighborhood as nh
import neighborhood_slurm as ns
import resolution_pilot as rp
import resolution_4097 as previous

HERE = Path(__file__).resolve().parent
PILOT = "N8193_M4_and_N4097_M8"
CONFIGS = {
    "N8193_M4": dict(N=8193, M=4, seeds=("seed25", "seed138"),
                     comparison_kind="N_resolution", baseline_N=4097, baseline_M=4),
    "N4097_M8": dict(N=4097, M=8, seeds=("seed101", "seed174"),
                     comparison_kind="M_robustness", baseline_N=4097, baseline_M=4),
}
SOURCES = (*previous.SOURCES, "resolution_nm.py", "test_resolution_nm.py",
           "submit_resolution_nm.sh", "run_resolution_nm.sbatch")
save_npz = rp.save_npz
read_rows = rp.read_rows


def frequency(n, m):
    if n <= 4*m*m or m < 1:
        raise ValueError("FAST requires positive M and N > 4*M*M")
    omega = (n-1)//(2*m)
    if omega//(2*m) < len(efast.ACTIVE)-1:
        raise ValueError("Insufficient distinct complementary frequencies for 14 parameters")
    return omega


def orders(y, n, m):
    """Use the same M-dependent frequency in canonical and recovery diagnostics."""
    from SALib.analyze.fast import compute_orders
    return compute_orders(y, n, m, frequency(n, m))


def task_rows(config):
    spec = CONFIGS[config]
    n, m = spec["N"], spec["M"]
    frequency(n, m)
    tasks, summaries = [], []
    for seed in spec["seeds"]:
        sid = len(summaries)+1
        summaries.append(dict(summary_id=sid, fit_seed=seed, N=n, M=m, config=config))
        for repeat in range(1, 6):
            for p, parameter in enumerate(efast.ACTIVE):
                tasks.append(dict(task_id=len(tasks)+1, summary_id=sid,
                    config=config, fit_seed=seed, fit_seed_number=int(seed[4:]), N=n, M=m,
                    replicate=repeat, parameter_index=p, focal_parameter=parameter,
                    phase_seed=202610080+int(seed[4:])*10+repeat-1,
                    parent_sample_id_start=p*n+1, parent_sample_id_end=(p+1)*n))
    return tasks, summaries


def write_rows(path, rows):
    if not rows:
        raise ValueError("Refuse an empty report table: "+str(path))
    efast.write_table(path, list(rows[0]), rows)


def prepare(args):
    root, baseline = Path(args.out_dir).resolve(), Path(args.baseline_dir).resolve()
    config, spec = args.config, CONFIGS[args.config]
    if root != baseline/"pilot"/PILOT/config:
        raise ValueError("Use the approved separate pilot/config directory")
    root.mkdir(parents=True, exist_ok=True)
    if (root/"slurm/submission_plan.json").exists():
        context(root, verify_sources=True)
        print("Retained verified "+config+" plan", flush=True)
        return
    old_root = baseline/"pilot/N4097"
    old = previous.context(old_root, verify_sources=True)
    complete = json.loads((old_root/"completion.json").read_text())
    if complete["status"] != "complete" or complete["n_tasks"] != 420:
        raise ValueError("Require complete six-seed N4097 baseline")
    if efast.salib_version() != "1.5.2" or efast.sha256(Path(old["sif"])) != ns.SIF_SHA256:
        raise ValueError("Require the approved SALib 1.5.2 image")
    if shutil.disk_usage(root).free < 100*1024**3:
        raise ValueError("Require 100 GiB available for raw outputs")
    selected = []
    endpoints = {r["fit_seed"]:r for r in efast.read_table(old_root/"seed_selection.tsv")}
    for i, seed in enumerate(spec["seeds"]):
        row = dict(endpoints[seed], pilot_order=i+1, config=config)
        row["selection_reason"] = ("follow-up to unresolved N2049-to-N4097 ploidy ST; "
            + ("increase N at M4" if spec["M"] == 4 else "compare M8 against M4 at fixed N4097"))
        for field, path in (("best_params", Path(old["fit_root"])/seed/"best_params.tsv"),
                            ("fit_config", Path(old["fit_root"])/seed/"fit_config.rds"),
                            ("parameter_table", Path(old["fit_root"])/seed/"parameter_table.csv"),
                            ("bounds", baseline/"bounds"/(seed+".tsv"))):
            if efast.sha256(path) != row[field+"_sha256"]:
                raise ValueError("Fitted endpoint inputs changed: "+str(path))
        selected.append(row)
    inputs = json.loads((baseline/"input_manifest.json").read_text())
    for name, field in (("seed_manifest.tsv", "seed_manifest_sha256"),
                        ("neighborhood_bounds.tsv", "neighborhood_bounds_sha256"),
                        ("operator_input_audit.json", "operator_audit_sha256")):
        if efast.sha256(baseline/name) != inputs[field]:
            raise ValueError("Baseline audit changed: "+name)
    if efast.sha256(Path(old["figure4_dir"])/"fixed_o2_dominant_ploidy_201grid.tsv") != inputs["reference_grid_sha256"]:
        raise ValueError("Reference oxygen grid changed")
    sources = {str(old_root/name):efast.sha256(old_root/name) for name in
               ("slurm/submission_plan.json", "completion.json", "seed_selection.tsv", "baseline_indices.npz")}
    values = {}
    with np.load(old_root/"baseline_indices.npz") as data:
        for seed in spec["seeds"]:
            i = list(data["fit_seeds"]).index(seed)
            j = list(data["resolutions"]).index(2049)
            values[(seed,2049)] = data["indices"][i,j].copy()
    for seed in spec["seeds"]:
        prefix = old_root/"summaries"/(seed+"_N4097")
        receipt_path = prefix.with_suffix(".receipt.json")
        receipt = json.loads(receipt_path.read_text())
        if receipt["status"] != "complete" or receipt["n_trajectories"] != 70:
            raise ValueError("Incomplete N4097 baseline summary")
        sources[str(receipt_path)] = efast.sha256(receipt_path)
        for path, field in ((prefix.with_suffix(".npz"), "indices_sha256"),
                            (prefix.with_name(prefix.name+"_repeat_diagnostics.tsv"), "diagnostics_sha256")):
            if efast.sha256(path) != receipt[field]:
                raise ValueError("Baseline summary hash differs")
            sources[str(path)] = efast.sha256(path)
        with np.load(prefix.with_suffix(".npz")) as data:
            values[(seed,4097)] = previous.checked_indices(data,seed,4097)
    save_npz(root/"baseline_indices.npz", indices=np.stack([
        np.stack([values[(s,n)] for n in (2049,4097)]) for s in spec["seeds"]]),
        fit_seeds=spec["seeds"], resolutions=(2049,4097), M=4,
        oxygen=np.arange(201)*.025, parameters=efast.ACTIVE, outputs=efast.OUTPUTS, index_names=nh.INDICES)
    write_rows(root/"seed_selection.tsv",selected)
    tasks, summaries = task_rows(config)
    directory = root/"slurm"
    directory.mkdir(exist_ok=True)
    write_rows(directory/"task_manifest.tsv",tasks)
    write_rows(directory/"summary_manifest.tsv",summaries)
    plan = copy.deepcopy(old)
    # Remove previous stage's comparison selectors to avoid ambiguous provenance.
    for key in ("primary_comparison", "comparison_resolutions", "baseline_selected_indices_sha256",
                "baseline_source_sha256"):
        plan.pop(key,None)
    plan.update(prepared_utc=ns.now(), git_commit=subprocess.check_output(
        ["git","-C",str(HERE),"rev-parse","HEAD"],text=True).strip(),
        result_root=str(root), previous_pilot_dir=str(old_root), config=config,
        resolutions=[spec["N"]], M=spec["M"], fit_seeds=list(spec["seeds"]),
        comparison_kind=spec["comparison_kind"],
        baseline_config=dict(N=4097,M=4), contextual_baseline_config=dict(N=2049,M=4),
        n_seeds=2, n_tasks=140, n_summaries=2, n_phases=5,
        total_operator_evaluations=sum(t["N"]*201 for t in tasks),
        comparison_source_sha256=sources,
        comparison_indices_sha256=efast.sha256(root/"baseline_indices.npz"),
        seed_selection_sha256=efast.sha256(root/"seed_selection.tsv"),
        task_manifest_sha256=efast.sha256(directory/"task_manifest.tsv"),
        summary_manifest_sha256=efast.sha256(directory/"summary_manifest.tsv"),
        code_sha256={name:efast.sha256(HERE/name) for name in SOURCES},
        scope="two selected diagnostic endpoints per configuration; all 14 parameters; no population estimate",
        convergence_policy="diagnostic_only", qos="xxlarge", time="12:00:00",
        specified_compute_node=None, array_concurrency_limit=None)
    nh.atomic_json(directory/"submission_plan.json",plan)
    print(json.dumps({k:plan[k] for k in ("config","fit_seeds","M","resolutions","n_tasks",
        "total_operator_evaluations","comparison_kind","git_commit")},indent=2),flush=True)


def context(root, verify_sources=False):
    plan = rp.context(root)
    config = plan["config"]
    spec = CONFIGS[config]
    if (root.name != config or plan["fit_seeds"] != list(spec["seeds"]) or
        plan["resolutions"] != [spec["N"]] or plan["M"] != spec["M"] or
        plan["comparison_kind"] != spec["comparison_kind"] or
        plan["n_tasks"] != 140 or plan["n_summaries"] != 2):
        raise ValueError("Plan differs from approved N/M configuration")
    if efast.sha256(root/"baseline_indices.npz") != plan["comparison_indices_sha256"]:
        raise ValueError("Selected baseline changed")
    if verify_sources:
        for name, digest in plan["comparison_source_sha256"].items():
            if efast.sha256(Path(name)) != digest:
                raise ValueError("Comparison source changed: "+name)
    return plan


def require_smoke(root):
    receipt = json.loads((root/"smoke/smoke_receipt.json").read_text())
    plan = context(root)
    if (receipt["status"] != "passed" or receipt["config"] != plan["config"] or
        receipt["submission_plan_sha256"] != efast.sha256(root/"slurm/submission_plan.json")):
        raise ValueError("Require passed smoke for this exact configuration/plan")

def ensure_design(root, seed, n, plan, smoke=False):
    m=plan["M"]
    target = root/("smoke/runs" if smoke else "runs")/seed/("N%d_M%d"%(n,m))
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
        signature = dict(bounds_sha256=endpoint["bounds_sha256"],N=n,M=m,smoke=smoke,
                         submission_plan_sha256=efast.sha256(root/"slurm/submission_plan.json"))
        if receipt.exists():
            if json.loads(receipt.read_text()) != signature:
                raise ValueError("Prepared design signature changed")
        else:
            efast.prepare(SimpleNamespace(fit_root=plan["fit_root"],figure4_dir=plan["figure4_dir"],
                out_dir=str(target),n=n,m=m,replicates=1 if smoke else 5,
                seed=202610080+int(seed[4:])*10,full_grid=not smoke,oxygen="0,2.5,5",
                fit_seed=seed,ranges_file=str(bounds)))
            nh.atomic_json(receipt,signature)
    return target

def verify_parent(parent, task, plan):
    meta=json.loads((parent/"metadata.json").read_text())
    endpoint=next(r for r in efast.read_table(Path(plan["result_root"])/"seed_selection.tsv")
                  if r["fit_seed"]==task["fit_seed"])
    expected=dict(fit_seed=task["fit_seed"],N=int(task["N"]),M=int(task["M"]),replicate=int(task["replicate"]),
        seed=int(task["phase_seed"]),n_samples=int(task["N"])*14,n_parameters=14,
        salib_version="1.5.2",bounds_sha256=endpoint["bounds_sha256"],
        best_params_sha256=endpoint["best_params_sha256"],fit_config_sha256=endpoint["fit_config_sha256"])
    if any(meta.get(k)!=v for k,v in expected.items()) or not np.allclose(meta["oxygen_pct"],np.arange(201)*.025,rtol=0,atol=1e-12):
        raise ValueError("Parent design identity/grid differs from submitted task")

def fallback(root, part):
    """Recover only negative-vector points; all other failures stop the task."""
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
                before[:,j,k]=orders(y,metadata["N"],metadata["M"])
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
        code_sha256={name:efast.sha256(HERE/name) for name in ("resolution_nm.py","neighborhood_slurm.py","repair_eigen_points_safe.R","perron_high_precision_safe.cpp")},
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
        cache_sha256=efast.sha256(part/"indices.npz"),analyzer_sha256=efast.sha256(HERE/"resolution_nm.py"),
        salib_version=efast.salib_version()))
    return result

def worker(root, task_id):
    require_smoke(root)
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
    if receipt["analyzer_sha256"]!=efast.sha256(HERE/"resolution_nm.py"):
        raise ValueError("Trajectory analyzer changed")
    with np.load(part/"indices.npz") as data:
        return {k:data[k] for k in data.files}

def seed_summary(root, summary_id):
    plan=context(root)
    summaries=efast.read_table(root/"slurm/summary_manifest.tsv")
    if not 1<=summary_id<=len(summaries):
        raise ValueError("Summary ID outside pilot manifest")
    summary=summaries[summary_id-1]
    seed,n,m=summary["fit_seed"],int(summary["N"]),int(summary["M"])
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
        path=destination/("%s_N%d_M%d.npz"%(seed,n,m))
        save_npz(path,indices=result,variance=variances,oxygen=oxygen,parameters=efast.ACTIVE,
                 outputs=efast.OUTPUTS,index_names=nh.INDICES,fit_seed=seed,N=n,M=m)
        rows=rp.repeat_diagnostics(seed,n,result,oxygen)
        for row in rows: row.update(M=m,config=plan["config"])
        table=destination/("%s_N%d_M%d_repeat_diagnostics.tsv"%(seed,n,m))
        efast.write_table(table,list(rows[0]),rows)
        nh.atomic_json(destination/("%s_N%d_M%d.receipt.json"%(seed,n,m)),dict(status="complete",
            fit_seed=seed,N=n,M=m,config=plan["config"],n_trajectories=len(tasks),indices_sha256=efast.sha256(path),
            diagnostics_sha256=efast.sha256(table),completed_utc=ns.now(),sources=sources,
            numerical_recovery=recovery_rows))
        print("Completed %s N%d summary from 70 validated trajectory caches"%(seed,n),flush=True)


def comparison_rows(seed, spec, low_values, new_values, oxygen):
    """A change in M is explicitly a setting comparison, never an N step."""
    rows = rp.resolution_diagnostics(seed,4097,spec['N'],low_values,new_values,oxygen)
    for row in rows:
        row.update(comparison_kind=spec['comparison_kind'], N_reference=row.pop('N_low'),
                   N_test=row.pop('N_high'), M_reference=4, M_test=spec['M'],
                   primary_comparison=True)
    return rows


def plot_comparison(root, spec, values, oxygen):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    destination=root/'figures'
    destination.mkdir(exist_ok=True)
    configs=((2049,4),(4097,4),(spec['N'],spec['M']))
    colors=('#999999','#E6550D','#087E8B')
    for output,k in (('ploidy_ST',0),('growth_ST',1)):
        fig,axes=plt.subplots(2,5,figsize=(17,6.5),sharex=True,sharey=True)
        for i,seed in enumerate(spec['seeds']):
            for j,parameter in enumerate(rp.KEY_PARAMETERS):
                ax=axes[i,j]
                for (n,m),color in zip(configs,colors):
                    curves=values[(seed,n,m)][:,1,efast.ACTIVE.index(parameter),:,k]
                    mean,_,span,count=rp.phase_stats(curves)
                    valid=count==5
                    low=np.min(np.where(np.isfinite(curves),curves,np.inf),axis=0)
                    high=np.max(np.where(np.isfinite(curves),curves,-np.inf),axis=0)
                    ax.plot(oxygen,np.where(valid,mean,np.nan),lw=1,color=color,label='N=%d, M=%d'%(n,m))
                    ax.fill_between(oxygen,low,high,where=valid,color=color,alpha=.12,lw=0)
                ax.axvspan(0,1,color='#CFE2F3',alpha=.18)
                ax.axvspan(3,5,color='#FCE5CD',alpha=.18)
                ax.set(xlim=(0,5),ylim=(0,1.02))
                if i==0: ax.set_title(parameter)
                if j==0: ax.set_ylabel(seed+'\nST')
                if i==1: ax.set_xlabel('Oxygen (%)')
        axes[0,0].legend(fontsize=7)
        kind='N resolution at fixed M=4' if spec['M']==4 else 'M setting robustness at fixed N=4097'
        fig.suptitle(output.replace('_',' ')+': '+kind+'\nFive-phase means and full phase ranges; bands are not confidence intervals',fontsize=12)
        fig.tight_layout(rect=(0,0,1,.91))
        for ext in ('pdf','png'):
            fig.savefig(destination/(output+'_comparison.'+ext),dpi=150)
        plt.close(fig)


def status(root):
    plan=context(root)
    states={'complete':[],'failed':[],'missing_or_running':[]}
    for i in range(1,141):
        path=root/'slurm/task_receipts'/('task_%05d.json'%i)
        state=json.loads(path.read_text()).get('status') if path.exists() else None
        states[state if state in ('complete','failed') else 'missing_or_running'].append(i)
    summaries=[]
    for path in (root/'summaries').glob('*.receipt.json'):
        receipt=json.loads(path.read_text())
        if receipt.get('status')=='complete': summaries.append(path.name)
    report=dict(checked_utc=ns.now(),config=plan['config'],total_tasks=140,
        completed_tasks=len(states['complete']),failed_tasks=len(states['failed']),
        missing_or_running_tasks=len(states['missing_or_running']),
        failed_task_ids=states['failed'],missing_task_ids=states['missing_or_running'],
        completed_seed_summaries=len(summaries),convergence_policy='diagnostic_only')
    nh.atomic_json(root/'slurm/status.json',report)
    print(json.dumps({k:v for k,v in report.items() if not k.endswith('_ids')},indent=2),flush=True)
    return report


def finalize(root):
    plan=context(root,verify_sources=True)
    spec=CONFIGS[plan['config']]
    with ns.locked(root/'slurm/finalize.lock'):
        state=status(root)
        if state['completed_tasks']!=140 or state['completed_seed_summaries']!=2:
            raise ValueError('Require all 140 trajectories and two validated summaries')
        values={}
        with np.load(root/'baseline_indices.npz') as data:
            oxygen=data['oxygen'].copy()
            for i,seed in enumerate(spec['seeds']):
                for j,n in enumerate((2049,4097)):
                    values[(seed,n,4)]=data['indices'][i,j].copy()
        for seed in spec['seeds']:
            prefix=root/'summaries'/('%s_N%d_M%d'%(seed,spec['N'],spec['M']))
            receipt=json.loads(prefix.with_suffix('.receipt.json').read_text())
            if (receipt['status']!='complete' or receipt['n_trajectories']!=70 or
                receipt['M']!=spec['M'] or receipt['N']!=spec['N']):
                raise ValueError('Incomplete or wrong-configuration summary')
            for path,field in ((prefix.with_suffix('.npz'),'indices_sha256'),
                (prefix.with_name(prefix.name+'_repeat_diagnostics.tsv'),'diagnostics_sha256')):
                if efast.sha256(path)!=receipt[field]: raise ValueError('Summary hash differs')
            with np.load(prefix.with_suffix('.npz')) as data:
                if int(data['M'])!=spec['M']: raise ValueError('Wrong summary M')
                values[(seed,spec['N'],spec['M'])]=previous.checked_indices(data,seed,spec['N'])
        repeats,deltas,bands,curves=[],[],[],[]
        configs=((2049,4),(4097,4),(spec['N'],spec['M']))
        for seed in spec['seeds']:
            for n,m in configs:
                v=values[(seed,n,m)]
                r=rp.repeat_diagnostics(seed,n,v,oxygen)
                b=rp.comparison_rows(seed,n,v,oxygen)
                for row in r+b: row.update(M=m,config='N%d_M%d'%(n,m))
                repeats.extend(r)
                bands.extend(b)
                mean,sd,span,count=rp.phase_stats(v)
                for t,index in enumerate(nh.INDICES):
                    for p,parameter in enumerate(efast.ACTIVE):
                        for j,o in enumerate(oxygen):
                            for k,output in enumerate(efast.OUTPUTS):
                                curves.append(dict(fit_seed=seed,N=n,M=m,index=index,
                                    parameter=parameter,O2_pct=float(o),output=output,
                                    phase_mean=float(mean[t,p,j,k]),phase_sd=float(sd[t,p,j,k]),
                                    phase_range=float(span[t,p,j,k]),n_defined_phases=int(count[t,p,j,k])))
            deltas.extend(comparison_rows(seed,spec,values[(seed,4097,4)],
                values[(seed,spec['N'],spec['M'])],oxygen))
        for name,rows in (('repeat_diagnostics.tsv',repeats),('configuration_comparison.tsv',deltas),
                          ('band_rank_mechanism_comparison.tsv',bands),('oxygen_curve_summary.tsv.gz',curves)):
            write_rows(root/name,rows)
        for name,rows in (('key_ploidy_ST_repeat_diagnostics.tsv',repeats),
                          ('key_ploidy_ST_configuration_comparison.tsv',deltas)):
            keep=[r for r in rows if r['key_parameter'] and r['index']=='ST' and r['output']==efast.OUTPUTS[0]]
            write_rows(root/name,keep)
        plot_comparison(root,spec,values,oxygen)
        nh.atomic_json(root/'completion.json',dict(status='complete',completed_utc=ns.now(),
            config=plan['config'],comparison_kind=spec['comparison_kind'],
            n_tasks=140,n_seed_summaries=2,diagnostic_only=True,
            interpretation='Diagnostic endpoints, not a population sample. M comparisons are not N convergence. Phase ranges are not confidence intervals.',
            artifact_sha256={str(p.relative_to(root)):efast.sha256(p) for p in sorted(root.rglob('*'))
                            if p.is_file() and (p.parent==root or p.parent==root/'figures') and p.name!='completion.json'}))
        print('Completed '+plan['config']+' diagnostic tables and comparative figures',flush=True)


def smoke(root):
    plan=context(root,verify_sources=True)
    spec=CONFIGS[plan['config']]
    subprocess.run([sys.executable,'-m','unittest','discover','-s',str(HERE),
                    '-p','test_resolution_nm.py','-v'],check=True)
    # Both selected endpoints get a full-N trajectory on three grid points.
    checks=[]
    for seed in spec['seeds']:
        p=efast.ACTIVE.index('buffer_beta' if spec['M']==4 else 'p_mis_base')
        target=ensure_design(root,seed,spec['N'],plan,smoke=True)
        parent=target/'runs'/('N%d_R1'%spec['N'])
        meta=json.loads((parent/'metadata.json').read_text())
        if meta['M']!=spec['M'] or meta['N']!=spec['N'] or meta['oxygen_pct']!=[0.,2.5,5.]:
            raise ValueError('Smoke metadata differs')
        part=ns.split_trajectory(parent,p)
        begin=time.monotonic()
        result=evaluate_part(root,part,plan)
        if result['indices'].shape!=(2,3,2) or not np.isfinite(result['indices']).all():
            raise ValueError('Smoke produced undefined or unexpected indices')
        checks.append(dict(fit_seed=seed,N=spec['N'],M=spec['M'],focal_parameter=efast.ACTIVE[p],
            phase_seed=meta['seed'],n_rows=spec['N']*3,oxygen_pct=meta['oxygen_pct'],
            elapsed_seconds=time.monotonic()-begin,output_sha256=efast.sha256(part/'outputs.tsv.gz')))
    subprocess.run(['Rscript','--vanilla',str(HERE/'resolution_pilot_smoke.R'),str(root),str(part)],check=True)
    proof=json.loads((root/'smoke/perron_validation.json').read_text())
    if proof['status']!='passed': raise ValueError('Smoke Perron validation failed')
    nh.atomic_json(root/'smoke/smoke_receipt.json',dict(status='passed',completed_utc=ns.now(),
        config=plan['config'],N=spec['N'],M=spec['M'],
        submission_plan_sha256=efast.sha256(root/'slurm/submission_plan.json'),
        unit_tests='test_resolution_nm.py',trajectory_checks=checks,perron_validation=proof,
        host=os.uname().nodename,slurm_job_id=os.environ.get('SLURM_JOB_ID')))
    print('Passed SIF tests, both endpoint trajectories and 50/100-digit recovery interface',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','worker','seed-summary','status','finalize','smoke'))
    parser.add_argument('--out-dir',required=True)
    parser.add_argument('--baseline-dir')
    parser.add_argument('--config',choices=tuple(CONFIGS))
    parser.add_argument('--task-id',type=int)
    parser.add_argument('--summary-id',type=int)
    args=parser.parse_args()
    root=Path(args.out_dir).resolve()
    if args.command=='prepare':
        if not args.config or not args.baseline_dir: parser.error('prepare requires config and baseline-dir')
        prepare(args)
    elif args.command=='worker': worker(root,args.task_id)
    elif args.command=='seed-summary': seed_summary(root,args.summary_id)
    elif args.command=='status': status(root)
    elif args.command=='smoke': smoke(root)
    else: finalize(root)


if __name__=='__main__': main()
