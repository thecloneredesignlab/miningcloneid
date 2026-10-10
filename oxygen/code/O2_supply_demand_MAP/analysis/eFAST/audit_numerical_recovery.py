#!/usr/bin/env python3
"""Collect compact, source-backed numerical-recovery diagnostics after completion."""
import argparse
import json
from pathlib import Path

import numpy as np
import efast
import neighborhood as nh
import neighborhood_slurm as ns


def audit(root):
    ids=set()
    for path in (root/"slurm").glob("numerical_recovery_plan*.json"):
        ids.update(json.loads(path.read_text())["task_ids"])
    points, impacts, n_patched=[],[],0
    for task_id in sorted(ids):
        task=ns.task_identity(task_id)
        part=root/"runs"/task["fit_seed"]/"runs"/("N513_R%d"%task["replicate"])/"trajectories"/("P%02d_%s"%(task["parameter_index"]+1,task["focal_parameter"]))
        path=part/"recovery/numerical_recovery.json"
        if not path.exists(): continue
        proof=json.loads(path.read_text())
        receipt=json.loads((part/"outputs.tsv.gz.receipt.json").read_text())
        assert receipt["recovery_manifest_sha256"]==efast.sha256(path)
        assert proof["status"]=="passed" and proof["no_matrix_regularization"] and proof["no_sampling_change"]
        corrections=part/"recovery/corrections.json"
        assert proof["artifact_sha256"]["corrections.json"]==efast.sha256(corrections)
        records=json.loads(corrections.read_text())["reports"]
        assert len(records)==proof["n_recomputed_points"]
        n_patched+=1
        for record in records:
            v50=np.asarray(record["vector50"]); v100=np.asarray(record["vector100"])
            assert np.all(v50>=0) and np.all(v100>=0)
            assert abs(v50.sum()-1)<1e-12 and np.sum(np.abs(v50-v100))<1e-10
            assert record["ploidy_precision_difference"]<1e-12 and record["growth_precision_difference"]<1e-12
            assert record["double_residual"]<1e-12 and 0<=record["leading_root_bound_width"]<1e-12
            original=record["original"]
            points.append(dict(task_id=task_id,fit_seed=task["fit_seed"],replicate=task["replicate"],
                focal_parameter=task["focal_parameter"],sample_id=record["sample_id"],O2_pct=record["O2_pct"],
                original_ploidy=original["dominant_mean_ploidy"],corrected_ploidy=record["dominant_mean_ploidy"],
                ploidy_change=record["ploidy_change"],original_growth=original["dominant_growth_rate"],
                corrected_growth=record["dominant_growth_rate"],growth_change=record["growth_change"],
                ploidy_precision_difference=record["ploidy_precision_difference"],
                growth_precision_difference=record["growth_precision_difference"],
                vector_precision_difference=record["vector_precision_difference"],
                double_residual=record["double_residual"],leading_root_bound_width=record["leading_root_bound_width"],
                original_spectral_gap=original["spectral_gap"],iterations50=record["iterations50"],iterations100=record["iterations100"]))
        for row in proof["index_impact"]:
            impacts.append(dict(task_id=task_id,fit_seed=task["fit_seed"],replicate=task["replicate"],
                focal_parameter=task["focal_parameter"],**row))
    assert points and impacts
    efast.write_table(root/"slurm/numerical_recovery_points.tsv",list(points[0]),points)
    efast.write_table(root/"slurm/numerical_recovery_index_impact.tsv",list(impacts[0]),impacts)
    status=json.loads((root/"slurm/status.json").read_text())
    report=dict(status="passed",checked_utc=ns.now(),n_patched_trajectories=n_patched,n_recomputed_points=len(points),
        n_recovery_task_ids=len(ids),max_absolute_ploidy_change=max(abs(r["ploidy_change"]) for r in points),
        max_absolute_growth_change=max(abs(r["growth_change"]) for r in points),
        max_ploidy_precision_difference=max(r["ploidy_precision_difference"] for r in points),
        max_growth_precision_difference=max(r["growth_precision_difference"] for r in points),
        max_double_residual=max(r["double_residual"] for r in points),
        max_leading_root_bound_width=max(r["leading_root_bound_width"] for r in points),
        n_index_undefined_pattern_changes=sum(r["undefined_pattern_changed"] for r in impacts),
        index_impact=[dict(index=t,output=o,max_absolute_index_difference=max(
            r["max_absolute_index_difference"] for r in impacts if r["index"]==t and r["output"]==o and r["max_absolute_index_difference"] is not None))
            for t in nh.INDICES for o in efast.OUTPUTS],
        completed_tasks=status["completed_tasks"],completed_seed_summaries=status["completed_seed_summaries"],
        model_and_samples_changed=False,spectral_gap_source="original canonical double-precision diagnostic")
    nh.atomic_json(root/"slurm/numerical_recovery_audit.json",report)
    print(json.dumps(report,indent=2),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir",required=True)
    audit(Path(parser.parse_args().out_dir).resolve())
