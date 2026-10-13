#!/usr/bin/env python3
"""Resolution, complete-trajectory, undefined-index and recovery checks."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import efast
import neighborhood as nh
import neighborhood_slurm as ns
import resolution_pilot as rp

HAS_SALIB = importlib.util.find_spec("SALib") is not None


class PilotChecks(unittest.TestCase):
    def test_manifest_covers_12_seeds_two_resolutions_five_phases_all_parameters(self):
        tasks,summaries=rp.task_rows([*rp.ORIGINAL_SEEDS,"seed174","seed494"])
        self.assertEqual(len(tasks),1680)
        self.assertEqual(len(summaries),24)
        self.assertEqual([t["task_id"] for t in tasks],list(range(1,1681)))
        self.assertEqual(sum(t["N"]*201 for t in tasks),519014160)
        self.assertEqual(len({(t["fit_seed"],t["N"],t["replicate"],t["focal_parameter"]) for t in tasks}),1680)
        self.assertEqual(len({t["phase_seed"] for t in tasks}),60)
        self.assertEqual(tasks[69]["N"],1025)
        self.assertEqual(tasks[70]["N"],2049)
        for summary in summaries:
            members=[t for t in tasks if t["summary_id"]==summary["summary_id"]]
            self.assertEqual(len(members),70)
            self.assertEqual({m["fit_seed"] for m in members},{summary["fit_seed"]})
            self.assertEqual({m["N"] for m in members},{summary["N"]})

    def test_repeat_diagnostics_detect_parameter_local_and_window_failures(self):
        oxygen=np.arange(201)*.025
        values=np.full((5,2,14,201,2),.2)
        values[0,1,efast.ACTIVE.index("buffer_beta"),0,0]=.8
        rows=rp.repeat_diagnostics("seed1",1025,values,oxygen)
        selected=[r for r in rows if r["parameter"]=="buffer_beta" and r["index"]=="ST" and r["output"]==efast.OUTPUTS[0]]
        by_window={r["oxygen_window"]:r for r in selected}
        self.assertTrue(by_window["overall_0_to_5pct"]["p90_pass"])
        self.assertEqual(by_window["overall_0_to_5pct"]["diagnostic_status"],"localized_review")
        self.assertAlmostEqual(by_window["low_0_to_1pct"]["phase_range_max"],.6)
        self.assertEqual(by_window["low_0_to_1pct"]["phase_range_max_O2_pct"],0.)
        self.assertEqual(by_window["high_3_to_5pct"]["diagnostic_status"],"passed_by_threshold")

    def test_resolution_primary_comparison_and_missing_phase_are_explicit(self):
        oxygen=np.arange(201)*.025
        a=np.full((5,2,14,201,2),.2);b=a.copy()
        b[:,1,0,:41,0]+=.08
        rows=rp.resolution_diagnostics("seed1",1025,2049,a,b,oxygen)
        selected=[r for r in rows if r["parameter"]=="lam_max" and r["index"]=="ST" and r["output"]==efast.OUTPUTS[0]]
        self.assertTrue(all(r["primary_comparison"] for r in selected))
        self.assertEqual(selected[1]["diagnostic_status"],"not_passed")
        self.assertAlmostEqual(selected[1]["band_delta"],.08)
        b[0,1,0,0,0]=np.nan
        rows=rp.resolution_diagnostics("seed1",1025,2049,a,b,oxygen)
        row=next(r for r in rows if r["parameter"]=="lam_max" and r["index"]=="ST" and r["output"]==efast.OUTPUTS[0] and r["oxygen_window"]=="low_0_to_1pct")
        self.assertEqual(row["diagnostic_status"],"undefined_or_incomplete")
        self.assertEqual(row["n_valid_oxygen"],40)
        self.assertIsNone(row["band_delta"])

    def test_undefined_constant_indices_remain_nan(self):
        x=np.full((5,3),np.nan)
        mean,sd,span,count=rp.phase_stats(x)
        self.assertTrue(np.isnan(mean).all())
        self.assertTrue(np.isnan(sd).all())
        self.assertTrue(np.isnan(span).all())
        self.assertTrue((count==0).all())

    def test_band_weights_respect_physical_oxygen_grid(self):
        oxygen=np.arange(201)*.025
        self.assertAlmostEqual(float(rp.band_scores(oxygen,oxygen,0,1)),.5)
        self.assertAlmostEqual(float(rp.band_scores(oxygen,oxygen,3,5)),4.)
        curve=oxygen.copy();curve[160]=np.nan
        self.assertTrue(np.isnan(rp.band_scores(curve,oxygen,3,5)))

    def test_mechanism_and_direction_checks_compare_endpoint_not_pooled_seeds(self):
        oxygen=np.arange(201)*.025
        values=np.full((5,2,14,201,2),.01)
        for p,param in enumerate(efast.ACTIVE):
            if efast.GROUP[param]=="death":values[:,:,p,:41,:]=.5
            if efast.GROUP[param]=="missegregation":values[:,:,p,120:,:]=.3
            if efast.GROUP[param]=="buffering":values[:,:,p,120:,:]=.2
        bands=[]
        for n in (513,1025,2049):
            copy=values.copy()
            if n==2049:
                for p,param in enumerate(efast.ACTIVE):
                    if efast.GROUP[param]=="missegregation":copy[:,:,p,120:,:]=.001
            bands.extend(rp.comparison_rows("seed1",n,copy,oxygen))
        comparisons=rp.stability_rows(bands)
        row=next(r for r in comparisons if r["N_high"]==2049 and r["parameter"]=="p_mis_base" and r["index"]=="ST" and r["output"]==efast.OUTPUTS[0])
        self.assertTrue(row["noemi_support_at_lower_N"])
        self.assertFalse(row["noemi_support_at_higher_N"])
        self.assertFalse(row["noemi_support_unchanged"])
        self.assertFalse(row["contrast_direction_unchanged"])

    @unittest.skipUnless(HAS_SALIB,"SALib checks run in the required SIF")
    def test_fft_matches_analytic_frequency_power_at_both_resolutions(self):
        from SALib.analyze.fast import compute_orders
        for n in rp.RESOLUTIONS:
            s=2*np.pi*np.arange(n)/n;omega=(n-1)//8
            y=np.sin(omega*s)+.5*np.sin(s)+.3*np.sin(2*omega*s)
            s1,st=compute_orders(y,n,4,omega)
            self.assertAlmostEqual(s1,1.09/1.34,places=12)
            self.assertAlmostEqual(st,1.09/1.34,places=12)

    @unittest.skipUnless(HAS_SALIB,"SALib checks run in the required SIF")
    def test_complete_sampling_blocks_match_additive_model_and_keep_14_columns(self):
        from SALib.sample.fast_sampler import sample
        from SALib.analyze.fast import compute_orders
        problem=dict(num_vars=14,names=list(efast.ACTIVE),bounds=[[0.,1.]]*14)
        for n in rp.RESOLUTIONS:
            x=sample(problem,n,M=4,seed=202610080+4940)
            self.assertEqual(x.shape,(n*14,14))
            y=x[:,0]+2*x[:,1]
            for p,expected in ((0,.2),(1,.8),(13,0.)):
                s1,st=compute_orders(y[p*n:(p+1)*n],n,4,(n-1)//8)
                self.assertAlmostEqual(s1,expected,delta=.01)
                self.assertAlmostEqual(st,expected,delta=.01)
            # Split preserves the actual SALib ordering and all simultaneous inputs.
            with tempfile.TemporaryDirectory() as tmp:
                parent=Path(tmp)
                efast.write_table(parent/"samples.tsv.gz",["sample_id",*efast.ACTIVE],
                    (dict(sample_id=i+1,**dict(zip(efast.ACTIVE,v))) for i,v in enumerate(x)))
                nh.atomic_json(parent/"metadata.json",dict(N=n,M=4,n_samples=n*14,
                    n_parameters=14,fit_seed="seed494",method="eFAST",oxygen_pct=[0.,5.],
                    samples_sha256=efast.sha256(parent/"samples.tsv.gz")))
                part=ns.split_trajectory(parent,13)
                block=efast.read_table(part/"samples.tsv.gz")
                self.assertEqual([int(r["sample_id"]) for r in block],list(range(1,n+1)))
                np.testing.assert_array_equal([[float(r[k]) for k in efast.ACTIVE] for r in block],x[13*n:])

    @unittest.skipUnless(HAS_SALIB,"SALib checks run in the required SIF")
    def test_streamed_recovery_preserves_good_rows_and_proof_detects_tampering(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);part=root/"part";part.mkdir()
            n=65;oxygen=[0.,2.5,5.]
            efast.write_table(part/"samples.tsv.gz",["sample_id","lam_max"],
                (dict(sample_id=i+1,lam_max=i) for i in range(n)))
            nh.atomic_json(part/"metadata.json",dict(N=n,M=4,n_samples=n,fit_seed="seed1",
                oxygen_pct=oxygen,samples_sha256=efast.sha256(part/"samples.tsv.gz")))
            rows=[dict(sample_id=i+1,O2_pct=o,dominant_mean_ploidy=2+.1*np.sin(i),
                dominant_growth_rate=.01+.001*np.cos(i),spectral_gap=.0001,
                eigenvector_nonnegative="FALSE" if i==3 and o==2.5 else "TRUE",status="ok")
                for i in range(n) for o in oxygen]
            original=part/"outputs.tsv.gz.part01.gz"
            efast.write_table(original,ns.OUTPUT_FIELDS,rows)
            before=efast.read_table(original)
            def fake_R(*args,**kwargs):
                recovery=part/"recovery"
                bad=json.loads((recovery/"failed_points.json").read_text())
                matrix=recovery/"point0001_matrix.rds";matrix.write_bytes(b"mock-original-matrix")
                reports=[dict(original=r,sample_id=int(r["sample_id"]),O2_pct=float(r["O2_pct"]),
                    dominant_mean_ploidy=float(r["dominant_mean_ploidy"]),
                    dominant_growth_rate=float(r["dominant_growth_rate"]),
                    matrix_file=matrix.name,matrix_sha256=efast.sha256(matrix)) for r in bad]
                nh.atomic_json(recovery/"corrections.json",dict(status="passed",n_points=len(bad),reports=reports))
            with patch.object(rp.subprocess,"run",side_effect=fake_R):rp.fallback(root,part)
            self.assertTrue(nh.verify_completion(part))
            after=efast.read_table(part/"outputs.tsv.gz")
            for x,y in zip(before,after):
                if x["eigenvector_nonnegative"]=="TRUE":self.assertEqual(x,y)
            self.assertEqual(sum(r["eigenvector_nonnegative"]=="TRUE" for r in after),len(after))
            (part/"recovery/point0001_matrix.rds").write_bytes(b"tampered")
            with self.assertRaisesRegex(ValueError,"source artifact changed"):nh.verify_completion(part)
            efast.write_table(original,ns.OUTPUT_FIELDS,rows[:-1])
            with self.assertRaisesRegex(ValueError,"complete"):rp.fallback(root,part)


if __name__=="__main__":unittest.main()
