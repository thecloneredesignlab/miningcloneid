#!/usr/bin/env python3
"""Small numerical/provenance checks; runnable in the production SALib SIF."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import efast
import neighborhood as nh


class NeighborhoodChecks(unittest.TestCase):
    def test_natural_span_before_log_transform_and_clipping(self):
        ranges = [dict(parameter="x", natural_lower=1., natural_upper=101., transform="log10"),
                  dict(parameter="y", natural_lower=0., natural_upper=1., transform="identity")]
        rows = nh.neighborhood_bounds(ranges, dict(x=2., y=.5))
        self.assertEqual((rows[0]["natural_lower"], rows[0]["natural_upper"]), (1., 12.))
        self.assertEqual(rows[0]["lower_clipped"], "TRUE")
        self.assertAlmostEqual(rows[0]["encoded_upper"], np.log10(12.))
        self.assertEqual((rows[1]["natural_lower"], rows[1]["natural_upper"]), (.4, .6))
        with self.assertRaises(ValueError):
            nh.neighborhood_bounds(ranges, dict(x=200., y=.5))

    def test_five_reproducible_independent_phase_designs(self):
        from SALib.sample import fast_sampler
        problem = dict(num_vars=14, names=list(efast.ACTIVE), bounds=[[0., 1.]] * 14)
        samples = [fast_sampler.sample(problem, 129, seed=1000 + i) for i in range(5)]
        self.assertTrue(all(x.shape == (14 * 129, 14) for x in samples))
        np.testing.assert_array_equal(samples[0], fast_sampler.sample(problem, 129, seed=1000))
        self.assertEqual(len({x.tobytes() for x in samples}), 5)

    def test_resume_preserves_designs_and_refuses_changed_seed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fit, figure, output = root / "fit", root / "figure", root / "out"
            fit.mkdir()
            figure.mkdir()
            (fit / "parameter_table.csv").write_text(
                "param_prototype,param_name,estimate,lower_bound,upper_bound\n" +
                "".join("%s,%s,TRUE,0,1\n" % (p, p) for p in efast.ACTIVE))
            efast.write_table(figure / "fixed_o2_dominant_ploidy_201grid.tsv", ["O2_pct"],
                              [dict(O2_pct=o) for o in np.linspace(0, 5, 201)])
            for name in ("continuous_ploidy_spearman_by_o2.tsv", "parameter_function_groups.tsv"):
                (figure / name).write_text("parameter\nlam_max\n")
            args = SimpleNamespace(fit_root=str(fit), figure4_dir=str(figure), out_dir=str(output),
                n=129, m=4, replicates=5, seed=123, full_grid=True, oxygen="", fit_seed="seed25")
            efast.prepare(args)
            hashes = [efast.sha256(p) for p in sorted(output.glob("runs/*/samples.tsv.gz"))]
            efast.prepare(args)
            self.assertEqual(hashes, [efast.sha256(p) for p in sorted(output.glob("runs/*/samples.tsv.gz"))])
            args.seed += 1
            with self.assertRaisesRegex(ValueError, "metadata differs"):
                efast.prepare(args)

    def test_streaming_indices_match_salib_and_reject_missing_rows(self):
        from SALib.sample import fast_sampler
        from SALib.analyze import fast
        problem = dict(num_vars=14, names=list(efast.ACTIVE), bounds=[[0., 1.]] * 14)
        n = 513
        x = fast_sampler.sample(problem, n, seed=123)
        y = x[:, 0] + 2 * x[:, 1] + 3 * x[:, 0] * x[:, 2]
        expected = fast.analyze(problem, y, seed=123)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = dict(N=n, M=4, n_samples=len(x), oxygen_pct=[0, 5], seed=123)
            (root / "metadata.json").write_text(json.dumps(metadata))
            rows = [dict(sample_id=i + 1, O2_pct=o, dominant_mean_ploidy=float(value),
                         dominant_growth_rate=1., spectral_gap=.01,
                         eigenvector_nonnegative="TRUE", status="ok")
                    for i, value in enumerate(y) for o in (0, 5)]
            efast.write_table(root / "outputs.tsv.gz", list(rows[0]), rows)
            result = nh.analyze_design(root)
            np.testing.assert_allclose(result["indices"][0, :, 0, 0], expected["S1"], atol=1e-14)
            np.testing.assert_allclose(result["indices"][1, :, 0, 0], expected["ST"], atol=1e-14)
            self.assertTrue(np.isnan(result["indices"][:, :, :, 1]).all())
            # Known additive contributions are recovered at adequate resolution.
            analytic = fast.analyze(problem, x[:, 0] + 2 * x[:, 1], seed=123)
            self.assertAlmostEqual(analytic["S1"][0], .2, delta=.025)
            self.assertAlmostEqual(analytic["S1"][1], .8, delta=.025)
            efast.write_table(root / "outputs.tsv.gz", list(rows[0]), rows[:-1])
            with self.assertRaisesRegex(ValueError, "Incomplete"):
                nh.analyze_design(root)

    def test_seed_summary_averages_five_phases_and_preserves_undefined(self):
        results = []
        for value in (.1, .2, .3, .4, .5):
            indices = np.full((2, 14, 2, 2), value)
            if value == .5:
                indices[0, 0, 0, 0] = np.nan
            results.append(dict(indices=indices, oxygen=np.array([0, 5]),
                gap_counts=np.array([[10, 0, 0, 0], [0, 0, 0, 0], [10, 0, 0, 0]]),
                gap_min=np.array([.01, np.inf, .01])))
        with tempfile.TemporaryDirectory() as temporary, patch.object(nh, "analyze_design", side_effect=results):
            rows = nh.seed_summary(Path(temporary), [513])
            affected = next(r for r in rows if r["output"] == "dominant_mean_ploidy" and
                            r["parameter"] == "lam_max" and r["O2_pct"] == 0)
            self.assertEqual(affected["S1_n_valid_repeats"], 4)
            self.assertTrue(np.isnan(affected["S1_mean"]))
            self.assertAlmostEqual(affected["ST_mean"], .3)
            self.assertAlmostEqual(affected["ST_range"], .4)

    def test_diagnostics_distinguish_failed_passed_and_untested_resolution(self):
        rows = [dict(N=n, output=output, parameter=parameter, O2_pct=o,
                     S1_range=.08, ST_range=.08, S1_delta_vs_max_N=.02, ST_delta_vs_max_N=.02)
                for n in (257, 513) for output in efast.OUTPUTS
                for parameter in efast.ACTIVE for o in (0, 5)]
        reports = nh.convergence_diagnostic_rows("seed1", rows, rows)
        self.assertTrue(all(r["convergence_status"] == "passed" for r in reports))
        untested = nh.convergence_diagnostic_rows("seed1", rows)
        self.assertTrue(all(r["convergence_status"] == "resolution_not_tested" for r in untested))
        self.assertTrue(all(np.isnan(r["resolution_delta_p90"]) for r in untested))
        for row in rows:
            row["ST_range"] = .5
            row["S1_delta_vs_max_N"] = .2
        reports = nh.convergence_diagnostic_rows("seed1", rows, rows)
        self.assertTrue(all(r["convergence_status"] == "not_passed" for r in reports))
        self.assertTrue(all(r["included_in_full_analysis"] == "TRUE" for r in reports))
        for row in rows:
            row["S1_range"] = np.nan
        reports = nh.convergence_diagnostic_rows("seed1", rows)
        self.assertTrue(all(r["repeat_status"] == "undefined" for r in reports if r["index"] == "S1"))

    def test_full_policy_accepts_failing_or_absent_pilot_but_checks_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "input_manifest.json").write_text('{}\n')
            nh.record_full_execution_policy(root)
            policy = json.loads((root / "execution_policy.json").read_text())
            self.assertFalse(policy["require_convergence_pass"])
            self.assertEqual(policy["pilot_status"], "not_available")
            report = root / "pilot" / "convergence_gate.json"
            nh.atomic_json(report, dict(status="needs_review",
                input_manifest_sha256=efast.sha256(root / "input_manifest.json")))
            nh.record_full_execution_policy(root)
            self.assertEqual(json.loads((root / "execution_policy.json").read_text())["pilot_status"], "needs_review")
            (root / "input_manifest.json").write_text('{"changed":true}\n')
            with self.assertRaisesRegex(ValueError, "inputs differ"):
                nh.record_full_execution_policy(root)

    def test_incremental_diagnostics_keep_all_completed_endpoints(self):
        rows = [dict(N=513, output=output, parameter=parameter, O2_pct=o,
                     S1_range=.08, ST_range=.5)
                for output in efast.OUTPUTS for parameter in efast.ACTIVE for o in (0, 5)]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for seed in ("seed10", "seed2"):
                nh.write_full_diagnostics(root, seed, rows)
            statuses = efast.read_table(root / "summaries" / "seed_convergence_status.tsv")
            self.assertEqual([r["fit_seed"] for r in statuses], ["seed2", "seed10"])
            self.assertTrue(all(r["convergence_status"] == "not_passed" and
                                r["included_in_full_analysis"] == "TRUE" for r in statuses))
            diagnostics = efast.read_table(root / "summaries" / "seed_convergence_diagnostics.tsv")
            self.assertEqual(len(diagnostics), 8)


if __name__ == "__main__":
    unittest.main()
