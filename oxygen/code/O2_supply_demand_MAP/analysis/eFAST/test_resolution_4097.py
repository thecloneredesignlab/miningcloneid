#!/usr/bin/env python3
"""N4097 manifest, Fourier estimator and cross-resolution diagnostic checks."""
import importlib.util
import unittest

import numpy as np
import efast
import neighborhood as nh
import resolution_pilot as rp
import resolution_4097 as run


class N4097Checks(unittest.TestCase):
    def test_complete_trajectory_manifest_and_inherited_phase_seeds(self):
        tasks, summaries = run.task_rows()
        self.assertEqual([r["fit_seed"] for r in summaries], list(run.SEEDS))
        self.assertEqual(len(tasks), 420)
        self.assertEqual([r["task_id"] for r in tasks], list(range(1, 421)))
        self.assertEqual(sum(r["N"]*201 for r in tasks), 345868740)
        self.assertEqual(len({(r["fit_seed"], r["replicate"], r["parameter_index"]) for r in tasks}), 420)
        old, _ = rp.task_rows(run.SEEDS)
        phases = {(r["fit_seed"], r["replicate"]): r["phase_seed"] for r in old}
        for r in tasks:
            self.assertEqual(r["phase_seed"], phases[(r["fit_seed"], r["replicate"])])
            self.assertEqual(r["parent_sample_id_end"]-r["parent_sample_id_start"]+1, 4097)
        for summary in summaries:
            rows = [r for r in tasks if r["summary_id"] == summary["summary_id"]]
            self.assertEqual(len(rows), 70)
            self.assertEqual({r["focal_parameter"] for r in rows}, set(efast.ACTIVE))
            self.assertEqual({r["replicate"] for r in rows}, {1, 2, 3, 4, 5})

    def test_primary_pair_and_localized_error_not_hidden_by_p90(self):
        oxygen = np.arange(201)*.025
        a = np.full((5, 2, 14, 201, 2), .2)
        b = a.copy()
        b[:, 1, 0, 77, 0] += .3
        rows = run.resolution_rows("seed138", 2049, 4097, a, b, oxygen)
        row = next(r for r in rows if r["parameter"] == efast.ACTIVE[0] and r["index"] == "ST"
                   and r["output"] == efast.OUTPUTS[0] and r["oxygen_window"] == rp.WINDOWS[0][0])
        self.assertTrue(row["primary_comparison"])
        self.assertTrue(row["p90_pass"])
        self.assertEqual(row["diagnostic_status"], "localized_review")
        self.assertAlmostEqual(row["absolute_delta_max_O2_pct"], 1.925)
        older = run.resolution_rows("seed138", 1025, 2049, a, b, oxygen)
        self.assertFalse(any(r["primary_comparison"] for r in older))
        b[0, 1, 0, 0, 0] = np.nan
        rows = run.resolution_rows("seed138", 2049, 4097, a, b, oxygen)
        row = next(r for r in rows if r["parameter"] == efast.ACTIVE[0] and r["index"] == "ST"
                   and r["output"] == efast.OUTPUTS[0] and r["oxygen_window"] == rp.WINDOWS[1][0])
        self.assertEqual(row["diagnostic_status"], "undefined_or_incomplete")
        self.assertFalse(row["p90_pass"])

    def test_baseline_endpoint_axes_and_phase_count_are_checked(self):
        data = dict(parameters=np.array(efast.ACTIVE), outputs=np.array(efast.OUTPUTS),
                    index_names=np.array(nh.INDICES), oxygen=np.arange(201)*.025,
                    fit_seed=np.array("seed25"), N=np.array(2049), indices=np.zeros((5, 2, 14, 201, 2)))
        self.assertEqual(run.checked_indices(data, "seed25", 2049).shape, (5, 2, 14, 201, 2))
        with self.assertRaisesRegex(ValueError, "endpoint/resolution"):
            run.checked_indices(data, "seed101", 2049)
        data["indices"] = data["indices"][:4]
        with self.assertRaisesRegex(ValueError, "five phases"):
            run.checked_indices(data, "seed25", 2049)

    def test_mechanism_stability_includes_new_pair_with_original_labels(self):
        oxygen = np.arange(201)*.025
        values = np.full((5, 2, 14, 201, 2), .1)
        bands = []
        for seed in run.SEEDS:
            for n in (*run.BASELINE_NS, run.N):
                v = values.copy()
                v[:, :, 0, 120:, :] = .3 if n < run.N else .01
                bands.extend(rp.comparison_rows(seed, n, v, oxygen))
        rows = run.stability_rows(bands)
        self.assertEqual(len(rows), 6*3*2*2*14)
        self.assertEqual({(r["N_low"], r["N_high"]) for r in rows}, set(run.PAIRS))
        row = next(r for r in rows if r["fit_seed"] == "seed25" and r["N_high"] == 4097
                   and r["parameter"] == efast.ACTIVE[0] and r["output"] == efast.OUTPUTS[0] and r["index"] == "ST")
        self.assertFalse(row["contrast_direction_unchanged"])
        self.assertLess(row["high_minus_low_at_higher_N"], 0)

    @unittest.skipUnless(importlib.util.find_spec("SALib"), "Run inside the specified SALib 1.5.2 SIF")
    def test_N4097_fft_and_full_design_match_known_variance_fractions(self):
        from SALib.sample.fast_sampler import sample
        from SALib.analyze.fast import compute_orders
        n, m = run.N, 4
        omega = (n-1)//(2*m)
        s = 2*np.pi*np.arange(n)/n
        y = np.sin(omega*s)+.5*np.sin(s)+.3*np.sin(2*omega*s)
        np.testing.assert_allclose(compute_orders(y, n, m, omega), [1.09/1.34]*2, rtol=0, atol=1e-12)
        problem = dict(num_vars=14, names=list(efast.ACTIVE), bounds=[[0., 1.]]*14)
        for repeat in range(5):
            x = sample(problem, n, M=m, seed=202610080+1380+repeat)
            self.assertEqual(x.shape, (14*n, 14))
            y = x[:, 0]+2*x[:, 1]
            for p, expected in ((0, .2), (1, .8), (13, 0.)):
                np.testing.assert_allclose(compute_orders(y[p*n:(p+1)*n], n, m, omega),
                                           [expected]*2, rtol=0, atol=.01)


if __name__ == "__main__": unittest.main()
