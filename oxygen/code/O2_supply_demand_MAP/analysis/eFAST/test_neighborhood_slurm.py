#!/usr/bin/env python3
"""Check trajectory completeness and equivalence without expensive model runs."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import efast
import neighborhood as nh
import neighborhood_slurm as ns


def receipt(directory, n_rows):
    metadata = json.loads((directory / "metadata.json").read_text())
    nh.atomic_json(directory / "outputs.tsv.gz.receipt.json", dict(fit_seed="seed1", n_rows=n_rows,
        metadata_sha256=efast.sha256(directory / "metadata.json"), samples_sha256=metadata["samples_sha256"],
        outputs_sha256=efast.sha256(directory / "outputs.tsv.gz"),
        evaluator_sha256=efast.sha256(ns.HERE / "evaluate_fixed_o2.R"), elapsed_seconds=1., workers=1))


class SlurmChecks(unittest.TestCase):
    def test_all_task_ids_cover_500_seeds_five_phases_fourteen_trajectories(self):
        rows = list(ns.task_rows())
        self.assertEqual(len(rows), 35000)
        self.assertEqual([r["task_id"] for r in rows], list(range(1, 35001)))
        for row in rows:
            self.assertEqual(row, ns.task_identity(row["task_id"]))
        self.assertEqual(len({(r["fit_seed"], r["replicate"], r["focal_parameter"]) for r in rows}), 35000)
        self.assertEqual(len({r["phase_seed"] for r in rows}), 2500)
        with self.assertRaises(ValueError):
            ns.task_identity(35001)

    def test_split_merge_matches_original_indices_variance_and_gap_qc(self):
        from SALib.sample import fast_sampler
        problem = dict(num_vars=14, names=list(efast.ACTIVE), bounds=[[0., 1.]] * 14)
        n = 129
        samples = fast_sampler.sample(problem, n, seed=123)
        y = samples[:, 0] + 2 * samples[:, 1] + 3 * samples[:, 0] * samples[:, 2]
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            sample_rows = [dict(sample_id=i + 1, **dict(zip(efast.ACTIVE, values))) for i, values in enumerate(samples)]
            efast.write_table(parent / "samples.tsv.gz", list(sample_rows[0]), sample_rows)
            metadata = dict(N=n, M=4, n_samples=len(samples), n_parameters=14,
                oxygen_pct=[0, 5], fit_seed="seed1", method="eFAST", samples_sha256=efast.sha256(parent / "samples.tsv.gz"))
            nh.atomic_json(parent / "metadata.json", metadata)
            outputs = [dict(sample_id=i + 1, O2_pct=o, dominant_mean_ploidy=value,
                dominant_growth_rate=1., spectral_gap=.01 if o==0 else .00001,
                eigenvector_nonnegative="TRUE", status="ok") for i, value in enumerate(y) for o in (0, 5)]
            efast.write_table(parent / "outputs.tsv.gz", ns.OUTPUT_FIELDS, outputs)
            receipt(parent, len(outputs))
            expected = nh.analyze_design(parent)
            (parent / "outputs.tsv.gz").unlink()
            (parent / "outputs.tsv.gz.receipt.json").unlink()
            for p in range(14):
                part = ns.split_trajectory(parent, p)
                block = efast.read_table(part / "samples.tsv.gz")
                np.testing.assert_allclose([[float(r[k]) for k in efast.ACTIVE] for r in block],
                    samples[p*n:(p+1)*n], rtol=0, atol=0)
                selected = [dict(row, sample_id=int(row["sample_id"])-p*n)
                    for row in outputs if p*n < int(row["sample_id"]) <= (p+1)*n]
                efast.write_table(part / "outputs.tsv.gz", ns.OUTPUT_FIELDS, selected)
                receipt(part, len(selected))
                result = ns.analyze_trajectory(part)
                np.testing.assert_allclose(result["indices"], expected["indices"][:, p], atol=1e-14, equal_nan=True)
                ns.write_trajectory_indices(part, result)
            ns.assemble_phase(parent)
            self.assertTrue(nh.verify_completion(parent))
            actual = nh.analyze_design(parent)
            for key in expected:
                np.testing.assert_allclose(actual[key], expected[key], atol=1e-14, equal_nan=True)
            restored = efast.read_table(parent / "outputs.tsv.gz")
            self.assertEqual([int(r["sample_id"]) for r in restored], [r["sample_id"] for r in outputs])
            # Cache corruption is detected before any raw assembly.
            (parent / "outputs.tsv.gz").unlink()
            (parent / "outputs.tsv.gz.receipt.json").unlink()
            first = parent / "trajectories" / "P01_lam_max" / "indices.npz"
            first.write_bytes(b"corrupted")
            with self.assertRaisesRegex(ValueError, "cache changed"):
                ns.assemble_phase(parent)

    def test_reject_incomplete_duplicate_and_nonvalid_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            part = Path(temporary)
            nh.atomic_json(part / "metadata.json", dict(N=65, M=4, oxygen_pct=[0, 5]))
            rows = [dict(sample_id=i+1, O2_pct=o, dominant_mean_ploidy=float(i),
                dominant_growth_rate=1., spectral_gap=.1, eigenvector_nonnegative="TRUE", status="ok")
                for i in range(65) for o in (0, 5)]
            efast.write_table(part / "outputs.tsv.gz", ns.OUTPUT_FIELDS, rows[:-1])
            with self.assertRaisesRegex(ValueError, "Incomplete"):
                ns.analyze_trajectory(part)
            efast.write_table(part / "outputs.tsv.gz", ns.OUTPUT_FIELDS, rows + rows[:1])
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                ns.analyze_trajectory(part)
            rows[0]["status"] = "failed"
            efast.write_table(part / "outputs.tsv.gz", ns.OUTPUT_FIELDS, rows)
            with self.assertRaisesRegex(ValueError, "Invalid trajectory"):
                ns.analyze_trajectory(part)


if __name__ == "__main__":
    unittest.main()
