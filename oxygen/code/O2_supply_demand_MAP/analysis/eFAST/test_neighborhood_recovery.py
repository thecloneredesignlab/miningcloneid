#!/usr/bin/env python3
"""Protect recovered-point provenance, targeting and unchanged scientific rows."""
import json
from pathlib import Path
import tempfile
import unittest

import efast
import neighborhood as nh
from recover_eigen import compressed_ids, invalid


class RecoveryChecks(unittest.TestCase):
    def test_target_ids_are_absolute_and_never_throttled(self):
        self.assertEqual(compressed_ids([27230,793,27228,794,795]),"793-795,27228,27230")
        self.assertNotIn("%",compressed_ids(range(1,35001)))

    def test_validity_requires_finite_outputs_and_nonnegative_vector(self):
        row=dict(status="ok",eigenvector_nonnegative="TRUE",dominant_mean_ploidy="2",dominant_growth_rate="-0.1")
        self.assertFalse(invalid(row))
        self.assertTrue(invalid(dict(row,eigenvector_nonnegative="FALSE")))
        self.assertTrue(invalid(dict(row,dominant_mean_ploidy="NA")))

    def test_recovery_receipt_rejects_modified_source_proof_or_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/"recovery").mkdir()
            nh.atomic_json(root/"metadata.json",dict(n_samples=1,oxygen_pct=[0],samples_sha256="samples"))
            output=root/"outputs.tsv.gz"; output.write_bytes(b"original output")
            source=root/"recovery/source.json"; source.write_text('{"original":"preserved"}')
            proof=dict(status="passed",outputs_sha256=efast.sha256(output),code_sha256={},
                artifact_sha256={"source.json":efast.sha256(source)})
            nh.atomic_json(root/"recovery/numerical_recovery.json",proof)
            receipt=dict(metadata_sha256=efast.sha256(root/"metadata.json"),samples_sha256="samples",
                outputs_sha256=efast.sha256(output),evaluator_sha256=efast.sha256(nh.HERE/"evaluate_fixed_o2.R"),
                n_rows=1,recovery_manifest_sha256=efast.sha256(root/"recovery/numerical_recovery.json"))
            nh.atomic_json(root/"outputs.tsv.gz.receipt.json",receipt)
            self.assertTrue(nh.verify_completion(root))
            source.write_text('{"original":"changed"}')
            with self.assertRaisesRegex(ValueError,"source artifact changed"):
                nh.verify_completion(root)
            source.write_text('{"original":"preserved"}')
            (root/"recovery/numerical_recovery.json").write_text('{}')
            with self.assertRaisesRegex(ValueError,"proof changed"):
                nh.verify_completion(root)
            nh.atomic_json(root/"recovery/numerical_recovery.json",proof)
            output.write_bytes(b"changed output")
            with self.assertRaisesRegex(ValueError,"receipt mismatch"):
                nh.verify_completion(root)


if __name__=="__main__": unittest.main()
