#!/usr/bin/env python3
"""Fault-injection and allocation-accounting tests; no models or real scheduler."""
import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("audit", Path(__file__).with_name("audit_culture_pilot_artifacts.py"))
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class ArtifactAudit(unittest.TestCase):
    def fixture(self, root):
        audit.write(root / "state.json", dict(status="pilot_fitting_finished", tasks={"fit_001": {"status": "accepted"}}))
        folder = root / "fits/fit_001"
        folder.mkdir(parents=True)
        (folder / "accepted_fit.rds").write_bytes(b"known immutable payload")
        audit.write(folder / "artifact_validation_v05.json", dict(pass_audit=True, **{"pass": True},
            target_constraint_tolerance=0, objective=3.))
        audit.write(folder / "validated_result.json", dict(accepted=True, launch_completed=True, objective=3.,
            artifact_validation_file="artifact_validation_v05.json",
            artifact_validation_sha256=audit.sha(folder / "artifact_validation_v05.json"),
            files_sha256={"accepted_fit.rds": audit.sha(folder / "accepted_fit.rds")}))
        return folder

    def test_valid_terminal_payload(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); self.fixture(root)
            result = audit.verify(root, True)
            self.assertEqual(result["payload_files"], 1)
            self.assertEqual(result["state_sha256"], audit.sha(root / "manifests/accepted_artifact_state_snapshot.json"))

    def test_modified_fit_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); folder = self.fixture(root)
            (folder / "accepted_fit.rds").write_bytes(b"modified payload")
            with self.assertRaisesRegex(RuntimeError, "payload changed"):
                audit.verify(root, True)

    def test_active_fit_prevents_final_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); self.fixture(root)
            audit.write(root / "state.json", dict(status="profile_running", tasks={"fit_001": {"status": "submitted"}}))
            with self.assertRaisesRegex(RuntimeError, "terminal study"):
                audit.verify(root, True)

    def test_changed_guard_is_rejected_even_when_score_is_unchanged(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); folder = self.fixture(root)
            guard = audit.read(folder / "artifact_validation_v05.json")
            guard["numeric_params_exact"] = {"p_misseg": "0.123"}
            audit.write(folder / "artifact_validation_v05.json", guard)
            with self.assertRaisesRegex(RuntimeError, "guard changed"):
                audit.verify(root, True)

    def test_unfingerprinted_guard_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); folder = self.fixture(root)
            result = audit.read(folder / "validated_result.json")
            del result["artifact_validation_sha256"]
            audit.write(folder / "validated_result.json", result)
            with self.assertRaisesRegex(RuntimeError, "Missing full-precision guard fingerprint"):
                audit.verify(root, True)

    def test_historical_guard_and_full_precision_supplement_are_both_verified(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); folder = self.fixture(root)
            audit.write(folder / "artifact_validation.json", {"pass": True, "objective": 3.})
            result = audit.read(folder / "validated_result.json")
            del result["artifact_validation_file"]
            result["artifact_validation_sha256"] = audit.sha(folder / "artifact_validation.json")
            audit.write(folder / "validated_result.json", result)
            audit.write(folder / "validation_supplement_v05.json", dict(
                original_validated_result_sha256=audit.sha(folder / "validated_result.json"),
                artifact_validation_sha256=audit.sha(folder / "artifact_validation_v05.json")))
            self.assertEqual(audit.verify(root, True)["accepted_fits"], 1)
            audit.write(folder / "artifact_validation.json", {"pass": True, "objective": 3., "changed": True})
            with self.assertRaisesRegex(RuntimeError, "guard changed"):
                audit.verify(root, True)

    def test_changed_validation_supplement_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); folder = self.fixture(root)
            audit.write(folder / "validation_supplement_v05.json", dict(original_validated_result_sha256="stale",
                artifact_validation_sha256=audit.sha(folder / "artifact_validation_v05.json")))
            with self.assertRaisesRegex(RuntimeError, "supplement changed"):
                audit.verify(root)

    def test_costs_count_only_study_parent_allocations(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            audit.write(root / "tasks/fit_001/submission.json", {"job_id": "101"})
            audit.write(root / "controllers/v06/fit_001/submission.json", {"job_id": "102"})
            text = ("JobID|JobName|State|ExitCode|ElapsedRaw|AllocCPUS|MaxRSS|TotalCPU|Start|End\n"
                    "101|fit|COMPLETED|0:0|3600|8||01:00:00|a|b\n"
                    "101.batch|batch|COMPLETED|0:0|3600|8|1G|01:00:00|a|b\n"
                    "102|controller|COMPLETED|0:0|60|1||00:01:00|a|b\n"
                    "999|unrelated|COMPLETED|0:0|7200|100||01:00:00|a|b\n")
            with patch.object(audit.subprocess, "run", return_value=SimpleNamespace(returncode=0, stdout=text, stderr="")):
                result = audit.collect_accounting(root, Path("/unused"))
            self.assertEqual(result["accounted_parent_jobs"], 2)
            self.assertAlmostEqual(result["total_allocated_core_hours"], 8 + 1/60)
            self.assertFalse(result["missing_job_ids"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
