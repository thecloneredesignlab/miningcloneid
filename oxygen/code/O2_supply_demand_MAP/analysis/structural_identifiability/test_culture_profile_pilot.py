#!/usr/bin/env python3
"""Analytic orchestration oracles; no scheduler calls or biological fitting."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location("pilot", Path(__file__).with_name("culture_profile_pilot.py"))
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)


class PilotPolicies(unittest.TestCase):
    def test_grid_keeps_exact_null_optimum_and_upper_without_duplicates(self):
        for target, maximum in [("p_misseg", .4), ("buffer_beta", 10.)]:
            for optimum in [0., maximum, maximum / 3]:
                grid = pilot.profile_grid(target, optimum)
                self.assertEqual(len(grid), 10)
                self.assertEqual(len(set(grid)), 10)
                self.assertEqual(grid[0], 0.)
                self.assertEqual(grid[-1], maximum)
                self.assertIn(optimum, grid)
                self.assertEqual(grid, sorted(grid))

    def test_lower_envelope_is_minimum_and_preserves_negative_delta(self):
        result = pilot.envelope([{"task_id": "a", "objective": 7}, {"task_id": "b", "objective": 2}], 3)
        self.assertEqual(result["winner"], "b")
        self.assertEqual(result["objective"], 2)
        self.assertEqual(result["delta"], -1)
        self.assertFalse(result["agreement"])
        self.assertIsNone(pilot.envelope([], 0)["objective"])

    def test_baseline_needs_different_starts_and_stronger_search(self):
        a = dict(objective=3., primary_family="seed144", role="baseline_initial", ordinal=1)
        b = dict(objective=3.0005, primary_family="independent_1", role="baseline_initial", ordinal=2)
        self.assertFalse(pilot.baseline_gate([a, b], 3.1)[0])
        strong = dict(objective=2.9999, primary_family="independent_1", role="baseline_strengthen", ordinal=3, origin_best_q=3.)
        self.assertTrue(pilot.baseline_gate([a, b, strong], 3.1)[0])
        same = dict(a, role="baseline_strengthen", ordinal=4, origin_best_q=3.)
        self.assertFalse(pilot.baseline_gate([a, same], 3.1)[0])
        warm_copy = dict(strong, objective=3., agreement_eligible=False)
        far = dict(b, objective=3.5)
        self.assertFalse(pilot.baseline_gate([a, far, warm_copy], 3.6)[0])
        improving = dict(strong, objective=2.9)
        self.assertFalse(pilot.baseline_gate([a, b, improving], 3.1)[0])
        self.assertFalse(pilot.baseline_gate([a, b, strong], 2.9)[0])

    def test_xxlarge_has_per_job_resources_and_no_array_throttle(self):
        for reserve in [False, True]:
            args = pilot.resource_args({"reserve": reserve})
            self.assertIn("--qos=xxlarge", args)
            self.assertIn("--cpus-per-task=8", args)
            self.assertIn("--time=12:00:00" if reserve else "--time=08:00:00", args)
            self.assertFalse(any("%" in arg for arg in args))

    def test_rejected_or_failed_result_cannot_win(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            pilot.write(root / "study.json", {"master_seed": 1})
            p = pilot.Pilot(root, root)
            p.state["tasks"] = {
                "a": dict(status="accepted", stage="baseline", primary_family="a"),
                "b": dict(status="rejected", stage="baseline", primary_family="b"),
                "c": dict(status="failed", stage="baseline", primary_family="c")}
            for key, q in [("a", 3.), ("b", -100.), ("c", -200.)]:
                pilot.write(root / "fits" / key / "validated_result.json", {"task_id": key, "objective": q})
            self.assertEqual(p.baseline()["task_id"], "a")

    def test_budget_is_total_work_not_a_concurrency_gate(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            pilot.write(root / "study.json", {"master_seed": 1})
            p = pilot.Pilot(root, root)
            p.state["tasks"] = {str(i): {"reserve": True, "status": "submitted"} for i in range(20)}
            self.assertEqual(p.remaining_reserve(), 0)
            self.assertIsNone(p.create_task("profile", "repair", {}, "random", reserve=True))
            p.state["tasks"] = {str(i): {"reserve": False, "status": "submitted"} for i in range(125)}
            self.assertIsNone(p.create_task("baseline", "extra", {}, "random"))

    def test_completed_tasks_do_not_spawn_replacement_controllers(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            pilot.write(root / "study.json", {"master_seed": 1})
            p = pilot.Pilot(root, root)
            for status in ("accepted", "failed", "rejected"):
                p.ensure_controller({"task_id": "closed", "status": status})
            self.assertFalse((root / "controllers").exists())

    def test_atomic_state_replacement_and_full_precision(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "state.json"
            pilot.write(path, {"zero": 0., "target": .0102907636326605})
            self.assertEqual(pilot.read(path)["target"], .0102907636326605)
            pilot.write(path, {"zero": 0.})
            self.assertEqual(pilot.read(path), {"zero": 0.})
            self.assertEqual(len(list(Path(folder).glob(".atomic_*"))), 0)

    def test_json_transport_does_not_relax_the_exact_zero_null(self):
        self.assertTrue(pilot.target_transport_agrees(.00349735724318028, .003497357243180285))
        self.assertFalse(pilot.target_transport_agrees(.0034, .003497357243180285))
        self.assertTrue(pilot.target_transport_agrees(0., 0.))
        self.assertFalse(pilot.target_transport_agrees(1e-15, 0.))


if __name__ == "__main__":
    unittest.main(verbosity=2)
