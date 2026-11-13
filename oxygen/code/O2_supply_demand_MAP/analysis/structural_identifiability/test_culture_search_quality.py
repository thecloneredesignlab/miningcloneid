#!/usr/bin/env python3
"""Check that optimization support respects actual anchor provenance."""
import importlib.util
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location("findings", Path(__file__).with_name("summarize_culture_profile_findings.py"))
findings = importlib.util.module_from_spec(spec)
spec.loader.exec_module(findings)


class SearchQuality(unittest.TestCase):
    def fixture(self, definitions):
        state = {"points": {"pm_00": {"target": "p_misseg", "value": 0., "task_ids": []}}, "tasks": {}}
        results = {}
        for tid, labels, role, score, status in definitions:
            state["points"]["pm_00"]["task_ids"].append(tid)
            state["tasks"][tid] = {"starts_labels": labels, "role": role, "status": status}
            if status == "accepted":
                results[tid] = {"objective": score, "optimizer_trace": {"local_convergence": 0, "iter_completed": 26},
                                "budget_limited": False, "winner_phase": "native"}
        return state, results

    def test_warm_repetition_does_not_count_as_independent_recovery(self):
        state, results = self.fixture([
            ("a", ["unrestricted_best"], "point_start_1", 3., "accepted"),
            ("b", ["neighbor_pm_01"], "point_start_5", 3.0001, "accepted"),
            ("c", ["independent_1"], "point_start_3", 4., "accepted")])
        _, summary = findings.search_quality(state, results, .001)
        self.assertTrue(summary[0]["repeated_search_agreement"])
        self.assertEqual(summary[0]["independent_anchor_families_near_minimum"], 0)
        self.assertEqual(summary[0]["independent_best_gap"], 1.)

    def test_two_distinct_independent_anchors_are_reported(self):
        state, results = self.fixture([
            ("a", ["independent_1"], "point_start_3", 3., "accepted"),
            ("b", ["independent_2"], "point_start_4", 3.0005, "accepted"),
            ("c", ["null_independent"], "null_independent", 4., "accepted")])
        endpoints, summary = findings.search_quality(state, results, .001)
        self.assertEqual(len(endpoints), 3)
        self.assertTrue(summary[0]["two_independent_anchors_near_minimum"])
        self.assertEqual(summary[0]["independent_anchor_families_near_minimum"], 2)

    def test_return_from_null_is_continuation_not_independent_recovery(self):
        state, results = self.fixture([
            ("a", ["unrestricted_best"], "point_start_1", 3., "accepted"),
            ("b", ["reverse_from_null_pm_00"], "point_reverse", 3.0001, "accepted")])
        endpoints, summary = findings.search_quality(state, results, .001)
        self.assertEqual(endpoints[1]["start_class"], "neighbor_or_reverse")
        self.assertEqual(summary[0]["independent_anchor_families_near_minimum"], 0)

    def test_mixed_strengthening_is_not_attributed_to_its_independent_member(self):
        state, results = self.fixture([
            ("a", ["unrestricted_best"], "point_start_1", 3., "accepted"),
            ("b", ["current_best", "independent_5"], "point_strengthen", 2.9, "accepted"),
            ("c", ["independent_2"], "point_start_4", None, "failed"),
            ("d", ["independent_3"], "point_start_5", None, "submitted")])
        _, summary = findings.search_quality(state, results, .001)
        row = summary[0]
        self.assertEqual(row["winner_start_class"], "mixed_warm_and_independent")
        self.assertEqual(row["independent_anchor_searches"], 0)
        self.assertIsNone(row["independent_best_gap"])
        self.assertEqual(row["active"], 1)
        self.assertEqual(row["failed_or_rejected"], 1)
        self.assertAlmostEqual(row["endpoint_score_spread"], .1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
