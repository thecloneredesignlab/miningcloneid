#!/usr/bin/env python3
"""Check observation-scope labeling, denominator units and duplicate protection."""
import importlib.util
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location("findings", Path(__file__).with_name("summarize_culture_profile_findings.py"))
findings = importlib.util.module_from_spec(spec)
spec.loader.exec_module(findings)


class FluxScope(unittest.TestCase):
    def observations(self):
        base = dict(cohort="2N", segment_id="recorded", observation_id="a", recorded_duration_days="3",
            measured_initial_cells="100", measured_final_cells="400", growth_scored="TRUE",
            karyotype_scored="FALSE", flow_scored="FALSE")
        return [base, dict(base, segment_id="prehistory", observation_id="b", recorded_duration_days="NA",
            measured_initial_cells="NA", measured_final_cells="NA", growth_scored="FALSE", karyotype_scored="TRUE")]

    def flux(self):
        row = dict(task_id="fit", cohort="2N", segment_id="recorded", division="10")
        for field in ["ordinary_generated_loss", "ordinary_retained_loss", "ordinary_retained_gain", "ordinary_unchanged",
                      "ordinary_biological_nonviability", "numerical_tail_pruning", "ordinary_boundary_loss",
                      "ordinary_boundary_gain", "wgd_retained", "wgd_boundary", "direct_death"]:
            row[field] = "0"
        row["ordinary_retained_loss"] = "4.5"
        return row

    def test_distribution_scored_prehistory_remains_unrecorded_for_flux(self):
        scope = findings.recording_scope(self.observations())
        self.assertEqual(scope[("2N", "recorded")]["recording_scope"], "counts_and_duration_recorded")
        self.assertEqual(scope[("2N", "prehistory")]["recording_scope"], "counts_or_duration_unrecorded")
        self.assertEqual(scope[("2N", "prehistory")]["chromosome_observations"], 1)

    def test_ordinary_division_denominator_excludes_wgd(self):
        _, totals = findings.flux_summaries([self.flux()], findings.recording_scope(self.observations()), {"fit": {"p_wgd": .1}})
        self.assertEqual(totals[0]["ordinary_divisions"], 9.)
        self.assertEqual(totals[0]["ordinary_retained_loss_per_ordinary_division"], .5)

    def test_duplicate_interval_cannot_be_counted_twice(self):
        with self.assertRaisesRegex(ValueError, "double count"):
            findings.flux_summaries([self.flux(), self.flux()], findings.recording_scope(self.observations()), {"fit": {"p_wgd": .1}})


if __name__ == "__main__":
    unittest.main(verbosity=2)
