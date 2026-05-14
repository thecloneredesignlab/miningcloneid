import importlib.util
import io
import math
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares


MODULE_PATH = Path(__file__).resolve().parents[2] / "code" / "invitro_fitting.py"
SPEC = importlib.util.spec_from_file_location("invitro_fitting", MODULE_PATH)
invitro_fitting = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(invitro_fitting)

LEGACY_MODULE_PATH = Path(__file__).resolve().parents[2] / "code" / "legacy_gemcitabine_pk.py"
LEGACY_SPEC = importlib.util.spec_from_file_location("legacy_gemcitabine_pk", LEGACY_MODULE_PATH)
legacy_gemcitabine_pk = importlib.util.module_from_spec(LEGACY_SPEC)
assert LEGACY_SPEC.loader is not None
LEGACY_SPEC.loader.exec_module(legacy_gemcitabine_pk)


class ExtracellularDecayTests(unittest.TestCase):
    def make_constant_surface(self, value_uM: float = 2.0) -> invitro_fitting.DfdctpSignalSurface:
        profile = invitro_fitting.DfdctpDoseProfile(
            sheet_name="synthetic",
            dose_uM=1.0,
            analyte_column="dFdCTP (ng/mL)",
            time_days=np.array([0.0, 1.0]),
            raw_signal_uM_values=np.array([value_uM, value_uM]),
            induced_signal_uM_values=np.array([value_uM, value_uM]),
            peak_signal_uM=value_uM,
            peak_time_days=0.0,
            tail_half_life_days=1.0,
            tail_decay_per_day=math.log(2.0),
        )
        return invitro_fitting.DfdctpSignalSurface(
            source_ploidy="synthetic",
            analyte_column="dFdCTP (ng/mL)",
            calibration_profiles_by_dose={1.0: profile},
            calibration_sheet_names=("synthetic",),
        )

    def make_paths_with_platemap(self, platemap_df: pd.DataFrame) -> invitro_fitting.ExperimentPaths:
        tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(tmpdir.cleanup)
        root = Path(tmpdir.name)
        exp_base = root / "data" / "InVitroData_Gemcitabine"
        exp_base.mkdir(parents=True, exist_ok=True)
        platemap_path = exp_base / "Gemcitabine_PlateMap_20240111.xlsx"
        platemap_df.to_excel(platemap_path, index=False)
        return invitro_fitting.ExperimentPaths(
            project_root=root,
            exp_base=exp_base,
            counts_raw=exp_base / "processed" / "counts_by_well_time.parquet",
            counts_agg=exp_base / "processed" / "counts_by_well_time_wellAggregated.parquet",
            platemap=platemap_path,
            pkpd_constants=exp_base / "drugKinetics" / "GemcitabineExposure_PKPD.xlsx",
            output_dir=root / "outputs",
        )

    def test_get_pk_reference_dose_uM_uses_documented_sheet_mapping(self):
        self.assertEqual(invitro_fitting.get_pk_reference_dose_uM("2N"), 1.0)
        self.assertEqual(invitro_fitting.get_pk_reference_dose_uM("4N"), 1.0)
        self.assertEqual(invitro_fitting.get_pk_reference_dose_uM("2N_lowInitialGemcitabine"), 0.1)
        self.assertEqual(invitro_fitting.get_pk_reference_dose_uM("4N_lowInitialGemcitabine"), 0.1)
        with self.assertRaises(KeyError):
            invitro_fitting.get_pk_reference_dose_uM("unknown")

    def test_validate_model_variant(self):
        for variant in [
            "delayed_death_only",
            "immediate_cytostasis_delayed_death",
            "delayed_cytostasis_delayed_death",
        ]:
            self.assertEqual(invitro_fitting.validate_model_variant(variant), variant)
        with self.assertRaises(ValueError):
            invitro_fitting.validate_model_variant("bad_variant")

    def test_cytostasis_multiplier(self):
        self.assertEqual(invitro_fitting.cytostasis_multiplier(0.0, 10.0), 1.0)
        mid = invitro_fitting.cytostasis_multiplier(1.0, 2.0)
        high = invitro_fitting.cytostasis_multiplier(2.0, 2.0)
        self.assertGreaterEqual(mid, 0.0)
        self.assertLessEqual(mid, 1.0)
        self.assertGreater(mid, high)
        with self.assertRaises(ValueError):
            invitro_fitting.cytostasis_multiplier(1.0, -1.0)

    def test_model2_reduces_growth_immediately(self):
        surface = self.make_constant_surface(0.1)
        t = np.linspace(0.0, 2.0, 5)
        legacy_cfg = invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False)
        model2_cfg = invitro_fitting.JointFitConfig(model_variant="immediate_cytostasis_delayed_death", use_hill_dose_gate=False, use_confluence_death=False)
        alive1, _ = invitro_fitting.simulate_joint_ext(
            t=t, params=[1.0, 0.0, 0.0, 1.0, 0.0], N0=100.0, D0=0.0, r=0.8, K=5000.0,
            dose_muM=0.1, dfdctp_signal_curve=surface, n_tr=2, model_variant="delayed_death_only",
            fit_config=legacy_cfg,
        )
        alive2, _ = invitro_fitting.simulate_joint_ext(
            t=t, params=[1.0, 0.0, 0.0, 100.0, 1.0, 0.0], N0=100.0, D0=0.0, r=0.8, K=5000.0,
            dose_muM=0.1, dfdctp_signal_curve=surface, n_tr=2, model_variant="immediate_cytostasis_delayed_death",
            fit_config=model2_cfg,
        )
        self.assertLess(alive2[-1], alive1[-1])

    def test_model3_has_delayed_cytostasis(self):
        surface = self.make_constant_surface(0.1)
        t = np.linspace(0.0, 2.0, 9)
        legacy_cfg = invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False)
        model2_cfg = invitro_fitting.JointFitConfig(model_variant="immediate_cytostasis_delayed_death", use_hill_dose_gate=False, use_confluence_death=False)
        model3_cfg = invitro_fitting.JointFitConfig(model_variant="delayed_cytostasis_delayed_death", use_hill_dose_gate=False, use_confluence_death=False)
        alive1, _ = invitro_fitting.simulate_joint_ext(
            t=t, params=[1.0, 0.0, 0.0, 1.0, 0.0], N0=100.0, D0=0.0, r=0.8, K=5000.0,
            dose_muM=0.1, dfdctp_signal_curve=surface, n_tr=2, model_variant="delayed_death_only",
            fit_config=legacy_cfg,
        )
        alive2, _ = invitro_fitting.simulate_joint_ext(
            t=t, params=[1.0, 0.0, 0.0, 100.0, 1.0, 0.0], N0=100.0, D0=0.0, r=0.8, K=5000.0,
            dose_muM=0.1, dfdctp_signal_curve=surface, n_tr=2, model_variant="immediate_cytostasis_delayed_death",
            fit_config=model2_cfg,
        )
        alive3, _ = invitro_fitting.simulate_joint_ext(
            t=t, params=[1.0, 0.0, 0.0, 100.0, 1.0, 0.0], N0=100.0, D0=0.0, r=0.8, K=5000.0,
            dose_muM=0.1, dfdctp_signal_curve=surface, n_tr=2, model_variant="delayed_cytostasis_delayed_death",
            fit_config=model3_cfg,
        )
        self.assertLess(abs(alive3[1] - alive1[1]), abs(alive2[1] - alive1[1]))
        self.assertLess(alive3[-1], alive1[-1])

    def test_model2_k_cyto_near_zero_matches_model1(self):
        surface = self.make_constant_surface(0.1)
        t = np.linspace(0.0, 2.0, 7)
        legacy_cfg = invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False)
        model2_cfg = invitro_fitting.JointFitConfig(model_variant="immediate_cytostasis_delayed_death", use_hill_dose_gate=False, use_confluence_death=False)
        alive1, dead1 = invitro_fitting.simulate_joint_ext(
            t=t, params=[1.0, 0.5, 0.2, 1.0, 0.0], N0=100.0, D0=0.0, r=0.8, K=5000.0,
            dose_muM=0.1, dfdctp_signal_curve=surface, n_tr=2, model_variant="delayed_death_only",
            fit_config=legacy_cfg,
        )
        alive2, dead2 = invitro_fitting.simulate_joint_ext(
            t=t, params=[1.0, 0.5, 0.2, 1e-8, 1.0, 0.0], N0=100.0, D0=0.0, r=0.8, K=5000.0,
            dose_muM=0.1, dfdctp_signal_curve=surface, n_tr=2, model_variant="immediate_cytostasis_delayed_death",
            fit_config=model2_cfg,
        )
        self.assertTrue(np.allclose(alive1, alive2, rtol=1e-8, atol=1e-8))
        self.assertTrue(np.allclose(dead1, dead2, rtol=1e-8, atol=1e-8))

    def test_import_has_no_stdout_or_path_warnings(self):
        stdout = io.StringIO()
        fresh_spec = importlib.util.spec_from_file_location("invitro_fitting_fresh", MODULE_PATH)
        fresh_module = importlib.util.module_from_spec(fresh_spec)
        assert fresh_spec.loader is not None
        with redirect_stdout(stdout):
            fresh_spec.loader.exec_module(fresh_module)
        self.assertEqual(stdout.getvalue(), "")

    def test_load_platemap_wide_format(self):
        platemap = pd.DataFrame(
            {
                "Row": ["A", "B"],
                "1": ["2N 0 nM", "4N 12.5 nM"],
                "2": ["2N 6.25 nM", "4N 25 nM"],
            }
        )
        paths = self.make_paths_with_platemap(platemap)
        parsed = invitro_fitting.load_platemap(paths=paths)
        self.assertEqual(set(parsed.columns), {"plate_row", "plate_col", "ploidy", "gem", "condition_raw"})
        self.assertEqual(len(parsed), 4)
        self.assertIn("2N", parsed["ploidy"].tolist())
        self.assertIn("12.5 nM", parsed["gem"].tolist())

    def test_load_platemap_long_well_condition_format(self):
        platemap = pd.DataFrame(
            {
                "well": ["A01", "B2"],
                "condition": ["2N 0 nM", "4N 25 nM"],
            }
        )
        paths = self.make_paths_with_platemap(platemap)
        parsed = invitro_fitting.load_platemap(paths=paths)
        self.assertEqual(parsed.loc[0, "plate_row"], "A")
        self.assertEqual(int(parsed.loc[0, "plate_col"]), 1)
        self.assertEqual(parsed.loc[1, "ploidy"], "4N")
        self.assertEqual(parsed.loc[1, "gem"], "25 nM")

    def test_load_platemap_long_row_col_separate_fields_format(self):
        platemap = pd.DataFrame(
            {
                "row": ["A", "H"],
                "column": [1, 12],
                "ploidy": ["2N", "4N"],
                "gem": ["0 nM", "50 nM"],
            }
        )
        paths = self.make_paths_with_platemap(platemap)
        parsed = invitro_fitting.load_platemap(paths=paths)
        self.assertEqual(parsed.loc[0, "condition_raw"], "2N | 0 nM")
        self.assertEqual(parsed.loc[1, "plate_row"], "H")
        self.assertEqual(int(parsed.loc[1, "plate_col"]), 12)

    def test_load_platemap_unsupported_format_raises(self):
        platemap = pd.DataFrame({"foo": ["A01"], "bar": ["2N 0 nM"]})
        paths = self.make_paths_with_platemap(platemap)
        with self.assertRaises(ValueError):
            invitro_fitting.load_platemap(paths=paths)

    def test_pk_tail_diagnostics_records_tail_points(self):
        df = pd.DataFrame(
            {
                "Timepoint": [0.0, 24.0, 48.0, 72.0],
                "dFdCTP (ng/mL)": [10.0, 25.0, 20.0, 15.0],
            }
        )
        profile = invitro_fitting.build_dfdctp_profile_from_sheet(df, sheet_name="2N", reference_dose_uM=1.0)
        self.assertIsNotNone(profile)
        diagnostics = profile.tail_fit_diagnostics
        self.assertIsNotNone(diagnostics)
        self.assertGreaterEqual(diagnostics.tail_n_points, 2)
        self.assertEqual(len(diagnostics.tail_time_days), diagnostics.tail_n_points)
        self.assertEqual(len(diagnostics.tail_signal_uM), diagnostics.tail_n_points)

    def test_pk_tail_diagnostics_records_fallback(self):
        df = pd.DataFrame(
            {
                "Timepoint": [0.0],
                "dFdCTP (ng/mL)": [10.0],
            }
        )
        profile = invitro_fitting.build_dfdctp_profile_from_sheet(df, sheet_name="2N", reference_dose_uM=1.0)
        self.assertIsNotNone(profile)
        diagnostics = profile.tail_fit_diagnostics
        self.assertTrue(diagnostics.fallback_used)
        self.assertFalse(diagnostics.tail_fit_used)

    def test_parse_pk_concentration_value_handles_censored_and_numeric_inputs(self):
        value, was_censored, upper = invitro_fitting.parse_pk_concentration_value("BDL", censored_strategy="nan")
        self.assertTrue(np.isnan(value))
        self.assertTrue(was_censored)
        self.assertTrue(np.isnan(upper))

        value, was_censored, upper = invitro_fitting.parse_pk_concentration_value("BQL", censored_strategy="nan")
        self.assertTrue(np.isnan(value))
        self.assertTrue(was_censored)
        self.assertTrue(np.isnan(upper))

        value, was_censored, upper = invitro_fitting.parse_pk_concentration_value("<0.5", censored_strategy="half_lod")
        self.assertTrue(math.isclose(value, 0.25, rel_tol=1e-12, abs_tol=1e-12))
        self.assertTrue(was_censored)
        self.assertTrue(math.isclose(upper, 0.5, rel_tol=1e-12, abs_tol=1e-12))

        value, was_censored, upper = invitro_fitting.parse_pk_concentration_value("<0.5", censored_strategy="nan")
        self.assertTrue(np.isnan(value))
        self.assertTrue(was_censored)
        self.assertTrue(math.isclose(upper, 0.5, rel_tol=1e-12, abs_tol=1e-12))

        value, was_censored, upper = invitro_fitting.parse_pk_concentration_value("1.23", censored_strategy="nan")
        self.assertTrue(math.isclose(value, 1.23, rel_tol=1e-12, abs_tol=1e-12))
        self.assertFalse(was_censored)
        self.assertTrue(np.isnan(upper))

        value, was_censored, upper = invitro_fitting.parse_pk_concentration_value(np.nan, censored_strategy="nan")
        self.assertTrue(np.isnan(value))
        self.assertFalse(was_censored)
        self.assertTrue(np.isnan(upper))

    def test_parse_pk_concentration_value_zero_strategy_is_explicit(self):
        value, was_censored, upper = invitro_fitting.parse_pk_concentration_value("BDL", censored_strategy="zero")
        self.assertEqual(value, 0.0)
        self.assertTrue(was_censored)
        self.assertTrue(np.isnan(upper))

        value, was_censored, upper = invitro_fitting.parse_pk_concentration_value("<0.5", censored_strategy="zero")
        self.assertEqual(value, 0.0)
        self.assertTrue(was_censored)
        self.assertTrue(math.isclose(upper, 0.5, rel_tol=1e-12, abs_tol=1e-12))

    def test_converted_extracellular_gem_signal_matches_legacy_expression(self):
        dose_muM = 12.5
        self.assertTrue(
            math.isclose(
                legacy_gemcitabine_pk.converted_extracellular_gem_signal(dose_muM),
                0.00971 * dose_muM,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )

    def test_gemcitabine_ng_per_ml_to_uM_uses_physical_conversion(self):
        concentration_ng_per_ml = 26.32
        self.assertTrue(
            math.isclose(
                legacy_gemcitabine_pk.gemcitabine_ng_per_ml_to_uM(concentration_ng_per_ml),
                concentration_ng_per_ml / legacy_gemcitabine_pk.GEMCITABINE_MOLECULAR_WEIGHT_NG_PER_NMOL,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )

    def test_dfdctp_ng_per_ml_to_uM_uses_molecular_weight_conversion(self):
        concentration_ng_per_ml = 100.62
        self.assertTrue(
            math.isclose(
                invitro_fitting.dfdctp_ng_per_ml_to_uM(concentration_ng_per_ml),
                concentration_ng_per_ml / 503.1,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )

    def test_dfdctp_signal_surface_preserves_calibrated_profiles(self):
        high = invitro_fitting.DfdctpDoseProfile(
            sheet_name="2N",
            dose_uM=1.0,
            analyte_column="dFdCTP (ng/mL)",
            time_days=np.array([0.0, 1.0, 2.0]),
            raw_signal_uM_values=np.array([1.0, 5.0, 3.0]),
            induced_signal_uM_values=np.array([0.0, 4.0, 2.0]),
            peak_signal_uM=4.0,
            peak_time_days=1.0,
            tail_half_life_days=1.0,
            tail_decay_per_day=math.log(2.0),
        )
        low = invitro_fitting.DfdctpDoseProfile(
            sheet_name="2N_lowInitialGemcitabine",
            dose_uM=0.1,
            analyte_column="dFdCTP (ng/mL)",
            time_days=np.array([0.0, 2.0, 3.0]),
            raw_signal_uM_values=np.array([1.0, 2.5, 2.0]),
            induced_signal_uM_values=np.array([0.0, 1.5, 1.0]),
            peak_signal_uM=1.5,
            peak_time_days=2.0,
            tail_half_life_days=1.0,
            tail_decay_per_day=math.log(2.0),
        )
        surface = invitro_fitting.DfdctpSignalSurface(
            source_ploidy="2N",
            analyte_column="dFdCTP (ng/mL)",
            calibration_profiles_by_dose={0.1: low, 1.0: high},
            calibration_sheet_names=("2N", "2N_lowInitialGemcitabine"),
        )
        self.assertTrue(math.isclose(surface(1.0, 1.0), 4.0, rel_tol=1e-12, abs_tol=1e-12))
        self.assertTrue(math.isclose(surface(2.0, 0.1), 1.5, rel_tol=1e-12, abs_tol=1e-12))
        self.assertGreaterEqual(surface(2.5, 0.5), 0.0)
        self.assertEqual(surface.get_matching_calibrated_dose(0.1000000000001), 0.1)
        self.assertEqual(surface.get_matching_calibrated_dose(0.9999999999999), 1.0)

    def test_baseline_subtract_treatment_induced_signal(self):
        raw = np.array([10.0, 25.0, 20.0])
        expected = np.array([0.0, 15.0, 10.0])
        result = invitro_fitting.baseline_subtract_treatment_induced_signal(raw)
        self.assertTrue(np.array_equal(result, expected))

    def test_build_dfdctp_profile_from_sheet_converts_and_normalizes(self):
        df = pd.DataFrame(
            {
                "Timepoint": [0.0, 0.0, 24.0, 24.0, 48.0, 48.0],
                "dFdCTP (ng/mL)": [10.0, 10.0, 25.0, 25.0, 20.0, 20.0],
            }
        )
        profile = invitro_fitting.build_dfdctp_profile_from_sheet(
            df,
            sheet_name="2N",
            reference_dose_uM=1.0,
        )
        self.assertIsNotNone(profile)
        assert profile is not None
        expected_signal = np.array([0.0, 15.0, 10.0]) / invitro_fitting.DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL
        self.assertTrue(np.allclose(profile.induced_signal_uM_values, expected_signal))
        self.assertTrue(math.isclose(profile.peak_signal_uM, expected_signal[1], rel_tol=1e-12, abs_tol=1e-12))
        self.assertGreaterEqual(np.min(profile.evaluate(np.array([0.0, 1.0, 2.0]))), 0.0)

    def test_dfdctp_signal_surface_is_nonnegative_and_monotone_in_dose(self):
        pk_sheets = {
            "2N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0], "dFdCTP (ng/mL)": [10.0, 25.0, 20.0]}),
            "2N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0], "dFdCTP (ng/mL)": [10.0, 17.5, 15.0]}),
        }
        curve = invitro_fitting.build_dfdctp_signal_curve_for_ploidy(pk_sheets, "2N")
        t_grid = np.linspace(0.0, 5.0, 11)
        low = curve(t_grid, 0.1)
        high = curve(t_grid, 1.0)
        self.assertTrue(np.all(low >= 0.0))
        self.assertTrue(np.all(high >= 0.0))
        self.assertTrue(np.any(high > low))

    def test_build_dfdctp_signal_curve_for_ploidy_does_not_need_gemcitabine_column(self):
        pk_sheets = {
            "2N": pd.DataFrame(
                {
                    "Timepoint": [0.0, 24.0, 48.0],
                    "dFdCTP (ng/mL)": [10.0, 25.0, 20.0],
                }
            ),
            "2N_lowInitialGemcitabine": pd.DataFrame(
                {
                    "Timepoint": [0.0, 24.0, 48.0],
                    "dFdCTP (ng/mL)": [10.0, 17.5, 15.0],
                }
            ),
        }
        curve = invitro_fitting.build_dfdctp_signal_curve_for_ploidy(pk_sheets, "2N")
        self.assertEqual(curve.analyte_column, "dFdCTP (ng/mL)")
        self.assertEqual(curve.calibration_sheet_names[0], "2N")
        self.assertGreater(curve.calibration_profiles_by_dose[1.0].peak_signal_uM, 0.0)

    def test_low_dose_peak_timing_is_not_forced_to_high_dose_peak_timing(self):
        pk_sheets = {
            "2N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0], "dFdCTP (ng/mL)": [10.0, 25.0, 20.0]}),
            "2N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0], "dFdCTP (ng/mL)": [10.0, 17.5, 15.0]}),
        }
        surface = invitro_fitting.build_dfdctp_signal_curve_for_ploidy(pk_sheets, "2N")
        low_profile = surface.calibration_profiles_by_dose[0.1]
        high_profile = surface.calibration_profiles_by_dose[1.0]
        self.assertNotEqual(low_profile.peak_time_days, high_profile.peak_time_days)
        self.assertTrue(math.isclose(surface(low_profile.peak_time_days, 0.1), low_profile.peak_signal_uM, rel_tol=1e-12, abs_tol=1e-12))
        self.assertTrue(math.isclose(surface(high_profile.peak_time_days, 1.0), high_profile.peak_signal_uM, rel_tol=1e-12, abs_tol=1e-12))

    def test_intermediate_dose_interpolates_between_dose_specific_profiles(self):
        pk_sheets = {
            "2N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0], "dFdCTP (ng/mL)": [10.0, 25.0, 20.0]}),
            "2N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0], "dFdCTP (ng/mL)": [10.0, 17.5, 15.0]}),
        }
        surface = invitro_fitting.build_dfdctp_signal_curve_for_ploidy(pk_sheets, "2N")
        t_eval = 2.0
        low_value = surface(t_eval, 0.1)
        mid_value = surface(t_eval, 0.316227766)
        high_value = surface(t_eval, 1.0)
        self.assertGreaterEqual(mid_value, min(low_value, high_value))
        self.assertLessEqual(mid_value, max(low_value, high_value))

    def test_zero_dose_returns_zero_signal_everywhere(self):
        pk_sheets = {
            "2N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0], "dFdCTP (ng/mL)": [10.0, 25.0, 20.0]}),
            "2N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0], "dFdCTP (ng/mL)": [10.0, 17.5, 15.0]}),
        }
        surface = invitro_fitting.build_dfdctp_signal_curve_for_ploidy(pk_sheets, "2N")
        values = surface(np.linspace(0.0, 5.0, 8), 0.0)
        self.assertTrue(np.array_equal(values, np.zeros_like(values)))

    def test_below_range_dose_uses_proportional_scaling_of_minimum_profile(self):
        pk_sheets = {
            "2N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0], "dFdCTP (ng/mL)": [10.0, 25.0, 20.0]}),
            "2N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0], "dFdCTP (ng/mL)": [10.0, 17.5, 15.0]}),
        }
        surface = invitro_fitting.build_dfdctp_signal_curve_for_ploidy(pk_sheets, "2N")
        t_eval = np.linspace(0.0, 5.0, 11)
        min_values = surface(t_eval, 0.1)
        below_values = surface(t_eval, 0.01)
        self.assertTrue(np.allclose(below_values, min_values * 0.1, rtol=1e-12, atol=1e-12))
        self.assertTrue(np.all((below_values > 0.0) == (min_values > 0.0)))

    def test_below_range_dose_preserves_minimum_profile_peak_timing(self):
        pk_sheets = {
            "2N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0], "dFdCTP (ng/mL)": [10.0, 25.0, 20.0]}),
            "2N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0], "dFdCTP (ng/mL)": [10.0, 17.5, 15.0]}),
        }
        surface = invitro_fitting.build_dfdctp_signal_curve_for_ploidy(pk_sheets, "2N")
        min_profile = surface.calibration_profiles_by_dose[0.1]
        peak_time = min_profile.peak_time_days
        peak_value_below = surface(peak_time, 0.01)
        nearby_values = surface(np.array([peak_time - 0.5, peak_time, peak_time + 0.5]), 0.01)
        self.assertGreaterEqual(peak_value_below, nearby_values[0])
        self.assertGreaterEqual(peak_value_below, nearby_values[-1])

    def test_doses_below_calibrated_minimum_trigger_extrapolation_status(self):
        pk_sheets = {
            "2N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0], "dFdCTP (ng/mL)": [10.0, 25.0, 20.0]}),
            "2N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0], "dFdCTP (ng/mL)": [10.0, 17.5, 15.0]}),
        }
        surface = invitro_fitting.build_dfdctp_signal_curve_for_ploidy(pk_sheets, "2N")
        status = surface.calibration_status_for_doses([0.003125, 0.1, 1.0])
        self.assertIn("below", status)
        self.assertIn("proportional scaling", status)

    def test_live_dead_driver_runs_without_gemcitabine_column(self):
        pk_sheets = {
            "2N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0], "dFdCTP (ng/mL)": [10.0, 25.0, 20.0]}),
            "2N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0], "dFdCTP (ng/mL)": [10.0, 17.5, 15.0]}),
        }
        surface = invitro_fitting.build_dfdctp_signal_curve_for_ploidy(pk_sheets, "2N")
        alive, dead = invitro_fitting.simulate_joint_ext(
            t=np.linspace(0.0, 5.0, 6),
            params=[0.5, 1.0, 0.2, 1.0, 0.0],
            N0=1000.0,
            D0=0.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.05,
            dfdctp_signal_curve=surface,
            n_tr=2,
            model_variant="delayed_death_only",
            fit_config=invitro_fitting.JointFitConfig(
                model_variant="delayed_death_only",
                use_hill_dose_gate=False,
                use_confluence_death=False,
            ),
        )
        self.assertTrue(np.all(np.isfinite(alive)))
        self.assertTrue(np.all(np.isfinite(dead)))

    def test_preferred_surface_builder_runs_without_gemcitabine_column(self):
        pk_sheets = {
            "2N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0], "dFdCTP (ng/mL)": [10.0, 25.0, 20.0]}),
            "2N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0], "dFdCTP (ng/mL)": [10.0, 17.5, 15.0]}),
            "4N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0], "dFdCTP (ng/mL)": [11.0, 31.0, 21.0]}),
            "4N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0], "dFdCTP (ng/mL)": [11.0, 18.0, 16.0]}),
        }
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            surfaces = invitro_fitting.build_preferred_dfdctp_signal_surfaces(pk_sheets)
        self.assertEqual(set(surfaces.keys()), {"2N", "4N"})
        self.assertIsInstance(surfaces["2N"], invitro_fitting.DfdctpSignalSurface)
        self.assertIn("dFdCTP", surfaces["2N"].analyte_column)

    def test_preferred_surface_builder_requires_dfdctp_column(self):
        pk_sheets = {
            "2N": pd.DataFrame({"Timepoint": [0.0, 24.0, 48.0]}),
            "2N_lowInitialGemcitabine": pd.DataFrame({"Timepoint": [0.0, 48.0, 72.0]}),
        }
        with self.assertRaises(ValueError):
            invitro_fitting.build_dfdctp_signal_curve_for_ploidy(pk_sheets, "2N")

    def test_assert_preferred_live_dead_driver_rejects_legacy_or_generic_callables(self):
        with self.assertRaises(TypeError):
            invitro_fitting.assert_preferred_live_dead_driver(lambda t, dose: 0.0)
        with self.assertRaises(TypeError):
            invitro_fitting.assert_preferred_live_dead_driver(object())

    def test_main_module_does_not_expose_legacy_parent_gemcitabine_symbols(self):
        legacy_only_symbols = [
            "MODEL_GEM_EXPOSURE_SIGNAL_PER_DOSE_UM",
            "GEMCITABINE_MOLECULAR_WEIGHT_NG_PER_NMOL",
            "ExtracellularExposureCurve",
            "gemcitabine_ng_per_ml_to_uM",
            "build_extracellular_exposure_curve_from_pk_sheet",
            "fit_constrained_extracellular_gem_decay",
            "simulate_intracellular_dfdctp_signal",
            "fit_intracellular_pk_parameters",
            "get_r_eta_parameters",
            "extract_eta",
            "extract_half_life",
            "estimate_eta_per_day",
        ]
        for symbol_name in legacy_only_symbols:
            self.assertFalse(hasattr(invitro_fitting, symbol_name), symbol_name)
            self.assertTrue(hasattr(legacy_gemcitabine_pk, symbol_name), symbol_name)

    def test_estimate_eta_per_day_uses_physical_uM_formula(self):
        delta_dfdctp_uM_per_hour = 0.25
        dose_muM = 1.0
        expected_eta = (delta_dfdctp_uM_per_hour / dose_muM) * 24.0
        refactored_eta = legacy_gemcitabine_pk.estimate_eta_per_day(
            delta_dfdctp_uM_per_hour=delta_dfdctp_uM_per_hour,
            dose_muM=dose_muM,
        )
        self.assertTrue(
            math.isclose(
                refactored_eta,
                expected_eta,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )

    def test_decay_rate_from_half_life_days_matches_ln2(self):
        self.assertTrue(
            math.isclose(
                invitro_fitting.decay_rate_from_half_life_days(1.0),
                math.log(2.0),
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )

    def test_decay_rate_from_half_life_days_validates_positive_input(self):
        with self.assertRaises(ValueError):
            invitro_fitting.decay_rate_from_half_life_days(0.0)

    def test_constrained_bolus_curve_uses_time_zero_bolus_and_dose_scaling(self):
        curve = legacy_gemcitabine_pk.build_extracellular_exposure_curve_from_profile(
            time_days=np.array([0.0, 1.0, 2.0]),
            gem_concentration_ng_per_ml=np.array([100.0, 50.0, 25.0]),
            gem_was_censored=np.array([False, False, False]),
            gem_censor_upper_bound_ng_per_ml=np.array([np.nan, np.nan, np.nan]),
            reference_dose_muM=1.0,
            fallback_half_life_days=1.0,
            analyte_column="Gemcitabine (ng/mL)",
            source_ploidy="synthetic",
            preferred_mode="constrained_bolus_exponential",
        )
        self.assertEqual(curve.mode, "constrained_bolus_exponential")
        self.assertTrue(math.isclose(curve(0.0, 1.0), 1.0, rel_tol=1e-12, abs_tol=1e-12))
        self.assertLess(curve(0.5, 1.0), curve(0.0, 1.0))
        self.assertTrue(math.isclose(curve(0.5, 0.5), 0.5 * curve(0.5, 1.0), rel_tol=1e-12, abs_tol=1e-12))

    def test_fit_exponential_pk_decay_model_recovers_half_life(self):
        t_days = np.array([0.0, 1.0, 2.0, 3.0])
        conc = 100.0 * np.exp(-math.log(2.0) * t_days)
        fit = invitro_fitting.fit_exponential_pk_decay_model(t_days, conc)
        self.assertIsNotNone(fit)
        assert fit is not None
        self.assertTrue(math.isclose(fit["k_ext_decay_per_day"], math.log(2.0), rel_tol=1e-6, abs_tol=1e-6))
        self.assertTrue(math.isclose(fit["half_life_days"], 1.0, rel_tol=1e-6, abs_tol=1e-6))

    def test_half_life_fallback_curve_halves_at_one_day(self):
        curve = legacy_gemcitabine_pk.build_extracellular_exposure_curve_from_profile(
            time_days=np.array([], dtype=float),
            gem_concentration_ng_per_ml=np.array([], dtype=float),
            gem_was_censored=np.array([], dtype=bool),
            gem_censor_upper_bound_ng_per_ml=np.array([], dtype=float),
            reference_dose_muM=1.0,
            fallback_half_life_days=1.0,
            analyte_column=None,
            source_ploidy="fallback",
            preferred_mode="half_life_fallback",
        )
        self.assertEqual(curve.mode, "legacy_half_life_fallback")
        self.assertTrue(math.isclose(curve(1.0, 1.0), 0.5 * curve(0.0, 1.0), rel_tol=1e-12, abs_tol=1e-12))

    def test_parse_pk_concentration_value_validates_strategy(self):
        with self.assertRaises(ValueError):
            invitro_fitting.parse_pk_concentration_value("BDL", censored_strategy="bad")

    def test_converted_extracellular_gem_signal_validates_inputs(self):
        with self.assertRaises(ValueError):
            legacy_gemcitabine_pk.converted_extracellular_gem_signal(0.0)
        with self.assertRaises(ValueError):
            legacy_gemcitabine_pk.converted_extracellular_gem_signal(1.0, conversion_factor=0.0)

    def test_estimate_eta_per_day_validates_inputs(self):
        with self.assertRaises(ValueError):
            legacy_gemcitabine_pk.estimate_eta_per_day(1.0, 0.0)

    def test_extracellular_concentration_starts_at_initial_value(self):
        dose_muM = 12.5
        k_ext = invitro_fitting.decay_rate_from_half_life_days(1.0)
        self.assertTrue(
            math.isclose(
                legacy_gemcitabine_pk.extracellular_gemcitabine_concentration(0.0, dose_muM, k_ext),
                0.00971 * dose_muM,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )

    def test_extracellular_concentration_halves_at_half_life(self):
        dose_muM = 25.0
        half_life_days = 1.5
        k_ext = invitro_fitting.decay_rate_from_half_life_days(half_life_days)
        c0 = legacy_gemcitabine_pk.extracellular_gemcitabine_concentration(0.0, dose_muM, k_ext)
        c_half = legacy_gemcitabine_pk.extracellular_gemcitabine_concentration(half_life_days, dose_muM, k_ext)
        self.assertTrue(math.isclose(c_half, 0.5 * c0, rel_tol=1e-12, abs_tol=1e-12))

    def test_long_half_life_recovers_constant_exposure_limit_over_five_days(self):
        curve = legacy_gemcitabine_pk.build_extracellular_exposure_curve_from_profile(
            time_days=np.array([], dtype=float),
            gem_concentration_ng_per_ml=np.array([], dtype=float),
            gem_was_censored=np.array([], dtype=bool),
            gem_censor_upper_bound_ng_per_ml=np.array([], dtype=float),
            reference_dose_muM=1.0,
            fallback_half_life_days=1e9,
            analyte_column=None,
            source_ploidy="fallback",
            preferred_mode="half_life_fallback",
        )
        dose_muM = 50.0
        t = np.linspace(0.0, 5.0, 11)
        conc = curve(t, dose_muM)
        old_constant = np.full_like(conc, legacy_gemcitabine_pk.converted_extracellular_gem_signal(dose_muM), dtype=float)
        max_relative_error = np.max(np.abs(conc - old_constant) / old_constant)
        self.assertLess(max_relative_error, 1e-8)

    def test_constrained_decay_penalizes_numeric_censor_upper_bounds(self):
        fit = legacy_gemcitabine_pk.fit_constrained_extracellular_gem_decay(
            time_days=np.array([0.0, 1.0]),
            observed_gem_uM=np.array([np.nan, np.nan]),
            censored_mask=np.array([False, True]),
            censor_upper_bound_uM=np.array([np.nan, 0.1]),
            reference_dose_muM=1.0,
        )
        self.assertIsNotNone(fit)
        assert fit is not None
        predicted = fit["C0_uM"] * math.exp(-fit["k_ext_decay_per_day"] * 1.0)
        self.assertLessEqual(predicted, 0.1 + 1e-8)
        self.assertEqual(fit["n_censoring_violations"], 0)

    def test_get_aligned_live_dead_data_inner_joins_misaligned_times(self):
        df = pd.DataFrame(
            [
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Alive", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 10},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Alive", "time_days": 1.0, "plate_row": "A", "plate_col": 1, "count": 11},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Alive", "time_days": 2.0, "plate_row": "A", "plate_col": 1, "count": 12},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Alive", "time_days": 3.0, "plate_row": "A", "plate_col": 1, "count": 13},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Dead", "time_days": 4.0, "plate_row": "A", "plate_col": 1, "count": 24},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Dead", "time_days": 2.0, "plate_row": "A", "plate_col": 1, "count": 22},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Dead", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 20},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Dead", "time_days": 3.0, "plate_row": "A", "plate_col": 1, "count": 23},
            ]
        )
        aligned = invitro_fitting.get_aligned_live_dead_data(df, gem_dose="10 nM", ploidy="2N")
        self.assertTrue(np.array_equal(aligned["t"], np.array([0.0, 2.0, 3.0])))
        self.assertTrue(np.array_equal(aligned["y_alive"][:, 0], np.array([10.0, 12.0, 13.0])))
        self.assertTrue(np.array_equal(aligned["y_dead"][:, 0], np.array([20.0, 22.0, 23.0])))
        self.assertEqual(aligned["dropped_timepoints"], 2)

    def test_get_aligned_live_dead_data_counts_transitional_toward_dead_by_default(self):
        df = pd.DataFrame(
            [
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Alive", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 10},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Transitional", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 2},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Dead", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 1},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Alive", "time_days": 1.0, "plate_row": "A", "plate_col": 1, "count": 12},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Transitional", "time_days": 1.0, "plate_row": "A", "plate_col": 1, "count": 3},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Dead", "time_days": 1.0, "plate_row": "A", "plate_col": 1, "count": 4},
            ]
        )
        aligned = invitro_fitting.get_aligned_live_dead_data(df, gem_dose="10 nM", ploidy="2N")
        self.assertTrue(np.array_equal(aligned["y_alive"][:, 0], np.array([10.0, 12.0])))
        self.assertTrue(np.array_equal(aligned["y_dead"][:, 0], np.array([3.0, 7.0])))
        self.assertFalse(aligned["count_transitional_as_alive"])

    def test_get_aligned_live_dead_data_can_count_transitional_toward_alive(self):
        df = pd.DataFrame(
            [
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Alive", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 10},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Transitional", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 2},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Dead", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 1},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Alive", "time_days": 1.0, "plate_row": "A", "plate_col": 1, "count": 12},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Transitional", "time_days": 1.0, "plate_row": "A", "plate_col": 1, "count": 3},
                {"gem": "10 nM", "ploidy": "2N", "phenotype": "Dead", "time_days": 1.0, "plate_row": "A", "plate_col": 1, "count": 4},
            ]
        )
        aligned = invitro_fitting.get_aligned_live_dead_data(
            df,
            gem_dose="10 nM",
            ploidy="2N",
            count_transitional_as_alive=True,
        )
        self.assertTrue(np.array_equal(aligned["y_alive"][:, 0], np.array([12.0, 15.0])))
        self.assertTrue(np.array_equal(aligned["y_dead"][:, 0], np.array([1.0, 4.0])))
        self.assertTrue(aligned["count_transitional_as_alive"])

    def test_plot_extracellular_exposure_curve_closes_saved_figure(self):
        curve = legacy_gemcitabine_pk.build_extracellular_exposure_curve_from_profile(
            time_days=np.array([0.0, 1.0]),
            gem_concentration_ng_per_ml=np.array([100.0, 50.0]),
            gem_was_censored=np.array([False, False]),
            gem_censor_upper_bound_ng_per_ml=np.array([np.nan, np.nan]),
            reference_dose_muM=1.0,
            fallback_half_life_days=1.0,
            analyte_column="Gemcitabine (ng/mL)",
            source_ploidy="synthetic",
            preferred_mode="constrained_bolus_exponential",
        )
        before = set(plt.get_fignums())
        with tempfile.TemporaryDirectory() as tmpdir:
            legacy_gemcitabine_pk.plot_extracellular_exposure_curve(
                "synthetic",
                curve,
                output_dir=Path(tmpdir),
                close_fig=True,
            )
        after = set(plt.get_fignums())
        self.assertEqual(before, after)

    def test_simulate_joint_ext_runs_for_synthetic_input(self):
        t = np.linspace(0.0, 5.0, 6)
        params_3 = [0.5, 1.0, 0.2, 1.0, 0.0]
        dfdctp_signal_curve = self.make_constant_surface(0.05)
        fit_config = invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False)
        alive, dead = invitro_fitting.simulate_joint_ext(
            t=t,
            params=params_3,
            N0=1000.0,
            D0=0.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.05,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=2,
            model_variant="delayed_death_only",
            fit_config=fit_config,
        )
        self.assertEqual(alive.shape, t.shape)
        self.assertEqual(dead.shape, t.shape)
        self.assertTrue(np.issubdtype(alive.dtype, np.floating))
        self.assertTrue(np.issubdtype(dead.dtype, np.floating))
        self.assertTrue(np.all(np.isfinite(alive)))
        self.assertTrue(np.all(np.isfinite(dead)))

    def test_simulate_joint_dfdctp_safe_valid_result(self):
        surface = self.make_constant_surface(0.05)
        result = invitro_fitting.simulate_joint_dfdctp_safe(
            t=np.linspace(0.0, 2.0, 5),
            params=[0.5, 1.0, 0.2, 1.0, 1e-8],
            N0=100.0,
            D0=0.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.05,
            dfdctp_signal_curve=surface,
            n_tr=2,
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
            model_variant="delayed_death_only",
        )
        self.assertTrue(result.success)
        self.assertTrue(np.all(np.isfinite(result.alive)))

    def test_simulate_joint_dfdctp_safe_invalid_params_fail(self):
        surface = self.make_constant_surface(0.05)
        result = invitro_fitting.simulate_joint_dfdctp_safe(
            t=np.linspace(0.0, 2.0, 5),
            params=[-1.0, 1.0, 0.2, 1.0, 1e-8],
            N0=100.0,
            D0=0.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.05,
            dfdctp_signal_curve=surface,
            n_tr=2,
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
            model_variant="delayed_death_only",
        )
        self.assertFalse(result.success)

    def test_zero_dose_baseline_death_creates_dead_cells(self):
        surface = self.make_constant_surface(0.05)
        t = np.linspace(0.0, 2.0, 5)
        fit_config = invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False)
        no_baseline = invitro_fitting.simulate_joint_dfdctp_safe(
            t=t,
            params=[0.5, 10.0, 1e-8, 1.0, 1e-8],
            N0=100.0,
            D0=0.0,
            r=1e-8,
            K=1e9,
            dose_muM=0.0,
            dfdctp_signal_curve=surface,
            n_tr=2,
            fit_config=fit_config,
            model_variant="delayed_death_only",
        )
        baseline = invitro_fitting.simulate_joint_dfdctp_safe(
            t=t,
            params=[0.5, 10.0, 1e-8, 1.0, 0.05],
            N0=100.0,
            D0=0.0,
            r=1e-8,
            K=1e9,
            dose_muM=0.0,
            dfdctp_signal_curve=surface,
            n_tr=2,
            fit_config=fit_config,
            model_variant="delayed_death_only",
        )
        self.assertTrue(no_baseline.success)
        self.assertTrue(baseline.success)
        self.assertEqual(float(surface(1.0, 0.0)), 0.0)
        self.assertLess(no_baseline.dead[-1], 1e-5)
        self.assertGreater(baseline.dead[-1], baseline.dead[0])
        self.assertLess(baseline.alive[-1], no_baseline.alive[-1])

    def test_treated_wells_include_baseline_and_drug_death(self):
        surface = self.make_constant_surface(0.05)
        t = np.linspace(0.0, 2.0, 5)
        fit_config = invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False)
        no_baseline = invitro_fitting.simulate_joint_dfdctp_safe(
            t=t,
            params=[0.5, 10.0, 1e-8, 1.0, 1e-8],
            N0=100.0,
            D0=0.0,
            r=1e-8,
            K=1e9,
            dose_muM=1.0,
            dfdctp_signal_curve=surface,
            n_tr=2,
            fit_config=fit_config,
            model_variant="delayed_death_only",
        )
        baseline = invitro_fitting.simulate_joint_dfdctp_safe(
            t=t,
            params=[0.5, 10.0, 1e-8, 1.0, 0.05],
            N0=100.0,
            D0=0.0,
            r=1e-8,
            K=1e9,
            dose_muM=1.0,
            dfdctp_signal_curve=surface,
            n_tr=2,
            fit_config=fit_config,
            model_variant="delayed_death_only",
        )
        self.assertTrue(no_baseline.success)
        self.assertTrue(baseline.success)
        self.assertGreaterEqual(baseline.dead[-1], no_baseline.dead[-1])
        self.assertLessEqual(baseline.alive[-1], no_baseline.alive[-1])

    def test_safe_objective_returns_penalty_on_exception(self):
        penalty = 123.0
        value = invitro_fitting.safe_objective_value(lambda x: (_ for _ in ()).throw(RuntimeError("boom")), np.array([1.0]), penalty)
        self.assertEqual(value, penalty)

    def test_residuals_global_joint_use_raw_count_units(self):
        dfdctp_signal_curve = self.make_constant_surface(0.0)
        t = np.array([0.0, 1.0])
        dose_data_list = [
            {
                "t": t,
                "y_alive": np.array([10.0, np.nan]),
                "y_dead": np.array([0.0, 1.0]),
                "N0": 10.0,
                "D0": 0.0,
                "dose_muM": 0.01,
            }
        ]
        params_3 = [1e-4, 1e-4, 1e-8, 1.0, 1e-8]
        residuals = invitro_fitting.residuals_global_joint(
            treatment_params=params_3,
            dose_data_list=dose_data_list,
            r_fixed=1e-8,
            K_fixed=100.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=2,
            fit_means_only=True,
            high_dose_weight=1.0,
            observation_channels="alive_dead",
            model_variant="delayed_death_only",
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
        )
        expected = np.array([0.0, 0.0, 1.0])
        self.assertTrue(np.allclose(residuals, expected, rtol=0.0, atol=1e-6))

    def test_poisson_nll_returns_finite_value(self):
        y = np.array([0.0, 3.0, 10.0])
        mu = np.array([0.2, 2.5, 8.0])
        nll = invitro_fitting.poisson_nll(y, mu)
        self.assertTrue(np.isfinite(nll))
        self.assertGreaterEqual(nll, 0.0)

    def test_negative_binomial_nll_returns_finite_value(self):
        y = np.array([0.0, 3.0, 10.0])
        mu = np.array([0.2, 2.5, 8.0])
        nll = invitro_fitting.negative_binomial_nll(y, mu, theta=25.0)
        self.assertTrue(np.isfinite(nll))
        self.assertGreaterEqual(nll, 0.0)

    def test_negative_binomial_nll_rejects_invalid_theta(self):
        with self.assertRaises(ValueError):
            invitro_fitting.negative_binomial_nll(np.array([1.0]), np.array([1.0]), theta=0.0)

    def test_negative_binomial_nll_approaches_poisson_for_large_theta(self):
        y = np.array([1.0, 4.0, 9.0])
        mu = np.array([1.2, 3.8, 8.7])
        poisson = invitro_fitting.poisson_nll(y, mu)
        negbin = invitro_fitting.negative_binomial_nll(y, mu, theta=1e9)
        self.assertTrue(math.isclose(negbin, poisson, rel_tol=1e-6, abs_tol=1e-6))

    def test_likelihood_objective_uses_alive_and_dead_observations(self):
        dfdctp_signal_curve = self.make_constant_surface(0.0)
        dose_data_list = [
            {
                "t": np.array([0.0, 1.0]),
                "y_alive": np.array([10.0, 9.0]),
                "y_dead": np.array([0.0, 1.0]),
                "N0": 10.0,
                "D0": 0.0,
                "dose_muM": 0.01,
            }
        ]
        params_3 = [0.1, 0.2, 0.3, 1.0, 1e-8]
        baseline_nll = invitro_fitting.live_dead_objective_nll(
            treatment_params=params_3,
            dose_data_list=dose_data_list,
            r_fixed=1e-8,
            K_fixed=100.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=2,
            objective="negative_binomial",
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
            observation_channels="alive_dead",
            fit_means_only=True,
            theta_alive=20.0,
            theta_dead=20.0,
        )
        modified = [dict(dose_data_list[0])]
        modified[0]["y_dead"] = np.array([0.0, 7.0])
        changed_nll = invitro_fitting.live_dead_objective_nll(
            treatment_params=params_3,
            dose_data_list=modified,
            r_fixed=1e-8,
            K_fixed=100.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=2,
            objective="negative_binomial",
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
            observation_channels="alive_dead",
            fit_means_only=True,
            theta_alive=20.0,
            theta_dead=20.0,
        )
        self.assertNotEqual(baseline_nll, changed_nll)

    def test_alive_only_likelihood_ignores_dead_observations(self):
        dfdctp_signal_curve = self.make_constant_surface(0.0)
        dose_data_list = [
            {
                "t": np.array([0.0, 1.0]),
                "y_alive": np.array([10.0, 9.0]),
                "y_dead": np.array([0.0, 1.0]),
                "N0": 10.0,
                "D0": 0.0,
                "dose_muM": 0.01,
            }
        ]
        params_3 = [0.1, 0.2, 0.3, 1.0, 1e-8]
        baseline_nll = invitro_fitting.live_dead_objective_nll(
            treatment_params=params_3,
            dose_data_list=dose_data_list,
            r_fixed=1e-8,
            K_fixed=100.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=2,
            objective="negative_binomial",
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
            observation_channels="alive_only",
            fit_means_only=True,
            theta_alive=20.0,
        )
        modified = [dict(dose_data_list[0])]
        modified[0]["y_dead"] = np.array([0.0, 7.0])
        changed_nll = invitro_fitting.live_dead_objective_nll(
            treatment_params=params_3,
            dose_data_list=modified,
            r_fixed=1e-8,
            K_fixed=100.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=2,
            objective="negative_binomial",
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
            observation_channels="alive_only",
            fit_means_only=True,
            theta_alive=20.0,
        )
        self.assertEqual(baseline_nll, changed_nll)

    def test_alive_only_residuals_ignore_dead_observations(self):
        dfdctp_signal_curve = self.make_constant_surface(0.0)
        t = np.array([0.0, 1.0])
        dose_data_list = [
            {
                "t": t,
                "y_alive": np.array([10.0, np.nan]),
                "y_dead": np.array([0.0, 9.0]),
                "N0": 10.0,
                "D0": 0.0,
                "dose_muM": 0.01,
            }
        ]
        residuals = invitro_fitting.residuals_global_joint(
            treatment_params=[1e-4, 1e-4, 1e-8, 1.0, 1e-8],
            dose_data_list=dose_data_list,
            r_fixed=1e-8,
            K_fixed=100.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=2,
            fit_means_only=True,
            observation_channels="alive_only",
            model_variant="delayed_death_only",
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
        )
        self.assertTrue(np.array_equal(residuals, np.array([0.0])))

    def test_likelihood_helpers_floor_zero_predictions(self):
        y = np.array([0.0, 1.0, 2.0])
        mu = np.array([0.0, 0.0, 0.0])
        self.assertTrue(np.isfinite(invitro_fitting.poisson_nll(y, mu)))
        self.assertTrue(np.isfinite(invitro_fitting.negative_binomial_nll(y, mu, theta=10.0)))

    def test_fit_live_dead_model_supports_negative_binomial_and_least_squares(self):
        dfdctp_signal_curve = self.make_constant_surface(0.02)
        dose_data_list = [
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 95.0, 90.0]),
                "y_dead": np.array([0.0, 4.0, 9.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.01,
            },
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 92.0, 82.0]),
                "y_dead": np.array([0.0, 7.0, 16.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.4,
            },
        ]
        nb_fit = invitro_fitting.fit_live_dead_model(
            dose_data_list=dose_data_list,
            r_opt=0.2,
            K_opt=1000.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=2,
            objective="negative_binomial",
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
            fit_means_only=True,
            max_nfev=200,
        )
        ls_fit = invitro_fitting.fit_live_dead_model(
            dose_data_list=dose_data_list,
            r_opt=0.2,
            K_opt=1000.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=2,
            objective="least_squares",
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
            fit_means_only=True,
            max_nfev=200,
        )
        self.assertIsNotNone(nb_fit)
        self.assertIsNotNone(ls_fit)
        self.assertEqual(nb_fit["objective"], "negative_binomial")
        self.assertEqual(ls_fit["objective"], "least_squares")
        self.assertTrue(np.isfinite(nb_fit["summary"]["nll"]))
        self.assertTrue(np.isfinite(nb_fit["summary"]["aic"]))
        self.assertTrue(np.isfinite(nb_fit["summary"]["bic"]))
        self.assertIsNone(ls_fit["summary"]["nll"])
        self.assertTrue(np.isfinite(ls_fit["summary"]["objective_value"]))

    def test_alive_only_negative_binomial_reports_one_dispersion_parameter(self):
        dfdctp_signal_curve = self.make_constant_surface(0.02)
        dose_data_list = [
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 95.0, 90.0]),
                "y_dead": np.array([0.0, 4.0, 9.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.01,
            }
        ]
        nb_fit = invitro_fitting.fit_live_dead_model(
            dose_data_list=dose_data_list,
            r_opt=0.2,
            K_opt=1000.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=2,
            objective="negative_binomial",
            observation_channels="alive_only",
            fit_config=invitro_fitting.JointFitConfig(model_variant="delayed_death_only", use_hill_dose_gate=False, use_confluence_death=False),
            fit_means_only=True,
            max_nfev=3,
        )
        self.assertEqual(nb_fit["summary"]["observation_channels"], "alive_only")
        self.assertIsNotNone(nb_fit["summary"]["theta_alive"])
        self.assertIsNone(nb_fit["summary"]["theta_dead"])
        self.assertEqual(nb_fit["summary"]["n_parameters"], 6)

    def test_legacy_fit_live_dead_model_supports_fixed_hill(self):
        dfdctp_signal_curve = self.make_constant_surface(0.02)
        dose_data_list = [
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 95.0, 89.0]),
                "y_dead": np.array([0.0, 4.0, 10.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.01,
            },
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 92.0, 82.0]),
                "y_dead": np.array([0.0, 7.0, 16.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.05,
            },
        ]
        fit = invitro_fitting.fit_live_dead_model(
            dose_data_list=dose_data_list,
            r_opt=0.2,
            K_opt=1000.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=2,
            objective="negative_binomial",
            fit_config=invitro_fitting.JointFitConfig(
                model_variant="delayed_death_only",
                use_hill_dose_gate=True,
                fit_hill_dose_gate=False,
                fixed_dose_gate_ec50_uM=0.0125,
                fixed_dose_gate_hill=2.0,
                use_confluence_death=False,
            ),
            fit_means_only=True,
            max_nfev=3,
        )
        self.assertIsNotNone(fit)
        self.assertTrue(np.isfinite(fit["summary"]["dose_gate_ec50_uM"]))
        self.assertTrue(np.isfinite(fit["summary"]["dose_gate_hill"]))
        self.assertFalse(fit["summary"]["fit_hill_dose_gate"])

    def test_legacy_fit_live_dead_model_supports_fitted_hill(self):
        dfdctp_signal_curve = self.make_constant_surface(0.02)
        dose_data_list = [
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 95.0, 89.0]),
                "y_dead": np.array([0.0, 4.0, 10.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.01,
            },
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 92.0, 82.0]),
                "y_dead": np.array([0.0, 7.0, 16.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.05,
            },
        ]
        fit = invitro_fitting.fit_live_dead_model(
            dose_data_list=dose_data_list,
            r_opt=0.2,
            K_opt=1000.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=2,
            objective="negative_binomial",
            fit_config=invitro_fitting.JointFitConfig(
                model_variant="delayed_death_only",
                use_hill_dose_gate=True,
                fit_hill_dose_gate=True,
                use_confluence_death=False,
            ),
            fit_means_only=True,
            max_nfev=3,
        )
        self.assertIsNotNone(fit)
        self.assertTrue(np.isfinite(fit["summary"]["dose_gate_ec50_uM"]))
        self.assertTrue(np.isfinite(fit["summary"]["dose_gate_hill"]))
        self.assertTrue(fit["summary"]["fit_hill_dose_gate"])

    def test_legacy_fit_live_dead_model_supports_hill_only_mode(self):
        dfdctp_signal_curve = self.make_constant_surface(0.02)
        dose_data_list = [
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 95.0, 89.0]),
                "y_dead": np.array([0.0, 4.0, 10.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.01,
            },
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 92.0, 82.0]),
                "y_dead": np.array([0.0, 7.0, 16.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.05,
            },
        ]
        fit = invitro_fitting.fit_live_dead_model(
            dose_data_list=dose_data_list,
            r_opt=0.2,
            K_opt=1000.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=2,
            objective="negative_binomial",
            fit_config=invitro_fitting.JointFitConfig(
                model_variant="delayed_death_only",
                use_hill_dose_gate=True,
                fit_hill_dose_gate=True,
                fit_beta_dose=False,
                fixed_beta_dose=1.0,
                use_confluence_death=False,
            ),
            fit_means_only=True,
            max_nfev=3,
        )
        self.assertIsNotNone(fit)
        self.assertAlmostEqual(fit["summary"]["beta_dose"], 1.0)
        self.assertTrue(np.isfinite(fit["summary"]["dose_gate_ec50_uM"]))
        self.assertTrue(np.isfinite(fit["summary"]["dose_gate_hill"]))

    def test_legacy_fit_live_dead_model_hill_and_confluence_are_compatible(self):
        dfdctp_signal_curve = self.make_constant_surface(0.02)
        dose_data_list = [
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 95.0, 89.0]),
                "y_dead": np.array([0.0, 4.0, 10.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.01,
            },
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([100.0, 92.0, 82.0]),
                "y_dead": np.array([0.0, 7.0, 16.0]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.05,
            },
        ]
        fit = invitro_fitting.fit_live_dead_model(
            dose_data_list=dose_data_list,
            r_opt=0.2,
            K_opt=1000.0,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=2,
            objective="negative_binomial",
            fit_config=invitro_fitting.JointFitConfig(
                model_variant="delayed_death_only",
                use_hill_dose_gate=True,
                fit_hill_dose_gate=False,
                use_confluence_death=True,
                fit_mu_confluence_death=True,
            ),
            fit_means_only=True,
            max_nfev=100,
        )
        self.assertIsNotNone(fit)
        self.assertTrue(fit["summary"]["use_hill_dose_gate"])
        self.assertTrue(fit["summary"]["use_confluence_death"])
        self.assertIn("mu_confluence_death", fit["summary"])

    def test_fit_joint_one_replicate_accepts_current_default_config(self):
        surface = self.make_constant_surface(0.02)
        dose_data_list = [
            {
                "t": np.array([0.0, 1.0, 2.0]),
                "y_alive": np.array([[100.0], [95.0], [89.0]]),
                "y_dead": np.array([[0.0], [4.0], [10.0]]),
                "N0": 100.0,
                "D0": 0.0,
                "dose_muM": 0.01,
                "dose_label": "10 nM",
            }
        ]
        fake_summary = {
            "objective_value": 1.0,
            "nll": 1.0,
            "aic": 2.0,
            "bic": 3.0,
            "rmse": 0.1,
            "theta_alive": 20.0,
            "theta_dead": 20.0,
            "use_hill_dose_gate": True,
            "fit_hill_dose_gate": True,
            "dose_gate_ec50_uM": 0.0125,
            "dose_gate_hill": 2.0,
            "beta_dose": 1.0,
            "use_confluence_death": True,
            "mu_confluence_death": 0.05,
            "confluence_death_exponent": 4.0,
        }
        fake_result = {
            "summary": fake_summary,
            "treatment_params": np.array([0.5, 25.0, 0.5, 1.0, 1.0, 0.02], dtype=float),
            "objective": "negative_binomial",
            "success": True,
        }
        with patch.object(invitro_fitting, "fit_live_dead_model", return_value=fake_result):
            result = invitro_fitting.fit_joint_one_replicate(
                dose_data_list=dose_data_list,
                r_opt=0.2,
                K_opt=1000.0,
                dfdctp_signal_curve=surface,
                ploidy="2N",
                rep_idx="A",
                objective="negative_binomial",
                observation_channels="alive_dead",
                fit_config=invitro_fitting.JointFitConfig(),
            )
        self.assertIsNotNone(result)

    def test_build_joint_fit_trajectories_uses_replicate_n0_d0(self):
        df = pd.DataFrame(
            [
                {"gem": "0 nM", "ploidy": "2N", "phenotype": "Alive", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 10},
                {"gem": "0 nM", "ploidy": "2N", "phenotype": "Transitional", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 2},
                {"gem": "0 nM", "ploidy": "2N", "phenotype": "Alive", "time_days": 1.0, "plate_row": "A", "plate_col": 1, "count": 11},
                {"gem": "0 nM", "ploidy": "2N", "phenotype": "Transitional", "time_days": 1.0, "plate_row": "A", "plate_col": 1, "count": 3},
                {"gem": "0 nM", "ploidy": "2N", "phenotype": "Dead", "time_days": 0.0, "plate_row": "A", "plate_col": 1, "count": 1},
                {"gem": "0 nM", "ploidy": "2N", "phenotype": "Dead", "time_days": 1.0, "plate_row": "A", "plate_col": 1, "count": 2},
            ]
        )
        trajectories = invitro_fitting.build_joint_fit_trajectories(df, ploidies=["2N"], gem_doses=["0 nM"], fit_t_max=5.0)
        self.assertEqual(len(trajectories), 1)
        self.assertEqual(trajectories[0].N0, 10.0)
        self.assertEqual(trajectories[0].D0, 3.0)

    def test_trim_finite_live_dead_observations_drops_nonfinite_pairs(self):
        t, alive, dead, dropped = invitro_fitting.trim_finite_live_dead_observations(
            t=np.array([0.0, 1.0, 2.0]),
            alive=np.array([10.0, 11.0, np.nan]),
            dead=np.array([1.0, 2.0, 3.0]),
        )
        self.assertTrue(np.array_equal(t, np.array([0.0, 1.0])))
        self.assertTrue(np.array_equal(alive, np.array([10.0, 11.0])))
        self.assertTrue(np.array_equal(dead, np.array([1.0, 2.0])))
        self.assertEqual(dropped, 1)

    def test_negative_binomial_rejects_fit_means_only_for_partial_pooling(self):
        surface = self.make_constant_surface(0.05)
        trajectories = [
            invitro_fitting.ReplicateTrajectory(
                ploidy="2N",
                dose_label="0 nM",
                dose_uM=0.0,
                replicate_id=("A", 1),
                t=np.array([0.0, 1.0]),
                alive=np.array([10.0, 11.0]),
                dead=np.array([0.0, 1.0]),
                N0=10.0,
                D0=0.0,
            )
        ]
        bad_config = invitro_fitting.JointFitConfig(fit_means_only=True)
        with self.assertRaises(ValueError):
            invitro_fitting.fit_joint_partial_pooling_model(
                trajectories,
                dfdctp_signal_curve_by_ploidy={"2N": surface},
                fit_config=bad_config,
                n_tr=2,
            )

    def test_joint_partial_pooling_returns_population_and_ploidy_parameters(self):
        surface = self.make_constant_surface(0.05)
        trajectories = [
            invitro_fitting.ReplicateTrajectory(
                ploidy="2N",
                dose_label="0 nM",
                dose_uM=0.0,
                replicate_id=("A", 1),
                t=np.array([0.0, 1.0]),
                alive=np.array([10.0, 11.0]),
                dead=np.array([0.0, 0.5]),
                N0=10.0,
                D0=0.0,
            ),
            invitro_fitting.ReplicateTrajectory(
                ploidy="4N",
                dose_label="0 nM",
                dose_uM=0.0,
                replicate_id=("E", 1),
                t=np.array([0.0, 1.0]),
                alive=np.array([12.0, 13.0]),
                dead=np.array([0.0, 0.3]),
                N0=12.0,
                D0=0.0,
            ),
        ]
        config = invitro_fitting.JointFitConfig(max_nfev=20, use_hill_dose_gate=False, use_confluence_death=False)
        result = invitro_fitting.fit_joint_partial_pooling_model(
            trajectories,
            dfdctp_signal_curve_by_ploidy={"2N": surface, "4N": surface},
            fit_config=config,
            n_tr=2,
        )
        self.assertIn("r", result["population_parameters"])
        self.assertIn("beta_dose", result["population_parameters"])
        self.assertIn("mu_base_death", result["population_parameters"])
        self.assertIn("2N", result["ploidy_parameters"])
        self.assertIn("4N", result["ploidy_parameters"])
        self.assertIn("beta_dose", result["ploidy_parameters"]["2N"])
        self.assertIn("beta_dose", result["ploidy_parameters"]["4N"])
        self.assertIn("mu_base_death", result["ploidy_parameters"]["2N"])
        self.assertIn("mu_base_death", result["ploidy_parameters"]["4N"])

    def test_partial_pooling_outputs_k_cyto_for_model2(self):
        surface = self.make_constant_surface(0.05)
        trajectories = [
            invitro_fitting.ReplicateTrajectory(
                ploidy="2N",
                dose_label="0 nM",
                dose_uM=0.0,
                replicate_id=("A", 1),
                t=np.array([0.0, 1.0]),
                alive=np.array([10.0, 11.0]),
                dead=np.array([0.0, 0.5]),
                N0=10.0,
                D0=0.0,
            ),
            invitro_fitting.ReplicateTrajectory(
                ploidy="4N",
                dose_label="0 nM",
                dose_uM=0.0,
                replicate_id=("E", 1),
                t=np.array([0.0, 1.0]),
                alive=np.array([12.0, 13.0]),
                dead=np.array([0.0, 0.3]),
                N0=12.0,
                D0=0.0,
            ),
        ]
        config = invitro_fitting.JointFitConfig(
            max_nfev=50,
            objective="least_squares",
            model_variant="immediate_cytostasis_delayed_death",
            use_hill_dose_gate=False,
            use_confluence_death=False,
        )
        result = invitro_fitting.fit_joint_partial_pooling_model(
            trajectories,
            dfdctp_signal_curve_by_ploidy={"2N": surface, "4N": surface},
            fit_config=config,
            n_tr=2,
        )
        self.assertIn("2N", result["ploidy_parameters"])
        self.assertIn("4N", result["ploidy_parameters"])
        self.assertIn("k_cyto", result["ploidy_parameters"]["2N"])
        self.assertIn("k_cyto", result["ploidy_parameters"]["4N"])
        self.assertIn("beta_dose", result["ploidy_parameters"]["2N"])
        self.assertIn("beta_dose", result["ploidy_parameters"]["4N"])
        self.assertIn("mu_base_death", result["ploidy_parameters"]["2N"])
        self.assertIn("mu_base_death", result["ploidy_parameters"]["4N"])

    def test_joint_partial_pooling_uses_four_starts_for_default_model2(self):
        surface = self.make_constant_surface(0.05)
        trajectories = [
            invitro_fitting.ReplicateTrajectory(
                ploidy="2N",
                dose_label="0 nM",
                dose_uM=0.0,
                replicate_id=("A", 1),
                t=np.array([0.0, 1.0]),
                alive=np.array([10.0, 11.0]),
                dead=np.array([0.0, 0.5]),
                N0=10.0,
                D0=0.0,
            ),
            invitro_fitting.ReplicateTrajectory(
                ploidy="4N",
                dose_label="0 nM",
                dose_uM=0.0,
                replicate_id=("E", 1),
                t=np.array([0.0, 1.0]),
                alive=np.array([12.0, 13.0]),
                dead=np.array([0.0, 0.3]),
                N0=12.0,
                D0=0.0,
            ),
        ]
        config = invitro_fitting.JointFitConfig(
            max_nfev=1,
            objective="least_squares",
            model_variant="immediate_cytostasis_delayed_death",
            use_hill_dose_gate=False,
        )
        result = invitro_fitting.fit_joint_partial_pooling_model(
            trajectories,
            dfdctp_signal_curve_by_ploidy={"2N": surface, "4N": surface},
            fit_config=config,
            n_tr=2,
        )
        self.assertEqual(len(result["optimizer_attempts"]), 4)

    def test_apply_dose_power_correction_beta_one_preserves_signal(self):
        for dose_uM in (0.01, 0.1, 1.0):
            corrected = invitro_fitting.apply_dose_power_correction(
                signal_uM=0.25,
                dose_uM=dose_uM,
                reference_dose_uM=0.1,
                beta_dose=1.0,
            )
            self.assertAlmostEqual(corrected, 0.25)

    def test_apply_dose_power_correction_sublinear_and_supralinear_behavior(self):
        ref = 0.1
        signal = 1.0

        below_ref_sublinear = invitro_fitting.apply_dose_power_correction(signal, 0.01, ref, 0.5)
        below_ref_linear = invitro_fitting.apply_dose_power_correction(signal, 0.01, ref, 1.0)
        below_ref_supralinear = invitro_fitting.apply_dose_power_correction(signal, 0.01, ref, 2.0)
        self.assertGreater(below_ref_sublinear, below_ref_linear)
        self.assertLess(below_ref_supralinear, below_ref_linear)

        above_ref_sublinear = invitro_fitting.apply_dose_power_correction(signal, 1.0, ref, 0.5)
        above_ref_linear = invitro_fitting.apply_dose_power_correction(signal, 1.0, ref, 1.0)
        above_ref_supralinear = invitro_fitting.apply_dose_power_correction(signal, 1.0, ref, 2.0)
        self.assertLess(above_ref_sublinear, above_ref_linear)
        self.assertGreater(above_ref_supralinear, above_ref_linear)

    def test_apply_dose_power_correction_zero_dose_returns_zero(self):
        corrected = invitro_fitting.apply_dose_power_correction(
            signal_uM=0.25,
            dose_uM=0.0,
            reference_dose_uM=0.1,
            beta_dose=1.5,
        )
        self.assertEqual(corrected, 0.0)

    def test_joint_fit_config_from_preset_expected_values(self):
        minimal = invitro_fitting.joint_fit_config_from_preset("minimal")
        self.assertFalse(minimal.use_hill_dose_gate)
        self.assertFalse(minimal.fit_beta_dose)
        self.assertFalse(minimal.use_confluence_death)

        beta_only = invitro_fitting.joint_fit_config_from_preset("beta_only")
        self.assertFalse(beta_only.use_hill_dose_gate)
        self.assertTrue(beta_only.fit_beta_dose)
        self.assertFalse(beta_only.use_confluence_death)

        hill_only = invitro_fitting.joint_fit_config_from_preset("hill_only")
        self.assertTrue(hill_only.use_hill_dose_gate)
        self.assertTrue(hill_only.fit_hill_dose_gate)
        self.assertFalse(hill_only.fit_beta_dose)
        self.assertEqual(hill_only.fixed_beta_dose, 1.0)
        self.assertFalse(hill_only.use_confluence_death)

        beta_hill_conf = invitro_fitting.joint_fit_config_from_preset("beta_hill_baseline_confluence")
        self.assertTrue(beta_hill_conf.use_hill_dose_gate)
        self.assertTrue(beta_hill_conf.fit_beta_dose)
        self.assertTrue(beta_hill_conf.use_confluence_death)
        self.assertTrue(beta_hill_conf.fit_mu_confluence_death)

    def test_joint_fit_config_preset_manual_override_works(self):
        cfg = invitro_fitting.joint_fit_config_from_preset(
            "hill_only",
            fixed_dose_gate_ec50_uM=0.025,
        )
        self.assertTrue(cfg.use_hill_dose_gate)
        self.assertFalse(cfg.fit_beta_dose)
        self.assertAlmostEqual(cfg.fixed_dose_gate_ec50_uM, 0.025)

    def test_default_config_is_not_most_complex_model(self):
        cfg = invitro_fitting.JointFitConfig()
        self.assertFalse(cfg.use_hill_dose_gate)
        self.assertFalse(cfg.use_confluence_death)

    def test_resolve_beta_dose_for_fixed_mode(self):
        cfg = invitro_fitting.JointFitConfig(
            fit_beta_dose=False,
            fixed_beta_dose=1.0,
        )
        self.assertNotIn("beta_dose", invitro_fitting.get_parameter_names_for_config(cfg))
        beta = invitro_fitting.resolve_beta_dose_for_params({}, cfg)
        self.assertEqual(beta, 1.0)

    def test_resolve_mu_confluence_death_disabled_and_fixed_modes(self):
        cfg_off = invitro_fitting.JointFitConfig(use_confluence_death=False)
        self.assertEqual(invitro_fitting.resolve_mu_confluence_death_for_params({}, cfg_off), 0.0)
        self.assertNotIn("mu_confluence_death", invitro_fitting.get_parameter_names_for_config(cfg_off))

        cfg_fixed = invitro_fitting.JointFitConfig(
            use_confluence_death=True,
            fit_mu_confluence_death=False,
            fixed_mu_confluence_death=0.05,
        )
        self.assertNotIn("mu_confluence_death", invitro_fitting.get_parameter_names_for_config(cfg_fixed))
        self.assertAlmostEqual(
            invitro_fitting.resolve_mu_confluence_death_for_params({}, cfg_fixed),
            0.05,
        )

    def test_resolve_hill_dose_gate_params_disabled_and_fixed_modes(self):
        unpacked = {"dose_gate_params": {"dose_gate_ec50_uM": 0.02, "dose_gate_hill": 3.0}}
        cfg_off = invitro_fitting.JointFitConfig(use_hill_dose_gate=False, fit_beta_dose=False)
        resolved_off = invitro_fitting.resolve_hill_dose_gate_params_for_unpacked(unpacked, cfg_off)
        self.assertFalse(resolved_off["use_hill_dose_gate"])

        cfg_fixed = invitro_fitting.JointFitConfig(
            use_hill_dose_gate=True,
            fit_hill_dose_gate=False,
            fixed_dose_gate_ec50_uM=0.0125,
            fixed_dose_gate_hill=2.0,
        )
        resolved_fixed = invitro_fitting.resolve_hill_dose_gate_params_for_unpacked(unpacked, cfg_fixed)
        self.assertTrue(resolved_fixed["use_hill_dose_gate"])
        self.assertAlmostEqual(resolved_fixed["dose_gate_ec50_uM"], 0.0125)
        self.assertAlmostEqual(resolved_fixed["dose_gate_hill"], 2.0)

    def test_simulate_joint_dfdctp_wrapper_matches_safe_for_fitted_hill_values(self):
        surface = self.make_constant_surface(0.05)
        cfg = invitro_fitting.JointFitConfig(
            model_variant="delayed_death_only",
            use_hill_dose_gate=False,
            use_confluence_death=False,
        )
        params = [0.3, 20.0, 0.2, 1.0, 0.02]
        t = np.linspace(0.0, 2.0, 5)
        safe = invitro_fitting.simulate_joint_dfdctp_safe(
            t=t,
            params=params,
            N0=100.0,
            D0=0.0,
            r=0.4,
            K=1000.0,
            dose_muM=0.05,
            dfdctp_signal_curve=surface,
            n_tr=2,
            fit_config=cfg,
            use_hill_dose_gate=True,
            dose_gate_ec50_uM=0.02,
            dose_gate_hill=3.0,
        )
        alive, dead = invitro_fitting.simulate_joint_dfdctp(
            t=t,
            params=params,
            N0=100.0,
            D0=0.0,
            r=0.4,
            K=1000.0,
            dose_muM=0.05,
            dfdctp_signal_curve=surface,
            n_tr=2,
            model_variant="delayed_death_only",
            fit_config=cfg,
            use_hill_dose_gate=True,
            dose_gate_ec50_uM=0.02,
            dose_gate_hill=3.0,
        )
        self.assertTrue(safe.success)
        np.testing.assert_allclose(alive, safe.alive, atol=1e-5, rtol=1e-6)
        np.testing.assert_allclose(dead, safe.dead, atol=1e-5, rtol=1e-6)

    def test_simulate_joint_dfdctp_wrapper_hill_values_change_predictions(self):
        surface = self.make_constant_surface(0.05)
        cfg = invitro_fitting.JointFitConfig(
            model_variant="delayed_death_only",
            use_hill_dose_gate=False,
            use_confluence_death=False,
        )
        params = [0.3, 20.0, 0.2, 1.0, 0.02]
        t = np.linspace(0.0, 2.0, 5)
        alive_1, dead_1 = invitro_fitting.simulate_joint_dfdctp(
            t=t, params=params, N0=100.0, D0=0.0, r=0.4, K=1000.0, dose_muM=0.05,
            dfdctp_signal_curve=surface, n_tr=2, fit_config=cfg,
            model_variant="delayed_death_only",
            use_hill_dose_gate=True, dose_gate_ec50_uM=0.02, dose_gate_hill=3.0,
        )
        alive_2, dead_2 = invitro_fitting.simulate_joint_dfdctp(
            t=t, params=params, N0=100.0, D0=0.0, r=0.4, K=1000.0, dose_muM=0.05,
            dfdctp_signal_curve=surface, n_tr=2, fit_config=cfg,
            model_variant="delayed_death_only",
            use_hill_dose_gate=True, dose_gate_ec50_uM=0.005, dose_gate_hill=1.0,
        )
        self.assertFalse(np.allclose(alive_1, alive_2))
        self.assertFalse(np.allclose(dead_1, dead_2))

    def test_simulate_joint_dfdctp_wrapper_hill_disabled_ignores_hill_values(self):
        surface = self.make_constant_surface(0.05)
        cfg = invitro_fitting.JointFitConfig(
            model_variant="delayed_death_only",
            use_hill_dose_gate=False,
            use_confluence_death=False,
        )
        params = [0.3, 20.0, 0.2, 1.0, 0.02]
        t = np.linspace(0.0, 2.0, 5)
        alive_1, dead_1 = invitro_fitting.simulate_joint_dfdctp(
            t=t, params=params, N0=100.0, D0=0.0, r=0.4, K=1000.0, dose_muM=0.05,
            dfdctp_signal_curve=surface, n_tr=2, fit_config=cfg,
            model_variant="delayed_death_only",
            use_hill_dose_gate=False, dose_gate_ec50_uM=0.02, dose_gate_hill=3.0,
        )
        alive_2, dead_2 = invitro_fitting.simulate_joint_dfdctp(
            t=t, params=params, N0=100.0, D0=0.0, r=0.4, K=1000.0, dose_muM=0.05,
            dfdctp_signal_curve=surface, n_tr=2, fit_config=cfg,
            model_variant="delayed_death_only",
            use_hill_dose_gate=False, dose_gate_ec50_uM=0.005, dose_gate_hill=1.0,
        )
        np.testing.assert_allclose(alive_1, alive_2)
        np.testing.assert_allclose(dead_1, dead_2)

    def test_simulate_joint_dfdctp_safe_accepts_beta_dose_for_model2(self):
        surface = self.make_constant_surface(0.05)
        result = invitro_fitting.simulate_joint_dfdctp_safe(
            t=np.linspace(0.0, 2.0, 5),
            params=[0.5, 1.0, 0.2, 10.0, 1.0, 1e-8],
            N0=100.0,
            D0=0.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.1,
            dfdctp_signal_curve=surface,
            n_tr=2,
            model_variant="immediate_cytostasis_delayed_death",
            fit_config=invitro_fitting.JointFitConfig(model_variant="immediate_cytostasis_delayed_death", use_hill_dose_gate=False, use_confluence_death=False),
        )
        self.assertTrue(result.success)
        self.assertEqual(result.alive.shape, (5,))
        self.assertEqual(result.dead.shape, (5,))

    def test_confluence_death_hazard_basic_behavior(self):
        self.assertEqual(
            invitro_fitting.confluence_death_hazard(0.0, 100.0, 0.2, 4.0),
            0.0,
        )
        self.assertEqual(
            invitro_fitting.confluence_death_hazard(50.0, 100.0, 0.0, 4.0),
            0.0,
        )
        low = invitro_fitting.confluence_death_hazard(25.0, 100.0, 0.2, 4.0)
        high = invitro_fitting.confluence_death_hazard(100.0, 100.0, 0.2, 4.0)
        self.assertLess(low, high)
        self.assertAlmostEqual(high, 0.2)
        with self.assertRaises(ValueError):
            invitro_fitting.confluence_death_hazard(10.0, 0.0, 0.2, 4.0)

    def test_low_density_confluence_death_is_small(self):
        low = invitro_fitting.confluence_death_hazard(25.0, 100.0, 0.2, 4.0)
        confluent = invitro_fitting.confluence_death_hazard(100.0, 100.0, 0.2, 4.0)
        self.assertLess(low, confluent * 0.01)

    def test_zero_dose_confluence_death_creates_late_dead_cells(self):
        surface = self.make_constant_surface(0.0)
        config = invitro_fitting.JointFitConfig(
            model_variant="immediate_cytostasis_delayed_death",
            use_hill_dose_gate=False,
            use_confluence_death=True,
            fit_mu_confluence_death=False,
            fixed_mu_confluence_death=0.4,
            confluence_death_exponent=4.0,
        )
        t = np.linspace(0.0, 6.0, 13)
        result = invitro_fitting.simulate_joint_dfdctp_safe(
            t=t,
            params=[0.5, 1.0, 0.2, 1.0, 1.0, 0.0],
            N0=50.0,
            D0=0.0,
            r=1.0,
            K=100.0,
            dose_muM=0.0,
            dfdctp_signal_curve=surface,
            n_tr=2,
            model_variant="immediate_cytostasis_delayed_death",
            fit_config=config,
        )
        self.assertTrue(result.success)
        dead_increments = np.diff(result.dead)
        self.assertGreater(dead_increments[-1], dead_increments[0])

    def test_disabling_confluence_death_recovers_previous_behavior(self):
        surface = self.make_constant_surface(0.0)
        t = np.linspace(0.0, 3.0, 7)
        config_disabled = invitro_fitting.JointFitConfig(
            model_variant="immediate_cytostasis_delayed_death",
            use_hill_dose_gate=False,
            use_confluence_death=False,
        )
        config_disabled_nonzero = invitro_fitting.JointFitConfig(
            model_variant="immediate_cytostasis_delayed_death",
            use_hill_dose_gate=False,
            use_confluence_death=False,
            fixed_mu_confluence_death=0.5,
        )
        params = [0.5, 1.0, 0.2, 10.0, 1.0, 0.02]
        sim_a = invitro_fitting.simulate_joint_dfdctp_safe(
            t=t, params=params, N0=50.0, D0=0.0, r=0.8, K=100.0, dose_muM=0.0,
            dfdctp_signal_curve=surface, n_tr=2,
            model_variant="immediate_cytostasis_delayed_death", fit_config=config_disabled
        )
        sim_b = invitro_fitting.simulate_joint_dfdctp_safe(
            t=t, params=params, N0=50.0, D0=0.0, r=0.8, K=100.0, dose_muM=0.0,
            dfdctp_signal_curve=surface, n_tr=2,
            model_variant="immediate_cytostasis_delayed_death", fit_config=config_disabled_nonzero
        )
        self.assertTrue(sim_a.success)
        self.assertTrue(sim_b.success)
        np.testing.assert_allclose(sim_a.alive, sim_b.alive)
        np.testing.assert_allclose(sim_a.dead, sim_b.dead)

    def test_mu_confluence_death_parameter_name_appears_only_when_enabled_and_fitted(self):
        cfg_off = invitro_fitting.JointFitConfig(use_confluence_death=False)
        cfg_fixed = invitro_fitting.JointFitConfig(use_confluence_death=True, fit_mu_confluence_death=False)
        cfg_fit = invitro_fitting.JointFitConfig(use_confluence_death=True, fit_mu_confluence_death=True)
        self.assertNotIn(
            invitro_fitting.CONFLUENCE_DEATH_PARAMETER_NAME,
            invitro_fitting.get_parameter_names_for_config(cfg_off),
        )
        self.assertNotIn(
            invitro_fitting.CONFLUENCE_DEATH_PARAMETER_NAME,
            invitro_fitting.get_parameter_names_for_config(cfg_fixed),
        )
        self.assertIn(
            invitro_fitting.CONFLUENCE_DEATH_PARAMETER_NAME,
            invitro_fitting.get_parameter_names_for_config(cfg_fit),
        )

    def test_simulate_joint_dfdctp_safe_accepts_confluence_death_for_model2(self):
        surface = self.make_constant_surface(0.05)
        config = invitro_fitting.JointFitConfig(
            model_variant="immediate_cytostasis_delayed_death",
            use_hill_dose_gate=False,
            use_confluence_death=True,
            fit_mu_confluence_death=True,
        )
        result = invitro_fitting.simulate_joint_dfdctp_safe(
            t=np.linspace(0.0, 2.0, 5),
            params=[0.5, 1.0, 0.2, 10.0, 1.0, 1e-8, 0.05],
            N0=100.0,
            D0=0.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.1,
            dfdctp_signal_curve=surface,
            n_tr=2,
            model_variant="immediate_cytostasis_delayed_death",
            fit_config=config,
        )
        self.assertTrue(result.success)

    def test_partial_pooling_outputs_mu_confluence_death_when_enabled(self):
        surface = self.make_constant_surface(0.0)
        trajectories = [
            invitro_fitting.ReplicateTrajectory(
                ploidy="2N", dose_label="0 nM", dose_uM=0.0, replicate_id=("A", 1),
                t=np.array([0.0, 1.0]), alive=np.array([10.0, 11.0]), dead=np.array([0.0, 0.5]),
                N0=10.0, D0=0.0,
            ),
            invitro_fitting.ReplicateTrajectory(
                ploidy="4N", dose_label="0 nM", dose_uM=0.0, replicate_id=("E", 1),
                t=np.array([0.0, 1.0]), alive=np.array([12.0, 13.0]), dead=np.array([0.0, 0.3]),
                N0=12.0, D0=0.0,
            ),
        ]
        config = invitro_fitting.JointFitConfig(
            max_nfev=50,
            objective="negative_binomial",
            model_variant="immediate_cytostasis_delayed_death",
            use_hill_dose_gate=False,
            use_confluence_death=True,
            fit_mu_confluence_death=True,
        )
        result = invitro_fitting.fit_joint_partial_pooling_model(
            trajectories,
            dfdctp_signal_curve_by_ploidy={"2N": surface, "4N": surface},
            fit_config=config,
            n_tr=2,
            start_indices=(0,),
        )
        self.assertIn(invitro_fitting.CONFLUENCE_DEATH_PARAMETER_NAME, result["population_parameters"])
        self.assertIn(invitro_fitting.CONFLUENCE_DEATH_PARAMETER_NAME, result["ploidy_parameters"]["2N"])
        self.assertIn(invitro_fitting.CONFLUENCE_DEATH_PARAMETER_NAME, result["ploidy_parameters"]["4N"])

    def test_joint_partial_pooling_objective_includes_control_and_treatment(self):
        surface = self.make_constant_surface(0.05)
        trajectories = [
            invitro_fitting.ReplicateTrajectory(
                ploidy="2N",
                dose_label="0 nM",
                dose_uM=0.0,
                replicate_id=("A", 1),
                t=np.array([0.0, 1.0]),
                alive=np.array([10.0, 11.0]),
                dead=np.array([0.0, 0.5]),
                N0=10.0,
                D0=0.0,
            ),
            invitro_fitting.ReplicateTrajectory(
                ploidy="2N",
                dose_label="12.5 nM",
                dose_uM=0.0125,
                replicate_id=("A", 2),
                t=np.array([0.0, 1.0]),
                alive=np.array([10.0, 9.0]),
                dead=np.array([0.0, 1.5]),
                N0=10.0,
                D0=0.0,
            ),
        ]
        config = invitro_fitting.JointFitConfig(max_nfev=10, use_hill_dose_gate=False, use_confluence_death=False)
        result = invitro_fitting.fit_joint_partial_pooling_model(
            trajectories,
            dfdctp_signal_curve_by_ploidy={"2N": surface},
            fit_config=config,
            n_tr=2,
        )
        self.assertEqual(result["n_observations"], 8)

    def test_joint_partial_pooling_ignores_nonfinite_observations_in_nll(self):
        surface = self.make_constant_surface(0.05)
        trajectories = [
            invitro_fitting.ReplicateTrajectory(
                ploidy="2N",
                dose_label="0 nM",
                dose_uM=0.0,
                replicate_id=("A", 1),
                t=np.array([0.0, 1.0, 2.0]),
                alive=np.array([10.0, 11.0, np.nan]),
                dead=np.array([0.0, 0.5, np.nan]),
                N0=10.0,
                D0=0.0,
            ),
            invitro_fitting.ReplicateTrajectory(
                ploidy="4N",
                dose_label="0 nM",
                dose_uM=0.0,
                replicate_id=("E", 1),
                t=np.array([0.0, 1.0]),
                alive=np.array([12.0, 13.0]),
                dead=np.array([0.0, 0.3]),
                N0=12.0,
                D0=0.0,
            ),
        ]
        config = invitro_fitting.JointFitConfig(max_nfev=10, use_hill_dose_gate=False, use_confluence_death=False)
        result = invitro_fitting.fit_joint_partial_pooling_model(
            trajectories,
            dfdctp_signal_curve_by_ploidy={"2N": surface, "4N": surface},
            fit_config=config,
            n_tr=2,
        )
        self.assertTrue(np.isfinite(result["posterior_objective"]) or not result["success"])
        self.assertEqual(result["n_observations"], 8)
        if not result["optimizer_attempts"].empty:
            self.assertTrue((result["optimizer_attempts"]["final_objective"] < config.large_objective_penalty).any())

    def test_simulate_intracellular_dfdctp_signal_runs_for_synthetic_input(self):
        t = np.linspace(0.0, 3.0, 7)
        exposure_curve = lambda t_days, dose_muM: legacy_gemcitabine_pk.converted_extracellular_gem_signal(dose_muM) * np.exp(-0.3 * np.asarray(t_days))
        signal = legacy_gemcitabine_pk.simulate_intracellular_dfdctp_signal(
            t,
            exposure_curve=exposure_curve,
            dose_muM=1.0,
            eta_per_day=12.0,
            k_decay_per_day=0.8,
            initial_signal=0.05,
        )
        self.assertEqual(signal.shape, t.shape)
        self.assertTrue(np.issubdtype(signal.dtype, np.floating))
        self.assertTrue(np.all(np.isfinite(signal)))
        self.assertGreaterEqual(np.min(signal), 0.0)

    def test_fit_intracellular_pk_parameters_recovers_synthetic_eta_and_k_decay(self):
        exposure_curve = legacy_gemcitabine_pk.build_extracellular_exposure_curve_from_profile(
            time_days=np.array([0.0, 1.0, 2.0, 3.0]),
            gem_concentration_ng_per_ml=np.array([100.0, 60.0, 36.0, 21.6]),
            gem_was_censored=np.array([False, False, False, False]),
            gem_censor_upper_bound_ng_per_ml=np.array([np.nan, np.nan, np.nan, np.nan]),
            reference_dose_muM=1.0,
            fallback_half_life_days=1.0,
            analyte_column="Gemcitabine (ng/mL)",
            source_ploidy="synthetic",
            preferred_mode="constrained_bolus_exponential",
        )
        t_days = np.array([0.0, 0.25, 0.5, 1.0, 2.0, 3.0])
        eta_true = 8.0
        k_decay_true = 1.2
        initial_signal = 0.03
        observed_signal = legacy_gemcitabine_pk.simulate_intracellular_dfdctp_signal(
            t_days,
            exposure_curve=exposure_curve,
            dose_muM=1.0,
            eta_per_day=eta_true,
            k_decay_per_day=k_decay_true,
            initial_signal=initial_signal,
        )
        df = pd.DataFrame(
            {
                "Timepoint": t_days * 24.0,
                "dFdCTP (ng/mL)": observed_signal * invitro_fitting.DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
            }
        )
        fit = legacy_gemcitabine_pk.fit_intracellular_pk_parameters(
            df,
            exposure_curve=exposure_curve,
            reference_dose_muM=1.0,
        )
        self.assertIsNotNone(fit)
        assert fit is not None
        self.assertTrue(math.isclose(fit["eta"], eta_true, rel_tol=1e-3, abs_tol=1e-3))
        self.assertTrue(math.isclose(fit["k_decay"], k_decay_true, rel_tol=1e-3, abs_tol=1e-3))
        self.assertGreater(fit["r2_confidence"], 0.999)

    def test_single_clone_signal_ode_joint_returns_finite_derivatives(self):
        y = [1000.0, 0.1, 0.0]
        dfdctp_signal_curve = self.make_constant_surface(0.05)
        deriv = invitro_fitting.single_clone_signal_ode_joint(
            y=y,
            t=1.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.05,
            dfdctp_signal_curve=dfdctp_signal_curve,
            k_tr=0.5,
            k_kill=2.0,
            k_clear=0.2,
            n_tr=1,
            model_variant="delayed_death_only",
        )
        self.assertTrue(np.all(np.isfinite(deriv)))
        kappa = 2.0 * 0.1
        self.assertTrue(np.isfinite(kappa))

    def test_single_clone_dfdctp_ode_joint_uses_expected_state_dimension(self):
        y = [1000.0, 0.1, 0.2, 0.3, 5.0]
        curve = self.make_constant_surface(2.0)
        deriv = invitro_fitting.single_clone_dfdctp_ode_joint(
            y=y,
            t=1.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.05,
            dfdctp_signal_curve=curve,
            k_tr=0.5,
            k_kill=2.0,
            k_clear=0.2,
            n_tr=3,
            model_variant="delayed_death_only",
        )
        self.assertEqual(len(deriv), 5)
        self.assertTrue(np.all(np.isfinite(deriv)))

    def test_fit_baseline_shared_rk_replicate_n0_recovers_shared_parameters(self):
        t = np.array([0.0, 1.0, 2.0, 3.0])
        n0_by_rep = np.array([800.0, 1000.0, 1200.0])
        true_r = 0.7
        true_K = 6000.0
        y = invitro_fitting.logistic_growth_replicate_n0(t, n0_by_rep, true_r, true_K)

        result = invitro_fitting.fit_baseline_shared_rk_replicate_n0(
            t_data=t,
            y_alive=y,
            y_dead=np.zeros_like(y),
            ploidy_label="synthetic",
            t_max=None,
            output_dir=None,
        )

        self.assertTrue(math.isclose(result["r"], true_r, rel_tol=1e-4, abs_tol=1e-4))
        self.assertTrue(math.isclose(result["K"], true_K, rel_tol=1e-4, abs_tol=1e-4))
        self.assertTrue(np.allclose(result["N0_by_replicate"], n0_by_rep))

    def test_fit_baseline_shared_rk_replicate_n0_handles_missing_values(self):
        t = np.array([0.0, 1.0, 2.0, 3.0])
        n0_by_rep = np.array([900.0, 1100.0, 1300.0])
        y = invitro_fitting.logistic_growth_replicate_n0(t, n0_by_rep, 0.6, 5500.0)
        y[2, 1] = np.nan

        result = invitro_fitting.fit_baseline_shared_rk_replicate_n0(
            t_data=t,
            y_alive=y,
            y_dead=np.zeros_like(y),
            ploidy_label="synthetic_missing",
            t_max=None,
            output_dir=None,
        )

        self.assertTrue(np.isfinite(result["r"]))
        self.assertTrue(np.isfinite(result["K"]))
        self.assertEqual(result["n_replicates"], 3)

    def test_replicate_specific_n0_fit_beats_legacy_mean_n0_when_initial_counts_differ(self):
        t = np.array([0.0, 1.0, 2.0, 3.0])
        n0_by_rep = np.array([700.0, 1000.0, 1400.0])
        true_r = 0.55
        true_K = 5000.0
        y = invitro_fitting.logistic_growth_replicate_n0(t, n0_by_rep, true_r, true_K)

        shared_result = invitro_fitting.fit_baseline_shared_rk_replicate_n0(
            t_data=t,
            y_alive=y,
            y_dead=np.zeros_like(y),
            ploidy_label="synthetic_compare",
            t_max=None,
            output_dir=None,
        )
        shared_residual = invitro_fitting.residuals_shared_rk_replicate_n0(
            [shared_result["r"], shared_result["K"]],
            t,
            y,
            n0_by_rep,
            None,
        )
        shared_sse = float(np.sum(shared_residual ** 2))

        n0_mean = float(np.mean(n0_by_rep))
        legacy = least_squares(
            invitro_fitting.residuals_fixed_N0,
            [0.5, np.nanmax(y)],
            args=(np.repeat(t, y.shape[1]), y.flatten(), n0_mean),
            bounds=([0, 0], [np.inf, np.inf]),
        )
        legacy_sse = float(np.sum(invitro_fitting.residuals_fixed_N0(legacy.x, np.repeat(t, y.shape[1]), y.flatten(), n0_mean) ** 2))

        self.assertLess(shared_sse, legacy_sse)

    def test_replicate_specific_n0_agrees_with_legacy_when_initial_counts_match(self):
        t = np.array([0.0, 1.0, 2.0, 3.0])
        n0_by_rep = np.array([1000.0, 1000.0, 1000.0])
        true_r = 0.65
        true_K = 5200.0
        y = invitro_fitting.logistic_growth_replicate_n0(t, n0_by_rep, true_r, true_K)

        shared_result = invitro_fitting.fit_baseline_shared_rk_replicate_n0(
            t_data=t,
            y_alive=y,
            y_dead=np.zeros_like(y),
            ploidy_label="synthetic_equal_n0",
            t_max=None,
            output_dir=None,
        )

        n0_mean = float(np.mean(n0_by_rep))
        legacy = least_squares(
            invitro_fitting.residuals_fixed_N0,
            [0.5, np.nanmax(y)],
            args=(np.repeat(t, y.shape[1]), y.flatten(), n0_mean),
            bounds=([0, 0], [np.inf, np.inf]),
        )

        self.assertTrue(math.isclose(shared_result["r"], legacy.x[0], rel_tol=1e-6, abs_tol=1e-6))
        self.assertTrue(math.isclose(shared_result["K"], legacy.x[1], rel_tol=1e-6, abs_tol=1e-6))

    def test_pk_workbook_smoke_builds_2n_and_4n_curves(self):
        if not invitro_fitting.PATHS["pkpd_constants"].exists():
            self.skipTest("PK workbook not available")
        sheets = invitro_fitting.import_and_clean_pkpd(censored_strategy="nan")
        for ploidy in ["2N", "4N"]:
            reference_dose = invitro_fitting.get_pk_reference_dose_uM(ploidy)
            curve = legacy_gemcitabine_pk.build_extracellular_exposure_curve_from_pk_sheet(
                sheets[ploidy],
                ploidy_label=ploidy,
                reference_dose_muM=reference_dose,
                fallback_half_life_days=1.0,
                preferred_mode="constrained_bolus_exponential",
            )
            vals = np.array([curve(0.0, reference_dose), curve(1.0, reference_dose), curve(5.0, reference_dose)], dtype=float)
            self.assertTrue(np.all(np.isfinite(vals)))
            self.assertTrue(math.isclose(curve(0.0, reference_dose), reference_dose, rel_tol=1e-12, abs_tol=1e-12))


if __name__ == "__main__":
    unittest.main()
