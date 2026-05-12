import importlib.util
import math
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares


MODULE_PATH = Path(__file__).resolve().parents[2] / "code" / "invitro_fitting.py"
SPEC = importlib.util.spec_from_file_location("invitro_fitting", MODULE_PATH)
invitro_fitting = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(invitro_fitting)


class ExtracellularDecayTests(unittest.TestCase):
    def test_get_pk_reference_dose_uM_uses_documented_sheet_mapping(self):
        self.assertEqual(invitro_fitting.get_pk_reference_dose_uM("2N"), 1000.0)
        self.assertEqual(invitro_fitting.get_pk_reference_dose_uM("4N"), 1000.0)
        self.assertEqual(invitro_fitting.get_pk_reference_dose_uM("2N_lowInitialGemcitabine"), 100.0)
        self.assertEqual(invitro_fitting.get_pk_reference_dose_uM("4N_lowInitialGemcitabine"), 100.0)
        with self.assertRaises(KeyError):
            invitro_fitting.get_pk_reference_dose_uM("unknown")

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
                invitro_fitting.converted_extracellular_gem_signal(dose_muM),
                0.00971 * dose_muM,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )

    def test_gemcitabine_ng_per_ml_to_uM_uses_physical_conversion(self):
        concentration_ng_per_ml = 26.32
        self.assertTrue(
            math.isclose(
                invitro_fitting.gemcitabine_ng_per_ml_to_uM(concentration_ng_per_ml),
                concentration_ng_per_ml / invitro_fitting.GEMCITABINE_MOLECULAR_WEIGHT_NG_PER_NMOL,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )

    def test_estimate_eta_per_day_uses_physical_uM_formula(self):
        delta_dfdctp_uM_per_hour = 0.25
        dose_muM = 1000.0
        expected_eta = (delta_dfdctp_uM_per_hour / dose_muM) * 24.0
        refactored_eta = invitro_fitting.estimate_eta_per_day(
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
        curve = invitro_fitting.build_extracellular_exposure_curve_from_profile(
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
        curve = invitro_fitting.build_extracellular_exposure_curve_from_profile(
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
            invitro_fitting.converted_extracellular_gem_signal(0.0)
        with self.assertRaises(ValueError):
            invitro_fitting.converted_extracellular_gem_signal(1.0, conversion_factor=0.0)

    def test_estimate_eta_per_day_validates_inputs(self):
        with self.assertRaises(ValueError):
            invitro_fitting.estimate_eta_per_day(1.0, 0.0)

    def test_extracellular_concentration_starts_at_initial_value(self):
        dose_muM = 12.5
        k_ext = invitro_fitting.decay_rate_from_half_life_days(1.0)
        self.assertTrue(
            math.isclose(
                invitro_fitting.extracellular_gemcitabine_concentration(0.0, dose_muM, k_ext),
                0.00971 * dose_muM,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )

    def test_extracellular_concentration_halves_at_half_life(self):
        dose_muM = 25.0
        half_life_days = 1.5
        k_ext = invitro_fitting.decay_rate_from_half_life_days(half_life_days)
        c0 = invitro_fitting.extracellular_gemcitabine_concentration(0.0, dose_muM, k_ext)
        c_half = invitro_fitting.extracellular_gemcitabine_concentration(half_life_days, dose_muM, k_ext)
        self.assertTrue(math.isclose(c_half, 0.5 * c0, rel_tol=1e-12, abs_tol=1e-12))

    def test_long_half_life_recovers_constant_exposure_limit_over_five_days(self):
        curve = invitro_fitting.build_extracellular_exposure_curve_from_profile(
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
        old_constant = np.full_like(conc, invitro_fitting.converted_extracellular_gem_signal(dose_muM), dtype=float)
        max_relative_error = np.max(np.abs(conc - old_constant) / old_constant)
        self.assertLess(max_relative_error, 1e-8)

    def test_constrained_decay_penalizes_numeric_censor_upper_bounds(self):
        fit = invitro_fitting.fit_constrained_extracellular_gem_decay(
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

    def test_simulate_joint_ext_runs_for_synthetic_input(self):
        t = np.linspace(0.0, 5.0, 6)
        params_3 = [0.5, 1.0, 0.2]
        exposure_curve = lambda t_days, dose_muM: invitro_fitting.converted_extracellular_gem_signal(dose_muM) * np.exp(-0.4 * np.asarray(t_days))
        alive, dead = invitro_fitting.simulate_joint_ext(
            t=t,
            params_3=params_3,
            N0=1000.0,
            D0=0.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.05,
            eta_fixed=10.0,
            k_decay_fixed=0.7,
            exposure_curve=exposure_curve,
            n_tr=2,
        )
        self.assertEqual(alive.shape, t.shape)
        self.assertEqual(dead.shape, t.shape)
        self.assertTrue(np.issubdtype(alive.dtype, np.floating))
        self.assertTrue(np.issubdtype(dead.dtype, np.floating))
        self.assertTrue(np.all(np.isfinite(alive)))
        self.assertTrue(np.all(np.isfinite(dead)))

    def test_simulate_intracellular_dfdctp_signal_runs_for_synthetic_input(self):
        t = np.linspace(0.0, 3.0, 7)
        exposure_curve = lambda t_days, dose_muM: invitro_fitting.converted_extracellular_gem_signal(dose_muM) * np.exp(-0.3 * np.asarray(t_days))
        signal = invitro_fitting.simulate_intracellular_dfdctp_signal(
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
        exposure_curve = invitro_fitting.build_extracellular_exposure_curve_from_profile(
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
        observed_signal = invitro_fitting.simulate_intracellular_dfdctp_signal(
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
        fit = invitro_fitting.fit_intracellular_pk_parameters(
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
        y = [1000.0, 0.25, 0.1, 0.0]
        exposure_curve = lambda t_days, dose_muM: invitro_fitting.converted_extracellular_gem_signal(dose_muM) * np.exp(-0.4 * np.asarray(t_days))
        deriv = invitro_fitting.single_clone_signal_ode_joint(
            y=y,
            t=1.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.05,
            eta=10.0,
            k_decay=0.7,
            exposure_curve=exposure_curve,
            k_tr=0.5,
            k_kill=2.0,
            k_clear=0.2,
            n_tr=1,
        )
        self.assertTrue(np.all(np.isfinite(deriv)))
        kappa = 2.0 * 0.1
        self.assertTrue(np.isfinite(kappa))

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
            curve = invitro_fitting.build_extracellular_exposure_curve_from_pk_sheet(
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
