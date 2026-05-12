import importlib.util
import math
import unittest
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


MODULE_PATH = Path(__file__).resolve().parents[2] / "code" / "invitro_fitting.py"
SPEC = importlib.util.spec_from_file_location("invitro_fitting", MODULE_PATH)
invitro_fitting = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(invitro_fitting)


class ExtracellularDecayTests(unittest.TestCase):
    def test_parse_pk_concentration_value_handles_censored_and_numeric_inputs(self):
        value, was_censored = invitro_fitting.parse_pk_concentration_value("BDL", censored_strategy="nan")
        self.assertTrue(np.isnan(value))
        self.assertTrue(was_censored)

        value, was_censored = invitro_fitting.parse_pk_concentration_value("BQL", censored_strategy="nan")
        self.assertTrue(np.isnan(value))
        self.assertTrue(was_censored)

        value, was_censored = invitro_fitting.parse_pk_concentration_value("<0.5", censored_strategy="half_lod")
        self.assertTrue(math.isclose(value, 0.25, rel_tol=1e-12, abs_tol=1e-12))
        self.assertTrue(was_censored)

        value, was_censored = invitro_fitting.parse_pk_concentration_value("<0.5", censored_strategy="nan")
        self.assertTrue(np.isnan(value))
        self.assertTrue(was_censored)

        value, was_censored = invitro_fitting.parse_pk_concentration_value("1.23", censored_strategy="nan")
        self.assertTrue(math.isclose(value, 1.23, rel_tol=1e-12, abs_tol=1e-12))
        self.assertFalse(was_censored)

        value, was_censored = invitro_fitting.parse_pk_concentration_value(np.nan, censored_strategy="nan")
        self.assertTrue(np.isnan(value))
        self.assertFalse(was_censored)

    def test_parse_pk_concentration_value_zero_strategy_is_explicit(self):
        value, was_censored = invitro_fitting.parse_pk_concentration_value("BDL", censored_strategy="zero")
        self.assertEqual(value, 0.0)
        self.assertTrue(was_censored)

        value, was_censored = invitro_fitting.parse_pk_concentration_value("<0.5", censored_strategy="zero")
        self.assertEqual(value, 0.0)
        self.assertTrue(was_censored)

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

    def test_estimate_eta_per_day_matches_legacy_formula(self):
        delta_c_hourly = 12.34
        dose_muM = 1.0
        legacy_eta = (
            (delta_c_hourly / (0.00971 * dose_muM))
            * 24.0
            / invitro_fitting.DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL
        )
        refactored_eta = invitro_fitting.estimate_eta_per_day(
            delta_dfdctp_signal_per_hour=delta_c_hourly,
            dose_muM=dose_muM,
        )
        self.assertTrue(
            math.isclose(
                refactored_eta,
                legacy_eta,
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
        with self.assertRaises(ValueError):
            invitro_fitting.estimate_eta_per_day(1.0, 1.0, conversion_factor=0.0)
        with self.assertRaises(ValueError):
            invitro_fitting.estimate_eta_per_day(
                1.0,
                1.0,
                molecular_weight_ng_per_nmol=0.0,
            )

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
        dose_muM = 50.0
        k_ext = invitro_fitting.decay_rate_from_half_life_days(1e9)
        t = np.linspace(0.0, 5.0, 11)
        conc = invitro_fitting.extracellular_gemcitabine_concentration(t, dose_muM, k_ext)
        old_constant = np.full_like(conc, 0.00971 * dose_muM, dtype=float)
        max_relative_error = np.max(np.abs(conc - old_constant) / old_constant)
        self.assertLess(max_relative_error, 1e-8)

    def test_simulate_joint_ext_runs_for_synthetic_input(self):
        t = np.linspace(0.0, 5.0, 6)
        params_3 = [0.5, 1.0, 0.2]
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
            k_ext_decay_per_day=invitro_fitting.decay_rate_from_half_life_days(1.0),
            n_tr=2,
        )
        self.assertEqual(alive.shape, t.shape)
        self.assertEqual(dead.shape, t.shape)
        self.assertTrue(np.issubdtype(alive.dtype, np.floating))
        self.assertTrue(np.issubdtype(dead.dtype, np.floating))
        self.assertTrue(np.all(np.isfinite(alive)))
        self.assertTrue(np.all(np.isfinite(dead)))

    def test_single_clone_signal_ode_joint_returns_finite_derivatives(self):
        y = [1000.0, 0.25, 0.1, 0.0]
        deriv = invitro_fitting.single_clone_signal_ode_joint(
            y=y,
            t=1.0,
            r=0.8,
            K=5000.0,
            dose_muM=0.05,
            eta=10.0,
            k_decay=0.7,
            k_ext_decay_per_day=invitro_fitting.decay_rate_from_half_life_days(1.0),
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


if __name__ == "__main__":
    unittest.main()
