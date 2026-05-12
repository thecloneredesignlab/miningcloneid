import importlib.util
import math
import unittest
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[2] / "code" / "invitro_fitting.py"
SPEC = importlib.util.spec_from_file_location("invitro_fitting", MODULE_PATH)
invitro_fitting = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(invitro_fitting)


class ExtracellularDecayTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
