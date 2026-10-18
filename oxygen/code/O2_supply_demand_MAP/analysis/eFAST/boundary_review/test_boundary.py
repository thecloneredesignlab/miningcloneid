"""Protect the scientific scope and comparison of unequal chromosome domains."""
import unittest
from boundary import selected, compare

THRESHOLDS = dict(ploidy_abs=.01, growth_abs=1e-6, growth_rel=1e-3,
                  distribution_TV=1e-3, upper_22_mass=1e-4)


class BoundaryChecks(unittest.TestCase):
    def test_only_approved_seeds_and_roles(self):
        rows = [dict(fit_seed=seed, roles=role) for seed in ('seed138', 'seed101', 'seed174')
                for role in ('abnormal', 'ordinary_control', 'oxygen_neighbor')]
        actual = selected(rows)
        self.assertEqual(len(actual), 4)
        self.assertFalse(any(r['fit_seed'] == 'seed174' or r['roles'] == 'oxygen_neighbor' for r in actual))

    def test_escaped_mass_is_not_renormalized_away(self):
        a = dict(vector_decimal=['1', '0'], ploidy='1', **{'lambda': '-1'})
        b = dict(vector_decimal=['0.5', '0', '0.5'], ploidy='2', **{'lambda': '-0.5'})
        r = compare(a, b, THRESHOLDS)
        self.assertEqual(float(r['distribution_TV']), .5)
        self.assertFalse(r['diagnostic_stable'])

    def test_boundary_mass_prevents_false_stability(self):
        a = dict(vector_decimal=['0'] * 29 + ['1'], ploidy='2', **{'lambda': '-1'})
        b = dict(vector_decimal=['0'] * 29 + ['1', '0'], ploidy='2', **{'lambda': '-1'})
        r = compare(a, b, THRESHOLDS)
        self.assertEqual(float(r['distribution_TV']), 0)
        self.assertFalse(r['diagnostic_stable'])


if __name__ == '__main__': unittest.main()
