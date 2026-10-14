"""Checks that protect the interpretation against false convergence claims."""
import io
import unittest
import zipfile
import numpy as np
import collect_resolution_results as collect
import report_resolution_pilot as report


class ReportChecks(unittest.TestCase):
    def test_local_spike_does_not_become_a_full_curve_pass(self):
        oxygen=np.arange(201)*.025
        difference=np.full(201,.01);difference[77]=.4
        q,maximum,at,count,passed,status=report.diagnostic(difference,oxygen,0,5,.05)
        self.assertTrue(passed)
        self.assertEqual(status,'localized_review')
        self.assertAlmostEqual(q,.01)
        self.assertAlmostEqual(maximum,.4)
        self.assertAlmostEqual(at,1.925)
        self.assertEqual(count,201)
        difference[1]=np.nan
        self.assertEqual(report.diagnostic(difference,oxygen,0,5,.05)[-1],'undefined_or_incomplete')
        self.assertFalse(report.diagnostic(difference,oxygen,0,5,.05)[-2])

    def test_incomplete_phase_spread_is_not_replaced_with_zero(self):
        values=np.ones((5,3));values[0,1]=np.nan;values[:,2]=np.nan
        mean,sd,span,count=report.stats(values)
        np.testing.assert_array_equal(count,[5,4,0])
        self.assertEqual(span[0],0.)
        self.assertTrue(np.isnan(span[1:]).all())
        self.assertTrue(np.isnan(sd[1:]).all())
        self.assertTrue(np.isnan(mean[2]))

    def test_stdlib_gap_reader_matches_numpy_and_rejects_object_payload(self):
        stream=io.BytesIO()
        a=np.arange(12,dtype=np.int64).reshape(3,4)
        b=np.array([1e-10,2e-5,.01],dtype=np.float64)
        np.savez_compressed(stream,gap_counts=a,gap_min=b,unsafe=np.array([object()],dtype=object))
        stream.seek(0)
        with zipfile.ZipFile(stream) as archive:
            shape,data=collect.numeric_npy(archive,'gap_counts')
            self.assertEqual(shape,a.shape)
            np.testing.assert_array_equal(np.array(data).reshape(shape),a)
            shape,data=collect.numeric_npy(archive,'gap_min')
            np.testing.assert_array_equal(data,b)
            with self.assertRaises(AssertionError):collect.numeric_npy(archive,'unsafe')

    def test_validator_rejects_fabricated_zero_for_undefined_or_wrong_values(self):
        verifier=report.Verifier()
        verifier.check('nan',np.nan,'undefined preserved')
        verifier.check('0.2',.2,'same numeric')
        with self.assertRaises(AssertionError):verifier.check('0',np.nan,'undefined became zero')
        with self.assertRaises(AssertionError):verifier.check('0.1',.2,'wrong table value')


if __name__=='__main__':unittest.main()
