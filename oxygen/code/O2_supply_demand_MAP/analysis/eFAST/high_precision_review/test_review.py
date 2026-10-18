"""Protect decimal precision verification and scientific point classification."""
from decimal import Decimal
import json
from pathlib import Path
import tempfile
import unittest
import review


class ReviewChecks(unittest.TestCase):
    def fixture(self,root):
        ident=dict(point_id='1',fit_seed='seed101',N='4097',M='8',replicate='1',parameter='p_mis_base',sample_id='1',O2_pct='4.025',roles='abnormal',original_ploidy='7',original_growth='0')
        report=dict(status='evaluated',identity=ident,double=dict(normalized_ploidy=7,lambda_real=0,negative_l1_mass=0,
            clipped_ploidy=7,signed_ploidy=7,raw_relative_residual=1e-15,normalized_relative_residual=1e-15),matrix=dict(n=2,strongly_connected=True),files_sha256={})
        for digits in (50,100):
            side=dict(converged=True,vector_decimal=['0.5','0.5'],lambda_='0',ploidy='1.5',relative_residual_inf='1e-90',
                      relative_tolerance='1e-32' if digits==50 else '1e-65',relative_upper_minus_lambda='1e-90',vector_change_l1='1e-90',vector_tolerance='1e-30',
                      max_relative_component_change='1e-90',component_tolerance='1e-20')
            side['lambda']=side.pop('lambda_')
            p=root/('precision_%d.json'%digits);review.dump(p,dict(right=side,left=side));report['files_sha256'][p.name]=review.sha(p)
        (root/'matrix.rds').write_bytes(b'unit fixture');report['files_sha256']['matrix.rds']=review.sha(root/'matrix.rds')
        review.dump(root/'result.json',report)

    def test_corrected_response_is_not_mistaken_for_failed_precision(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p);r=review.audit_point(p)
            self.assertTrue(r['precision_verified']);self.assertEqual(r['classification'],'double_response_changed')
            self.assertEqual(r['ploidy_change'],-5.5)

    def test_tiny_components_remain_nonzero_in_decimal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p)
            for digits in (50,100):
                z=json.loads((p/('precision_%d.json'%digits)).read_text())
                for side in ('right','left'):z[side]['vector_decimal']=['1','1e-600']
                f=p/('precision_%d.json'%digits);review.dump(f,z)
                report=json.loads((p/'result.json').read_text());report['files_sha256'][f.name]=review.sha(f);review.dump(p/'result.json',report)
            r=review.audit_point(p);self.assertTrue(r['precision_verified'])
            self.assertGreater(Decimal(r['left_right_overlap']),0)

    def test_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p);(p/'matrix.rds').write_bytes(b'changed')
            with self.assertRaises(AssertionError):review.audit_point(p)


if __name__=='__main__':unittest.main()
