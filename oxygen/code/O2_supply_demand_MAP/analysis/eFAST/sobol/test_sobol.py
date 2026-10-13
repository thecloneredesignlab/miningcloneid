#!/usr/bin/env python3
"""Numerical and sample-order checks for second-order Sobol and response surfaces."""
import tempfile
from pathlib import Path
import unittest

import numpy as np
import core as c
import report
import refine


def unit_ranges():
    return [dict(parameter=p,encoded_lower=0.,encoded_upper=1.,natural_lower=0.,natural_upper=1.,transform='identity') for p in c.ACTIVE]


class SobolChecks(unittest.TestCase):
    def test_full_natural_span_neighborhood_and_log_distribution(self):
        rr=[dict(parameter='x',natural_lower=1.,natural_upper=101.,encoded_lower=0.,encoded_upper=np.log10(101),transform='log10')]
        local=c.bounds_for(rr,dict(x=2.))
        self.assertEqual((local[0]['natural_lower'],local[0]['natural_upper']),(1.,12.))
        self.assertAlmostEqual(local[0]['encoded_upper'],np.log10(12.))
        self.assertEqual(c.bounds_for(rr),rr)

    def test_five_independent_scrambles_and_nested_prefixes(self):
        rr=unit_ranges(); rr[0].update(encoded_lower=-6.,encoded_upper=-1.,natural_lower=1e-6,natural_upper=.1,transform='log10')
        xx=[c.samples(rr,1024,100+i) for i in range(5)]
        self.assertEqual(len({x.tobytes() for x in xx}),5)
        for n in (256,512): np.testing.assert_array_equal(c.prefix(xx[0],n),c.samples(rr,n,100))
        self.assertEqual(xx[0].shape,(30*1024,14))
        self.assertTrue(np.all((xx[0][:,0]>=1e-6)&(xx[0][:,0]<=.1)))
        larger=c.samples(rr,4096,100)
        for n in (1024,2048): np.testing.assert_array_equal(c.prefix(larger,n),c.samples(rr,n,100))

    def test_split_and_reassemble_matches_interleaved_salib_design(self):
        x=c.samples(unit_ranges(),128,19)
        payload=np.stack([x[:,0],x[:,1]],axis=1)[:,None,:]
        rebuilt=c.assemble([payload[b::c.STEP] for b in range(c.STEP)])
        np.testing.assert_array_equal(rebuilt,payload)
        with self.assertRaises(ValueError): c.assemble([payload[0::c.STEP]])

    def test_additive_and_known_product_second_order(self):
        x=c.samples(unit_ranges(),4096,29)
        si,_=c.analyze(unit_ranges(),x[:,0]+2*x[:,1])
        np.testing.assert_allclose(si['S1'][:2],[.2,.8],atol=.01)
        self.assertLess(max(abs(si['S2'][i,j]) for i,j in c.PAIRS),.02)
        product,_=c.analyze(unit_ranges(),x[:,0]*x[:,1])
        self.assertAlmostEqual(product['S2'][0,1],1/7,delta=.01)
        np.testing.assert_allclose(product['S1'][:2],[3/7,3/7],atol=.01)

    def test_constant_outputs_and_missing_or_duplicate_model_rows(self):
        si,var=c.analyze(unit_ranges(),np.ones(30*256))
        self.assertEqual(var,0); self.assertTrue(np.isnan(si['S2']).all())
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'out.tsv.gz'
            rows=[dict(sample_id=i,O2_pct=o,dominant_mean_ploidy=2+i,dominant_growth_rate=.1,
                spectral_gap=.01,eigen_residual=1e-14,eigenvector_nonnegative='TRUE',status='ok',solver='double') for i in (1,2) for o in (0,5)]
            c.table(path,rows); values,qc=c.output_array(path,2,[0,5])
            self.assertEqual(values.shape,(2,2,2)); self.assertEqual(qc['n_points'],4)
            c.table(path,rows[:-1])
            with self.assertRaisesRegex(ValueError,'Incomplete'): c.output_array(path,2,[0,5])
            c.table(path,rows+[rows[0]])
            with self.assertRaisesRegex(ValueError,'Duplicate'): c.output_array(path,2,[0,5])

    def test_negative_estimates_and_failed_checks_are_retained(self):
        rr=[]
        for n in (256,512,1024):
            for repeat in range(1,6):
                rr.append(dict(scope='global',case_id='global',fit_seed='seed25',N=n,O2_pct=0,output=c.OUTPUTS[0],
                    index='S2',parameter_i=c.ACTIVE[0],parameter_j=c.ACTIVE[1],repeat=repeat,
                    value=-.02,bootstrap_half_width=.1,bootstrap_lower=-.12,output_variance=1.))
        rows=report.summarize_rows(rr,dict(repeat_range_tolerance=.05,resolution_delta_tolerance=.02))
        self.assertTrue(all(r['value_mean']<0 and r['retained']=='TRUE' for r in rows))
        rr[-1]['value']=.8
        rows=report.summarize_rows(rr,dict(repeat_range_tolerance=.05,resolution_delta_tolerance=.02))
        self.assertEqual(rows[-1]['convergence_status'],'not_passed')
        self.assertEqual(rows[-1]['n_valid_repeats'],5)

    def test_exact_anchor_and_nonadditive_surface(self):
        gi=report.anchored_grid(.01,10,.173,True)
        gj=report.anchored_grid(0,1,.381,False)
        ai=int(np.flatnonzero(gi==.173)[0]); aj=int(np.flatnonzero(gj==.381)[0])
        self.assertEqual(len(gi),41); self.assertEqual((gi[0],gi[-1]),(.01,10.))
        near_boundary=report.anchored_grid(.005,.2999999999999997,.3,True)
        self.assertEqual(near_boundary[-1],.3)
        with self.assertRaises(ValueError): report.anchored_grid(.005,.3,.31,True)
        additive=gi[:,None]**2+np.cos(gj[None,:])
        np.testing.assert_allclose(report.interaction_contrast(additive,ai,aj),0,atol=2e-14)
        product=gi[:,None]*gj[None,:]
        np.testing.assert_allclose(report.interaction_contrast(product,ai,aj),
            (gi-.173)[:,None]*(gj-.381)[None,:],atol=2e-14)

    def test_refined_resolution_diagnostics_use_actual_ladder(self):
        rr=[]
        for n,value in ((1024,.9),(2048,.1),(4096,.11)):
            for repeat in range(1,6):
                rr.append(dict(scope='global',case_id='global',fit_seed='seed25',N=n,O2_pct=0,output=c.OUTPUTS[0],
                    index='S2',parameter_i=c.ACTIVE[0],parameter_j=c.ACTIVE[1],repeat=repeat,
                    value=value,bootstrap_half_width=.01,bootstrap_lower=value-.01,output_variance=1.))
        plan=dict(resolutions=[1024,2048,4096],repeat_range_tolerance=.05,resolution_delta_tolerance=.02)
        rows=report.summarize_rows(rr,plan)
        self.assertEqual(rows[-1]['resolution_low'],2048)
        self.assertEqual(rows[-1]['resolution_high'],4096)
        self.assertAlmostEqual(rows[-1]['delta_N2048_N4096'],.01)
        self.assertEqual(rows[-1]['convergence_status'],'passed')
        self.assertNotIn('delta_N512_N1024',rows[-1])

    def test_refinement_prefix_fusion_restores_order_and_rejects_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); paths=[root/name for name in ('old.tsv.gz','tail.tsv.gz','combined.tsv.gz')]
            def row(i,o,value):
                return dict(sample_id=i,O2_pct=o,dominant_mean_ploidy=value,dominant_growth_rate=.1,
                    spectral_gap=.01,eigen_residual=1e-14,eigenvector_nonnegative='TRUE',status='ok',solver='double')
            old=[row(i,o,i) for i in (1,2) for o in (0,5)]
            tail=[row(i,o,i+2) for i in (1,2,3,4) for o in (0,5)]
            c.table(paths[0],old); c.table(paths[1],tail)
            refine.fuse_outputs(*paths,2,4,[0,5])
            values,qc=c.output_array(paths[2],6,[0,5])
            np.testing.assert_array_equal(values[:,0,0],np.arange(1,7))
            self.assertEqual(qc['n_points'],12)
            c.table(paths[1],tail[:-1])
            with self.assertRaisesRegex(ValueError,'Incomplete'):
                refine.fuse_outputs(*paths,2,4,[0,5])
            unchanged,_=c.output_array(paths[2],6,[0,5])
            np.testing.assert_array_equal(values,unchanged)


if __name__=='__main__': unittest.main()
