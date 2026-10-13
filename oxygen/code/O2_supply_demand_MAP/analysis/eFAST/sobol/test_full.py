#!/usr/bin/env python3
"""Mathematical ordering, lossless storage and full-run diagnostic checks."""
import tempfile
from pathlib import Path
import unittest
import numpy as np
import core as c
import full_core as f
import full_report
from unittest.mock import patch


class FullChecks(unittest.TestCase):
    def test_full_task_bijection_and_boundaries(self):
        self.assertEqual(f.N_TASKS,75150)
        for task in range(1,f.N_TASKS+1):
            ci,repeat,block=f.decode(task)
            self.assertEqual((ci-1)*150+(repeat-1)*30+block+1,task)
        self.assertEqual(f.decode(1),(1,1,0)); self.assertEqual(f.decode(150),(1,5,29))
        self.assertEqual(f.decode(151),(2,1,0)); self.assertEqual(f.decode(75150),(501,5,29))
        with self.assertRaises(ValueError): f.decode(75151)

    def test_transport_matches_direct_salib_and_nested_prefix(self):
        unit=[dict(parameter=p,encoded_lower=0.,encoded_upper=1.,natural_lower=0.,natural_upper=1.,transform='identity') for p in c.ACTIVE]
        bounds=[dict(r,encoded_lower=-2.,encoded_upper=3.,natural_lower=-2.,natural_upper=3.) for r in unit]
        bounds[0].update(encoded_lower=-6.,encoded_upper=-1.,natural_lower=1e-6,natural_upper=.1,transform='log10')
        u=c.samples(unit,256,99); x=f.scale_unit(u,bounds)
        np.testing.assert_array_equal(x,c.samples(bounds,256,99))
        np.testing.assert_array_equal(x[:128*30],c.samples(bounds,128,99))
        np.testing.assert_array_equal(f.scale_unit(u[17::30],bounds),x[17::30])

    def test_compact_codec_preserves_values_and_rejects_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); rows=[dict(sample_id=i,O2_pct=o,dominant_mean_ploidy=1.2345678901234567+i,
                dominant_growth_rate=-.01234567890123456,spectral_gap=1e-11,eigen_residual=4e-14,
                eigenvector_nonnegative='TRUE',solver='double' if i==1 else 'decimal50_100',status='ok') for i in [1,2] for o in [0.,.025,5.]]
            c.table(p/'out.tsv.gz',rows)
            values,flags,gaps,residual=f.read_model(p/'out.tsv.gz',2,[0.,.025,5.])
            baseline,_=c.output_array(p/'out.tsv.gz',2,[0.,.025,5.])
            np.testing.assert_array_equal(values,baseline)
            f.atomic_npz(p/'binary.npz',outputs=values,solver_flags=flags)
            with np.load(p/'binary.npz') as z:
                self.assertEqual(z['outputs'].dtype,np.float64); np.testing.assert_array_equal(z['outputs'],values)
                self.assertEqual(int((z['solver_flags']==2).sum()),3)
            c.table(p/'out.tsv.gz',rows[:-1])
            with self.assertRaisesRegex(ValueError,'Incomplete'): f.read_model(p/'out.tsv.gz',2,[0.,.025,5.])
            c.table(p/'out.tsv.gz',rows+[rows[0]])
            with self.assertRaisesRegex(ValueError,'Duplicate'): f.read_model(p/'out.tsv.gz',2,[0.,.025,5.])

    def test_checkpoint_resume_and_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'checkpoint.npz'; oxygen=np.array([0.,.025,5.])
            v,flags,done,gaps,residual,elapsed=f.checkpoint(p,'key',2,oxygen)
            v[:,0]=[[1,.2],[2,.3]]; flags[:,0]=1; done[0]=True; gaps[0]=1e-8
            f.atomic_npz(p,identity=np.array('key'),oxygen=oxygen,outputs=v,solver_flags=flags,done=done,
                         minimum_gap=gaps,maximum_residual=residual,elapsed_seconds=np.array(12.))
            loaded=f.checkpoint(p,'key',2,oxygen); np.testing.assert_array_equal(loaded[0],v)
            self.assertEqual(loaded[-1],12.)
            with self.assertRaisesRegex(ValueError,'identity'): f.checkpoint(p,'other',2,oxygen)

    def test_nonconverged_and_undefined_indices_retained(self):
        values=np.ones((2,5,3,2,119))*.02
        values[-1,4,0,0,0]=.5; values[-1,:,1,0,1]=np.nan; values[-1,:,2,0,2]=-.01
        mean,spread,delta,status=f.diagnostics(values)
        self.assertEqual(status[0,0,0],2); self.assertEqual(status[1,0,1],0)
        self.assertTrue(np.isnan(mean[1,0,1])); self.assertEqual(mean[2,0,2],-.01)
        self.assertEqual(mean.shape,(3,2,119))

    def test_whole_case_reconstruction_and_index_export(self):
        old_n,old_ns,old_oxygen=f.N,f.NS,c.OXYGEN
        try:
            f.N=128; f.NS=[64,128]; c.OXYGEN=[0.,.025,5.]
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); (root/'slurm').mkdir(); (root/'inputs').mkdir()
                bounds=[dict(parameter=p,encoded_lower=0.,encoded_upper=1.,natural_lower=0.,natural_upper=1.,transform='identity') for p in c.ACTIVE]
                case=dict(scope='global',case_id='global',fit_seed='seed25',best_params_sha256='fit',fit_config_sha256='cfg')
                directory=f.case_dir(root,case); directory.mkdir(parents=True)
                c.table(directory/'parameter_bounds.tsv',bounds); case['bounds_sha256']=c.sha(directory/'parameter_bounds.tsv')
                p=dict(unit_sha256={},git_commit='test',engine_sha256='engine',sif=dict(sha256='sif'),oxygen_pct=c.OXYGEN)
                for repeat in range(1,6):
                    x=c.samples(bounds,128,2026101300+repeat)
                    unit=root/'inputs'/('unit_R%02d.npz'%repeat); f.atomic_npz(unit,unit=x)
                    p['unit_sha256'][str(repeat)]=c.sha(unit)
                c.save_json(root/'slurm/submission_plan.json',p)
                for repeat in range(1,6):
                    for block,name in enumerate(c.BLOCKS):
                        x=f.samples(root,case,repeat,block)
                        y=np.stack([x[:,0]+2*x[:,1],x[:,0]*x[:,1]],axis=1)
                        output=np.repeat(y[:,None,:],3,axis=1)
                        path=directory/'runs'/('R%02d'%repeat)/'blocks'/name/'outputs.npz'
                        f.atomic_npz(path,outputs=output,oxygen=np.asarray(c.OXYGEN))
                        c.save_json(str(path)+'.receipt.json',dict(status='complete',identity=f.identity(p,case,repeat,block,x),sha256=c.sha(path),
                            n_points=384,n_recovered=0,minimum_spectral_gap=1e-8,maximum_residual=1e-14))
                with patch.object(full_report,'plot_indices'):
                    full_report.summarize_case(root,p,case)
                with np.load(directory/'indices_summary.npz') as z:
                    np.testing.assert_allclose(z['mean'][:,0,:2],np.tile([.2,.8],(3,1)),atol=.04)
                    self.assertAlmostEqual(z['mean'][0,1,28],1/7,delta=.04)
                    self.assertEqual(z['mean'].shape,(3,2,119))
                rows=c.read(directory/'indices_summary.tsv.gz')
                self.assertEqual(len(rows),3*2*119)
                self.assertTrue(all(r['included_in_full_analysis']=='TRUE' for r in rows))
                self.assertTrue((directory/'completion.json').exists())
        finally:
            f.N,f.NS,c.OXYGEN=old_n,old_ns,old_oxygen


if __name__=='__main__': unittest.main()
