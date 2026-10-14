#!/usr/bin/env python3
"""Independent global output isolation, interrupted resume and source export checks."""
from contextlib import contextmanager
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
import core as c
import full_core as f
import full_report
import global_postprocess as g


@contextmanager
def fixture():
    previous=(f.N,f.NS,c.OXYGEN)
    f.N=128; f.NS=[64,128]; c.OXYGEN=[0.,.025,5.]
    try:
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            case=dict(scope='global',case_id='global',fit_seed='seed25',
                      best_params_sha256='best',fit_config_sha256='cfg')
            source=f.case_dir(root,case); source.mkdir(parents=True)
            bounds=[dict(parameter=p,encoded_lower=0.,encoded_upper=1.,
                         natural_lower=0.,natural_upper=1.,transform='identity') for p in c.ACTIVE]
            c.table(source/'parameter_bounds.tsv',bounds)
            case['bounds_sha256']=c.sha(source/'parameter_bounds.tsv')
            p=dict(unit_sha256={},git_commit='model-version',engine_sha256='engine',
                   sif=dict(sha256='sif'),oxygen_pct=c.OXYGEN)
            for repeat in range(1,6):
                unit=root/'inputs'/('unit_R%02d.npz'%repeat)
                f.atomic_npz(unit,unit=c.samples(bounds,f.N,2026101300+repeat))
                p['unit_sha256'][str(repeat)]=c.sha(unit)
            c.save_json(root/'slurm/submission_plan.json',p)
            for repeat in range(1,6):
                for bi,block in enumerate(c.BLOCKS):
                    x=f.samples(root,case,repeat,bi)
                    y=np.stack([x[:,0]+2*x[:,1],x[:,0]*x[:,1]],axis=1)
                    path=source/'runs'/('R%02d'%repeat)/'blocks'/block/'outputs.npz'
                    f.atomic_npz(path,outputs=np.repeat(y[:,None,:],3,axis=1),oxygen=np.asarray(c.OXYGEN))
                    c.save_json(str(path)+'.receipt.json',dict(status='complete',
                                identity=f.identity(p,case,repeat,bi,x),sha256=c.sha(path),
                                n_points=f.N*3,n_recovered=0,minimum_spectral_gap=1e-8,
                                maximum_residual=1e-14))
            yield root,p,case,source
    finally:
        f.N,f.NS,c.OXYGEN=previous


class GlobalChecks(unittest.TestCase):
    def test_separate_outputs_idempotence_and_known_model(self):
        with fixture() as (root,p,case,source):
            before={str(x.relative_to(source)):c.sha(x) for x in source.rglob('*') if x.is_file()}
            target=root/'global/standalone_N4096'; tag={'postprocess':'version-one'}
            with patch.object(full_report,'plot_indices'):
                full_report.summarize_case(root,p,case,target_dir=target,analysis_metadata=tag)
            self.assertFalse((source/'completion.json').exists())
            after={str(x.relative_to(source)):c.sha(x) for x in source.rglob('*') if x.is_file()}
            self.assertEqual(before,after)
            with np.load(target/'indices_summary.npz') as z:
                np.testing.assert_allclose(z['mean'][:,0,:2],np.tile([.2,.8],(3,1)),atol=.04)
                self.assertAlmostEqual(z['mean'][0,1,28],1/7,delta=.04)
            with patch.object(c,'analyze',side_effect=AssertionError('Completed analysis was recomputed')):
                full_report.summarize_case(root,p,case,target_dir=target,analysis_metadata=tag)
            with self.assertRaisesRegex(ValueError,'identity changed'):
                full_report.summarize_case(root,p,case,target_dir=target,analysis_metadata={'postprocess':'version-two'})
            with open(target/'indices_summary.tsv.gz','ab') as handle: handle.write(b'changed')
            with self.assertRaisesRegex(ValueError,'Summary changed'):
                full_report.summarize_case(root,p,case,target_dir=target,analysis_metadata=tag)

    def test_interrupted_resume_retains_completed_repeat(self):
        with fixture() as (root,p,case,source):
            target=root/'global/standalone_N4096'; original=c.analyze; count=[0]
            def interrupt(*args,**kwargs):
                count[0]+=1
                if count[0]==13: raise RuntimeError('Deliberate interruption after one complete repeat')
                return original(*args,**kwargs)
            with patch.object(full_report,'plot_indices'),patch.object(c,'analyze',side_effect=interrupt):
                with self.assertRaisesRegex(RuntimeError,'Deliberate interruption'):
                    full_report.summarize_case(root,p,case,target_dir=target)
            with np.load(target/'indices_checkpoint.npz') as z:
                self.assertTrue(z['done'][0].all())
                self.assertFalse(z['done'][1:].any())
                saved=z['values'][:,0].copy()
            with patch.object(full_report,'plot_indices'),patch.object(c,'analyze',wraps=original) as calls:
                full_report.summarize_case(root,p,case,target_dir=target)
                self.assertEqual(calls.call_count,48)
            with np.load(target/'indices_repeats.npz') as z:
                np.testing.assert_array_equal(z['values'][:,0],saved)
            self.assertFalse((target/'indices_checkpoint.npz').exists())

    def test_repeat_export_retains_negative_undefined_and_intervals(self):
        shape=(2,5,3,2,len(f.METRICS)); values=np.zeros(shape); widths=np.ones(shape)*.02
        variance=np.ones(shape[:-1]); values[1,2,1,0,28]=-.03; values[:,:,2,1,:]=np.nan
        rows=list(g.repeat_rows(values,widths,variance,[0,.025,5],[2048,4096]))
        self.assertEqual(len(rows),2*5*3*2*119)
        r=next(r for r in rows if (r['N'],r['repeat'],r['O2_pct'],r['output'],r['index'],r['parameter_i'],r['parameter_j'])
               ==(4096,3,.025,c.OUTPUTS[0],'S2',c.ACTIVE[0],c.ACTIVE[1]))
        self.assertEqual(r['value'],-.03); self.assertAlmostEqual(r['bootstrap_confidence_low'],-.05)
        self.assertAlmostEqual(r['bootstrap_confidence_high'],-.01)
        self.assertTrue(any(np.isnan(r['value']) for r in rows))
        self.assertTrue(all(r['included']=='TRUE' for r in rows))


if __name__=='__main__': unittest.main()
