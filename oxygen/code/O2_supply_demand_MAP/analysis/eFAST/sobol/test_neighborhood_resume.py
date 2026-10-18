#!/usr/bin/env python3
"""Legacy provenance, cache ordering, interruption and no-recomputation tests."""
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
import core as c
import full_core as f
import full_report
import recover_neighborhood as r
from test_global import fixture


class ResumeChecks(unittest.TestCase):
    def test_legacy_identity_temporary_and_zero_variance(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp); values,widths,variance,done=r.blank()
            done[0,0,:]=True; variance[:,0,0,:]=0
            # Undefined indices with certified zero variance are valid and retained.
            f.atomic_npz(source/'indices_checkpoint.npz',identity=np.array('B-key'),
                         values=values,widths=widths,variance=variance,done=done)
            got,record=r.load_legacy(source,{'B-key'})
            self.assertEqual(record['n_done'],2)
            done[0,1,:]=True; variance[:,0,1,:]=0
            f.atomic_npz(source/'indices_checkpoint.npz.incomplete',identity=np.array('B-key'),
                         values=values,widths=widths,variance=variance,done=done)
            got,record=r.load_legacy(source,{'B-key'})
            self.assertEqual(record['n_done'],4)
            self.assertTrue(record['path'].endswith('.incomplete'))
            with self.assertRaisesRegex(ValueError,'Legacy identity'):
                r.load_legacy(source,{'different-inputs'})
            (source/'indices_checkpoint.npz.incomplete').write_bytes(b'broken atomic temporary')
            got,record=r.load_legacy(source,{'B-key'})
            self.assertEqual(record['n_done'],2)
            self.assertEqual(len(record['rejected_temporary']),1)
            variance[:,0,0,:]=np.nan
            with self.assertRaisesRegex(ValueError,'certified'):
                r.validate_arrays(values,widths,variance,done)

    def test_cache_matches_original_interleaving_and_checks_source(self):
        with fixture() as (root,p,case,source):
            case['case_index']=1; p['cases']=[case]
            recovery=root/'recovery'; task=dict(task_id=1,case_index=1,case_id='global',repeat=1,
                oxygen_indices=[0,2],needed_outputs=[[0,1],[0,1]])
            rp=dict(tasks=[task]); c.save_json(recovery/'plan.json',rp)
            path,oxygen,receipt=r.source_cache(root,recovery,p,rp,task)
            arrays=[]
            for block in c.BLOCKS:
                with np.load(source/'runs/R01/blocks'/block/'outputs.npz') as z: arrays.append(z['outputs'])
            expected=c.assemble(arrays)[:,[0,2],:].transpose(1,0,2)
            np.testing.assert_array_equal(np.load(path),expected)
            path2,oxygen2,_=r.source_cache(root,recovery,p,rp,task)
            self.assertEqual(path,path2); self.assertEqual(oxygen,oxygen2)
            path.chmod(0o644)
            mm=np.load(path,mmap_mode='r+'); mm[0,0,0]+=1; mm.flush(); del mm
            with self.assertRaisesRegex(ValueError,'cache slice changed'):
                r.shard(root,recovery,p,rp,1)

    def test_partial_shards_resume_and_keep_all_legacy_values(self):
        with fixture() as (root,p,case,source):
            case['case_index']=1; p['cases']=[case]; recovery=root/'recovery'
            bounds=c.read(source/'parameter_bounds.tsv'); expected=r.blank()
            ordered={}
            for repeat in range(1,6):
                arrays=[]
                for block in c.BLOCKS:
                    with np.load(source/'runs'/('R%02d'%repeat)/'blocks'/block/'outputs.npz') as z: arrays.append(z['outputs'])
                ordered[repeat]=c.assemble(arrays).transpose(1,0,2)
                for oi in range(3):
                    for k in range(2):
                        for ni,n in enumerate(f.NS):
                            si,var=c.analyze(bounds,ordered[repeat][oi,:n*c.STEP,k],seed=2026101300+repeat+10000+n+oi*10+k)
                            expected[0][ni,repeat-1,oi,k],expected[1][ni,repeat-1,oi,k]=f.metric_arrays(si)
                            expected[2][ni,repeat-1,oi,k]=var
                        expected[3][repeat-1,oi,k]=True
            legacy=list(r.blank()); legacy[3][0]=True; legacy[3][1,0]=True
            for a,truth in zip(legacy[:3],expected[:3]): a[:,legacy[3]]=truth[:,legacy[3]]
            checkpoint=recovery/'checkpoints/global.npz'
            f.atomic_npz(checkpoint,values=legacy[0],widths=legacy[1],variance=legacy[2],done=legacy[3])
            tasks=[]
            for repeat in range(1,6):
                oxygen=np.flatnonzero(~legacy[3][repeat-1].all(axis=1)).tolist()
                for start in range(0,len(oxygen),2):
                    oo=oxygen[start:start+2]
                    tasks.append(dict(task_id=len(tasks)+1,case_id='global',case_index=1,repeat=repeat,
                        oxygen_indices=oo,needed_outputs=[[k for k in range(2) if not legacy[3][repeat-1,j,k]] for j in oo]))
            item=dict(case=case,checkpoint_sha256=c.sha(checkpoint),task_ids=list(range(1,len(tasks)+1)),legacy=dict(n_done=int(legacy[3].sum())))
            rp=dict(tasks=tasks,cases=[item]); c.save_json(recovery/'plan.json',rp)
            calls=[0]
            def analyze(y,seed):
                calls[0]+=1
                if calls[0]==3: raise RuntimeError('simulate interruption between analyzer cells')
                return c.analyze(bounds,y,seed=seed)
            task=tasks[0]
            with self.assertRaisesRegex(RuntimeError,'simulate interruption'):
                r.calculate_shard(recovery,task,ordered[task['repeat']],{i:i for i in range(3)},analyze)
            calls=[0]
            def resumed(y,seed):
                calls[0]+=1; return c.analyze(bounds,y,seed=seed)
            for task in tasks:
                r.calculate_shard(recovery,task,ordered[task['repeat']],{i:i for i in range(3)},resumed)
            self.assertEqual(calls[0],int((~legacy[3]).sum())*2-2)
            def forbidden(*args,**kwargs): raise AssertionError('Completed cells were recalculated')
            for task in tasks: r.calculate_shard(recovery,task,None,None,forbidden)
            merged=r.merge_case(recovery,rp,item)
            for actual,truth in zip(merged,expected[:3]): np.testing.assert_array_equal(actual,truth)
            with patch.object(full_report,'plot_indices'):
                r.finalize_case(root,recovery,p,rp,1)
            r.finalize_case(root,recovery,p,rp,1)
            with np.load(source/'indices_summary.npz') as z:
                np.testing.assert_allclose(z['mean'][:,0,:2],np.tile([.2,.8],(3,1)),atol=.04)
                self.assertAlmostEqual(z['mean'][0,1,28],1/7,delta=.04)
            self.assertTrue(c.load_json(source/'completion.json')['all_original_checkpoint_cells_preserved'])


if __name__=='__main__': unittest.main()
