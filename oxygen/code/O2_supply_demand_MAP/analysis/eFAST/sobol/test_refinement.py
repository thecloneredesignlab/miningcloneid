#!/usr/bin/env python3
"""Focused tests for selection, exact-tail reuse, ordering and submission scope."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import refinement_run as r
import submit_refinement as s


class PolicyTests(unittest.TestCase):
    def row(self, pair, oi=10, output=r.OUTPUTS[0], status='not_passed', value='-.1'):
        return dict(index='S2',parameter_i=pair[0],parameter_j=pair[1],O2_pct=str(oi*.025),
                    output=output,convergence_status=status,value_mean=value)

    def test_union_has_no_output_or_pair_duplicates(self):
        rows=[self.row(r.EXPECTED_PAIRS[0]),self.row(r.EXPECTED_PAIRS[1]),
              self.row(r.EXPECTED_PAIRS[0],output=r.OUTPUTS[1])]
        found=r.selected_rows(rows,set(r.EXPECTED_PAIRS))
        self.assertEqual(list(found),[10]); self.assertEqual(len(found[10]),3)

    def test_selection_keeps_negative_and_undefined_failures(self):
        self.assertEqual(sorted(r.selected_rows([self.row(r.EXPECTED_PAIRS[0],10),
            self.row(r.EXPECTED_PAIRS[0],11,status='undefined')],set(r.EXPECTED_PAIRS))),[10,11])

    def test_passed_and_unregistered_pairs_not_selected(self):
        self.assertEqual(r.selected_rows([self.row(r.EXPECTED_PAIRS[0],status='passed'),
                         self.row(('lam_max','p_wgd'))],set(r.EXPECTED_PAIRS)),{})

    def test_off_grid_fails(self):
        row=self.row(r.EXPECTED_PAIRS[0]); row['O2_pct']='0.123'
        with self.assertRaises(ValueError): r.selected_rows([row],set(r.EXPECTED_PAIRS))

    def test_global_tasks_cover_only_new_samples_once(self):
        seen=set()
        for n,count in [(8192,150),(16384,300),(32768,600)]:
            for task in range(1,count+1):
                cid,rep,block,begin,oxygen=r.decode_model('global',n,task,[])
                key=(rep,block,begin)
                self.assertNotIn(key,seen); seen.add(key)
                self.assertGreaterEqual(begin,4096); self.assertEqual(oxygen,list(range(201)))
        self.assertEqual(len(seen),1050)
        for rep in range(1,6):
            for block in range(30):
                self.assertEqual(sorted(x[2] for x in seen if x[:2]==(rep,block)),list(range(4096,32768,4096)))

    def test_neighborhood_targets_not_padded_or_extended(self):
        cases=[dict(case=dict(case_id='seed283'),oxygen_indices=[0,19,200])]
        for task in range(1,151):
            cid,rep,block,begin,ois=r.decode_model('neighborhood10pct',8192,task,cases)
            self.assertEqual((cid,begin,ois),('seed283',4096,[0,19,200]))
        with self.assertRaises(ValueError): r.decode_model('neighborhood10pct',8192,151,cases)
        with self.assertRaises(ValueError): r.decode_model('neighborhood10pct',16384,1,cases)

    def test_stage_dependencies_include_all_required_tails(self):
        stages=s.stage_specs(dict(cases=[{}]*497)); by={x[0]:x for x in stages}
        self.assertEqual(by['model_n8192'][4],74550)
        self.assertEqual(by['analyze_n8192'][4],2485)
        self.assertTrue({'model_g8192','model_g16384','model_g32768','analyze_g16384'}<=set(by['analyze_g32768'][-1]))
        self.assertEqual(sum(x[4] for x in stages if x[1]=='model'),75600)
        self.assertEqual(sum(x[4] for x in stages if x[1]=='analyze'),2875)

    def test_resources_no_node_or_concurrency_limit(self):
        stage=s.stage_specs(dict(cases=[{}]*497))[1]
        args=s.sbatch_args(Path('/tmp/frozen_refinement'),stage,['123'])
        self.assertIn('--qos=xxlarge',args); self.assertIn('--time=12:00:00',args)
        self.assertIn('--hold',args)
        self.assertFalse(any(x.startswith(('--nodelist','--nodes=')) for x in args))
        self.assertFalse(any('%' in x for x in args if x.startswith('--array=')))

    def test_checkpoint_and_final_paths_use_absolute_runtime_snapshot(self):
        shell=Path(__file__).with_name('run_refinement.sbatch').read_text()
        self.assertIn('source "$snapshot/retry_container_startup.sh"',shell)
        self.assertIn('"$snapshot/refinement_run.py"',shell)
        self.assertNotIn('dirname "$0"',shell)
        self.assertIn('trap',shell)


@unittest.skipUnless(importlib.util.find_spec('numpy') and importlib.util.find_spec('SALib'),
                     'Numerical design checks require original SIF')
class NumericalTests(unittest.TestCase):
    def test_five_scrambles_exact_prefix_and_interleaving(self):
        import numpy as np
        from SALib.sample import sobol
        problem=dict(num_vars=14,names=['p%d'%i for i in range(14)],bounds=[[0.,1.]]*14)
        for repeat in range(1,6):
            low=sobol.sample(problem,4096,calc_second_order=True,scramble=True,seed=2026101300+repeat)
            high=sobol.sample(problem,32768,calc_second_order=True,scramble=True,seed=2026101300+repeat)
            self.assertTrue(np.array_equal(low,high[:4096*30]))
            pieces=[low]+[high[a*30:(a+4096)*30] for a in range(4096,32768,4096)]
            self.assertTrue(np.array_equal(np.concatenate(pieces),high))
            # Matrix blocks must be interleaved within each base sample, not appended as blocks.
            tail=high[4096*30:8192*30]
            fused=np.stack([tail[b::30] for b in range(30)],axis=1).reshape((-1,14))
            self.assertTrue(np.array_equal(fused,tail))

    def test_combined_output_matches_one_shot_sobol_without_old_evaluations(self):
        import numpy as np
        from SALib.sample import sobol
        from SALib.analyze import sobol as analyzer
        problem=dict(num_vars=3,names=['a','b','c'],bounds=[[0.,1.]]*3)
        x=sobol.sample(problem,256,calc_second_order=True,scramble=True,seed=42)
        # Synthetic model calls are counted: only the unseen tail can be evaluated.
        f=lambda a:a[:,0]+a[:,1]*a[:,2]+a[:,0]*a[:,1]
        prefix=f(x[:128*8]); count=0
        def evaluate_new(a):
            nonlocal count
            count+=len(a); return f(a)
        y=np.concatenate([prefix,evaluate_new(x[128*8:])])
        self.assertEqual(count,128*8)
        a=analyzer.analyze(problem,y,calc_second_order=True,seed=77,num_resamples=20)
        b=analyzer.analyze(problem,f(x),calc_second_order=True,seed=77,num_resamples=20)
        for key in ('S1','ST','S2','S1_conf','ST_conf','S2_conf'):
            self.assertTrue(np.allclose(a[key],b[key],equal_nan=True))


if __name__=='__main__': unittest.main(verbosity=2)
