#!/usr/bin/env python3
"""Configuration, cache identity, FFT and recovery tests for the N/M follow-up."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import efast
import neighborhood as nh
import neighborhood_slurm as ns
import resolution_nm as run
import resolution_pilot as rp


class NMChecks(unittest.TestCase):
    def test_scope_and_phase_mapping(self):
        total=0
        for config,spec in run.CONFIGS.items():
            tasks,summaries=run.task_rows(config)
            self.assertEqual(len(tasks),140)
            self.assertEqual([s['fit_seed'] for s in summaries],list(spec['seeds']))
            self.assertEqual([t['task_id'] for t in tasks],list(range(1,141)))
            self.assertEqual(len({(t['fit_seed'],t['replicate'],t['parameter_index']) for t in tasks}),140)
            total+=sum(t['N']*201 for t in tasks)
            for t in tasks:
                self.assertEqual((t['N'],t['M']),(spec['N'],spec['M']))
                self.assertEqual(t['phase_seed'],202610080+int(t['fit_seed'][4:])*10+t['replicate']-1)
                self.assertEqual(t['focal_parameter'],efast.ACTIVE[t['parameter_index']])
                self.assertEqual(t['parent_sample_id_end']-t['parent_sample_id_start']+1,spec['N'])
            for s in summaries:
                rows=[t for t in tasks if t['summary_id']==s['summary_id']]
                self.assertEqual(len(rows),70)
                self.assertEqual({t['replicate'] for t in rows},set(range(1,6)))
                self.assertEqual({t['focal_parameter'] for t in rows},set(efast.ACTIVE))
        self.assertEqual(total,345840600)

    def test_frequency_depends_on_M(self):
        self.assertEqual(run.frequency(8193,4),1024)
        self.assertEqual(run.frequency(4097,8),256)
        self.assertEqual(run.frequency(4097,4),512)
        with self.assertRaises(ValueError): run.frequency(257,8)
        with self.assertRaises(ValueError): run.frequency(513,8)

    def test_M_comparison_never_labelled_N_convergence_and_local_max_retained(self):
        oxygen=np.arange(201)*.025
        a=np.full((5,2,14,201,2),.2)
        b=a.copy()
        b[:,1,0,77,0]+=.3
        for config,spec in run.CONFIGS.items():
            row=next(r for r in run.comparison_rows(spec['seeds'][0],spec,a,b,oxygen)
                     if r['parameter']==efast.ACTIVE[0] and r['index']=='ST'
                     and r['output']==efast.OUTPUTS[0] and r['oxygen_window']==rp.WINDOWS[0][0])
            self.assertEqual(row['comparison_kind'],spec['comparison_kind'])
            self.assertEqual((row['N_reference'],row['M_reference']),(4097,4))
            self.assertEqual((row['N_test'],row['M_test']),(spec['N'],spec['M']))
            self.assertNotIn('N_low',row)
            self.assertTrue(row['p90_pass'])
            self.assertEqual(row['diagnostic_status'],'localized_review')
            self.assertAlmostEqual(row['absolute_delta_max_O2_pct'],1.925)

    def test_prepared_signature_and_parent_identity_include_M(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            seed='seed101'
            fit=root/'fit'/seed
            fit.mkdir(parents=True)
            for name in ('best_params.tsv','fit_config.rds','parameter_table.csv'):
                (fit/name).write_text('test')
            bounds=root/'bounds'/(seed+'.tsv')
            bounds.parent.mkdir()
            bounds.write_text('test')
            (root/'slurm').mkdir()
            (root/'slurm/submission_plan.json').write_text('test')
            row=dict(fit_seed=seed,best_params_sha256=efast.sha256(fit/'best_params.tsv'),
                     fit_config_sha256=efast.sha256(fit/'fit_config.rds'),
                     parameter_table_sha256=efast.sha256(fit/'parameter_table.csv'),
                     bounds_sha256=efast.sha256(bounds))
            run.write_rows(root/'seed_selection.tsv',[row])
            plan=dict(M=8,baseline_dir=str(root),fit_root=str(root/'fit'),
                      figure4_dir='unused',result_root=str(root))
            with patch.object(efast,'prepare') as prepare:
                target=run.ensure_design(root,seed,4097,plan)
                self.assertEqual(prepare.call_args.args[0].m,8)
                self.assertIn('N4097_M8',str(target))
                self.assertEqual(json.loads((target/'prepared.json').read_text())['M'],8)
                run.ensure_design(root,seed,4097,plan)
                self.assertEqual(prepare.call_count,1)
                signature=json.loads((target/'prepared.json').read_text())
                signature['M']=4
                (target/'prepared.json').write_text(json.dumps(signature))
                with self.assertRaisesRegex(ValueError,'signature'): run.ensure_design(root,seed,4097,plan)
            parent=root/'parent'
            parent.mkdir()
            task=run.task_rows('N4097_M8')[0][0]
            meta=dict(fit_seed=seed,N=4097,M=8,replicate=1,seed=task['phase_seed'],
                      n_samples=4097*14,n_parameters=14,salib_version='1.5.2',
                      bounds_sha256=row['bounds_sha256'],best_params_sha256=row['best_params_sha256'],
                      fit_config_sha256=row['fit_config_sha256'],oxygen_pct=(np.arange(201)*.025).tolist())
            (parent/'metadata.json').write_text(json.dumps(meta))
            run.verify_parent(parent,task,plan)
            meta['M']=4
            (parent/'metadata.json').write_text(json.dumps(meta))
            with self.assertRaisesRegex(ValueError,'identity'): run.verify_parent(parent,task,plan)

    def test_submit_resources_dependencies_and_idempotency(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            binary=root/'bin'
            binary.mkdir()
            for name,body in {
                'apptainer':'#!/bin/bash\nexit 0\n',
                'flock':'#!/bin/bash\nexit 0\n',
                'sbatch':'''#!/bin/bash
printf '%s\\n' "$*" >> "$MOCK_SUBMISSIONS"
count=$(wc -l < "$MOCK_SUBMISSIONS")
echo $((9000+count))
''',
            }.items():
                path=binary/name
                path.write_text(body)
                path.chmod(0o755)
            record=root/'submitted.txt'
            env=dict(os.environ,PATH=str(binary)+os.pathsep+os.environ['PATH'],
                     EFAST_PILOT_ROOT=str(root/'pilot'),MOCK_SUBMISSIONS=str(record))
            script=run.HERE/'submit_resolution_nm.sh'
            for repeat in range(2):
                subprocess.run(['bash',str(script),'submit'],env=env,check=True,capture_output=True,text=True)
                calls=record.read_text().splitlines()
                self.assertEqual(len(calls),8)
            self.assertEqual(sum('--array=1-140 ' in c for c in calls),2)
            self.assertEqual(sum('--array=1-2 ' in c for c in calls),2)
            self.assertEqual(sum('--dependency=afterok:' in c for c in calls),2)
            self.assertEqual(sum('--dependency=afterany:' in c for c in calls),4)
            for call in calls:
                self.assertIn('--qos=xxlarge',call)
                self.assertIn('--time=12:00:00',call)
                self.assertNotIn('--nodelist',call)
                self.assertNotIn('--node=',call)
                array=[arg for arg in call.split() if arg.startswith('--array=')]
                self.assertTrue(all('%' not in arg for arg in array))

    @unittest.skipUnless(importlib.util.find_spec('SALib'),'Requires approved SALib 1.5.2 SIF')
    def test_targeted_fallback_M8_impact_has_no_spurious_frequency_change(self):
        # Mock only the R recovery handoff; the real 50/100-digit solver is
        # exercised on an actual canonical matrix in each compute-node smoke.
        n,m=4097,8
        s=2*np.pi*np.arange(n)/n
        y=np.sin(run.frequency(n,m)*s)+.5*np.sin(s)
        with tempfile.TemporaryDirectory() as directory:
            part=Path(directory)
            metadata=dict(N=n,M=m,oxygen_pct=[0.],fit_seed='seed101',samples_sha256='fixture')
            (part/'metadata.json').write_text(json.dumps(metadata))
            original=part/'outputs.tsv.gz.part01.gz'
            rows=(dict(sample_id=i+1,O2_pct=0.,dominant_mean_ploidy=2+y[i],
                       dominant_growth_rate=.1+.02*y[i],spectral_gap=.01,
                       eigenvector_nonnegative='FALSE' if i==0 else 'TRUE',status='ok')
                  for i in range(n))
            efast.write_table(original,ns.OUTPUT_FIELDS,rows)
            def recover(*args,**kwargs):
                recovery=part/'recovery'
                bad=json.loads((recovery/'failed_points.json').read_text())[0]
                (recovery/'matrix.rds').write_bytes(b'fixture')
                report=dict(sample_id=1,O2_pct=0.,original=bad,
                            dominant_mean_ploidy=float(bad[efast.OUTPUTS[0]]),
                            dominant_growth_rate=float(bad[efast.OUTPUTS[1]]),
                            matrix_file='matrix.rds',matrix_sha256=efast.sha256(recovery/'matrix.rds'))
                (recovery/'corrections.json').write_text(json.dumps(dict(status='passed',n_points=1,reports=[report])))
            with patch.object(run.subprocess,'run',side_effect=recover):
                run.fallback(part,part)
            proof=json.loads((part/'recovery/numerical_recovery.json').read_text())
            self.assertEqual(proof['n_recomputed_points'],1)
            self.assertEqual(proof['n_unchanged_good_rows'],n-1)
            for impact in proof['index_impact']:
                self.assertLessEqual(impact['max_absolute_index_difference'],1e-12)
                self.assertFalse(impact['undefined_pattern_changed'])

    @unittest.skipUnless(importlib.util.find_spec('SALib'),'Requires approved SALib 1.5.2 SIF')
    def test_Fourier_recovery_diagnostic_uses_M8_frequency(self):
        from SALib.analyze.fast import compute_orders
        for config,spec in run.CONFIGS.items():
            n,m=spec['N'],spec['M']
            omega=run.frequency(n,m)
            s=2*np.pi*np.arange(n)/n
            y=np.sin(omega*s)+.5*np.sin(s)+.3*np.sin(2*omega*s)
            np.testing.assert_allclose(run.orders(y,n,m),[1.09/1.34]*2,atol=1e-12,rtol=0)
            if m==8:
                try:
                    wrong=compute_orders(y,n,m,(n-1)//8)
                except IndexError:
                    # The old M4 frequency can address bins beyond Nyquist.
                    pass
                else:
                    self.assertGreater(float(np.max(np.abs(np.asarray(wrong)-run.orders(y,n,m)))),.5)

    @unittest.skipUnless(importlib.util.find_spec('SALib'),'Requires approved SALib 1.5.2 SIF')
    def test_full_sampler_and_SALib_analyzer_known_variance_fractions(self):
        from SALib.sample.fast_sampler import sample
        from SALib.analyze.fast import analyze
        problem=dict(num_vars=14,names=list(efast.ACTIVE),bounds=[[0.,1.]]*14)
        for config,spec in run.CONFIGS.items():
            n,m=spec['N'],spec['M']
            for repeat in range(5):
                x=sample(problem,n,M=m,seed=202610080+int(spec['seeds'][0][4:])*10+repeat)
                self.assertEqual(x.shape,(14*n,14))
                y=x[:,0]+2*x[:,1]
                result=analyze(problem,y,M=m,num_resamples=2,seed=13579)
                for p,expected in ((0,.2),(1,.8),(13,0.)):
                    actual=run.orders(y[p*n:(p+1)*n],n,m)
                    np.testing.assert_allclose(actual,[result['S1'][p],result['ST'][p]],atol=1e-12,rtol=0)
                    np.testing.assert_allclose(actual,[expected]*2,atol=.01,rtol=0)

    @unittest.skipUnless(importlib.util.find_spec('SALib'),'Requires approved SALib 1.5.2 SIF')
    def test_actual_trajectory_analyzer_matches_M8_orders_and_undefined_policy(self):
        n,m=4097,8
        s=2*np.pi*np.arange(n)/n
        y=np.sin(run.frequency(n,m)*s)+.5*np.sin(s)
        with tempfile.TemporaryDirectory() as directory:
            part=Path(directory)
            (part/'metadata.json').write_text(json.dumps(dict(N=n,M=m,oxygen_pct=[0.,2.5,5.])))
            rows=(dict(sample_id=i+1,O2_pct=o,dominant_mean_ploidy=2+y[i],
                       dominant_growth_rate=1. if o==5. else .1+.02*y[i],
                       spectral_gap=.01,eigenvector_nonnegative='TRUE',status='ok')
                  for i in range(n) for o in (0.,2.5,5.))
            efast.write_table(part/'outputs.tsv.gz',ns.OUTPUT_FIELDS,rows)
            result=ns.analyze_trajectory(part)
            np.testing.assert_allclose(result['indices'][:,0,0],run.orders(2+y,n,m),atol=1e-12,rtol=0)
            np.testing.assert_allclose(result['indices'][:,1,1],run.orders(.1+.02*y,n,m),atol=1e-12,rtol=0)
            self.assertTrue(np.isnan(result['indices'][:,2,1]).all())
            self.assertEqual(result['variance'][2,1],0.)


if __name__=='__main__': unittest.main()
