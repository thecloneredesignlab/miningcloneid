"""Contract checks for scheduling, complete trajectories and scientific summaries."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from common import (ACTIVE, OUTPUTS, INITIAL_NS, CHUNK, decode, model_count, next_n, valid_n,
                    resolutions,previous_resolution,efast_directory)
from submit import stages, sbatch_args


class DesignContracts(unittest.TestCase):
    def test_complete_nonoverlapping_trajectories(self):
        for n in INITIAL_NS+(32769,65537,131073):
            groups={}
            for task in range(1,model_count(n)+1):
                rep,pi,begin,end=decode(n,task)
                groups.setdefault((rep,pi),[]).append((begin,end))
            self.assertEqual(len(groups),70)
            for spans in groups.values():
                self.assertEqual(spans[0][0],0)
                self.assertEqual(spans[-1][1],n)
                self.assertTrue(all(left[1]==right[0] for left,right in zip(spans,spans[1:])))
                self.assertEqual(sum(end-begin for begin,end in spans),n)

    def test_no_upper_n_bound(self):
        for n in [32769,65537,131073,1048577]:
            self.assertEqual(valid_n(n),n)
            self.assertEqual(next_n(n),2*n-1)
            self.assertGreater(model_count(n),0)

    def test_reject_invalid_positions(self):
        for n,task in [(2049,0),(2049,model_count(2049)+1)]:
            with self.assertRaises(ValueError):decode(n,task)

    def test_refinement_bridges_previous_batch(self):
        p=dict(reference_resolutions=list(INITIAL_NS),resolutions=[32769,65537,131073])
        self.assertEqual(resolutions(p),list(INITIAL_NS)+p['resolutions'])
        self.assertEqual(previous_resolution(p,32769),16385)
        self.assertEqual(previous_resolution(p,65537),32769)
        self.assertEqual(previous_resolution(p,131073),65537)
        self.assertEqual(efast_directory('/new',p,16385),Path('/new/inputs/previous_efast/N16385'))
        self.assertEqual(efast_directory('/new',p,32769),Path('/new/efast/N32769'))
        with self.assertRaises(ValueError):resolutions(dict(reference_resolutions=[16385],resolutions=[16385]))

    def test_refinement_runs_all_three_levels_without_sobol_jobs(self):
        ns=[32769,65537,131073];plan=stages(ns,refinement=True);by_name={s[0]:s for s in plan}
        self.assertEqual(len(plan),12)
        self.assertEqual(sum(s[3] for s in plan if s[1]=='model'),15890)
        self.assertFalse(any('sobol' in s[1] for s in plan))
        self.assertEqual(by_name['model_32769'][-1],['smoke'])
        self.assertEqual(by_name['model_65537'][-1],['report_32769'])
        self.assertEqual(by_name['model_131073'][-1],['report_65537'])
        self.assertEqual(sum(ns)*14*5*201,3227362530)
        for s in plan:
            args=sbatch_args(Path('/scope/global/refinement'),s,[])
            self.assertIn('--qos=xxlarge',args);self.assertIn('--time=12:00:00',args)
            self.assertFalse(any(a.startswith(('--nodelist','--constraint')) for a in args))
            self.assertFalse(any('%' in a for a in args if a.startswith('--array=')))

    def test_oxygen_trend_synthesis_uses_new_values(self):
        from finish_report import oxygen_preference_comparison
        def values(high):
            return [dict(output=o,index='ST',parameter_i='gamma_mu',oxygen_band=b,band_mean=v)
                    for o in OUTPUTS for b,v in [('low',.3),('middle',.2),('high',high)]]
        self.assertTrue(all(r['same_trend'] for r in oxygen_preference_comparison(values(.1),values(.2))))
        self.assertFalse(any(r['same_trend'] for r in oxygen_preference_comparison(values(.1),values(.5))))

    def test_reference_cache_requires_identical_point_and_certificate(self):
        import workflow as w
        from common import save,sha,model_path
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)/'new';old=Path(temp)/'old'
            design=w.design_file(old,65,1);design.parent.mkdir(parents=True)
            x=np.arange(14,dtype=float)
            np.savez_compressed(design,natural=np.tile(x,(65*14,1)))
            save(root/'inputs/previous_batch/preparation.json',dict(design_sha256={str(design.relative_to(old)):sha(design)}))
            save(root/'inputs/previous_batch/plan.json',dict(scope='global',method='eFAST'))
            source=model_path(old,65,1,0,0);source.parent.mkdir(parents=True)
            values=np.arange(65*201*2,dtype=float).reshape(65,201,2)
            np.savez_compressed(source,outputs=values,solver_flags=np.ones((65,201),dtype=int),
                oxygen=np.arange(201)/40.,minimum_gap=np.ones(201),maximum_residual=np.full(201,1e-15))
            p=dict(reference_resolutions=[65],previous_root=str(old),oxygen_pct=[i/40. for i in range(201)])
            with patch.object(w,'np',np,create=True),patch.object(w,'f',create=True) as engine:
                engine.complete.return_value=True
                cached=w.reference_point(root,p,1,0,x)
                self.assertTrue(np.array_equal(cached[0],values[0]))
                self.assertTrue(cached[-1]['bitwise_parameter_equality_verified'])
                self.assertEqual(cached[-1]['reused_parameter_oxygen_points'],201)
                self.assertIsNone(w.reference_point(root,p,1,0,x+1e-10))
                engine.complete.return_value=False
                with self.assertRaises(ValueError):w.reference_point(root,p,1,0,x)

    def test_scheduling_contract(self):
        plan=stages(INITIAL_NS)
        self.assertEqual(sum(s[3] for s in plan if s[1]=='model'),2380)
        for stage in plan:
            args=sbatch_args(Path('/scope/global'),stage,[])
            self.assertIn('--qos=xxlarge',args)
            self.assertIn('--time=12:00:00',args)
            self.assertIn('--export=ALL',args)
            self.assertFalse(any(a.startswith(('--nodelist','--constraint')) for a in args))
            self.assertFalse(any('%' in a for a in args if a.startswith('--array=')))
            self.assertFalse(any(a.startswith('--export=ALL,') for a in args))
            if stage[1]=='model':self.assertEqual(stage[-1],['smoke'])
        self.assertTrue(all(s[1] not in ['neighborhood','sobol_model'] for s in plan))

    def test_r_model_body_preserved(self):
        here=Path(__file__).resolve().parent
        old=here.parent/'sobol/evaluate.R'
        if not old.exists():self.skipTest('Frozen canonical engine independently checked by preparation')
        expected=old.read_text().replace(
            '# Canonical fixed-O2 matrix, checked dominant eigenpair, unchanged-matrix recovery.',
            '# eFAST inputs; identical frozen Sobol fixed-O2 matrix and eigenpair recovery.')
        expected=expected.replace('here <- dirname(normalizePath(script))',
            "here <- normalizePath(Sys.getenv('EFAST_CANONICAL_HELPERS', unset='__missing_canonical_helpers__'))")
        expected=expected.replace("stopifnot(meta$method %in% c('Sobol_second_order', 'Sobol_response_surface'))",
                                  "stopifnot(identical(meta$method, 'eFAST'))")
        self.assertEqual((here/'evaluate_efast.R').read_text(),expected)

    def test_physical_oxygen_integration(self):
        import plot_helpers as gp
        oxygen=np.arange(201)/40.
        for lo,hi in [(0,1),(1,3),(3,5)]:
            w=gp.weights(oxygen,lo,hi)
            self.assertAlmostEqual(w.sum(),1.)
            self.assertAlmostEqual(float(w@oxygen),(lo+hi)/2)

    def test_monitor_scope_and_renewal(self):
        import time_extension as monitor
        args=monitor.sbatch_args('12345',predecessor='12346')
        self.assertIn('--qos=xxlarge',args)
        self.assertIn('--time=12:00:00',args)
        self.assertIn('--export=ALL',args)
        self.assertIn('--dependency=afterany:12346',args)
        self.assertEqual(monitor.HOURS,6)
        self.assertEqual(monitor.seconds('06:00:00'),21600)

    def test_frozen_report_ignores_canonical_name_collision(self):
        import sys
        from workflow import reporting
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);snapshot=root/'slurm/code_snapshot';snapshot.mkdir(parents=True)
            (snapshot/'report.py').write_text("SOURCE='correct_global_finalize'\n")
            wrong=root/'canonical';wrong.mkdir()
            (wrong/'report.py').write_text("raise RuntimeError('wrong canonical report imported')\n")
            sys.path.insert(0,str(wrong))
            try:self.assertEqual(reporting(root).SOURCE,'correct_global_finalize')
            finally:sys.path.remove(str(wrong))

    def test_scientific_group_and_pair_dimensions(self):
        from report import scientific_tables
        from common import table,read
        oxygen=np.arange(201)/40.
        metrics=[('S1',p,'') for p in ACTIVE]+[('ST',p,'') for p in ACTIVE]+[
            ('S2',ACTIVE[i],ACTIVE[j]) for i in range(14) for j in range(i+1,14)]
        groups=['growth','missegregation','missegregation','missegregation','buffering','buffering',
                'buffering','genome_doubling','growth','growth','death','death','shared_oxygen_stress','shared_oxygen_stress']
        v=np.full((5,201,2,119),.1)
        v[:,:,:,14:28]=.3
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            table(root/'inputs/model_structure/parameter_equation_map.tsv',
                  [dict(parameter=p,primary_group=g) for p,g in zip(ACTIVE,groups)])
            ranking,group,hyp=scientific_tables(root,root/'out',v,oxygen,metrics,'test',1)
            focal=read(root/'out/focal_36_pair_rankings.tsv')
            self.assertEqual(len(focal),36*2*4)
            self.assertEqual(len(hyp),2*2*3)
            self.assertTrue(all(abs(r['descriptive_parameter_mean']-(.1 if r['index']=='S1' else .3))<1e-12 for r in group))
            gaps=read(root/'out/ST_minus_S1.tsv.gz')
            self.assertEqual(len(gaps),201*2*14)
            self.assertTrue(all(abs(float(r['ST_minus_S1'])-.2)<1e-12 for r in gaps))


if __name__=='__main__':unittest.main()
