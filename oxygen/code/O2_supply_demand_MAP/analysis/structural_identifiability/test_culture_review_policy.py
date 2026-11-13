"""Tests for independent recovery, budgets, idempotency and launch policy."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import time
import os
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('review', Path(__file__).with_name('culture_review.py'))
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


class PolicyTests(unittest.TestCase):
    def test_warm_repetition_is_not_independent_recovery(self):
        values = [dict(task_id=str(i), objective=4., mode='global', pure_independent=False) for i in range(4)]
        values.append(dict(task_id='local', objective=4., tested_directional_stability=True, local_improvement=0))
        self.assertFalse(module.recovered(values)['pass_gate'])
        values[0]['pure_independent'] = values[1]['pure_independent'] = True
        self.assertTrue(module.recovered(values, True)['pass_gate'])

    def test_new_lower_candidate_resets_agreement(self):
        values = [dict(task_id=str(i), objective=4., mode='global', pure_independent=True) for i in range(2)]
        values += [dict(task_id='local', objective=3.9, tested_directional_stability=True, local_improvement=.1)]
        self.assertFalse(module.recovered(values)['pass_gate'])

    def test_baseline_requires_nonimproving_stronger_check(self):
        values = [dict(task_id=str(i), objective=4., mode='global', pure_independent=True) for i in range(2)]
        values += [dict(task_id='local', objective=4., tested_directional_stability=True, local_improvement=.02)]
        self.assertTrue(module.recovered(values)['pass_gate'])
        self.assertFalse(module.recovered(values, True)['pass_gate'])

    def test_budget_uses_worst_case_until_scheduler_closure(self):
        reservations = {'global':dict(kind='global',maximum_core_hours=128),
                        'screen':dict(kind='nonfit',maximum_core_hours=1,actual_core_hours=.02)}
        self.assertEqual(module.commitment(reservations),128.02)
        reservations['global']['actual_core_hours']=8
        self.assertEqual(module.commitment(reservations),8.02)
        self.assertEqual(module.commitment(reservations,'nonfit'),.02)

    def bare_review(self, folder):
        r = module.Review.__new__(module.Review)
        r.project = r.study = Path(folder)
        r.config = dict(global_attempt_cap=74,local_attempt_cap=68,reserve_global_cap=12,reserve_local_cap=20,master_seed=1)
        r.state = dict(tasks={},points={'pm_00':dict(domain='original',fixed={'p_misseg':0},priority='P1',task_ids=[])},
            global_attempts=0,local_attempts=0,reserve_global_attempts=0,reserve_local_attempts=0,next_ordinal=1,reservations={},events=[])
        r.save=lambda:None; r.launch=lambda t:None
        return r

    def test_task_slots_are_idempotent_after_partial_launch(self):
        with tempfile.TemporaryDirectory() as folder:
            r=self.bare_review(folder)
            first=r.create('pm_00','global',pure=True,slot='independent_1')
            again=r.create('pm_00','global',pure=True,slot='independent_1')
            self.assertEqual(first,again);self.assertEqual(r.state['global_attempts'],1)
            self.assertTrue(r.state['tasks'][first]['pure_independent'])
            self.assertIsNone(r.state['tasks'][first]['source_fit'])

    def test_attempt_and_reserve_caps_do_not_impose_concurrency(self):
        with tempfile.TemporaryDirectory() as folder:
            r=self.bare_review(folder)
            for i in range(74):self.assertIsNotNone(r.create('pm_00','global',pure=True))
            self.assertIsNone(r.create('pm_00','global',pure=True))
            self.assertEqual(len(r.state['tasks']),74)
            r.state['reserve_local_attempts']=20
            self.assertIsNone(r.create('pm_00','local',reserve=True))

    def test_global_resources_respect_verified_xxlarge_wall_limit(self):
        with tempfile.TemporaryDirectory()as folder:
            r=self.bare_review(folder);r.wrapper=Path(folder)/'wrapper.sh'
            r.controller=lambda *args:None;observed={}
            def submit(key,args,charge,kind):
                observed.update(arguments=args,charge=charge,kind=kind);return '123'
            r.submit=submit
            task=dict(task_id='test',mode='global',status='ready')
            module.Review.launch(r,task)
            self.assertIn('--time=12:00:00',observed['arguments'])
            self.assertEqual(observed['charge'],96)
            self.assertIn('--cpus-per-task=8',observed['arguments'])

    def test_screening_upgrade_is_validated_and_bound_to_each_future_task(self):
        with tempfile.TemporaryDirectory()as folder:
            r=self.bare_review(folder);root=Path(folder)
            (root/'screen.R').write_text('initialization only');(root/'screen.sh').write_text('wrapper')
            module.write(root/'screen_validation.json',{'pass':True,'script_sha256':module.sha(root/'screen.R')})
            module.write(root/'manifests/screening_versions.json',{'active':dict(script_path='screen.R',script_sha256=module.sha(root/'screen.R'),wrapper_path='screen.sh',wrapper_sha256=module.sha(root/'screen.sh'),validation_path='screen_validation.json')})
            tid=r.create('pm_00','global',pure=True)
            self.assertEqual(r.state['tasks'][tid]['screening_script_sha256'],module.sha(root/'screen.R'))
            self.assertTrue(r.state['tasks'][tid]['pure_independent'])
            self.assertIsNone(r.state['tasks'][tid]['source_fit'])
            (root/'screen.R').write_text('changed source')
            with self.assertRaisesRegex(RuntimeError,'Initialization source changed'):
                r.create('pm_00','global',pure=True)

    def test_total_and_nonfit_budgets_fail_before_submission(self):
        with tempfile.TemporaryDirectory() as folder:
            r=self.bare_review(folder)
            r.state['reservations']={'used':dict(kind='global',maximum_core_hours=9999)}
            with patch.object(module.subprocess,'run') as run:
                with self.assertRaisesRegex(RuntimeError,'Total allocation'):
                    r.submit('new',[],128,'global')
                run.assert_not_called()
            r.state['reservations']={'used':dict(kind='nonfit',maximum_core_hours=120)}
            with self.assertRaisesRegex(RuntimeError,'Screening/validation'):
                r.submit('new',[],1,'nonfit')

    def test_ambiguous_intent_never_resubmits(self):
        with tempfile.TemporaryDirectory() as folder:
            r=self.bare_review(folder)
            module.write(Path(folder)/'jobs/new/submission_intent.json',{})
            with patch.object(module.subprocess,'run')as run:
                with self.assertRaisesRegex(RuntimeError,'Ambiguous'):
                    r.submit('new',[],1,'nonfit')
                run.assert_not_called()

    def graph_review(self, folder, independent_gap=0):
        r=self.bare_review(folder);r.started=time.time();r.state['points']={};r.state['status']='validation_pending'
        r.verify=lambda:None;r.accounting=lambda:{};r.refresh=lambda a:False;r.report=lambda:None
        r.launch=lambda t:t.update(status='global_submitted' if t['mode']=='global' else 'local_submitted')
        source=Path(folder)/'source.rds';source.write_bytes(b'fixture')
        r.old={'points':{}}
        pm=[0,.0001,.000327,.001069,.003497,.010454,.011437,.037402,.122315,.4]
        bb=[0,.001,.003728,.013895,.051795,.19307,.719686,1.60162,2.682696,10]
        for prefix,target,grid in [('pm','p_misseg',pm),('bb','buffer_beta',bb)]:
            for i,v in enumerate(grid):r.old['points']['%s_%02d'%(prefix,i)]=dict(target=target,value=v)
        params=dict(p_misseg=.010454,buffer_beta=1.60162,p_mis_base=.0001,k_o_mis=.001)
        def old_values(pid):
            q=4 if pid=='base_original' else 4.2
            return [dict(task_id='old',objective=q,source_fit='source.rds',params=params),
                    dict(task_id='alt',objective=q+.02,source_fit='source.rds',params=dict(params,k_o_mis=.002),parent_independent_anchor=True)]
        r.old_values=old_values
        def values(pid):
            p=r.state['points'][pid];result=old_values(p['parent_point'])
            for tid in p['task_ids']:
                t=r.state['tasks'][tid]
                if t['status']!='accepted':continue
                result.append(dict(task_id=tid,objective=result[0]['objective']+(independent_gap if t['pure_independent']and pid=='base_original'else 0),
                    mode=t['mode'],pure_independent=t['pure_independent'],source_fit='source.rds',params=params,
                    tested_directional_stability=t['mode']=='local',local_improvement=0))
            return result
        r.values=values
        module.write(Path(folder)/'tests/validation_status.json',dict(integration_pass=True,code_manifest_sha256=module.sha(source)))
        module.write(Path(folder)/'tests/controller_policy_status.json',{'pass':True,'controller_sha256':module.sha(Path(__file__).with_name('culture_review.py'))})
        r.script=Path(__file__).with_name('culture_review.py')
        manifest=Path(folder)/'manifests/code_manifest.tsv';manifest.parent.mkdir(exist_ok=True);manifest.write_bytes(source.read_bytes())
        return r

    def test_priority_graph_releases_all_ready_jobs_then_p2(self):
        with tempfile.TemporaryDirectory()as folder:
            r=self.graph_review(folder);r.advance()
            self.assertEqual(r.state['global_attempts'],18);self.assertEqual(r.state['local_attempts'],16)
            self.assertEqual(len(r.state['tasks']),34)
            self.assertEqual(sum(t['status']not in module.TERMINAL for t in r.state['tasks'].values()),34)
            for t in r.state['tasks'].values():t['status']='accepted'
            r.advance()
            self.assertEqual(r.state['global_attempts'],44);self.assertEqual(r.state['local_attempts'],42)

    def test_unrecovered_baseline_prevents_p2_release(self):
        with tempfile.TemporaryDirectory()as folder:
            r=self.graph_review(folder,independent_gap=.02);r.advance()
            for t in r.state['tasks'].values():t['status']='accepted'
            r.advance()
            self.assertFalse(r.state['baseline_gate']['pass_gate'])
            self.assertFalse(any(t['priority']=='P2'for t in r.state['tasks'].values()))
            self.assertEqual(r.state['reserve_global_attempts'],2)

    def test_scientific_hold_prevents_new_release(self):
        with tempfile.TemporaryDirectory()as folder:
            r=self.graph_review(folder);r.state['scientific_hold']='injected validation defect';r.advance()
            self.assertEqual(len(r.state['tasks']),0);self.assertEqual(r.state['status'],'validation_hold')

    def test_verified_restricted_improvement_repairs_baseline_while_globals_run(self):
        with tempfile.TemporaryDirectory()as folder:
            r=self.graph_review(folder);r.advance()
            for tid in r.state['points']['pm_05']['task_ids']:
                if r.state['tasks'][tid]['mode']=='local':r.state['tasks'][tid]['status']='accepted';break
            original=r.values
            def values(pid):
                result=original(pid)
                if pid=='pm_05':result[-1]=dict(result[-1],objective=3.9)
                return result
            r.values=values;r.advance()
            self.assertEqual(r.state['reserve_local_attempts'],1)
            self.assertEqual(r.state['reserve_global_attempts'],0)
            repairs=[t for t in r.state['tasks'].values()if t.get('slot')=='early_baseline_repair']
            self.assertEqual(len(repairs),1);self.assertEqual(repairs[0]['source_objective'],3.9)

    def test_screen_failure_recovery_does_not_wait_for_unrelated_searches(self):
        with tempfile.TemporaryDirectory()as folder:
            r=self.graph_review(folder);r.advance();root=Path(folder)
            for tid in r.state['points']['pm_09']['task_ids']:
                if r.state['tasks'][tid]['mode']=='global':r.state['tasks'][tid]['status']='screen_failed'
            (root/'screen.R').write_text('independent proposal');(root/'screen.sh').write_text('wrapper')
            module.write(root/'screen_validation.json',{'pass':True,'script_sha256':module.sha(root/'screen.R')})
            module.write(root/'manifests/screening_versions.json',{'active':dict(script_path='screen.R',script_sha256=module.sha(root/'screen.R'),wrapper_path='screen.sh',wrapper_sha256=module.sha(root/'screen.sh'),validation_path='screen_validation.json')})
            r.advance()
            self.assertEqual(r.state['reserve_global_attempts'],2)
            self.assertTrue(r.state['points']['pm_09']['screen_recovery_registered'])
            self.assertTrue(any(t['status']=='local_submitted'for t in r.state['tasks'].values()))

    def test_hold_blocks_planned_task_launch_in_refresh(self):
        with tempfile.TemporaryDirectory()as folder:
            r=self.bare_review(folder);r.create('pm_00','global',pure=True)
            r.state['scientific_hold']='injected hold'
            r.launch=lambda t:self.fail('A held study launched a new task')
            r.refresh({})

    def test_bookkeeping_error_is_not_hidden_as_an_infeasible_evaluation(self):
        with tempfile.TemporaryDirectory()as folder:
            r=self.bare_review(folder);r.started=time.time()
            tid=r.create('pm_00','global',pure=True);t=r.state['tasks'][tid]
            t.update(status='global_submitted',execution_mode='global',job_id='1')
            out=Path(folder)/'fits'/tid;module.write(out/'candidate_result.json',{})
            log=out/'evaluation_logs/pid_1.tsv';log.parent.mkdir()
            log.write_text('error\ncheckpoint_write_failed: injected\n')
            with patch.dict(os.environ,{'SLURM_JOB_ID':'2'}),patch.object(module.subprocess,'run')as run:
                r.refresh({'1':dict(State='COMPLETED',ExitCode='0:0')})
                run.assert_not_called()
            self.assertEqual(t['status'],'rejected');self.assertIn('scientific_hold',r.state)


if __name__=='__main__':
    result=unittest.main(exit=False).result
    if os.getenv('O2REVIEW_TEST_STATUS_PATH'):
        module.write(Path(os.environ['O2REVIEW_TEST_STATUS_PATH']),dict(
            **{'pass':result.wasSuccessful()},tests_run=result.testsRun,
            controller_sha256=module.sha(Path(__file__).with_name('culture_review.py')),
            test_sha256=module.sha(__file__),completed_utc=module.now()))
    raise SystemExit(0 if result.wasSuccessful()else 1)
