"""Decision and reservation regressions; zero model calls or Slurm dispatch."""
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import policy as p


def row(task,cls='existing_local',key='root',phase='baseline',seconds=7200):
    return dict(task=task,phase=phase,slot_class=cls,slot_key=key,source_task='',
        decision_fit='',reservation_seconds=str(seconds),polish_escrow='0',
        check_escrow='900',consumed_parent='',policy_version=p.VERSION,recorded_utc='fixture')


class PolicyTests(unittest.TestCase):
    def test_thresholds_and_native_interrupt(self):
        for gain,expected in ((.000999,'resolution_stop'),(.001,'continue'),(.001001,'continue')):
            self.assertEqual(p.local_decision(gain,.002,True),expected)
        self.assertEqual(p.local_decision(0,.002001,True),'continue')
        self.assertEqual(p.local_decision(.1,0,False),'unresolved_check')
        self.assertEqual(p.local_decision(0,0,True,True),'accepted')
        self.assertEqual(p.gain_compare(14-13.998,.002),0)
    def test_cumulative_fixed_anchor(self):
        self.assertFalse(p.batch_eligible(14,13.998,13,13))
        self.assertTrue(p.batch_eligible(14,14-3*.0008,13,13))
        self.assertFalse(p.batch_eligible(14,14,13,13))
    def test_root_count_survives_component_release_failure(self):
        ledger=[row('inherit'),row('renamed_component'),row('new_release')]
        with self.assertRaisesRegex(ValueError,'allowance'):p.check_capacity([],ledger,row('fourth'))
    def test_exact_phase_edge_pending_dependents_and_closure(self):
        a=row('global','baseline_global','production:1:hybrid',seconds=57600);a['polish_escrow']='7200'
        base=row('used','audit','audit',seconds=p.CAPS['baseline']-p.CLOSURE-57600-7200-900)
        base['check_escrow']='0'
        self.assertEqual(p.check_capacity([], [base], a)['committed'],p.CAPS['baseline'])
        over=deepcopy(a);over['reservation_seconds']=str(57601)
        with self.assertRaisesRegex(ValueError,'commitment'):p.check_capacity([], [base],over)
    def test_failed_cost_and_escrow_consumed_once(self):
        a=row('global','baseline_global','production:1:hybrid',seconds=57600);a['polish_escrow']='7200'
        b=row('polish','baseline_polish',a['slot_key']);b['consumed_parent']='global/polish'
        jobs=[dict(task='global',state='FAILED',job_id='1',actual_core_seconds='800')]
        z=p.positions(jobs,[a,b])['baseline']
        self.assertEqual(z['actual'],800);self.assertEqual(z['live'],7200)
        self.assertEqual(z['escrow'],1800)
    def test_phases_do_not_borrow(self):
        a=row('base','audit','x',seconds=p.CAPS['baseline']-p.CLOSURE);a['check_escrow']='0'
        z=p.positions([], [a]);self.assertEqual(z['baseline']['committed'],p.CAPS['baseline'])
        self.assertEqual(z['profile']['committed'],p.CLOSURE)
    def test_profile_case_and_repair_do_not_reset_with_epoch(self):
        for cls,key,phase,limit in [('profile_global','production:4','profile',2),('baseline_repair','production','baseline',2),('profile_repair','data_only','profile',1)]:
            ledger=[row('old_epoch_'+str(i),cls,key,phase) for i in range(limit)]
            with self.assertRaisesRegex(ValueError,'allowance'):p.check_capacity([],ledger,row('new_epoch',cls,key,phase))
    def test_closure_escrow_is_spendable(self):
        base=row('used','audit','a',seconds=p.CAPS['baseline']-p.CLOSURE);base['check_escrow']='0'
        close=row('close','closure','c',seconds=3600);close['check_escrow']='0'
        z=p.check_capacity([], [base], close)
        self.assertEqual(z['committed'],p.CAPS['baseline']);self.assertEqual(z['closure'],3*3600)
    def test_real_reservation_rejects_relabelled_and_reused_roots(self):
        with TemporaryDirectory() as name:
            s=Path(name);(s/'tables').mkdir()
            p.w.typed(s/'tables/resolution_policy.tsv',dict(policy_version=p.VERSION,manifest='frozen',manifest_sha256='fixed'))
            spec=dict(phase='baseline',slot_class='baseline_global',slot_key='production:1:uniform',
                cores=8,wall_seconds=7200,max_calls=200000,policy='production',root_id='same',root_family='uniform',
                batch=1,branch='uniform',init_mode='uniform',source_task='')
            own=dict(task='old',root_id='same',root_family='archive',kind='baseline',policy='production',state='COMPLETED')
            with patch.object(p.w,'sha',return_value='fixed'),patch.object(p.w,'verify_manifest'),patch.object(p.w,'gate',return_value={'schema':'2'}):
                with self.assertRaisesRegex(ValueError,'lineage family'):p.reserve(s,[own],'new',spec)
                own['root_family']='uniform'
                with self.assertRaisesRegex(ValueError,'genuinely new'):p.reserve(s,[own],'new',spec)
                spec.update(branch='hybrid',init_mode='hybrid',slot_key='production:1:hybrid')
                with self.assertRaisesRegex(ValueError,'existing own uniform'):p.reserve(s,[],'new',spec)
            self.assertFalse((s/'tables/resolution_assignments.tsv').exists())
    def test_finite_terminal_requires_closed_parents_and_explicit_outcomes(self):
        import driver
        with TemporaryDirectory() as name:
            s=Path(name);(s/'tables').mkdir();path=s/'tables/strategy_review.tsv'
            def outcomes(a,b):p.w.table(path,[dict(policy='production',decision=a),dict(policy='data_only',decision=b)],['policy','decision'])
            outcomes('accepted','continue_local');self.assertFalse(driver.baseline_terminal(s,[]))
            outcomes('accepted','resolution_unresolved')
            self.assertFalse(driver.baseline_terminal(s,[dict(state='RUNNING')]))
            self.assertTrue(driver.baseline_terminal(s,[dict(state='COMPLETED')]))
            self.assertEqual(p.w.kv(s/'tables/resolution_baseline_terminal.tsv')['next_stage'],'profiles_not_released')
            outcomes('accepted','accepted');self.assertTrue(driver.baseline_terminal(s,[]))
            self.assertEqual(p.w.kv(s/'tables/resolution_baseline_terminal.tsv')['next_stage'],'conditional_profile_required')
    def test_free_reference_reserves_its_required_full_check(self):
        with TemporaryDirectory() as name:
            s=Path(name);(s/'tables').mkdir()
            p.w.typed(s/'tables/resolution_policy.tsv',dict(policy_version=p.VERSION,manifest='frozen',manifest_sha256='fixed'))
            p.w.typed(s/'tables/resolution_review_gate.tsv',dict(evidence_manifest='proof'))
            p.w.table(s/'tables/resolution_target_replays.tsv',[dict(policy='production',selection='target',root_id='informed')],['policy','selection','root_id'])
            spec=dict(phase='baseline',slot_class='reference',slot_key='target:production:hash',kind='baseline',
                cores=1,wall_seconds=900,max_calls=1,policy='production',root_id='informed',root_family='informed',
                reference_purpose='target_replay',selection_name='target',source_task='')
            with patch.object(p.w,'sha',return_value='fixed'),patch.object(p.w,'verify_manifest'):
                p.reserve(s,[],'q_reference',spec)
            a=p.assignments(s);self.assertEqual(a[0]['check_escrow'],'900')
            self.assertEqual(p.positions([],a)['baseline']['committed'],p.CLOSURE+1800)
    def test_pending_commitment_holds_without_premature_budget_endpoint(self):
        import driver
        with TemporaryDirectory() as name:
            s=Path(name);(s/'tables').mkdir()
            p.w.table(s/'tables/resolution_assignments.tsv',[row('active')],p.FIELDS)
            self.assertEqual(driver.resource_hold(s,[dict(task='active',state='RUNNING')],
                'production','fixture','baseline','batch'),'wait')
            self.assertFalse((s/'tables/resolution_budget_stops.tsv').exists())
            self.assertEqual(driver.resource_hold(s,[dict(task='active',state='COMPLETED')],
                'production','fixture','baseline','batch'),'stop')
            self.assertEqual(p.w.read_table(s/'tables/resolution_budget_stops.tsv')[0]['phase'],'baseline')
    def test_stale_optional_check_defers_to_same_root_checked_descendant(self):
        import driver
        with TemporaryDirectory() as name:
            s=Path(name);(s/'tests'/'new_check').mkdir(parents=True)
            p.w.typed(s/'tests/new_check/result.tsv',{'pass':True,'scope':'fixed_center_all_four_scales_axis_coupled_observation','complete':True,'representative_gap':.00005,'directional_improvement':.00005,'best_objective':13.9})
            old=row('old');new=row('new');ledger=[old,new]
            by={t:dict(task=t,root_id='root',policy='production',density_index='0',state='COMPLETED',exit_code='0:0',mode='local' if t!='new_check' else 'resolution_check') for t in ('old','new','new_check')}
            proof={'new':[by['new_check']]}
            self.assertTrue(driver.defer_old_representative(s,old,ledger,by,proof,14))
            self.assertFalse(driver.defer_old_representative(s,old,ledger,by,proof,13.8))
            by['new']['root_id']='other'
            self.assertFalse(driver.defer_old_representative(s,old,ledger,by,proof,14))
    def test_unreviewed_failure_holds_even_without_a_new_submission(self):
        import driver
        with TemporaryDirectory() as name:
            s=Path(name);(s/'tables').mkdir()
            jobs=[dict(task='bad',state='FAILED',exit_code='1:0')]
            with self.assertRaisesRegex(ValueError,'Unreviewed closed'):driver.assert_reviewed(s,jobs)
            p.w.table(s/'tables/reviewed_failures.tsv',[dict(task='bad')],['task'])
            driver.assert_reviewed(s,jobs)
    def test_pending_complete_check_blocks_stale_baseline_strategy(self):
        import driver
        with TemporaryDirectory() as name:
            # No current review packet is available. In particular a stale
            # finite-repair decision must not be consulted during this gap.
            jobs=[dict(task='new_check',mode='resolution_check',state='RUNNING')]
            with patch.object(driver,'submit') as dispatch:
                driver.baseline(Path(name),jobs)
                dispatch.assert_not_called()
    def test_pending_full_checks_coalesce_the_allocated_audit_wave(self):
        import driver
        with TemporaryDirectory() as name:
            jobs=[dict(task='pending_complete_check',mode='resolution_check',state='RUNNING')]
            with patch.object(driver.w,'submit') as dispatch:
                self.assertFalse(driver.audit(Path(name),jobs))
                dispatch.assert_not_called()
    def test_later_repair_cannot_retroactively_earn_second_global_batch(self):
        with TemporaryDirectory() as name:
            s=Path(name);(s/'tables').mkdir();collection=s/'collections/c_fixture';collection.mkdir(parents=True)
            p.w.table(collection/'jobs.tsv',[dict(task='closed_batch',state='COMPLETED')],['task','state'])
            p.w.table(collection/'manifest.tsv',[dict(path='jobs.tsv',sha256=p.w.sha(collection/'jobs.tsv'))],['path','sha256'])
            p.w.table(s/'tables/resolution_batches.tsv',[dict(policy='production',batch=1,before_uniform=14,before_feasible=13)],['policy','batch','before_uniform','before_feasible'])
            p.w.table(s/'tables/frozen_evidence.tsv',[dict(path='collections/c_fixture/manifest.tsv',sha256=p.w.sha(collection/'manifest.tsv'))],['path','sha256'])
            frozen=dict(policy='production',batch=1,before_uniform=14,before_feasible=13,after_uniform=13.9981,after_feasible=13,
                protected_gain=.0019,eligible='false',collection='c_fixture',collection_manifest_sha256=p.w.sha(collection/'manifest.tsv'),evidence_manifest='tables/frozen_evidence.tsv')
            path=s/'tables/resolution_batch_outcomes.tsv';p.w.table(path,[frozen],list(frozen))
            self.assertFalse(p.completed_batch_eligible(s,'production'))
            p.w.table(s/'tables/baseline_policy_gates.tsv',[dict(policy='production',best_objective=12.99)],['policy','best_objective'])
            self.assertFalse(p.completed_batch_eligible(s,'production'))
            frozen['eligible']='true';p.w.table(path,[frozen],list(frozen))
            with self.assertRaisesRegex(ValueError,'Inconsistent frozen'):p.completed_batch_eligible(s,'production')
    def test_second_repair_requires_complete_material_progress(self):
        with TemporaryDirectory() as name:
            s=Path(name);(s/'tables').mkdir();old=[row('repair1','baseline_repair','production')]
            jobs=[dict(task='repair1',state='COMPLETED',exit_code='0:0')]
            evidence=dict(task='repair1',complete='true',decision='resolution_stop',window_gain=.0009,full_check_gain=0)
            path=s/'tables/resolution_local_decisions.tsv'
            p.w.table(path,[evidence],list(evidence));self.assertFalse(p.repair_continuation_allowed(s,jobs,old))
            evidence.update(decision='continue',window_gain=.001)
            p.w.table(path,[evidence],list(evidence));self.assertTrue(p.repair_continuation_allowed(s,jobs,old))
            evidence.update(complete='false');p.w.table(path,[evidence],list(evidence))
            self.assertFalse(p.repair_continuation_allowed(s,jobs,old))
            jobs[0].update(state='FAILED',exit_code='1:0');self.assertFalse(p.repair_continuation_allowed(s,jobs,old))


if __name__=='__main__':unittest.main()
