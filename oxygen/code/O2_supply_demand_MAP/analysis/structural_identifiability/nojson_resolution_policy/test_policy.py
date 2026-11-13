#!/usr/bin/env python3
"""Decision and migration tests using metadata; no Slurm or model execution."""
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.dont_write_bytecode=True
import policy as P
import migrate


def job(tid,role='controller',mode=None,state='RUNNING',job_id='123',phase='control',reserved='12'):
    return dict(slot=tid,task=tid,mode=mode or role,role=role,phase=phase,job_id=job_id,token='owned_'+tid,
                state=state,exit_code='0:0',cpus='1',wall_seconds='43200',reserved_ch=reserved,
                actual_ch='0.1' if state in {'COMPLETED','CANCELLED','FAILED'} else '',submitted_utc=P.L.now())


class Tests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.project=Path(self.tmp.name)
        self.study=self.project/'oxygen/results/structural_identifiability/study'
        for name in ['state','tasks','tables','logs','submissions','policy/requests','reserve/1','reserve/2']:
            (self.study/name).mkdir(parents=True)
        self.release=self.study/'policy/releases/v1';self.release.mkdir(parents=True)
        P.L.write_table(self.release/'source_manifest.tsv',[],['path','sha256'])
        self.release_patch=patch.object(P,'RELEASE',self.release);self.release_patch.start()
        self.state=P.L.State(self.study);self.data=self.state.load()

    def tearDown(self):
        self.release_patch.stop();self.tmp.cleanup()

    def test_strict_resolution_boundary(self):
        self.assertFalse(P.repair_trigger(.001,0.))
        self.assertTrue(P.repair_trigger(.0010001,0.))
        self.assertFalse(P.repair_trigger(.0009999,0.))
        self.assertFalse(P.repair_trigger(3.,3.-1e-8))
        with self.assertRaises(ValueError):P.repair_trigger(float('nan'),3.)

    def test_small_records_do_not_advance_anchor(self):
        bg=dict(pass_gate=True,**{'pass':True},repair_trigger=False)
        with patch.object(P,'gates',return_value=(bg,[],[])),patch.object(P,'save_anchor') as save:
            for _ in range(3):self.assertTrue(P.repair_baseline(self.state,self.data,self.project))
            save.assert_not_called()
        self.assertFalse((self.study/'reserve/1/decision.tsv').exists())

    def test_cumulative_gains_are_not_pairwise(self):
        fixed=3.
        gains=[fixed-q for q in [2.9996,2.9992,2.9988]]
        self.assertEqual([P.repair_trigger(fixed,fixed-g) for g in gains],[False,False,True])
        self.assertTrue(all(a-b<.001 for a,b in zip([fixed,2.9996,2.9992],[2.9996,2.9992,2.9988])))

    def test_partial_repair_does_not_install_anchor(self):
        (self.study/'reserve/1/CLOSED').write_text('consumed')
        self.data['tasks']=[dict(task_id='base_R2_B_global',mode='global',phase='reserve_2',point='base')]
        P.L.write_kv(self.study/'tasks/base_R2_B_global.tsv',self.data['tasks'][0])
        bg=dict(**{'pass':False},repair_trigger=True)
        with patch.object(P,'gates',return_value=(bg,[],[])),patch.object(P.W,'open_round'),\
             patch.object(P.W,'chain',return_value=False),patch.object(P,'save_anchor') as save:
            self.assertIsNone(P.repair_baseline(self.state,self.data,self.project));save.assert_not_called()
        self.assertFalse((self.study/'reserve/2/CLOSED').exists())

    def test_early_success_closes_consumed_round_once(self):
        (self.study/'reserve/1/CLOSED').write_text('consumed')
        P.L.write_kv(self.study/'reserve/2/decision.tsv',dict(reason='baseline_priority'))
        self.state.save(self.data)
        bg=dict(**{'pass':True},repair_trigger=False)
        with patch.object(P,'gates',return_value=(bg,[],[])),patch.object(P,'save_anchor') as save:
            self.assertTrue(P.repair_baseline(self.state,self.data,self.project))
            self.assertTrue(P.repair_baseline(self.state,self.data,self.project));self.assertEqual(save.call_count,1)
        self.assertIsNone(P.W.round_number(self.study))

    def test_exhausted_rounds_stop_deficient_reference(self):
        for n in (1,2):(self.study/'reserve'/str(n)/'CLOSED').write_text('consumed')
        with patch.object(P,'gates',return_value=(dict(**{'pass':False},repair_trigger=False),[],[])):
            self.assertFalse(P.repair_baseline(self.state,self.data,self.project))

    def test_existing_B_origin_cannot_be_recreated(self):
        (self.study/'reserve/1/CLOSED').write_text('consumed')
        row=dict(task_id='base_R2_B_global',mode='global',phase='reserve_2',point='base')
        self.data['tasks']=[row];P.L.write_kv(self.study/'tasks/base_R2_B_global.tsv',row)
        bg=dict(**{'pass':False},repair_trigger=False,repeat_pass=True,representative_pass=True,families='A')
        with patch.object(P,'gates',return_value=(bg,[],[])),patch.object(P.W,'open_round'),\
             patch.object(P.W,'chain',return_value=True),patch.object(P.W,'add') as add:
            P.repair_baseline(self.state,self.data,self.project);add.assert_not_called()
        self.assertTrue((self.study/'reserve/2/CLOSED').exists())

    def test_representative_does_not_inherit_or_overwrite_candidate_flag(self):
        candidate=dict(task_id='a',point='pm_05',origin='a',family='A',mode='local',local_stable=False,objective=3.)
        ev=dict(action='endpoint',local_evidence=True,task_id='a',objective=3.00000002,
                verified_fit='endpoint.rds',evidence_job_id='456')
        with patch.object(P,'evidence_results',return_value=[ev]):
            r=P.witnesses(self.data,self.study,[candidate])
        self.assertFalse(candidate['local_stable']);self.assertEqual(len(r),1)
        self.assertGreater(r[0]['witness_objective'],candidate['objective'])

    def test_recovery_uses_actual_witness_without_added_tolerance(self):
        r=dict(task_id='a',point='pm_05',origin='root',family='A',witness_objective=3.00105)
        with patch.object(P.L,'accepted',return_value={'objective':4.}):
            self.assertEqual(P.root_witnesses(self.data,self.study,[r],3.,'pm_05'),[])
            r['witness_objective']=3.00099
            self.assertEqual(len(P.root_witnesses(self.data,self.study,[r],3.,'pm_05')),1)

    def test_own_global_chain_is_required(self):
        r=dict(task_id='warm',point='pm_05',origin='warm',family='A',witness_objective=3.)
        with patch.object(P.L,'accepted',return_value=None):
            self.assertEqual(P.root_witnesses(self.data,self.study,[r],3.,'pm_05'),[])

    def test_descendants_count_one_origin(self):
        a=dict(task_id='a',point='base',origin='root',family='A',witness_objective=3.)
        b=dict(a,task_id='a_repeat',witness_objective=3.00001)
        roots=P.root_witnesses(self.data,self.study,[a,b],3.,'base')
        self.assertEqual(len(roots),1);self.assertFalse(P.has_families(roots))

    def test_original_frozen_gate_is_restored_after_read(self):
        sentinel=lambda *x:[]
        with patch.object(P.W,'candidates',sentinel),patch.object(P,'OLD_GATES',side_effect=RuntimeError('read failure')):
            with self.assertRaises(RuntimeError):P.legacy_gates(self.data,self.study,self.project)
            self.assertIs(P.W.candidates,sentinel)

    def endpoint_folder(self):
        folder=self.study/'fits/local';folder.mkdir(parents=True)
        fields=['mesh','complete','after','maximum_tested_improvement','improvement']
        rows=[dict(mesh=h,complete='TRUE',after=3.00000002,maximum_tested_improvement=1e-6,improvement=1e-6) for h in [.001,.0001]]
        P.L.write_table(folder/'local_sweeps.tsv',rows,fields)
        P.L.write_kv(folder/'method_sequence.tsv',dict(total_calls=10294))
        return folder,rows,fields,dict(local_stable=False,budget_limited=False,objective=3.)

    def test_endpoint_subprecision_gap_can_be_screened(self):
        folder,_,_,result=self.endpoint_folder()
        self.assertAlmostEqual(P.endpoint_evidence(folder,result),3.00000002)
        result['objective']=2.9998;self.assertIsNone(P.endpoint_evidence(folder,result))

    def test_missing_sweep_or_overall_budget_is_not_waived(self):
        folder,rows,fields,result=self.endpoint_folder()
        result['budget_limited']=True;self.assertIsNone(P.endpoint_evidence(folder,result))
        result['budget_limited']=False;P.L.write_table(folder/'local_sweeps.tsv',rows[:1],fields)
        self.assertIsNone(P.endpoint_evidence(folder,result))

    def test_incomplete_or_improving_sweep_is_not_waived(self):
        folder,rows,fields,result=self.endpoint_folder()
        rows[-1]['complete']='FALSE';P.L.write_table(folder/'local_sweeps.tsv',rows,fields)
        self.assertIsNone(P.endpoint_evidence(folder,result))
        rows[-1]['complete']='TRUE';rows[-1]['maximum_tested_improvement']=.00011
        P.L.write_table(folder/'local_sweeps.tsv',rows,fields)
        self.assertIsNone(P.endpoint_evidence(folder,result))

    def test_boundary_improvement_cannot_be_hidden(self):
        folder,_,_,result=self.endpoint_folder()
        P.L.write_table(folder/'boundary_probes.tsv',[dict(endpoint=3.00000002,probe=2.9999)],['endpoint','probe'])
        self.assertIsNone(P.endpoint_evidence(folder,result))

    def test_unknown_uid_or_failed_evidence_is_not_success(self):
        self.data['jobs']=[job('e','check',mode='policy_evidence',state='FAILED',phase='policy_evidence',reserved='.5')]
        with self.assertRaises(ValueError):P.prepare_evidence(self.state,self.data,self.project)

    def test_policy_submission_keeps_xxlarge_no_throttle_and_same_budget(self):
        self.data['jobs']=[job('B','global',job_id='999',reserved='96')]
        with patch.object(P.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'456\n','')) as call:
            r=P.policy_submit(self.state,self.data,self.project,'e','policy_validate','check','policy_validation')
        command=call.call_args[0][0]
        self.assertIn('--qos=xxlarge',command);self.assertIn('--no-requeue',command)
        self.assertFalse(any('--array' in a for a in command));self.assertEqual(r['job_id'],'456')
        self.assertEqual(P.W.counts(self.data),(1,0));self.assertEqual(P.L.committed(self.data),96.5)

    def test_six_hour_evidence_cap_and_original_total_cap(self):
        self.data['jobs']=[job('e','check',phase='policy_evidence',reserved='6')]
        with self.assertRaises(ValueError):P.policy_submit(self.state,self.data,self.project,'e2','policy_evidence','check','policy_evidence')
        self.data['jobs']=[job('large','global',reserved='4000')]
        with self.assertRaises(ValueError):P.policy_submit(self.state,self.data,self.project,'e','policy_evidence','check','policy_evidence')

    def test_ambiguous_submission_is_not_duplicated(self):
        with patch.object(P.subprocess,'run',side_effect=subprocess.TimeoutExpired('sbatch',55)) as call:
            with self.assertRaises(subprocess.TimeoutExpired):P.policy_submit(self.state,self.data,self.project,'e','policy_validate','check','policy_validation')
            P.policy_submit(self.state,self.data,self.project,'e','policy_validate','check','policy_validation')
            self.assertEqual(call.call_count,1)
        self.assertEqual(self.data['jobs'][0]['state'],'AMBIGUOUS')

    def owner_events(self):
        self.state.event(self.data,'resolution_policy_activated',str(self.release.relative_to(self.study)))
        self.state.event(self.data,'resolution_policy_fingerprint',P.L.sha(self.release/'source_manifest.tsv'))

    def test_single_writer_and_manifest_identity(self):
        self.owner_events();self.data['jobs']=[job('c')]
        with patch.dict(os.environ,SLURM_JOB_ID='123'):
            self.assertTrue(P.assert_owner(self.data,self.study))
            self.data['jobs'].append(job('other',job_id='456'))
            with self.assertRaises(ValueError):P.assert_owner(self.data,self.study)
        self.data['jobs'].pop();(self.release/'source_manifest.tsv').write_text('changed')
        with self.assertRaises(ValueError):P.assert_owner(self.data,self.study)

    def test_renewal_waits_for_predecessor_closure(self):
        self.owner_events();self.data['jobs']=[job('parent',job_id='111'),job('child',job_id='222')]
        self.state.event(self.data,'resolution_controller_renewal','111 222')
        with patch.dict(os.environ,SLURM_JOB_ID='222'):
            self.assertFalse(P.assert_owner(self.data,self.study))
            self.data['jobs'][0]['state']='COMPLETED'
            self.assertTrue(P.assert_owner(self.data,self.study))

    def test_handover_only_cancels_receipt_selected_controllers(self):
        self.data['jobs']=[job('controller',job_id='111'),job('B','global',job_id='999',reserved='96')]
        before=P.W.counts(self.data)
        with patch.object(migrate.subprocess,'run',return_value=subprocess.CompletedProcess([],0)) as call:
            ids=migrate.cancel_dispatchers(self.state,self.data)
        self.assertEqual(ids,['111']);self.assertEqual(call.call_args[0][0],[str(P.L.SLURM/'scancel'),'111'])
        self.assertEqual(P.W.counts(self.data),before);self.assertEqual(self.data['jobs'][1]['state'],'RUNNING')

    def test_historical_generation_corruption_is_rejected(self):
        self.state.save(self.data);g=(self.study/'state/CURRENT').read_text().strip()
        (self.study/'state'/g/'events.tsv').write_text('corrupt')
        with self.assertRaises(ValueError):migrate.read_generation(self.study/'state'/g)

    def gate_fixture(self):
        P.L.write_table(self.study/'cases.tsv',[dict(point='pm_05',target='p_misseg',value=.01)],['point','target','value'])
        def row(tid,point,origin,family,q,stable,alternative=False):
            return dict(task_id=tid,point=point,origin=origin,family=family,objective=q,
                        mode='local',local_stable=stable,alternative=alternative,phase='repeat',
                        start_objective=q+.00001,source_fit=tid+'.rds')
        cs=[row('base_A','base','A','A',3.,False,True),row('base_B','base','B','B',3.0005,True),
            row('point_A','pm_05','pA','A',4.,True,True),row('point_B','pm_05','pB','B',4.0005,True)]
        ev=dict(action='endpoint',local_evidence=True,task_id='base_A',objective=3.00000002,
                verified_fit='endpoint.rds',evidence_job_id='456')
        return cs,ev

    def test_full_gate_uses_representative_and_preserves_lowest_status(self):
        cs,ev=self.gate_fixture()
        with patch.object(P,'candidates',return_value=cs),patch.object(P,'evidence_results',return_value=[ev]),\
             patch.object(P,'anchor',return_value=dict(objective=3.,task_id='anchor')),\
             patch.object(P.L,'accepted',return_value={'objective':4.}):
            bg,points,_=P.gates(self.data,self.study,self.project)
        self.assertTrue(bg['pass']);self.assertTrue(points[0]['pass_gate'])
        self.assertFalse(bg['lowest_candidate_local_stable']);self.assertFalse(cs[0]['local_stable'])
        self.assertGreater(bg['representative_gap'],0.)

    def test_full_gate_rejects_outside_representative_band(self):
        cs,ev=self.gate_fixture();ev['objective']=3.0001001
        with patch.object(P,'candidates',return_value=cs),patch.object(P,'evidence_results',return_value=[ev]),\
             patch.object(P,'anchor',return_value=dict(objective=3.,task_id='anchor')),\
             patch.object(P.L,'accepted',return_value={'objective':4.}):
            bg,points,_=P.gates(self.data,self.study,self.project)
        self.assertFalse(bg['pass']);self.assertFalse(points[0]['reference_pass'])

    def quality_case(self,objective=3.0005):
        (self.study/'reserve/1/CLOSED').write_text('consumed')
        row=dict(task_id='base_R2_B_global',mode='global',phase='reserve_2',point='base')
        self.data['tasks']=[row];P.L.write_kv(self.study/'tasks/base_R2_B_global.tsv',row)
        bg=dict(**{'pass':False},repair_trigger=False,repeat_pass=True,representative_pass=True,
                families='A',objective=3.)
        candidate=dict(task_id='base_R2_B_local',point='base',mode='local',family='B',origin='base_R2_B',
                       objective=objective,local_stable=False,source_fit='source.rds')
        return bg,candidate

    def test_missing_family_near_local_evidence_is_repaired_first(self):
        bg,r=self.quality_case()
        with patch.object(P,'gates',return_value=(bg,[],[r])),patch.object(P.W,'open_round'),\
             patch.object(P.W,'chain',return_value=True),patch.object(P.W,'task',return_value={}) as task,\
             patch.object(P.W,'add') as add,patch.object(P.W,'ready',return_value=False):
            self.assertIsNone(P.repair_baseline(self.state,self.data,self.project))
            self.assertEqual(task.call_args[0][2:7],('local','reserve_2','base','base_R2_B','B'))
            self.assertEqual(add.call_count,1)
        self.assertFalse((self.study/'reserve/2/CLOSED').exists())

    def test_outside_band_does_not_trigger_quality_polishing(self):
        bg,r=self.quality_case(3.002)
        with patch.object(P,'gates',return_value=(bg,[],[r])),patch.object(P.W,'open_round'),\
             patch.object(P.W,'chain',return_value=True),patch.object(P.W,'add') as add:
            P.repair_baseline(self.state,self.data,self.project);add.assert_not_called()
        self.assertTrue((self.study/'reserve/2/CLOSED').exists())

    def test_quality_repair_cannot_exceed_shared_local_cap(self):
        bg,r=self.quality_case()
        self.data['jobs']=[job(str(i),'local',phase='reserve_2',state='COMPLETED',job_id=str(i+100),reserved='6') for i in range(8)]
        with patch.object(P,'gates',return_value=(bg,[],[r])),patch.object(P.W,'open_round'),\
             patch.object(P.W,'chain',return_value=True),patch.object(P.W,'add') as add:
            P.repair_baseline(self.state,self.data,self.project);add.assert_not_called()
        self.assertTrue((self.study/'reserve/2/CLOSED').exists());self.assertEqual(P.W.counts(self.data)[1],8)

    def test_compatible_evidence_requires_identical_R_and_controls(self):
        old=self.study/'policy/releases/old';old.mkdir()
        for name in ['evidence.R','evidence.sh','controls.tsv']:
            (old/name).write_text('unchanged');(self.release/name).write_text('unchanged')
        P.L.write_table(old/'source_manifest.tsv',[dict(path=n,sha256=P.L.sha(old/n)) for n in ['evidence.R','evidence.sh','controls.tsv']],['path','sha256'])
        old_sha=P.L.sha(old/'source_manifest.tsv')
        P.L.write_table(self.release/'evidence_compatibility.tsv',[dict(release_path=str(old.relative_to(self.study)),sha256=old_sha)],['release_path','sha256'])
        self.assertIn(old_sha,P.compatible_evidence_manifests(self.study))
        (self.release/'evidence.R').write_text('changed')
        with self.assertRaises(ValueError):P.compatible_evidence_manifests(self.study)


if __name__=='__main__':
    P.initialize(Path(sys.argv.pop(1)))
    unittest.main()
