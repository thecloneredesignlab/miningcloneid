import unittest,io,os,shutil,subprocess
from tempfile import TemporaryDirectory
from pathlib import Path
from contextlib import redirect_stdout
from unittest.mock import patch
import workflow as w
import driver

def job(task,phase,seconds,state='RUNNING',actual=''):
    return dict(task=task,phase=phase,state=state,reserved_core_seconds=str(seconds),actual_core_seconds=actual)

class BudgetTests(unittest.TestCase):
    def test_declared_package_exact_integer_caps(self):
        self.assertEqual(sum(w.CAPS.values()),128*3600)
        planned=[job(t,p,n*wall) for t,(p,n,wall) in w.PLAN.items()]
        pos=w.position(planned)
        for phase in ('engineering','fit_A','fit_B','verification'):self.assertEqual(pos[phase]['remaining'],0)
        self.assertEqual(pos['buffer']['remaining'],7200)
    def test_phase_cannot_borrow_even_with_spare_other_phases(self):
        with self.assertRaisesRegex(ValueError,'Phase cap'):w.position([job('too_large','engineering',28801)])
    def test_failed_allocations_remain_charged(self):
        x=w.position([job('bad','engineering',3600,'FAILED','17'),job('live','engineering',3600)])['engineering']
        self.assertEqual((x['actual'],x['live'],x['remaining']),(17,3600,28800-3617))
    def test_forward_slots_never_reset_for_new_task_or_failed_attempt(self):
        with TemporaryDirectory() as name,redirect_stdout(io.StringIO()):
            s=Path(name)
            for i in range(64):w.reserve_forward(s,'new_task_'+str(i),'999','123')
            with self.assertRaisesRegex(ValueError,'64 full-forward'):w.reserve_forward(s,'retry_release2','1000','124')
            self.assertEqual(len(w.read(s/'accounting/engineering_forward_slots.tsv')),64)
    def test_undeclared_dispatch_is_rejected_without_slurm(self):
        with TemporaryDirectory() as name:
            with self.assertRaisesRegex(ValueError,'No declared allowance'):w.submit(Path(name),'global_C')
    def test_initial_record_counter_is_atomic_and_bounded(self):
        with TemporaryDirectory() as name,redirect_stdout(io.StringIO()):
            s=Path(name)
            for i in range(260):w.initial_ordinal(s,'global_A')
            self.assertEqual(w.kv(s/'fits/global_A/initial_counter.tsv')['count'],'256')
            self.assertTrue((s/'fits/global_A/initial_counter_COMPLETE').exists())
    def test_batch_uses_dispatch_release_when_slurm_relocates_script(self):
        # Mock launchers only: this fixture performs no container/model call.
        with TemporaryDirectory() as name:
            s=Path(name)/'study';release=s/'releases/r99';release.mkdir(parents=True)
            (release/'workflow.py').write_text('import sys\nassert sys.argv[1]=="verify"\n')
            w.typed(s/'dispatch/task/spec.tsv',dict(mode='prepare',release='r99'))
            spool=Path(name)/'spool';spool.mkdir();shutil.copy2(Path(__file__).with_name('batch.sh'),spool/'batch')
            bin=Path(name)/'mockbin';bin.mkdir()
            (bin/'sha256sum').write_text('#!/bin/sh\nprintf "'+w.IMAGE_SHA+'  fixture_only\\n"\n')
            (bin/'apptainer').write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "$D0_MOCK_RECORD"\n')
            for p in bin.iterdir():p.chmod(0o755)
            env=dict(os.environ,PATH=str(bin)+':'+os.environ['PATH'],SLURM_JOB_ID='fixture',SLURM_CPUS_PER_TASK='1',D0_MOCK_RECORD=str(Path(name)/'record'))
            subprocess.run(['bash',str(spool/'batch'),str(s),'task'],env=env,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            record=(Path(name)/'record').read_text();self.assertIn(str(release/'engineering.R'),record)
            self.assertEqual((s/'runtime/task/launch_exit_code.txt').read_text().strip(),'0')
    def test_extra_attempt_requires_authorization_and_preserves_original_costs(self):
        with TemporaryDirectory() as name:
            s=Path(name)
            with self.assertRaisesRegex(ValueError,'explicit extra'):w.extra_attempt_config(s,[],'prepare')
            d=s/'authorizations/attempt_extension_01';d.mkdir(parents=True)
            w.typed(d/'approval.tsv',dict(authorized=True,task='engineering_extra_01',cores=1,wall_seconds=3600,total_core_hour_cap=128,engineering_core_hour_cap=8,forward_cap=64))
            w.table(d/'manifest.tsv',[dict(path='approval.tsv',sha256=w.sha(d/'approval.tsv'))],['path','sha256'])
            rows=[dict(job('prepare','engineering',3600,'FAILED','9'),mode='prepare'),dict(job('engineering_retry','buffer',7200,'FAILED','14'),mode='prepare')]
            self.assertEqual(w.extra_attempt_config(s,rows,'prepare'),('buffer',1,3600))
            rows[0]['actual_core_seconds']='0'
            with self.assertRaisesRegex(ValueError,'costs'):w.extra_attempt_config(s,rows,'prepare')
    def test_extra_attempt_cannot_raise_forward_limit(self):
        with TemporaryDirectory() as name:
            s=Path(name);d=s/'authorizations/attempt_extension_01';d.mkdir(parents=True)
            w.typed(d/'approval.tsv',dict(authorized=True,task='engineering_extra_01',cores=1,wall_seconds=3600,total_core_hour_cap=128,engineering_core_hour_cap=8,forward_cap=65))
            w.table(d/'manifest.tsv',[dict(path='approval.tsv',sha256=w.sha(d/'approval.tsv'))],['path','sha256'])
            with self.assertRaisesRegex(ValueError,'alter approved caps'):w.extra_attempt_config(s,[],'prepare')
    def test_local_continuation_requires_new_consent_and_retains_cost_and_calls(self):
        with TemporaryDirectory() as name:
            s=Path(name)
            with self.assertRaisesRegex(ValueError,'explicit local continuation'):w.local_resume_config(s,[],'local_A')
            d=s/'authorizations/local_A_resume_01';d.mkdir(parents=True)
            values=dict(authorized=True,task='local_A_resume_01',cores=1,wall_seconds=7178,remaining_callbacks=11181,
                remaining_soft_seconds=6878,original_stage_core_seconds=7200,original_stage_callbacks=11183,
                total_core_hour_cap=128,engineering_core_hour_cap=8,forward_cap=64)
            def seal():
                w.typed(d/'approval.tsv',values)
                w.table(d/'manifest.tsv',[dict(path='approval.tsv',sha256=w.sha(d/'approval.tsv'))],['path','sha256'])
            seal();old=s/'fits/local_A';old.mkdir(parents=True)
            w.table(old/'evaluations_1.tsv',[dict(call=1),dict(call=2)],['call'])
            w.typed(old/'start_identity.tsv',dict(difference=0))
            rows=[dict(task='local_A',job_id='21902503',state='FAILED',actual_core_seconds='22')]
            self.assertEqual(w.local_resume_config(s,rows,'local_A'),('fit_A',1,7178))
            values['remaining_callbacks']=11182;seal()
            with self.assertRaisesRegex(ValueError,'cannot expand'):w.local_resume_config(s,rows,'local_A')
            values['remaining_callbacks']=11181;seal();rows[0]['actual_core_seconds']='0'
            with self.assertRaisesRegex(ValueError,'cost must remain'):w.local_resume_config(s,rows,'local_A')
    def test_missing_A_continuation_consent_holds_A_but_allows_original_B(self):
        with TemporaryDirectory() as name:
            s=Path(name);w.typed(s/'study.tsv',dict(project=name))
            w.typed(s/'tables/engineering_gate.tsv',dict(pass_gate=True,**{'pass':'true','cost_fits':'true'}))
            w.table(s/'tables/reviewed_failures.tsv',[dict(task='local_A')],['task'])
            rows=[dict(task=t,mode=t,state='COMPLETED',exit_code='0:0') for t in ('prepare','parallel','engineering_check','global_A','global_B')]
            rows.append(dict(task='local_A',mode='local_A',state='FAILED',exit_code='1:0'))
            with patch.object(w,'verify'),patch.object(w,'sha',return_value='metadata_fixture'),patch.object(w,'reconcile',return_value=rows),patch.object(driver,'accepted'),patch.object(w,'submit') as submit:
                driver.advance(s)
            submit.assert_called_once_with(s,'local_B',mode='local_B')
    def test_reviewed_plot_failure_never_becomes_successful_allocation(self):
        with TemporaryDirectory() as name:
            s=Path(name);p=s/'checks/check_A';p.mkdir(parents=True)
            result=dict(**{'pass':'true','prediction_recovery':'false','origin':'A','model_forward_calls':'2','data_sha256':'fixture','objective':'1'})
            w.typed(p/'result.tsv',result);w.typed(p/'FAILED.tsv',dict(message='unable to start device PNG'))
            w.typed(s/'tables/engineering_gate.tsv',dict(data_sha256='fixture'))
            w.table(p/'candidate_checks.tsv',[dict(source_task=t,saved=q,fresh=q,difference=0,independent=q) for t,q in [('global_A',2),('local_A_resume_01',1)]],['source_task','saved','fresh','difference','independent'])
            w.table(p/'parameter_recovery.tsv',[dict(parameter='fixture_'+str(i),estimated='TRUE' if i<14 else 'FALSE') for i in range(20)],['parameter','estimated'])
            w.table(p/'function_recovery.tsv',[dict(function_name='fixture')],['function_name'])
            w.text(p/'verified_candidate.rds','model-free fixture; not a real RDS')
            rel=s/'releases/r05';w.text(rel/'fit.R','fixture source binding')
            w.table(rel/'manifest.tsv',[dict(path='fit.R',sha256=w.sha(rel/'fit.R'))],['path','sha256'])
            spec=s/'dispatch/check_A/spec.tsv';w.typed(spec,dict(release_manifest_sha256=w.sha(rel/'manifest.tsv')))
            w.text(s/'runtime/check_A/image_sha256.txt',w.IMAGE_SHA);w.text(s/'runtime/check_A/launch_exit_code.txt','1')
            w.text(s/'logs/check_A_21904146.log','fixture X11 presentation failure')
            j=dict(task='check_A',phase='verification',job_id='21904146',state='FAILED',exit_code='1:0',actual_core_seconds='22',release='r05',dispatch_sha256=w.sha(spec))
            driver.numeric_check_evidence(s,j)
            evidence=w.kv(s/'tables/check_A_numerical_evidence.tsv')
            self.assertEqual(evidence['allocation_success'],'false')
            self.assertEqual(evidence['end_environment_guard_completed'],'false')
            with self.assertRaisesRegex(ValueError,'Allocation did not succeed'):driver.accepted(s,j)
            w.typed(p/'FAILED.tsv',dict(message='different numerical error'))
            with self.assertRaisesRegex(ValueError,'Only known'):driver.numeric_check_evidence(s,j)
            w.typed(p/'FAILED.tsv',dict(message='unable to start device PNG'));result['prediction_recovery']='true';w.typed(p/'result.tsv',result)
            with self.assertRaisesRegex(ValueError,'Only known'):driver.numeric_check_evidence(s,j)
if __name__=='__main__':unittest.main()
