#!/usr/bin/env python3
"""Formal budget, wall-time, coherent-state and receipt safety checks."""
import tempfile
import unittest
import shutil
from pathlib import Path
from unittest import mock
import workflow as w

class FormalWorkflowTests(unittest.TestCase):
    def fixture_study(self,root):
        project=Path(__file__).resolve().parents[8]
        pilot=project/'oxygen/results/structural_identifiability/20261111_invivo_identifiability_202133Z'
        # Locate the real immutable debit proof without assuming a machine root.
        if not (pilot/'prior').exists() and not (pilot/'tables/jobs.tsv').exists():
            for parent in Path(__file__).resolve().parents:
                candidate=parent/'oxygen/results/structural_identifiability/20261111_invivo_identifiability_202133Z'
                if (candidate/'tables/jobs.tsv').exists(): pilot=candidate;break
        for name in ('manifests','prior','logs','tables','releases/test'):
            (root/name).mkdir(parents=True,exist_ok=True)
        for name in ('jobs.tsv','raw_accounting.txt'): shutil.copy2(pilot/'tables'/name,root/'prior'/name)
        w.table(root/'manifests/budget.tsv',[dict(role='engineering_serial',attempts='2',cores='1',wall_seconds='900',cap_core_seconds='1800')],['role','attempts','cores','wall_seconds','cap_core_seconds'])
        w.table(root/'manifests/task_plan.tsv',[dict(task='fixture',role='engineering_serial',mode='engineering_serial',policy='production',source_task='',required_check='',selector='seed')],['task','role','mode','policy','source_task','required_check','selector'])
        (root/'manifests/task_specs.rds').write_bytes(b'opaque metadata only')
        paths=list((root/'manifests').iterdir())+list((root/'prior').iterdir())
        w.table(root/'manifests/input_manifest.tsv',[dict(path=str(p.relative_to(root)),sha256=w.sha(p)) for p in paths],['path','sha256'])
        (root/'releases/test/batch.sh').write_text('exit 0\n')
        w.table(root/'releases/test/manifest.tsv',[dict(path='batch.sh',sha256=w.sha(root/'releases/test/batch.sh'))],['path','sha256'])
    def test_submission_callback_is_idempotent(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.fixture_study(root);calls=[]
            def fake(cmd,**kwargs):
                calls.append(cmd)
                if cmd[0].endswith('sbatch'): return mock.Mock(returncode=0,stdout='12345\n',stderr='')
                return mock.Mock(returncode=0,stdout='12345|RUNNING|0:0|1|1||\n',stderr='')
            with mock.patch.object(w.subprocess,'run',side_effect=fake):
                w.submit(root,'fixture','test');w.submit(root,'fixture','test')
            sbatch=[c for c in calls if c[0].endswith('sbatch')]
            self.assertEqual(len(sbatch),1);self.assertIn('--time=00:15:00',sbatch[0])
            self.assertEqual(len(w.State(root).load()),1)
    def test_lost_receipt_keeps_reservation_and_blocks_retry(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.fixture_study(root)
            with mock.patch.object(w.subprocess,'run',return_value=mock.Mock(returncode=1,stdout='',stderr='lost receipt')) as run:
                with self.assertRaises(RuntimeError): w.submit(root,'fixture','test')
                with self.assertRaises(ValueError): w.submit(root,'fixture','test')
                self.assertEqual(run.call_count,1)
            rows=w.State(root).load();self.assertEqual(len(rows),1);self.assertEqual(w.charge(rows[0]),900)
    def test_minute_resolution(self):
        self.assertEqual(w.slurm_time(900),'00:15:00')
        self.assertEqual(w.slurm_time(1800),'00:30:00')
        self.assertEqual(w.slurm_time(3600),'01:00:00')
        with self.assertRaises(ValueError): w.slurm_time(30)
    def test_integer_dual_caps_and_failed_attempts(self):
        b=[dict(role='global',attempts='8',cores='8',wall_seconds='3600',cap_core_seconds='230400')]
        self.assertEqual(w.check_budget([],b,'global')[1],28800)
        failed=dict(role='global',state='FAILED',reserved_core_seconds='28800',actual_core_seconds='5')
        self.assertEqual(w.charge(failed),5)
        with self.assertRaises(ValueError): w.check_budget([failed]*8,b,'global')
        other=dict(role='other',state='RUNNING',reserved_core_seconds='620000',actual_core_seconds='')
        with self.assertRaises(ValueError): w.check_budget([other],b,'global')
        other['reserved_core_seconds']='610000'
        with self.assertRaises(ValueError): w.check_budget([other],b,'global',prior=60000)
        closed=dict(role='other',state='COMPLETED',reserved_core_seconds='28800',actual_core_seconds='')
        self.assertEqual(w.charge(closed),28800)
    def test_role_ceiling_is_separate(self):
        b=[dict(role='check',attempts='2',cores='1',wall_seconds='900',cap_core_seconds='900')]
        j=dict(role='check',state='COMPLETED',reserved_core_seconds='900',actual_core_seconds='1')
        with self.assertRaises(ValueError): w.check_budget([j],b,'check')
    def test_coherent_state_and_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            s=w.State(Path(d))
            with s.locked(): s.save([]);self.assertEqual(s.load(),[])
            current=(Path(d)/'state/CURRENT').read_text().strip()
            (Path(d)/'state'/current/'jobs.tsv').write_text('corrupt')
            with self.assertRaises(ValueError): s.load()
    def test_uncertain_receipt_is_not_retried(self):
        rows=[dict(job_id='',comment='owned',state='INTENT')]
        with self.assertRaises(ValueError): w.resolve_matches(rows,[('12','owned'),('13','owned')])
        with self.assertRaises(ValueError): w.resolve_matches(rows,[('12','other')])
        w.resolve_matches(rows,[('12','owned'),('12','owned')]);self.assertEqual(rows[0]['job_id'],'12')
    def test_format_and_payload_rejection(self):
        item=mock.Mock();item.is_file.return_value=True;item.suffix='.json'
        with mock.patch.object(Path,'rglob',return_value=[item]):
            with self.assertRaises(ValueError): w.audit_formats(Path('/unused'))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            with self.assertRaises(ValueError): w.verify_payload(p)
            (p/'record.tsv').write_text('value\n1\n')
            w.table(p/'payload_manifest.tsv',[dict(path='record.tsv',sha256=w.sha(p/'record.tsv'))],['path','sha256'])
            (p/'COMPLETE').write_text(w.sha(p/'payload_manifest.tsv')+'\n');w.verify_payload(p)
            (p/'record.tsv').write_text('value\n2\n')
            with self.assertRaises(ValueError): w.verify_payload(p)

if __name__=='__main__': unittest.main()
