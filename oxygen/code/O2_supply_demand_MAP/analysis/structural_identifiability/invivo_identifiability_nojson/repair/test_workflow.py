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
        study=Path(__file__).resolve().parents[0]
        for parent in Path(__file__).resolve().parents:
            candidate=parent/'oxygen/results/structural_identifiability/20261111_invivo_baseline_repair_045318Z'
            if (candidate/'prior/formal/jobs.tsv').exists(): study=candidate;break
        if not (study/'prior/formal/jobs.tsv').exists():
            candidate=Path(__file__).resolve().parents[1]
            if (candidate/'prior/formal/jobs.tsv').exists(): study=candidate
        for name in ('manifests','logs','tables','releases/test'):
            (root/name).mkdir(parents=True,exist_ok=True)
        shutil.copytree(study/'prior',root/'prior')
        w.table(root/'manifests/budget.tsv',[dict(role='engineering_serial',attempts='2',cores='1',wall_seconds='900',cap_core_seconds='1800',phase='repair')],['role','attempts','cores','wall_seconds','cap_core_seconds','phase'])
        w.table(root/'manifests/task_plan.tsv',[dict(task='fixture',role='engineering_serial',mode='engineering_serial',policy='production',source_task='',required_check='',selector='reference',root_id='',root_family='reference',stage='engineering')],['task','role','mode','policy','source_task','required_check','selector','root_id','root_family','stage'])
        (root/'manifests/task_specs.rds').write_bytes(b'opaque metadata only')
        paths=list((root/'manifests').iterdir())+[p for p in (root/'prior').rglob('*') if p.is_file()]
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
    def test_integer_phase_aggregate_caps_and_failed_attempts(self):
        b=[dict(role='global',attempts='8',cores='8',wall_seconds='3600',cap_core_seconds='230400',phase='repair'),
           dict(role='other',attempts='20',cores='1',wall_seconds='900',cap_core_seconds='700000',phase='profile')]
        self.assertEqual(w.check_budget([],b,'global')[1],28800)
        failed=dict(role='global',state='FAILED',reserved_core_seconds='28800',actual_core_seconds='5')
        self.assertEqual(w.charge(failed),5)
        with self.assertRaises(ValueError): w.check_budget([failed]*8,b,'global')
        with self.assertRaises(ValueError): w.check_budget([dict(failed,actual_core_seconds='125000')],b,'global')
        other=dict(role='other',state='RUNNING',reserved_core_seconds='450000',actual_core_seconds='')
        with self.assertRaises(ValueError): w.check_budget([other],b,'global')
        with self.assertRaises(ValueError): w.check_budget([],b,'global',prior=0)
        self.assertEqual(w.charge(dict(other,state='COMPLETED')),450000)
    def test_role_ceiling_is_separate(self):
        b=[dict(role='check',attempts='2',cores='1',wall_seconds='900',cap_core_seconds='900',phase='repair')]
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

    def test_profile_and_extended_search_require_specific_gates(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.fixture_study(root)
            w.table(root/'tables/engineering_gate.tsv',[
                dict(key='pass',value='true'),dict(key='release',value='test'),
                dict(key='input_manifest_sha256',value=w.sha(root/'manifests/input_manifest.tsv')),
                dict(key='release_manifest_sha256',value=w.sha(root/'releases/test/manifest.tsv'))],['key','value'])
            row=dict(mode='global',stage='profile',policy='production',source_task='',required_check='')
            with self.assertRaises(FileNotFoundError): w.gate(root,row,'test',[])
            row['stage']='extended'
            with self.assertRaises(FileNotFoundError): w.gate(root,row,'test',[])
            self.assertFalse((root/'state/CURRENT').exists())

if __name__=='__main__': unittest.main()
