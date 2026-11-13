"""Pure controller fixtures: mocked Slurm, never an actual submission."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
import workflow as w

class Controls(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='ivc_controller_fixture_');self.s=Path(self.tmp.name)
        for name in ('tables','manifests','releases/cont01','dispatch','logs','state'): (self.s/name).mkdir(parents=True,exist_ok=True)
        for name in ('manifests/input_manifest.tsv','releases/cont01/manifest.tsv'):(self.s/name).write_text('path\tsha256\n')
    def tearDown(self):self.tmp.cleanup()
    def row(self,task,**values):
        j={k:'' for k in w.FIELDS};j.update(task=task,mode='local',role='baseline_local',kind='baseline',policy='production',root_id='uniform_one',root_family='uniform',release='cont01',job_id='123',state='COMPLETED',exit_code='0:0',cores='1',wall_seconds='1000000',reserved_core_seconds='1000000',actual_core_seconds='1000000',density_index='0');j.update(values);return j
    def run_submit(self,rows,**options):
        w.State(self.s).save(rows)
        with patch.object(w,'reconcile'),patch.object(w,'gate',return_value={'release':'cont01'}),patch.object(w,'verify_payload'),patch.object(w.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout='900001\n',stderr='')) as fake:
            w.submit(self.s,'next',root_id='uniform_one',root_family='uniform',**options)
            return fake
    def test_no_total_or_role_attempt_cap(self):
        self.run_submit([self.row('old_'+str(i),job_id=str(1000+i)) for i in range(100)])
        self.assertEqual(len(w.State(self.s).load()),101)
    def test_cross_root_rejected(self):
        with self.assertRaisesRegex(ValueError,'Cross-root'):self.run_submit([self.row('parent',root_id='another')],selector='task',source_task='parent')
    def test_active_baseline_blocks_profile(self):
        with self.assertRaisesRegex(ValueError,'Baseline accounting'):self.run_submit([self.row('running',state='RUNNING',actual_core_seconds='')],kind='profile')
    def test_uncertain_receipt_rejected(self):
        with self.assertRaisesRegex(ValueError,'Uncertain'):self.run_submit([self.row('unknown',job_id='',state='INTENT')])
    def test_existing_task_rejected(self):
        with self.assertRaisesRegex(ValueError,'Existing task'):self.run_submit([self.row('next')])

if __name__=='__main__':unittest.main()
