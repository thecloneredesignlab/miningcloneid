"""Mock Slurm/site commands; never submit or extend a real job in these tests."""
import os
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import time_extension as monitor


class MonitorTest(unittest.TestCase):
    def setUp(self):
        patcher = patch.object(monitor, 'TARGETS', ('21126822', '21129865'))
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_slurm_durations(self):
        self.assertEqual(monitor.seconds('30-12:00:00'), 2635200)
        self.assertEqual(monitor.seconds('5:59:59'), 21599)
        self.assertEqual(monitor.seconds('29:48'), 1788)

    def test_frozen_target_registry_and_invalid_scope(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(monitor, 'ROOT', Path(folder)):
            path = Path(folder)/'targets.json'
            path.write_text(json.dumps({'target_jobs': ['1001', '1002']}))
            monitor.configure()
            self.assertEqual(monitor.TARGETS, ('1001', '1002'))
            path.write_text(json.dumps({'target_jobs': ['1001', '1001']}))
            with self.assertRaises(ValueError): monitor.configure()

    def test_submission_resources_and_handover(self):
        args = monitor.sbatch_args('21126822', '1234')
        for item in ['--qos=small', '--mem=1G', '--cpus-per-task=1',
                     '--time=30-12:00:00', '--dependency=afterany:1234',
                     '--signal=B:USR1@600']:
            self.assertIn(item, args)
        self.assertFalse(any(x.startswith('--nodelist') or x.startswith('--array') for x in args))

    def test_existing_monitor_is_reused(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(monitor, 'ROOT', Path(folder)), \
                patch.object(monitor, 'queue', return_value=[['1234', 'RUNNING']]), \
                patch.object(monitor, 'command') as command, patch.object(monitor, 'event'):
            self.assertEqual(monitor.submit('21126822'), '1234')
            command.assert_not_called()

    def test_parent_numeric_id_selects_exact_internal_task(self):
        result = type('Result', (), dict(returncode=0, stderr='', stdout=(
            'JobId=21126822 ArrayJobId=21126822 ArrayTaskId=25318 JobState=RUNNING\n'
            'JobId=21126825 ArrayJobId=21126822 ArrayTaskId=27551 JobState=RUNNING\n')))()
        with patch.object(monitor, 'command', return_value=result):
            self.assertEqual(monitor.job_details('21126822')['ArrayTaskId'], '25318')

    def test_renewal_ignores_current_but_reuses_successor(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(monitor, 'ROOT', Path(folder)), \
                patch.object(monitor, 'queue', return_value=[['1234', 'RUNNING'], ['1235', 'PENDING']]), \
                patch.object(monitor, 'command') as command, patch.object(monitor, 'event'):
            self.assertEqual(monitor.submit('21126822', predecessor='1234'), '1235')
            command.assert_not_called()

    def test_extension_checks_membership_and_readback(self):
        before = dict(JobId='999', ArrayJobId='21126822', JobState='RUNNING',
                      TimeLimit='12:00:00', RunTime='6:01:00')
        after = dict(before, TimeLimit='18:00:00')
        result = type('Result', (), dict(returncode=0, stdout='accepted', stderr=''))()
        with patch.dict(os.environ, RECOVERY_TIME_TARGET='21126822'), \
                patch.object(monitor, 'job_details', side_effect=[before, after]), \
                patch.object(monitor, 'command', return_value=result) as command, \
                patch.object(monitor, 'event') as event:
            self.assertEqual(monitor.extend('999', '6'), 0)
            command.assert_called_once_with(['/app/bin/s-add-job-hours', '999', '6'])
            self.assertTrue(event.call_args[1]['verified'])

    def test_ordinary_job_is_supported(self):
        before = dict(JobId='21129865', JobState='RUNNING', TimeLimit='12:00:00', RunTime='7:00:00')
        after = dict(before, TimeLimit='18:00:00')
        result = type('Result', (), dict(returncode=0, stdout='', stderr=''))()
        with patch.dict(os.environ, RECOVERY_TIME_TARGET='21129865'), \
                patch.object(monitor, 'job_details', side_effect=[before, after]), \
                patch.object(monitor, 'command', return_value=result), patch.object(monitor, 'event'):
            self.assertEqual(monitor.extend('21129865', '6'), 0)

    def test_reject_unrelated_job(self):
        before = dict(JobId='999', ArrayJobId='12345')
        with patch.dict(os.environ, RECOVERY_TIME_TARGET='21126822'), \
                patch.object(monitor, 'job_details', return_value=before), \
                patch.object(monitor, 'command') as command:
            with self.assertRaises(ValueError):
                monitor.extend('999', '6')
            command.assert_not_called()

    def test_skip_if_threshold_has_already_been_restored(self):
        before = dict(JobId='999', ArrayJobId='21126822', JobState='RUNNING',
                      TimeLimit='18:00:00', RunTime='6:01:00')
        with patch.dict(os.environ, RECOVERY_TIME_TARGET='21126822'), \
                patch.object(monitor, 'job_details', return_value=before), \
                patch.object(monitor, 'command') as command, patch.object(monitor, 'event'):
            self.assertEqual(monitor.extend('999', '6'), 0)
            command.assert_not_called()

    def test_monitor_near_deadline_submits_successor(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(monitor, 'ROOT', Path(folder)), \
                patch.dict(os.environ, SLURM_JOB_ID='1234', SLURM_JOB_END_TIME='1'), \
                patch.object(monitor.signal, 'signal'), patch.object(monitor, 'event'), \
                patch.object(monitor, 'submit', return_value='1235') as submit:
            self.assertEqual(monitor.watch('21126822'), 0)
            submit.assert_called_once_with('21126822', predecessor='1234')

    def test_monitor_exits_when_original_loop_and_target_finish(self):
        child = type('Child', (), dict(poll=lambda self: 0, returncode=0))()
        with tempfile.TemporaryDirectory() as folder, patch.object(monitor, 'ROOT', Path(folder)), \
                patch.dict(os.environ, SLURM_JOB_ID='1234', SLURM_JOB_END_TIME='9999999999'), \
                patch.object(monitor.signal, 'signal'), patch.object(monitor, 'event'), \
                patch.object(monitor.time, 'sleep'), patch.object(monitor, 'queue', return_value=[]), \
                patch.object(monitor.subprocess, 'Popen', return_value=child), \
                patch.object(monitor, 'submit') as submit:
            self.assertEqual(monitor.watch('21126822'), 0)
            submit.assert_not_called()


if __name__ == '__main__':
    unittest.main()
