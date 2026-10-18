"""Exercise the observed container startup failure without running any model."""
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import submit_retry5734 as driver


class StartupRetryTest(unittest.TestCase):
    def run_case(self, message, failures, expected_calls, expected_code):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); (root/'count').write_text('0')
            fake = root/'fake.py'
            fake.write_text('''from pathlib import Path
import os,sys
p=Path(os.environ['TEST_COUNT']); n=int(p.read_text())+1; p.write_text(str(n))
if n<=int(os.environ['TEST_FAILURES']):
    print(os.environ['TEST_MESSAGE']); sys.exit(127)
print('payload completed')
''')
            env = dict(os.environ, TEST_COUNT=str(root/'count'), TEST_MESSAGE=message,
                       TEST_FAILURES=str(failures), SOBOL_STARTUP_RETRIES='3',
                       SOBOL_STARTUP_DELAY_SECONDS='0')
            helper = Path(__file__).parent/'retry_container_startup.sh'
            result = subprocess.run(['bash', '-c',
                'set -euo pipefail; scratch="$1"; source "$2"; launch "$3" "$4"',
                'test', str(root), str(helper), sys.executable, str(fake)],
                env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            self.assertEqual(result.returncode, expected_code, result.stdout+result.stderr)
            self.assertEqual(int((root/'count').read_text()), expected_calls)

    def test_short_user_error_recovers(self):
        self.run_case('FATAL:   user: unknown userid 107865', 2, 3, 0)

    def test_long_user_error_recovers(self):
        self.run_case("FATAL: Couldn't determine user account information: user: unknown userid 107865", 1, 2, 0)

    def test_unrelated_failure_is_not_retried(self):
        self.run_case('FATAL: image not found', 20, 1, 127)

    def test_analysis_failure_is_not_retried(self):
        self.run_case('Traceback: invalid analysis data', 20, 1, 127)

    def test_identity_failure_has_bounded_retries(self):
        self.run_case('FATAL: user: unknown userid 107865', 20, 3, 127)

    def test_initial_success(self):
        self.run_case('unused', 0, 1, 0)

    def test_minimal_tasks_and_scheduler_options(self):
        self.assertEqual([s[1] for s in driver.STAGES], ['5734','173','','','1-54',''])
        with patch.object(driver, 'sha', return_value='original-plan'):
            args = driver.arguments(Path('/tmp/repo'), driver.STAGES[0], hold=True)
            self.assertIn('--hold', args)
            for stage in driver.STAGES:
                args = driver.arguments(Path('/tmp/repo'), stage, predecessor='123')
                self.assertIn('--qos=xxlarge', args)
                self.assertIn('--time=12:00:00', args)
                self.assertIn('--dependency=afterok:123', args)
                self.assertFalse(any(x.startswith(('--nodelist', '--exclude')) for x in args))
                self.assertFalse(any('%' in x for x in args if x.startswith('--array=')))


if __name__ == '__main__':
    unittest.main()
