import ast,io,unittest
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
import workflow as w

def row(task,phase,seconds,state='RUNNING',actual=''):
    return dict(task=task,phase=phase,reserved_core_seconds=str(seconds),state=state,actual_core_seconds=str(actual))

class WorkflowTests(unittest.TestCase):
    def test_exact_whole_and_complete_chain_caps(self):
        self.assertEqual(sum(w.CAPS.values()),800*3600)
        for origin in 'ABCDEF':
            self.assertEqual(sum(c*t for k,(p,c,t) in w.PLAN.items() if k.endswith('_'+origin)),105*3600)
        w.position([row(k,p,c*t) for k,(p,c,t) in w.PLAN.items()])
    def test_failed_attempts_are_charged_with_live_reservations(self):
        z=w.position([row('bad','engineering',3600,'FAILED',17),row('live','engineering',3600)])['engineering']
        self.assertEqual((z['actual'],z['live'],z['remaining']),(17,3600,24*3600-3617))
    def test_no_phase_borrowing(self):
        with self.assertRaisesRegex(ValueError,'Phase cap'):w.position([row('over','engineering',24*3600+1)])
    def test_unknown_accounting_phase_fails_closed(self):
        with self.assertRaisesRegex(ValueError,'Unknown accounting phase'):w.position([row('bad','misspelled_phase',1)])
    def test_origin_ceiling(self):
        with self.assertRaisesRegex(ValueError,'Origin cap'):w.position([row('one_A','primary',105*3600),row('extra_A','primary',1)])
    def test_fingerprint_fail_closed(self):
        with TemporaryDirectory() as t:
            s=Path(t);w.text(s/'input','registered');w.table(s/'manifest.tsv',[dict(path='input',sha256=w.sha(s/'input'))],['path','sha256']);w.verify(s,s/'manifest.tsv');w.text(s/'input','drift')
            with self.assertRaisesRegex(ValueError,'Fingerprint mismatch'):w.verify(s,s/'manifest.tsv')
    def test_unknown_task_never_reaches_scheduler(self):
        with TemporaryDirectory() as t:
            with self.assertRaisesRegex(ValueError,'Undeclared'):w.submit(Path(t),'tumor')
    def test_initial_inputs_have_finite_atomic_counter(self):
        with TemporaryDirectory() as t,redirect_stdout(io.StringIO()):
            s=Path(t)
            for i in range(260):w.initial_ordinal(s,'global_A')
            self.assertEqual(w.kv(s/'fits/global_A/initial_counter.tsv')['count'],'256')
            self.assertTrue((s/'fits/global_A/initial_counter_COMPLETE').exists())
    def test_forward_allowance_never_resets(self):
        with TemporaryDirectory() as t,redirect_stdout(io.StringIO()):
            s=Path(t)
            for i in range(128):w.reserve_forward(s,'attempt'+str(i),'fixture','fixture')
            with self.assertRaisesRegex(ValueError,'ceiling'):w.reserve_forward(s,'newrelease','fixture','fixture')
    def test_python_has_no_forbidden_runtime_imports(self):
        for path in Path(__file__).parent.glob('*.py'):
            tree=ast.parse(path.read_text())
            for node in ast.walk(tree):
                if isinstance(node,ast.Import):self.assertFalse(any(x.name.split('.')[0] in {'json','jsonlines','jpype','js2py'} for x in node.names))
                elif isinstance(node,ast.ImportFrom):self.assertNotIn((node.module or '').split('.')[0],{'json','jsonlines','jpype','js2py'})
if __name__=='__main__':unittest.main()
