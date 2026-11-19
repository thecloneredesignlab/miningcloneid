"""Safety invariants for the independent reduction audit; no scheduler calls."""
import concurrent.futures,importlib.util,tempfile,unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location('workflow',Path(__file__).with_name('workflow.py'))
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
class AuditTests(unittest.TestCase):
 def test_closed_and_live_allocation_accounting(self):
  rows=[dict(phase='time',state='FAILED',actual_core_seconds='900',reserved_core_seconds='1800'),dict(phase='time',state='RUNNING',actual_core_seconds='',reserved_core_seconds='1800')]
  x={r['phase']:r for r in w.position(rows)}['time'];self.assertEqual(x['actual']+x['committed'],2700)
 def test_phase_ceiling(self):
  with self.assertRaises(ValueError):w.position([dict(phase='time',state='SUBMITTED',actual_core_seconds='',reserved_core_seconds=str(w.CAPS['time']+1))])
 def test_caps_sum_and_call_caps(self):
  self.assertEqual(sum(w.CAPS.values()),96*3600);self.assertEqual(w.CALL_CAPS,dict(native=640,reference=2000,history=3))
 def test_concurrent_call_reservation_and_ceiling(self):
  with tempfile.TemporaryDirectory() as p:
   s=Path(p);(s/'state').mkdir()
   with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
    vals=list(pool.map(lambda n:w.reserve_call(s,'history',str(n),'test'),range(3)))
   self.assertEqual(sorted(vals),[1,2,3]);self.assertEqual(len(w.read(s/'accounting/call_ledger.tsv')),3)
   with self.assertRaises(ValueError):w.reserve_call(s,'history','over','test')
   with self.assertRaises(ValueError):w.reserve_call(s,'unknown','over','test')
 def test_closed_study_rejects_calls(self):
  with tempfile.TemporaryDirectory() as p:
   s=Path(p);(s/'COMPLETE').write_text('closed')
   with self.assertRaises(ValueError):w.reserve_call(s,'native','forbidden','test')
 def test_stale_payload_is_not_accepted(self):
  with tempfile.TemporaryDirectory() as p:
   s=Path(p);(s/'evidence.tsv').write_text('evidence')
   w.table(s/'payload_manifest.tsv',[dict(path='evidence.tsv',sha256=w.sha(s/'evidence.tsv'))],['path','sha256'])
   (s/'COMPLETE').write_text(w.sha(s/'payload_manifest.tsv'));self.assertTrue(w.payload_ok(s))
   (s/'evidence.tsv').write_text('changed')
   with self.assertRaises(ValueError):w.payload_ok(s)
 def test_time_gate_cannot_be_skipped(self):
  with tempfile.TemporaryDirectory() as p:
   with self.assertRaises((ValueError,FileNotFoundError)):w.admission(Path(p),'time')
if __name__=='__main__':unittest.main(verbosity=2)
