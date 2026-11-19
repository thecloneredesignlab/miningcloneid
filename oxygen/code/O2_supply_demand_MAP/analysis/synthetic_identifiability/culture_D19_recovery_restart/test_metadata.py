"""Static/metadata fixtures; no simulator, optimizer or response loading."""
from pathlib import Path
import fcntl,os,sys,tempfile
from common import *

def run():
 with tempfile.TemporaryDirectory(prefix='d19_metadata_') as d:
  p=Path(d);a=(p/'owner').open('a');b=(p/'owner').open('a');fcntl.flock(a,fcntl.LOCK_EX|fcntl.LOCK_NB)
  denied=False
  try:fcntl.flock(b,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:denied=True
  assert denied,'Second owner acquired lock'
  rows=[dict(job_id='submission_pending',task='D_A_bench')]
  assert any(not r['job_id'].isdigit() for r in rows),'Ambiguous submit was silently treated as absent'
  table(p/'rows.tsv',[dict(key='known_control',value='native RDS scalar')]);assert kv(p/'rows.tsv')['known_control']=='native RDS scalar'
  rows=[dict(path='payload',bytes=3,sha256=sha(p/'rows.tsv'))];(p/'payload').write_text('bad');table(p/'source_manifest.tsv',rows)
  rejected=False
  try:verify(p)
  except ValueError:rejected=True
  assert rejected,'Source corruption accepted'
 # No discretionary seconds cutoff, random screening, D/N dispatch or external backend.
 root=Path(__file__).resolve().parent
 for f in root.rglob('*'):
  if f.is_file():
   assert f.suffix.lower() not in {'.js','.java','.jar','.json','.jsonl'},str(f)
 print('PASS single_owner uncertain_dispatch atomic_tables stale_source excluded_artifacts')
if __name__=='__main__':run()
