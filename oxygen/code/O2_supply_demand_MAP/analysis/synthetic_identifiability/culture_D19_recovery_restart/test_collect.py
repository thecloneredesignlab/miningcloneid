"""Metadata-only thin-collection fixtures; no biological data or forwards."""
from pathlib import Path
import tempfile
from common import read,table,sha
from collect import select,validate

def main():
 with tempfile.TemporaryDirectory(prefix='d19_collect_fixture_') as tmp:
  repo=Path(tmp);s=repo/'study'
  for rel in ['inputs/public/fixture.txt','verification/D_A_replay2b/state.rds',
              'verification/D_A_replay2b/trace.123','verification/F_E_de2/state.rds',
              'verification/C1/audit_gate.tsv','verification/C1/phase_manifest.tsv',
              'logs/C1_123.txt','collect_release/collect.py','collect_release_v2/collect.py','controller/owner.lock','controller/SELECTION_LOCKED']:
   p=s/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('Non-scientific fixture\n')
  original=repo/'untouched.txt';original.write_text('Preserved fixture\n')
  table(s/'source_manifest.tsv',[dict(path='inputs/public/fixture.txt',bytes=(s/'inputs/public/fixture.txt').stat().st_size,sha256=sha(s/'inputs/public/fixture.txt'))])
  table(s/'accounting/preservation.tsv',[dict(path='untouched.txt',sha256=sha(original))])
  table(s/'reports/final_gate.tsv',[dict(active_owned_jobs=0)])
  inv=[]
  for rel in ['inputs/public/fixture.txt','verification/D_A_replay2b/state.rds','verification/D_A_replay2b/trace.123','verification/F_E_de2/state.rds']:
   p=s/rel;inv.append(dict(path=rel,bytes=p.stat().st_size,sha256=sha(p),authority='HPC'))
  table(s/'reports/authoritative_hpc_inventory.tsv',inv)
  try:select(s)
  except AssertionError:pass
  else:raise AssertionError('Incomplete fixture was accepted')
  (s/'controller/FINISHED').write_text('Fixture finished\n');select(s)
  records=read(s/'reports/local_collection_selection.tsv');names={r['path'] for r in records}
  assert {'verification/C1/audit_gate.tsv','verification/C1/phase_manifest.tsv','logs/C1_123.txt','collect_release/collect.py','collect_release_v2/collect.py','verification/D_A_replay2b/state.rds','controller/FINISHED','controller/SELECTION_LOCKED'}<=names
  assert not {'controller/owner.lock','verification/D_A_replay2b/trace.123','verification/F_E_de2/state.rds'}&names
  validate(s,repo)
  loc={r['path']:r['local_scope'] for r in read(s/'reports/artifact_locations.tsv')}
  assert loc['verification/C1/audit_gate.tsv']=='VERIFIED_LOCAL_COPY'
  assert loc['verification/D_A_replay2b/trace.123']=='HPC_ONLY'
  original.write_text('Changed preservation fixture\n')
  try:validate(s,repo)
  except AssertionError:pass
  else:raise AssertionError('Changed preservation fixture was accepted')
  original.write_text('Preserved fixture\n')
  (s/'inputs/public/fixture.txt').write_text('Corrupt selected fixture\n')
  try:validate(s,repo)
  except AssertionError:pass
  else:raise AssertionError('Corrupt selected file was accepted')
 print('COLLECTION_FIXTURES_PASS: completion gate, C1/log/source inclusion, raw exclusions, checksums, preservation, verified locations')

if __name__=='__main__':main()
