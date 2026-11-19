"""Seal tested code and inputs before any allocation. Does not submit jobs."""
from pathlib import Path
import shutil
import sys
from common import *
s=Path(sys.argv[1]).resolve();p=Path(__file__).resolve().parent
assert not (s/'source_manifest.tsv').exists()
assert all(r['passed']=='True' for r in read(s/'reports/static_tests.tsv'))
for f in p.iterdir():
 if f.is_file():shutil.copy2(f,s/'release'/f.name)
files=[]
for folder in ('inputs','runtime','release'):
 for f in sorted((s/folder).rglob('*')):
  if f.is_file():files.append(dict(path=str(f.relative_to(s)),bytes=f.stat().st_size,sha256=sha(f)))
table(s/'source_manifest.tsv',files);verify(s)
table(s/'accounting/site_limits.tsv',[dict(qos='xxlarge',max_seconds=43200,checked='2026-11-16; live sacctmgr xxlarge MaxWall 12:00:00')])
table(s/'accounting/authority.tsv',read(s/'accounting/authorization.tsv'))
table(s/'reports/case_registry.tsv',[dict(case='D19',free_mechanism=14,free_nuisance=5,known_nuisance='sigma_growth control scale',primary='A-D',reserve='E-F; one round'),dict(case='N20',free_mechanism=14,free_nuisance=6,known_nuisance='none',primary='A-D',reserve='E-F; one round')])
(s/'reports/structural_limits.md').write_text('''# Scientific limits

This is conditional finite culture recovery, not proof of structural identifiability or global optimization. Sparse discrete/reset structural identification remains NOT_PROVED. The twelve direct order-2 cases remain unresolved. sigma_kary has verified binary64 flat behavior near this truth; it remains free and its raw uniqueness is not an acceptance target. Historical conditioning counts, selected-day growth and actual experimental OLS growth are different contracts. Flow/live-dead, profiles and tumor are outside this package.
''')
print('SEALED_SOURCE',len(files),sha(s/'source_manifest.tsv'))
