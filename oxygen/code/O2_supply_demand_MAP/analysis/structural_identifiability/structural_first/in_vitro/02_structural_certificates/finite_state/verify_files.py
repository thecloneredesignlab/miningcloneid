"""Verify an admitted release without importing a scientific library."""
from pathlib import Path
import sys
from protocol import read,sha

study=Path(sys.argv[1]);release=sys.argv[2]
rows=read(study/('source_manifest_'+release+'.tsv'))
for r in rows:
    p=study/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],r['path']
print('ADMITTED_SOURCE_IDENTITIES_PASSED',len(rows),flush=True)
