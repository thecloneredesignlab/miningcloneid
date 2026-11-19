from pathlib import Path
import sys
from common import verify,read,sha
s=Path(sys.argv[1]);verify(s)
for r in read(s/'repair_manifest.tsv'):
 p=s/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==int(r['bytes'])
print('BASE_AND_REPAIR_SOURCE_PASS')
