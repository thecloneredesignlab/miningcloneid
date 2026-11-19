from pathlib import Path
import sys
from common import table,sha
folder=Path(sys.argv[1])
table(folder/'phase_manifest.tsv',[dict(path=str(p.relative_to(folder)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(folder.rglob('*')) if p.is_file() and p.name not in ('phase_manifest.tsv','collection_console.txt')])
print('PHASE_COLLECTED')
