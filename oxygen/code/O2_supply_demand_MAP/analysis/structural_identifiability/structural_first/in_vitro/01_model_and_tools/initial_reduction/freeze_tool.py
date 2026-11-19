"""Preserve a verified node-local depot without shared-filesystem mode changes."""
from pathlib import Path
import csv,hashlib,sys,tarfile
s=Path(sys.argv[1]);depot=Path(sys.argv[2]);task=sys.argv[3]
forbidden={'.json','.jsonl','.ndjson','.js','.java','.jar','.class'}
removed=[]
for root in [depot,s/'toolchain']:
 for p in sorted(root.rglob('*')):
  if p.is_file() and p.suffix.lower() in forbidden:
   removed.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),reason='Unused upstream documentation or test fixture; prohibited artifact format'));p.unlink()
with (s/'results'/task/'removed_upstream_fixtures.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,['path','sha256','reason'],delimiter='\t');w.writeheader();w.writerows(removed)
archive=s/'toolchain/verified_depot.tar.gz';assert not archive.exists()
with tarfile.open(archive,'w:gz') as tar:
 tar.add(depot,arcname='depot',filter=lambda x:None if '/compiled/' in x.name or x.name.endswith('/compiled') else x)
with tarfile.open(archive) as tar:assert not any(Path(x.name).suffix.lower() in forbidden for x in tar.getmembers() if x.isfile())
digest=hashlib.sha256(archive.read_bytes()).hexdigest();(s/'toolchain/verified_depot.sha256').write_text(digest+'\n')
print('Frozen verified node-local depot',digest,'removed unused fixtures',len(removed),flush=True)
