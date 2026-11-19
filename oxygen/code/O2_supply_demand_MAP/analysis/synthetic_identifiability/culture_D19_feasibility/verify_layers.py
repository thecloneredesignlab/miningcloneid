"""Read-only integrity check for scientific, controller and closure code layers."""
from pathlib import Path
import csv,hashlib,sys
s=Path(sys.argv[1])
for name in ('source_manifest.tsv','control_manifest.tsv','closure_manifest.tsv'):
 with (s/name).open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 for r in rows:
  p=s/r['path'];assert p.stat().st_size==int(r['bytes'])
  h=hashlib.sha256()
  with p.open('rb') as f:
   for block in iter(lambda:f.read(1048576),b''):h.update(block)
  assert h.hexdigest()==r['sha256'],str(p)
print('ALL_CODE_LAYERS_PASS')
