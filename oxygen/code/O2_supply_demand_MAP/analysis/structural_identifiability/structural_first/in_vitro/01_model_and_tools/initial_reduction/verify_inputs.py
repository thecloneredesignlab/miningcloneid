from pathlib import Path
import sys
from workflow import read,sha
s=Path(sys.argv[1]);settings={r['key']:r['value'] for r in read(s/'study.tsv')};project=Path(settings['project'])
for row in read(s/'inputs/manifest.tsv'):
 assert sha(s/'inputs'/row['path'])==row['sha256'],'Input changed: '+row['path']
for row in read(s/'inputs/external_sources.tsv'):
 assert sha(project/row['path'])==row['sha256'],'External source changed: '+row['path']
print('Frozen input and read-only source checks passed',flush=True)
