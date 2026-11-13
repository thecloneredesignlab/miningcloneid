#!/usr/bin/env python3
"""Pull closed review artifacts and atomically install a captured HPC state."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--project', type=Path,
    default=Path('/Users/4482173/Documents/GitHub/StructuralAndPracticalIdentifiability'))
args = parser.parse_args()
relative = Path('oxygen/results/structural_identifiability/20261110_culture_v2_optimization_review')
study = args.project.resolve() / relative
remote = '/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/' + str(relative)
host = '4482173@red.moffitt.org'
ssh = ['ssh', '-o', 'BatchMode=yes', '-o', 'IdentitiesOnly=yes',
       '-i', '/Users/4482173/.ssh/id_ed25519', '-o', 'ServerAliveInterval=15',
       '-o', 'ServerAliveCountMax=3']
script = '''from pathlib import Path
import base64,hashlib,json
from datetime import datetime,timezone
s=Path(%r)
raw=(s/'state.json').read_bytes();state=json.loads(raw);files=[]
for tid,task in state['tasks'].items():
 if task['status'] in ('accepted','failed','rejected','screen_failed'):
  for name in ('fits','populations'):
   folder=s/name/tid
   if folder.exists():
    files.extend(str(f.relative_to(s))for f in folder.rglob('*')if f.is_file())
folder=s/'report_inputs'
if folder.exists():
 files.extend(str(f.relative_to(s))for f in folder.rglob('*')if f.is_file()and not f.is_symlink())
print(json.dumps(dict(raw_state_base64=base64.b64encode(raw).decode(),
 state_sha256=hashlib.sha256(raw).hexdigest(),files=files,
 captured_utc=datetime.now(timezone.utc).isoformat())))
''' % remote
captured = json.loads(subprocess.check_output(ssh + [host, 'python3 -'],
    input=script, text=True, timeout=120))
import base64
raw = base64.b64decode(captured['raw_state_base64']); state = json.loads(raw)
if hashlib.sha256(raw).hexdigest() != captured['state_sha256']:
    raise RuntimeError('Captured state fingerprint mismatch')
common = ['rsync', '-az', '-e', shlex.join(ssh)]
subprocess.check_call(common + ['--exclude=/runtime/', '--exclude=/fits/',
    '--exclude=/populations/', '--exclude=/state.json', '--exclude=/report_inputs/',
    host + ':' + remote + '/', str(study) + '/'])
with tempfile.TemporaryDirectory(prefix='culture_review_sync_') as temp:
    listing = Path(temp) / 'closed_files.txt'
    listing.write_text(''.join(name + '\n' for name in sorted(captured['files'])))
    subprocess.check_call(common + ['--files-from=' + str(listing),
        host + ':' + remote + '/', str(study) + '/'])
for tid, task in state['tasks'].items():
    if task['status'] != 'accepted':
        continue
    owner = study / 'fits' / tid
    result = json.loads((owner / 'validated_result.json').read_text())
    for name, expected in result['files_sha256'].items():
        if hashlib.sha256((owner / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Closed accepted payload mismatch: ' + tid + '/' + name)
pending = study / 'state.sync_pending.json'; pending.write_bytes(raw)
pending.replace(study / 'state.json')
record = dict(captured_utc=captured['captured_utc'],
    completed_utc=datetime.now(timezone.utc).isoformat(),
    state_sha256=captured['state_sha256'], tasks=len(state['tasks']),
    accepted=sum(t['status']=='accepted' for t in state['tasks'].values()),
    runtime_caches_copied=False, report_input_plain_metadata_copied=True,
    direction='HPC to local; no remote state overwritten')
(study / 'manifests/state_sync_snapshot.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
