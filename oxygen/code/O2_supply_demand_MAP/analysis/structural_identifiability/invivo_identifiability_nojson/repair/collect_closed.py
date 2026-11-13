#!/usr/bin/env python3
"""Collect only closed owned payloads and control evidence; no active fit copying."""
import argparse
import re
import shlex
import subprocess
from pathlib import Path
import workflow as w

def main():
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);a=p.parse_args();local=a.study.resolve()
    if not re.fullmatch('20261111_invivo_baseline_repair_045318Z',local.name): raise ValueError('Unexpected study identity')
    remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+local.name
    ssh=['ssh','-o','BatchMode=yes','-o','IdentitiesOnly=yes','-i','/Users/4482173/.ssh/id_ed25519','-o','ControlMaster=auto','-o','ControlPath=/private/tmp/ivi_repair_ssh_%C','-o','ControlPersist=600']
    target='4482173@red.moffitt.org'
    script="""from pathlib import Path
import sys
s=Path(REMOTE);sys.path.insert(0,str(s/'implementation'));import workflow
rows=workflow.State(s).load()
generation=(s/'state/CURRENT').read_text().strip()
files=['STATUS.md','tables/jobs.tsv','tables/raw_accounting.txt','state/'+generation+'/']
if (s/'tests/start_identity_snapshots').exists(): files.append('tests/start_identity_snapshots/')
for j in rows:
 if not workflow.closed(j): continue
 files.append('dispatch/'+j['task']+'/')
 if j['state']=='COMPLETED' and j['exit_code']=='0:0':
  scope='fits' if j['mode'] in ('global','local') else 'tests'
  folder=s/scope/j['task'];workflow.verify_payload(folder);files.append(scope+'/'+j['task']+'/')
  if j['mode']=='engineering_parallel':
   for child in sorted((s/'tests').glob(j['task']+'_*')):
    if child.is_dir():
     workflow.verify_payload(child);files.append('tests/'+child.name+'/')
 for name in (j['task']+'_submission.txt',j['task']+'_receipt.txt',j['task']+'_'+j['job_id']+'.log'):
  if (s/'logs'/name).exists(): files.append('logs/'+name)
 for name in ('image_sha256.txt','launch_exit_code.txt'):
  path='runtime/'+j['task']+'/'+name
  if (s/path).exists(): files.append(path)
print('GEN\\t'+generation)
print('\\n'.join(files))
""".replace('REMOTE',repr(remote))
    result=subprocess.run(ssh+[target,'python3','-'],input=script,check=True,universal_newlines=True,stdout=subprocess.PIPE)
    lines=result.stdout.splitlines()
    marker=lines.pop(0).split('\t')
    if len(marker)!=2 or marker[0]!='GEN' or not re.fullmatch('g_[0-9a-f]+',marker[1]): raise ValueError('Invalid generation')
    generation=marker[1];paths=lines
    if not paths or any(path.startswith('/') or '..' in Path(path).parts for path in paths): raise ValueError('Invalid collection path')
    listing=local/'tables/collection_files.txt';listing.write_text('\n'.join(paths)+'\n')
    subprocess.run(['rsync','-az','--recursive','--files-from='+str(listing),'-e',' '.join(shlex.quote(v) for v in ssh),target+':'+remote+'/',str(local)+'/'],check=True)
    state_root=local/'state'/generation
    w.verify_manifest(state_root,state_root/'manifest.tsv')
    if w.read_table(state_root/'jobs.tsv')!=w.read_table(local/'tables/jobs.tsv'): raise ValueError('Collection crossed control generations')
    w.atomic_text(local/'state/CURRENT',generation+'\n')
    count=0
    for path in paths:
        if path.endswith('/') and (local/path/'COMPLETE').exists():
            w.verify_payload(local/path);count+=1
    print('Collected and verified',count,'closed owned payload directories; active fits were not transferred.')

if __name__=='__main__': main()
