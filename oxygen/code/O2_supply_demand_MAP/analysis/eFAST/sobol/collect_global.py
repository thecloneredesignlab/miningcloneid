#!/usr/bin/env python3
"""Collect only completed independent global artifacts from an allowlisted inventory."""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import tempfile


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as handle:
        for block in iter(lambda:handle.read(1048576),b''): h.update(block)
    return h.hexdigest()


def remote_snapshot(host,root):
    code="""import csv,json,pathlib,hashlib,subprocess,sys
p=pathlib.Path(sys.argv[1])
s=json.loads((p/'status.json').read_text())
if s.get('status')!='complete' or s.get('scope')!='global': raise SystemExit('Independent global analysis is not complete')
with open(str(p/'slurm/submission_jobs.tsv')) as h:
    jobs=[r['job_id'] for r in csv.DictReader(h,delimiter='\\t') if r['role']=='global_postprocess']
if len(jobs)!=1: raise SystemExit('Expected one global post-processing job')
a=subprocess.check_output(['/cm/shared/apps/slurm/current/bin/sacct','-X','-n','-P','-j',jobs[0],'--format=State,ExitCode'],universal_newlines=True).strip()
if a.splitlines()[0].split('|')[:2]!=['COMPLETED','0:0']: raise SystemExit('Global job is not COMPLETED with exit 0:0: '+a)
print(json.dumps(dict(status=s,job_id=jobs[0],accounting=a,manifest_sha256=hashlib.sha256((p/'collection_manifest.tsv').read_bytes()).hexdigest())))
"""
    cmd='python3 -c '+shlex.quote(code)+' '+shlex.quote(root)
    return json.loads(subprocess.check_output(['ssh','-o','BatchMode=yes',host,cmd],universal_newlines=True))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host',default='4482173@red.moffitt.org')
    parser.add_argument('--remote-dir',default='/share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/sobol/full_N4096/global/standalone_N4096')
    parser.add_argument('--local-dir',default=str(Path(__file__).resolve().parents[6]/'oxygen/results/eFAST/sobol/full_N4096/global/standalone_N4096'))
    args=parser.parse_args(); before=remote_snapshot(args.host,args.remote_dir)
    target=Path(args.local_dir).resolve(); target.mkdir(parents=True,exist_ok=True)
    remote=args.host+':'+args.remote_dir.rstrip('/')+'/'
    subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','-e','ssh',
                    remote+'collection_manifest.tsv',str(target)+'/'],check=True)
    if sha(target/'collection_manifest.tsv')!=before['manifest_sha256']: raise ValueError('Inventory changed')
    with open(target/'collection_manifest.tsv') as handle: rows=list(csv.DictReader(handle,delimiter='\t'))
    paths=[]
    for r in rows:
        p=Path(r['path'])
        if p.is_absolute() or '..' in p.parts or 'runs' in p.parts or 'code_snapshot' in p.parts:
            raise ValueError('Unsafe collection path')
        paths.append(str(p))
    with tempfile.NamedTemporaryFile('wb') as handle:
        handle.write(('\0'.join(paths)+'\0').encode()); handle.flush()
        subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','--from0',
                        '--files-from='+handle.name,'-e','ssh',remote,str(target)+'/'],check=True)
    after=remote_snapshot(args.host,args.remote_dir)
    if after['manifest_sha256']!=before['manifest_sha256']: raise ValueError('Inventory changed during collection')
    for r in rows:
        p=target/r['path']
        if not p.is_file() or p.stat().st_size!=int(r['bytes']) or sha(p)!=r['sha256']:
            raise ValueError('Collected artifact changed: '+r['path'])
        p.chmod(0o644)
    receipt=dict(collected_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 source=args.host+':'+args.remote_dir,source_manifest_sha256=after['manifest_sha256'],
                 final_job_id=after['job_id'],final_accounting=after['accounting'],
                 n_verified_files=len(rows),verified_bytes=sum(int(r['bytes']) for r in rows),
                 status=after['status'],collector_sha256=sha(Path(__file__)))
    (target/'artifact_collection_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__': main()
