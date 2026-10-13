#!/usr/bin/env python3
"""Collect completed Sobol deliverables and verify every copied source checksum."""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import shlex
import subprocess


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as handle:
        for b in iter(lambda:handle.read(1048576),b''): h.update(b)
    return h.hexdigest()


def remote_snapshot(host,root):
    code='''import csv,hashlib,json,pathlib,subprocess,sys
p=pathlib.Path(sys.argv[1])
s=json.loads((p/'status.json').read_text())
if s.get('status')!='complete': raise SystemExit('Pilot is not complete')
with open(p/'slurm/submission_jobs.tsv') as f:
    job=[r['job_id'] for r in csv.DictReader(f,delimiter='\\t') if r['role'].startswith('finish')][-1]
account=subprocess.check_output(['/cm/shared/apps/slurm/current/bin/sacct','-X','-j',job,'-n','-P','--format=State,ExitCode'],universal_newlines=True).strip()
if account.splitlines()[0].split('|')[:2]!=['COMPLETED','0:0']: raise SystemExit('Final job is not COMPLETED with exit 0:0: '+account)
manifest=p/'collection_manifest.tsv'
print(json.dumps(dict(status=s,manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),final_job_id=job,final_accounting=account)))
'''
    command='python3 -c '+shlex.quote(code)+' '+shlex.quote(root)
    return json.loads(subprocess.check_output(['ssh','-o','BatchMode=yes',host,command],text=True))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host',default='4482173@red.moffitt.org')
    parser.add_argument('--remote-dir',default='/share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/sobol')
    parser.add_argument('--local-dir',default=str(Path(__file__).resolve().parents[6]/'oxygen/results/eFAST/sobol'))
    args=parser.parse_args(); before=remote_snapshot(args.host,args.remote_dir)
    target=Path(args.local_dir).resolve(); target.mkdir(parents=True,exist_ok=True)
    filters=['--exclude=**/runs/***','--exclude=**/recovery/***','--exclude=**/logs/***',
        '--exclude=slurm/runtime_validation/***','--exclude=inputs/fixed_o2_dominant_ploidy_201grid.tsv',
        '--exclude=**/samples.tsv.gz','--exclude=**/outputs.tsv.gz*','--exclude=**/*.lock',
        '--include=*/','--include=*.json','--include=*.tsv','--include=*.tsv.gz','--include=*.md',
        '--include=*.pdf','--include=*.png','--exclude=*']
    subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','--prune-empty-dirs','-e','ssh',*filters,
        args.host+':'+args.remote_dir.rstrip('/')+'/',str(target)+'/'],check=True)
    after=remote_snapshot(args.host,args.remote_dir)
    if before['manifest_sha256']!=after['manifest_sha256'] or sha(target/'collection_manifest.tsv')!=after['manifest_sha256']:
        raise ValueError('Source inventory changed during collection')
    verified=[]; retained_hpc=[]
    with open(target/'collection_manifest.tsv') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            relative=Path(row['path'])
            if relative.is_absolute() or '..' in relative.parts: raise ValueError('Invalid inventory path')
            p=target/relative
            if p.is_file():
                if p.stat().st_size!=int(row['bytes']) or sha(p)!=row['sha256']:
                    raise ValueError('Copied artifact checksum mismatch: '+row['path'])
                p.chmod(0o644)
                verified.append(row)
            else: retained_hpc.append(row['path'])
    receipt=dict(collected_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source=args.host+':'+args.remote_dir,source_manifest_sha256=after['manifest_sha256'],
        status=after['status'],final_job_id=after['final_job_id'],final_accounting=after['final_accounting'],
        collector_sha256=sha(Path(__file__)),
        n_verified_files=len(verified),verified_bytes=sum(int(r['bytes']) for r in verified),
        verified_files=verified,n_files_retained_on_hpc=len(retained_hpc),
        exclusions='raw sampling blocks, ordered large arrays, recovery matrices, logs and full external Figure 4 reference table remain on HPC')
    (target/'artifact_collection_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('n_verified_files','verified_bytes','n_files_retained_on_hpc','final_job_id')},indent=2))


if __name__=='__main__': main()
