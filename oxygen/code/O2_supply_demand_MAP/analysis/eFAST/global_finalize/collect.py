#!/usr/bin/env python3
"""Collect only completed global exports, with before/after SHA256 verification."""
import argparse
from pathlib import Path
import json
import os
import shlex
import shutil
import subprocess
import tempfile

from common import load,save,sha,now


def inventory(root):
    root=Path(root)
    if load(root/'completion.json')['status']!='complete':
        raise ValueError('Final global comparison has not completed')
    candidates=set()
    for name in ['plan.json','preparation.json','smoke.json','launch_receipt.json','completion.json',
                 'methods_results_en.md','meeting_update.tsv','export_manifest.json','slurm/jobs.json',
                 'time_extension/targets.json','time_extension/submission_jobs.tsv']:
        candidates.add(root/name)
    for name in ['hpc_sync_receipt.json','preparation_recovery_receipt.json','reporting_recovery_receipt.json','slurm/reporting_recovery_jobs.json',
                 'time_extension_reporting/targets.json','time_extension_reporting/submission_jobs.tsv']:
        if (root/name).is_file():candidates.add(root/name)
    for directory in (root/'slurm').glob('reporting_overlay_*'):
        candidates.update(p for p in directory.iterdir() if p.is_file())
    for name in ['sobol','comparison','validation','slurm/code_snapshot']:
        candidates.update(p for p in (root/name).rglob('*') if p.is_file())
    candidates.update(p for p in (root/'inputs').rglob('*') if p.is_file() and 'designs' not in p.relative_to(root/'inputs').parts)
    for directory in (root/'efast').glob('N*'):
        candidates.update(p for p in directory.iterdir() if p.is_file())
        candidates.update(p for p in (directory/'figures').rglob('*') if p.is_file())
    candidates.update(p for p in (root/'slurm/logs').glob('*') if p.is_file() and p.name.startswith(('prepare_','smoke_','sobol_','report_','final_')))
    rows=[]
    for path in sorted(candidates):
        rel=path.relative_to(root)
        if '__pycache__' in rel.parts or path.suffix=='.pyc':
            continue
        if any(x in rel.parts for x in ['blocks','analysis','neighborhood10pct','designs']):
            raise ValueError('Raw/cache/neighborhood path in export selection')
        rows.append(dict(path=str(rel),size=path.stat().st_size,sha256=sha(path)))
    return dict(source_root=str(root),files=rows,n_files=len(rows),total_bytes=sum(r['size'] for r in rows),
                model_calls=0,neighborhood_files_selected=0)


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--root',type=Path,required=True)
    a.add_argument('--inventory',action='store_true')
    a.add_argument('--remote',default='4482173@red.moffitt.org')
    a.add_argument('--remote-root',type=Path,default=Path('/share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/global_finalize_20261112'))
    a.add_argument('--remote-script',type=Path,help='Explicit newer read-only collector; does not change frozen computation')
    args=a.parse_args()
    if args.inventory:
        print(json.dumps(inventory(args.root)));return
    remote_script=args.remote_script or args.remote_root/'slurm/code_snapshot/collect.py'
    command='bash -lc '+shlex.quote('python3 '+shlex.quote(str(remote_script))+' --inventory --root '+shlex.quote(str(args.remote_root)))
    def remote_inventory():
        return json.loads(subprocess.check_output(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20',args.remote,command],universal_newlines=True))
    before=remote_inventory()
    args.root.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='global_efast_collect_') as temporary:
        stage=Path(temporary)/'stage';stage.mkdir()
        filelist=Path(temporary)/'files.txt'
        filelist.write_text(''.join(r['path']+'\n' for r in before['files']))
        subprocess.run(['rsync','-az','--no-perms','--no-owner','--no-group','--files-from='+str(filelist),
                        args.remote+':'+str(args.remote_root)+'/',str(stage)+'/'],check=True)
        after=remote_inventory()
        if before!=after:
            raise ValueError('Remote exports changed during collection; keep local originals')
        for row in before['files']:
            path=stage/row['path'];dest=args.root/row['path']
            if path.stat().st_size!=row['size'] or sha(path)!=row['sha256']:
                raise ValueError('Transferred file checksum differs')
            if dest.exists() and sha(dest)!=row['sha256']:
                raise ValueError('Refusing to replace a differing local export: '+str(dest))
        for row in before['files']:
            dest=args.root/row['path']
            if not dest.exists():
                dest.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(str(stage/row['path']),str(dest))
            if sha(dest)!=row['sha256']:raise ValueError('Local readback differs')
    save(args.root/'collection_receipt.json',dict(status='complete',n_files=before['n_files'],total_bytes=before['total_bytes'],
         all_file_hashes_verified=True,source_unchanged_during_transfer=True,model_calls=0,neighborhood_files_selected=0,
         source_root=str(args.remote_root),completed_utc=now()))
    save(args.root/'collection_manifest.json',before)
    print('Collected and verified',before['n_files'],'global exports')


if __name__=='__main__':main()
