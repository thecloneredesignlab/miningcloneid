#!/usr/bin/env python3
"""Manifest existing eFAST trajectories for read-only local spectrum diagnostics.

No model execution, new samples, cluster jobs or mutation of scientific inputs.
The standard-library inventory runs on HPC; the exact listed files are copied
into a local ignored raw cache and checked against their original receipts.
"""
import argparse
import csv
import datetime
import json
from pathlib import Path
import shutil

from collect_resolution_nm_results import PARAMETERS, PILOT, sha, read, write

KEYS=('buffer_beta','buffer_smax','p_mis_base','mu_hp','O2_crit')
FOCUS={'seed25':('p_mis_base','mu_hp','O2_crit','buffer_beta'),
       'seed138':KEYS,'seed101':KEYS,'seed174':KEYS}


def datasets():
    for seed in FOCUS:
        yield seed,4097,4,'N4097'
        if seed in ('seed25','seed138'):yield seed,8193,4,PILOT+'/N8193_M4'
        else:yield seed,4097,8,PILOT+'/N4097_M8'


def inventory(repo,out):
    base=repo/'oxygen/results/eFAST/neighborhood10pct/pilot'
    code=repo/'oxygen/code/O2_supply_demand_MAP/analysis/eFAST'
    out.mkdir(parents=True,exist_ok=True)
    tasks=[];files={};plans={}
    for seed,n,m,relative in datasets():
        root=base/relative
        if relative not in plans:
            plan=json.loads((root/'slurm/submission_plan.json').read_text())
            for name,digest in plan['code_sha256'].items():assert sha(code/name)==digest,name
            plans[relative]=dict(path=str((root/'slurm/submission_plan.json').relative_to(repo)),
                                 sha256=sha(root/'slurm/submission_plan.json'),code_sha256=plan['code_sha256'])
        design=root/'runs'/seed/('N4097' if relative=='N4097' else 'N%d_M%d'%(n,m))
        ranges=design/'parameter_ranges.tsv';assert ranges.is_file()
        spec=read(ranges);assert [r['parameter'] for r in spec]==list(PARAMETERS)
        name=str(ranges.relative_to(repo));files[name]=dict(path=name,bytes=ranges.stat().st_size,sha256=sha(ranges))
        for repeat in range(1,6):
            parent=design/'runs'/('N%d_R%d'%(n,repeat))
            for parameter in FOCUS[seed]:
                p=PARAMETERS.index(parameter)
                part=parent/'trajectories'/('P%02d_%s'%(p+1,parameter))
                meta=json.loads((part/'metadata.json').read_text())
                receipt=json.loads((part/'outputs.tsv.gz.receipt.json').read_text())
                index=json.loads((part/'indices.receipt.json').read_text())
                assert (meta['fit_seed'],meta['N'],meta['M'],meta['replicate'],meta['trajectory_parameter_index'])==(seed,n,m,repeat,p)
                assert meta['seed']==202610080+int(seed[4:])*10+repeat-1
                assert meta['n_samples']==n and receipt['n_rows']==n*201
                assert len(meta['oxygen_pct'])==201 and all(abs(o-j*.025)<1e-12 for j,o in enumerate(meta['oxygen_pct']))
                assert receipt['metadata_sha256']==sha(part/'metadata.json')==index['metadata_sha256']
                assert receipt['samples_sha256']==meta['samples_sha256']
                assert index['outputs_sha256']==receipt['outputs_sha256']
                expected={'outputs.tsv.gz':receipt['outputs_sha256'],'samples.tsv.gz':meta['samples_sha256'],
                          'metadata.json':receipt['metadata_sha256'],'indices.npz':index['cache_sha256']}
                for filename in ('outputs.tsv.gz','samples.tsv.gz','metadata.json','indices.npz',
                                 'outputs.tsv.gz.receipt.json','indices.receipt.json'):
                    path=part/filename;assert path.is_file()
                    if filename=='outputs.tsv.gz':digest=expected[filename]
                    else:
                        digest=sha(path)
                        if filename in expected:assert digest==expected[filename]
                    name=str(path.relative_to(repo));files[name]=dict(path=name,bytes=path.stat().st_size,sha256=digest)
                tasks.append(dict(trajectory_id=len(tasks)+1,fit_seed=seed,N=n,M=m,replicate=repeat,
                    parameter=parameter,parameter_index=p,part_path=str(part.relative_to(repo)),
                    ranges_path=str(ranges.relative_to(repo)),phase_seed=meta['seed'],
                    output_bytes=(part/'outputs.tsv.gz').stat().st_size,output_sha256=receipt['outputs_sha256'],
                    samples_sha256=meta['samples_sha256'],metadata_sha256=receipt['metadata_sha256'],
                    indices_sha256=index['cache_sha256'],has_certified_recovery='recovery_manifest_sha256' in receipt))
    assert len(tasks)==190
    write(out/'trajectory_manifest.tsv',tasks)
    write(out/'raw_file_manifest.tsv',[files[k] for k in sorted(files)])
    (out/'rsync_files.txt').write_text('\n'.join(sorted(files))+'\n')
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_repo=str(repo),
        n_trajectories=len(tasks),n_source_files=len(files),total_compressed_bytes=sum(r['bytes'] for r in files.values()),
        plans=plans,collector_sha256=sha(Path(__file__)),output_hash_scope='Original output receipt SHA256; complete local files must match after transfer',
        n_parameters_in_each_design=14,diagnostic_focal_parameters=FOCUS,new_jobs_submitted=False)
    (out/'inventory.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='plans'},indent=2),flush=True)


def verify(root,cache):
    rows=read(root/'raw_file_manifest.tsv')
    total=sum(int(r['bytes']) for r in rows)
    assert shutil.disk_usage(cache).free>512*1024**2
    for i,row in enumerate(rows):
        relative=Path(row['path']);assert not relative.is_absolute() and '..' not in relative.parts
        path=cache/relative
        assert path.is_file() and not path.is_symlink() and path.stat().st_size==int(row['bytes']),row['path']
        assert sha(path)==row['sha256'],row['path']
        if (i+1)%100==0:print('Verified %d/%d source files'%(i+1,len(rows)),flush=True)
    record=dict(status='passed',checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        n_source_files=len(rows),total_bytes=total,n_trajectories=190,
        raw_manifest_sha256=sha(root/'raw_file_manifest.tsv'),full_raw_output_SHA256_checked=True,
        new_jobs_submitted=False)
    (root/'raw_collection_receipt.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('inventory','verify-local'))
    p.add_argument('--repo',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--cache',type=Path)
    a=p.parse_args()
    if a.command=='inventory':inventory(a.repo.resolve(),a.out.resolve())
    else:verify(a.out.resolve(),a.cache.resolve())


if __name__=='__main__':main()
