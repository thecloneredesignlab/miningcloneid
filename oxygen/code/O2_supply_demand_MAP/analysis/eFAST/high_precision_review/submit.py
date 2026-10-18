#!/usr/bin/env python3
"""Preflight frozen inputs and submit the approved smoke or point array once."""
import argparse
import json
from pathlib import Path
import subprocess
from review import HERE,sha,read,dump


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('role',choices=('smoke','main'));p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();root=a.root.resolve();repo=HERE.parents[5]
    (root/'slurm').mkdir(exist_ok=True)
    selection=json.loads((root/'selection.json').read_text());marker=root/'slurm'/(a.role+'_submission.json')
    if marker.exists():print(marker.read_text());return
    assert subprocess.check_output(['git','branch','--show-current'],cwd=repo,universal_newlines=True).strip()=='soft_couping_eFAST'
    assert sha(root/'point_manifest.tsv')==selection['point_manifest_sha256']
    assert sha(root/'relationships.tsv')==selection['relationships_sha256']
    for name,digest in selection['code_sha256'].items():assert sha(HERE/name)==digest,name
    checked={}
    for plan in selection['source_plans'].values():
        for name,digest in plan['code_sha256'].items():
            path=HERE.parent/name
            if str(path) not in checked:checked[str(path)]=sha(path)
            assert checked[str(path)]==digest,name
        if plan['sif'] not in checked:checked[plan['sif']]=sha(plan['sif'])
        assert checked[plan['sif']]==plan['sif_sha256']
        archive=plan['archives']['rcpp_template.tar.gz'];assert sha(archive['path'])==archive['sha256']
    if a.role=='main':
        receipt=json.loads((root/'smoke_receipt.json').read_text())
        assert receipt['status']=='passed' and receipt['selection_sha256']==sha(root/'selection.json')
        assert (root/'compiled_review_template.tar.gz').stat().st_size>0
        assert json.loads((root/'cache_portability.json').read_text())['status']=='passed'
        smoke=json.loads((root/'smoke_points.json').read_text())
        ids=[int(r['point_id']) for r in read(root/'point_manifest.tsv') if int(r['point_id']) not in smoke]
    else:ids=[]
    cmd=['sbatch','--parsable','--qos=xxlarge','--time=12:00:00','--cpus-per-task=1','--mem=4G',
         '--job-name=efast_hp_'+a.role,'--export=ALL,REVIEW_ROLE='+a.role,
         '--output='+str(root/'slurm'/('%s_%%A_%%a.out'%a.role))]
    if ids:cmd+=['--array='+','.join(map(str,ids))]
    cmd+=[str(HERE/'run.sbatch')]
    job=subprocess.check_output(cmd,cwd=repo,universal_newlines=True).strip().split(';')[0];assert job.isdigit()
    record=dict(job_id=job,role=a.role,array_ids=ids,command=cmd,git_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,universal_newlines=True).strip(),
        selection_sha256=sha(root/'selection.json'),code_sha256=selection['code_sha256'],qos='xxlarge',time='12:00:00',
        specified_node=None,array_concurrency_limit=None,sif_sha256=next(iter(selection['source_plans'].values()))['sif_sha256'])
    dump(marker,record);print(json.dumps(record,indent=2))


if __name__=='__main__':main()
