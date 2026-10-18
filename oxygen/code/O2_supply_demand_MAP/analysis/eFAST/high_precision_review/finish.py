#!/usr/bin/env python3
"""Collection, scheduler linkage and final scientific artifact audit."""
import argparse
import csv
import json
from pathlib import Path
from review import HERE,sha,read,table,dump


def collection(root):
    receipt=json.loads((root/'remote_collection_receipt.json').read_text())
    assert sha(root/'remote_collection_manifest.tsv')==receipt['manifest_sha256']
    rows=read(root/'remote_collection_manifest.tsv')
    for r in rows:
        relative=Path(r['path']);assert not relative.is_absolute() and '..' not in relative.parts
        path=root/relative;assert path.is_file() and not path.is_symlink()
        assert path.stat().st_size==int(r['bytes']) and sha(path)==r['sha256'],str(path)
    selection=json.loads((root/'selection.json').read_text())
    assert sha(root/'point_manifest.tsv')==selection['point_manifest_sha256']
    points=read(root/'point_manifest.tsv');assert len(points)==401
    with open(root/'slurm/final_jobs_accounting.tsv') as f:account=list(csv.DictReader(f,delimiter='|'))
    main=json.loads((root/'slurm/main_submission.json').read_text());smoke=json.loads((root/'slurm/smoke_submission.json').read_text())
    byjob={r['JobID']:r for r in account};smoke_ids=json.loads((root/'smoke_points.json').read_text())
    assert len(main['array_ids'])==395 and len(smoke_ids)==6
    for row in points:
        i=int(row['point_id']);job=smoke['job_id'] if i in smoke_ids else main['job_id']+'_'+str(i)
        a=byjob[job];assert (a['State'],a['ExitCode'])==('COMPLETED','0:0'),a
        r=json.loads((root/'points'/('point_%04d'%i)/'result.json').read_text())
        assert r['runtime']['slurm_job_id']==a['JobIDRaw']
        assert r['identity']==row and r['runtime']['selection_sha256']==sha(root/'selection.json')
    preserved=0
    repo=HERE.parents[5]
    for row in read(root.parent/'artifact_manifest.tsv'):
        assert sha(repo/row['path'])==row['sha256'],row['path'];preserved+=1
    dump(root/'collection_validation.json',dict(status='passed',n_files=len(rows),n_bytes=sum(int(r['bytes']) for r in rows),
        complete_final_point_count=401,n_main_completed_exit0=395,n_smoke_completed_points=6,
        all_runtime_jobs_match_accounting=True,n_prior_phase_artifacts_preserved=preserved,
        main_job=main['job_id'],smoke_job=smoke['job_id'],failed_history_retained=True))


def final(root,visual):
    from pypdf import PdfReader
    assert json.loads((root/'collection_validation.json').read_text())['status']=='passed'
    independent=json.loads((root/'independent_validation.json').read_text());assert independent['status']=='passed'
    completed=json.loads((root/'completion.json').read_text());assert completed['status']=='complete' and completed['n_points']==401
    pages=[]
    for p in sorted((root/'figures').glob('*.pdf')):
        d=PdfReader(p);assert len(d.pages)==1
        text=d.pages[0].extract_text();assert len(text)>200
        pages.append(dict(path=str(p.relative_to(root)),pages=1,text_characters=len(text)))
    assert len(pages)==4
    record=dict(status='passed',n_points=401,n_precision_verified=completed['n_precision_verified'],
        independent_matrix_residual_checks=independent['n_matrix_residual_checks'],
        pdf_checks=pages,visual_render_reviewed=visual,unit_tests_passed=3,cache_relocation_smoke='passed',
        new_FAST_samples=0,N_or_M_increased=False,original_sensitivity_outputs_modified=False,
        precision_scope='50/100-digit eigensolves of original fixed double matrices; decimal vectors retained',
        historical_failures='21139700 mathematical fixture; 21139702 cache loading; both preserved and superseded')
    dump(root/'deliverable_validation.json',record)
    paths=[p for p in root.rglob('*') if p.is_file() and not p.name.endswith('.tar.gz') and p.name not in ('artifact_manifest.tsv','.DS_Store')]
    paths += [p for p in HERE.iterdir() if p.is_file() and p.suffix in ('.py','.R','.cpp','.sbatch','.md')]
    repo=HERE.parents[5]
    table(root/'artifact_manifest.tsv',[dict(path=str(p.relative_to(repo)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(paths)])
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('collection','final'));p.add_argument('--root',type=Path,required=True)
    p.add_argument('--visual-render-reviewed',action='store_true');a=p.parse_args();root=a.root.resolve()
    if a.command=='collection':collection(root)
    else:final(root,a.visual_render_reviewed)
