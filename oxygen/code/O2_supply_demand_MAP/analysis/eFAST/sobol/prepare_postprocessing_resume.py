#!/usr/bin/env python3
"""Archive the sampling plan and record the surface-only roundoff correction."""
import argparse
from pathlib import Path
import shutil
import subprocess

import core as c
import run


def main():
    p=argparse.ArgumentParser(); p.add_argument('--out-dir',required=True)
    root=Path(p.parse_args().out_dir).resolve()
    path=root/'slurm/submission_plan.json'; old=c.load_json(path)
    archive=root/'slurm/submission_plan_before_surface_roundoff.json'
    if archive.exists(): raise ValueError('Postprocessing amendment already exists; inspect before another attempt')
    for name,expected in old['input_sha256'].items():
        if c.sha(root/name)!=expected: raise ValueError('Sampling inputs changed: '+name)
    new_hashes=run.file_hashes()
    changed=[name for name,digest in old['code_sha256'].items() if new_hashes.get(name)!=digest]
    allowed={str((run.HERE/name).relative_to(run.REPO)) for name in ('report.py','test_sobol.py')}
    if set(changed)-allowed: raise ValueError('This amendment permits only the surface-grid and boundary regression-test changes')
    for case in old['cases']:
        target=c.case_directory(root,case); receipt=c.load_json(target/'completion.json')
        if receipt['status']!='complete' or c.sha(target/'indices_repeats.tsv.gz')!=receipt['indices_sha256'] or c.sha(target/'indices_summary.tsv')!=receipt['summary_sha256']:
            raise ValueError('Original case indices incomplete or changed')
    if any(root.glob('*/pilot/*/responses/*/outputs.tsv.gz')):
        raise ValueError('Existing response evaluations require a separate provenance audit')
    shutil.copyfile(path,archive)
    if (root/'status.json').exists(): shutil.copyfile(root/'status.json',root/'slurm/status_before_surface_roundoff.json')
    amended=dict(old)
    amended['sampling_code_sha256']=old['code_sha256']
    amended['code_sha256']=new_hashes
    amended['postprocessing_update']=dict(recorded_utc=c.now(),
        reason='Best fitted values differ from log-converted bound endpoints by machine roundoff; retain exact response anchors within 64 epsilons',
        changed_existing_files=changed,original_plan_sha256=c.sha(archive),
        postprocessing_git_commit=subprocess.check_output(['git','-C',str(run.REPO),'rev-parse','HEAD'],text=True).strip(),
        sampling_and_indices_recomputed=False,sobol_sampling_bounds_changed=False)
    c.save_json(path,amended)
    print('Archived original plan; retained all original samples and case indices; amended surface postprocessing provenance')


if __name__=='__main__': main()
