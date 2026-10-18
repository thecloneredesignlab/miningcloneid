#!/usr/bin/env python3
"""Submit the approved smoke or array once after frozen-input validation."""
import argparse
import json
from pathlib import Path
import subprocess
from boundary import HERE, REPO, preflight, sha, read, dump


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('role', choices=('smoke', 'main'))
    p.add_argument('--root', type=Path, required=True)
    a = p.parse_args(); root = a.root.resolve()
    (root / 'slurm').mkdir(exist_ok=True)
    marker = root / 'slurm' / (a.role + '_submission.json')
    if marker.exists(): print(marker.read_text()); return
    assert subprocess.check_output(['git', 'branch', '--show-current'], cwd=REPO, universal_newlines=True).strip() == 'soft_couping_eFAST'
    s = preflight(root, image=True)
    ids = []
    if a.role == 'main':
        smoke = json.loads((root / 'smoke_receipt.json').read_text())
        assert smoke['status'] == 'passed' and smoke['selection_sha256'] == sha(root / 'selection.json')
        assert json.loads((root / 'cache_portability.json').read_text())['status'] == 'passed'
        assert sha(root / 'compiled_review_template.tar.gz') == smoke['compiled_template_sha256']
        ids = [int(r['task_id']) for r in read(root / 'task_manifest.tsv') if int(r['task_id']) not in s['smoke_tasks']]
        assert len(ids) == 640
    cmd = ['sbatch', '--parsable', '--qos=xxlarge', '--time=12:00:00', '--cpus-per-task=1', '--mem=4G',
           '--job-name=efast_boundary_' + a.role, '--export=ALL,BOUNDARY_ROLE=' + a.role,
           '--output=' + str(root / 'slurm' / (a.role + '_%A_%a.out'))]
    if ids: cmd += ['--array=' + ','.join(map(str, ids))]
    cmd += [str(HERE / 'run.sbatch')]
    job = subprocess.check_output(cmd, cwd=REPO, universal_newlines=True).strip().split(';')[0]
    assert job.isdigit()
    dump(marker, dict(job_id=job, role=a.role, array_ids=ids, command=cmd,
         git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, universal_newlines=True).strip(),
         selection_sha256=sha(root / 'selection.json'), code_sha256=s['code_sha256'],
         sif_sha256=next(iter(s['source_plans'].values()))['sif_sha256'],
         qos='xxlarge', time='12:00:00', specified_node=None, array_concurrency_limit=None))
    print(marker.read_text())


if __name__ == '__main__': main()
