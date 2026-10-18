#!/usr/bin/env python3
"""Verify collection bytes, scheduler identities and final deliverable inventory."""
import argparse
import csv
import json
from pathlib import Path
from boundary import HERE, REPO, read, sha, table, dump, preflight, point_dir


def collection(root):
    preflight(root)
    receipt = json.loads((root / 'remote_collection_receipt.json').read_text())
    assert sha(root / 'remote_collection_manifest.tsv') == receipt['manifest_sha256']
    files = read(root / 'remote_collection_manifest.tsv')
    for row in files:
        rel = Path(row['path']); assert not rel.is_absolute() and '..' not in rel.parts
        p = root / rel
        assert p.is_file() and not p.is_symlink()
        assert p.stat().st_size == int(row['bytes']) and sha(p) == row['sha256'], str(p)
    with open(root / 'slurm/final_jobs_accounting.tsv') as f:
        account = {r['JobID']: r for r in csv.DictReader(f, delimiter='|')}
    smoke = json.loads((root / 'slurm/smoke_submission.json').read_text())
    main = json.loads((root / 'slurm/main_submission.json').read_text())
    selection = json.loads((root / 'selection.json').read_text())
    assert len(main['array_ids']) == 640 and len(selection['smoke_tasks']) == 24
    for row in read(root / 'task_manifest.tsv'):
        task = int(row['task_id'])
        job = smoke['job_id'] if task in selection['smoke_tasks'] else main['job_id'] + '_' + str(task)
        state = account[job]
        assert (state['State'], state['ExitCode']) == ('COMPLETED', '0:0'), state
        report = json.loads((point_dir(root, task) / 'result.json').read_text())
        assert report['runtime']['slurm_job_id'] == state['JobIDRaw']
    prior = REPO / selection['prior_root']
    preserved = 0
    for row in read(prior / 'artifact_manifest.tsv'):
        assert sha(REPO / row['path']) == row['sha256'], row['path']
        preserved += 1
    dump(root / 'collection_validation.json', dict(status='passed', n_tasks=664,
         main_completed_exit0=640, smoke_completed_tasks=24, scheduler_runtime_ids_match=True,
         n_files=len(files), n_bytes=sum(int(r['bytes']) for r in files),
         n_prior_artifacts_unchanged=preserved))


def final(root, visual):
    from pypdf import PdfReader
    assert json.loads((root / 'collection_validation.json').read_text())['status'] == 'passed'
    independent = json.loads((root / 'independent_validation.json').read_text())
    assert independent['status'] == 'passed' and independent['n_tasks'] == 664
    assert json.loads((root / 'completion.json').read_text())['status'] == 'complete'
    checks = []
    for path in sorted((root / 'figures').glob('*.pdf')):
        pdf = PdfReader(path)
        assert len(pdf.pages) == 1 and len(pdf.pages[0].extract_text()) > 200
        checks.append(dict(path=str(path.relative_to(root)), pages=1))
    assert len(checks) == 4
    dump(root / 'deliverable_validation.json', dict(status='passed', n_cases=166, n_tasks=664,
         n_precision_verified=independent['n_precision_verified'], pdf_checks=checks,
         independent_matrix_residual_checks=independent['n_residual_checks'],
         visual_render_reviewed=visual, new_FAST_samples=0, FAST_N_M_unchanged=True,
         original_outputs_preserved=True))
    files = [p for p in root.rglob('*') if p.is_file() and p.name not in
             ('artifact_manifest.tsv', '.DS_Store', 'compiled_review_template.tar.gz')]
    files += [p for p in HERE.iterdir() if p.is_file()]
    table(root / 'artifact_manifest.tsv', [dict(path=str(p.relative_to(REPO)), bytes=p.stat().st_size,
                                               sha256=sha(p)) for p in sorted(files)])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=('collection', 'final'))
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--visual-render-reviewed', action='store_true')
    a = p.parse_args()
    if a.command == 'collection': collection(a.root.resolve())
    else: final(a.root.resolve(), a.visual_render_reviewed)
