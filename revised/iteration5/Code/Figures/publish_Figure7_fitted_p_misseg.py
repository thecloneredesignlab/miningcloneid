#!/usr/bin/env python3
"""Allowlist-only first publication, after independent and visual acceptance."""
import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for data in iter(lambda: f.read(1048576), b''):
            h.update(data)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--workspace', type=Path, required=True)
    p.add_argument('--candidate-root', type=Path, required=True)
    p.add_argument('--render-root', type=Path, required=True)
    p.add_argument('--validation-root', type=Path, required=True)
    p.add_argument('--protected-before', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('--publish', choices=['yes'], required=True)
    a = p.parse_args()
    w = a.workspace.resolve()
    target = w / 'data/Figures/Figure7/fitted_p_misseg'
    assert not target.exists(), 'The new condition already exists; refuse an implicit replacement'
    done = json.loads((a.candidate_root / 'provenance/completion.json').read_text())
    accepted = json.loads((a.validation_root / 'result.json').read_text())
    visual = json.loads((a.validation_root / 'visual_review.json').read_text())
    assert done['status'] == 'validated_complete'
    assert accepted['status'] == 'passed' and accepted['fingerprint'] == done['fingerprint']
    assert visual['status'] == 'passed' and visual['pages_reviewed'] == 25
    protected = json.loads(a.protected_before.read_text())
    assert all(sha(w / r['path']) == r['sha256'] for r in protected)
    images = []
    for d in ['Figures', 'manuscript/Figures']:
        for stem in ['Figure7_fitted_p_misseg', 'supp_fig7_8_fitted_p_misseg']:
            for ext in ['pdf', 'png']:
                rel = Path(d) / (stem + '.' + ext)
                src = a.render_root / rel
                assert src.is_file() and not (w / rel).exists(), str(rel)
                images.append((rel, src, sha(src)))
    stage = Path(tempfile.mkdtemp(prefix='.fitted_p_misseg_', dir=target.parent))
    published = []
    try:
        for d in ['prepared', 'tables', 'postprocess', 'provenance']:
            shutil.copytree(a.candidate_root / d, stage / d)
        # The retained dataset contains no task, route, cache, log or audit root.
        for file in stage.rglob('*'):
            assert not file.is_symlink()
        manifest = stage / 'provenance/validated_files.tsv'
        import csv
        with manifest.open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        assert all(sha(stage / r['relative_path']) == r['sha256'] for r in rows)
        assert sha(stage / 'prepared/analysis_bundle.rds') == done['bundle_sha256']
        assert sha(manifest) == done['manifest_sha256']
        os.rename(stage, target)
        published.append(target)
        for rel, src, expected in images:
            dest = w / rel
            fd, tmp = tempfile.mkstemp(prefix='.fitted_', dir=dest.parent)
            os.close(fd)
            try:
                shutil.copy2(src, tmp)
                assert sha(Path(tmp)) == expected
                os.link(tmp, dest)  # exclusive publication; never replace an existing file
            finally:
                Path(tmp).unlink(missing_ok=True)
            published.append(dest)
        assert all(sha(w / r['path']) == r['sha256'] for r in protected)
        receipt = dict(status='published', fingerprint=done['fingerprint'],
                       new_data_directory=str(target), protected_original_files=len(protected),
                       protected_original_files_unchanged=True,
                       images=[dict(path=str(rel), sha256=s) for rel, src, s in images],
                       existing_products_replaced=False)
        a.receipt.write_text(json.dumps(receipt, indent=2))
        print(json.dumps(receipt))
    except BaseException:
        # Roll back only the exact new products created by this invocation.
        for path in reversed(published):
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
        raise
    finally:
        if stage.exists():
            shutil.rmtree(stage)


if __name__ == '__main__':
    main()
