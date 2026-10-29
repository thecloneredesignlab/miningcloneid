#!/usr/bin/env python3
"""Collect completed global Sobol exports read-only; never evaluate a model.

The same standard-library-only file runs on RED with --inventory. All copies
use an explicit allowlist, a staging directory, and SHA-256 verification.
Existing runtime/provenance files and archived neighborhoods are not replaced.
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import zipfile

RESOLUTIONS = (65536, 131072)
RESULT_FILES = (
    'completion.json', 'indices_repeats.npz', 'indices_summary.npz',
    'indices_summary.tsv.gz', 'not_passed.tsv.gz', 'classification.tsv',
    'row_order.tsv',
)
PROVENANCE_FILES = (
    'plan.json', 'preparation.json', 'source_audit.json',
    'source_solver_audit.tsv', 'status.json', 'smoke.json',
    'launch_receipt.json', 'slurm/jobs.json',
)


def now():
    return datetime.datetime.utcnow().isoformat() + 'Z'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for part in iter(lambda: handle.read(8 * 1024 ** 2), b''):
            h.update(part)
    return h.hexdigest()


def save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write('\n')


def inventory(root):
    root = Path(root).resolve()
    plan = json.loads((root / 'plan.json').read_text())
    assert plan['scope'] == 'global' and plan['resolutions'] == list(RESOLUTIONS)
    assert (plan['n_repeats'], plan['n_oxygen'], plan['n_parameters'], plan['n_pairs']) == (5, 201, 14, 91)
    chosen = []
    for n in RESOLUTIONS:
        directory = root / ('N%d' % n)
        completion = json.loads((directory / 'completion.json').read_text())
        assert completion['status'] == 'complete' and completion['N'] == n
        assert (completion['n_repeats'], completion['n_oxygen'], completion['n_indices']) == (5, 201, 47838)
        assert digest(directory / 'indices_summary.tsv.gz') == completion['summary_sha256']
        chosen.extend(directory / name for name in RESULT_FILES)
        for child in ('figures', 'resolution_comparison'):
            folder = directory / child
            assert folder.is_dir(), str(folder)
            chosen.extend(p for p in sorted(folder.rglob('*')) if p.is_file())
    chosen.extend(root / rel for rel in PROVENANCE_FILES)
    chosen.extend((root / 'slurm/code_snapshot' / name) for name in plan['controller_sha256'])
    chosen.extend(sorted((root / 'inputs').glob('*.tsv')))
    for n, job in ((65536, '21262950'), (131072, '21262953')):
        chosen.extend(sorted((root / 'slurm/logs').glob('report_%d_%s_*.*' % (n, job))))
    records = []
    for path in sorted(set(chosen)):
        assert path.is_file() and not path.is_symlink(), str(path)
        assert root in path.resolve().parents, str(path)
        rel = str(path.relative_to(root))
        assert not any(x in Path(rel).parts for x in ('shards', 'increments', 'units', 'neighborhood10pct'))
        records.append(dict(path=rel, size=path.stat().st_size, sha256=digest(path)))
    for name, expected in plan['controller_sha256'].items():
        assert digest(root / 'slurm/code_snapshot' / name) == expected, name
    bounds = Path(plan['source_initial_root']) / 'global/full/global/parameter_bounds.tsv'
    assert digest(bounds) == plan['source_files'][str(bounds)]
    return dict(created_utc=now(), source_root=str(root), files=records,
                n_files=len(records), total_bytes=sum(r['size'] for r in records),
                external_evidence=[dict(path=str(bounds), export_name='parameter_bounds.tsv',
                                        sha256=digest(bounds), text=bounds.read_text())],
                model_calls=0, neighborhood_files_selected=0)


def remote_inventory(remote, remote_root):
    command = 'bash -lc ' + shlex.quote('python3 - --inventory --root ' + shlex.quote(str(remote_root)))
    p = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20', remote, command],
                       input=Path(__file__).read_text(), stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True, timeout=300)
    if p.returncode:
        raise RuntimeError('Remote inventory failed: ' + p.stderr[-4000:])
    return json.loads(p.stdout)


def verify(path, record):
    path = Path(path)
    assert path.is_file() and not path.is_symlink(), str(path)
    assert path.stat().st_size == record['size'], str(path) + ': size mismatch'
    assert digest(path) == record['sha256'], str(path) + ': SHA-256 mismatch'
    if str(path).endswith('.gz'):
        with gzip.open(str(path), 'rb') as handle:
            for _ in iter(lambda: handle.read(8 * 1024 ** 2), b''):
                pass
    if path.suffix == '.npz':
        with zipfile.ZipFile(str(path)) as archive:
            assert archive.testzip() is None, str(path) + ': NPZ CRC failure'


def collect(args):
    root = Path(args.local_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    collection_id = args.collection_id or datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
    collection = root / ('collection_' + collection_id)
    collection.mkdir(exist_ok=False)
    stage = root / ('.incoming_' + collection_id)
    stage.mkdir(exist_ok=False)
    manifest = remote_inventory(args.remote, args.remote_root)
    save(collection / 'remote_manifest.json', manifest)
    print(json.dumps(dict(stage='inventory', collection=str(collection), files=manifest['n_files'],
                          bytes=manifest['total_bytes'])), flush=True)
    allowlist = collection / 'files_from.txt'
    allowlist.write_text(''.join(r['path'] + '\n' for r in manifest['files']))
    # Refuse differing existing result files; shared runtime files are archived
    # separately in this collection's provenance directory.
    for record in manifest['files']:
        if record['path'].split('/')[0] in ['N%d' % n for n in RESOLUTIONS]:
            existing = root / record['path']
            if existing.exists():
                verify(existing, record)
    command = ['rsync', '-a', '-z', '--no-perms', '--no-owner', '--no-group',
               '--files-from=' + str(allowlist), '-e', 'ssh -o BatchMode=yes -o ConnectTimeout=20',
               args.remote + ':' + str(args.remote_root).rstrip('/') + '/', str(stage) + '/']
    with (collection / 'rsync.stdout.log').open('x') as out, (collection / 'rsync.stderr.log').open('x') as err:
        p = subprocess.run(command, stdout=out, stderr=err)
    save(collection / 'transfer.json', dict(command=command, returncode=p.returncode, completed_utc=now()))
    if p.returncode:
        raise RuntimeError('Transfer failed; preserved staging directory: ' + str(stage))
    for record in manifest['files']:
        verify(stage / record['path'], record)
    final_inventory = remote_inventory(args.remote, args.remote_root)
    assert final_inventory['files'] == manifest['files'], 'Remote export changed during collection'
    save(collection / 'remote_manifest_after.json', final_inventory)
    for n in RESOLUTIONS:
        source, dest = stage / ('N%d' % n), root / ('N%d' % n)
        if dest.exists():
            # Earlier preflight verified every selected destination file.
            shutil.rmtree(str(source))
        else:
            os.rename(str(source), str(dest))
    provenance = collection / 'provenance'
    os.rename(str(stage), str(provenance))
    for record in manifest['external_evidence']:
        content = record['text'].encode('utf-8')
        assert hashlib.sha256(content).hexdigest() == record['sha256']
        path = provenance / 'reference' / record['export_name']
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as handle:
            handle.write(content)
    rows = []
    for record in manifest['files']:
        is_result = record['path'].split('/')[0] in ['N%d' % n for n in RESOLUTIONS]
        path = (root if is_result else provenance) / record['path']
        verify(path, record)
        rows.append(dict(record, local_path=str(path), sha256_verified=True, container_crc_verified=path.suffix in ('.gz', '.npz')))
    with (collection / 'verified_files.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    save(collection / 'collection_receipt.json', dict(status='complete', completed_utc=now(),
         source_root=args.remote_root, local_root=str(root), n_files=len(rows),
         total_bytes=manifest['total_bytes'], all_file_hashes_verified=True,
         source_unchanged_during_transfer=True, model_calls=0,
         neighborhood_files_selected=0, collector_sha256=digest(__file__),
         manifest_sha256=digest(collection / 'remote_manifest.json')))
    print(json.dumps(dict(stage='complete', collection=str(collection), files=len(rows), model_calls=0)), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', action='store_true')
    parser.add_argument('--root')
    parser.add_argument('--remote', default='4482173@red.moffitt.org')
    parser.add_argument('--remote-root')
    parser.add_argument('--local-root')
    parser.add_argument('--collection-id')
    args = parser.parse_args()
    if args.inventory:
        print(json.dumps(inventory(args.root), ensure_ascii=False))
    else:
        if not args.remote_root or not args.local_root:
            parser.error('--remote-root and --local-root are required')
        collect(args)


if __name__ == '__main__':
    main()
