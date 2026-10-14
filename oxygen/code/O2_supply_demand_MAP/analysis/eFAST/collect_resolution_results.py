#!/usr/bin/env python3
"""Read-only HPC receipts audit and allowlisted collection of the finished pilot.

Python 3.6 compatible. Run HPC audit in a login shell so sacct is on PATH.
No model evaluation, Slurm submission, source edit, or raw-output transfer.
"""
import argparse
import ast
import csv
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import struct
import sys
import zipfile


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def files(root):
    found = [p for p in root.iterdir() if p.is_file() and p.suffix in ('.json', '.tsv', '.gz', '.npz')]
    for directory, suffixes in [('summaries', ('.npz', '.tsv', '.json')),
                                 ('figures', ('.pdf', '.png')),
                                 ('slurm', ('.json', '.tsv')), ('smoke', ('.json',))]:
        found.extend(p for p in (root / directory).iterdir() if p.is_file() and p.suffix in suffixes)
    return sorted(set(found))


def manifest(root):
    print('path\tsha256\tbytes')
    for path in files(root):
        print('{}\t{}\t{}'.format(path.relative_to(root), sha(path), path.stat().st_size))


def numeric_npy(archive, name):
    """Read only fixed numeric gap arrays; no pickle, NumPy install or evaluation."""
    payload = archive.read(name + '.npy')
    assert payload[:6] == b'\x93NUMPY'
    major = payload[6]
    length_bytes = 2 if major == 1 else 4
    assert major in (1, 2, 3)
    header_length = int.from_bytes(payload[8:8+length_bytes], 'little')
    offset = 8 + length_bytes + header_length
    header = ast.literal_eval(payload[8+length_bytes:offset].decode('utf-8').strip())
    assert not header['fortran_order'] and header['descr'] in ('<i8', '<f8')
    count = 1
    for size in header['shape']:count *= size
    code = 'q' if header['descr'] == '<i8' else 'd'
    return header['shape'], struct.unpack('<' + str(count) + code, payload[offset:])


def gap_qc(repo):
    root = repo / 'oxygen/results/eFAST/neighborhood10pct/pilot/N1025_N2049'
    with (root / 'slurm/task_manifest.tsv').open() as handle:
        tasks = list(csv.DictReader(handle, delimiter='\t'))
    summaries = {}
    for task in tasks:
        part = root / 'runs' / task['fit_seed'] / ('N' + task['N']) / 'runs' / ('N%s_R%s' % (task['N'], task['replicate'])) / 'trajectories' / ('P%02d_%s' % (int(task['parameter_index'])+1, task['focal_parameter']))
        path = part / 'indices.npz'
        digest = sha(path)
        receipt = json.loads((root / 'slurm/task_receipts' / ('task_%05d.json' % int(task['task_id']))).read_text())
        assert digest == receipt['indices_sha256'] and receipt['status'] == 'complete'
        with zipfile.ZipFile(str(path)) as archive:
            shape, counts = numeric_npy(archive, 'gap_counts');assert shape == (3,4)
            shape, minimum = numeric_npy(archive, 'gap_min');assert shape == (3,)
        key = (task['fit_seed'], int(task['N']))
        if key not in summaries:summaries[key] = dict(counts=[0]*12, minimum=[float('inf')]*3, hashes=[], n=0)
        record = summaries[key]
        record['counts'] = [a+b for a,b in zip(record['counts'], counts)]
        record['minimum'] = [min(a,b) for a,b in zip(record['minimum'], minimum)]
        record['hashes'].append(digest);record['n'] += 1
    writer = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    writer.writerow(['fit_seed','N','oxygen_band','n_trajectories','n_evaluations','fraction_gap_lt_1e_6','fraction_gap_lt_1e_4','fraction_gap_lt_1e_3','minimum_spectral_gap','source_cache_hash_chain_sha256'])
    for (seed,n), record in sorted(summaries.items()):
        assert record['n'] == 70
        digest = hashlib.sha256(('\n'.join(record['hashes'])+'\n').encode()).hexdigest()
        for j, (label, oxygen_count) in enumerate([('low_0_to_1pct',41),('middle_1_to_3pct',79),('high_3_to_5pct',81)]):
            count = record['counts'][4*j:4*j+4]
            assert count[0] == 70*n*oxygen_count
            writer.writerow([seed,n,label,70,count[0],count[1]/count[0],count[2]/count[0],count[3]/count[0],record['minimum'][j],digest])


def audit(repo):
    root = repo / 'oxygen/results/eFAST/neighborhood10pct/pilot/N1025_N2049'
    code = repo / 'oxygen/code/O2_supply_demand_MAP/analysis/eFAST'
    plan = json.loads((root / 'slurm/submission_plan.json').read_text())
    for name, digest in plan['code_sha256'].items():
        assert sha(code / name) == digest, name
    with (root / 'slurm/task_manifest.tsv').open() as handle:
        tasks = list(csv.DictReader(handle, delimiter='\t'))
    assert len(tasks) == 1680
    proofs = []
    for task in tasks:
        receipt = json.loads((root / 'slurm/task_receipts' / ('task_%05d.json' % int(task['task_id']))).read_text())
        assert receipt['status'] == 'complete' and receipt['fit_seed'] == task['fit_seed'] and receipt['N'] == task['N'], task
        part = root / 'runs' / task['fit_seed'] / ('N' + task['N']) / 'runs' / ('N%s_R%s' % (task['N'], task['replicate'])) / 'trajectories' / ('P%02d_%s' % (int(task['parameter_index']) + 1, task['focal_parameter']))
        assert sha(part / 'indices.npz') == receipt['indices_sha256'], str(part)
        assert sha(part / 'outputs.tsv.gz.receipt.json') == receipt['trajectory_receipt_sha256'], str(part)
        proof_path = part / 'recovery/numerical_recovery.json'
        if proof_path.exists():
            proof = json.loads(proof_path.read_text())
            for name, digest in proof['artifact_sha256'].items():
                assert sha(part / 'recovery' / name) == digest, str(part / name)
            for name, digest in proof['code_sha256'].items():
                assert sha(code / name) == digest, name
            proofs.append(dict(task_id=int(task['task_id']), fit_seed=task['fit_seed'], N=int(task['N']),
                focal_parameter=task['focal_parameter'], replicate=int(task['replicate']),
                proof_sha256=sha(proof_path), proof=proof,
                corrections=json.loads((part / 'recovery/corrections.json').read_text())))
    raw = subprocess.check_output(['sacct', '-X', '-j', '20833838,20833839,20833840',
        '--format=JobID,State,ExitCode,ElapsedRaw,Start,End', '-n', '-P'], universal_newlines=True)
    entries = {}
    for line in raw.splitlines():
        parts = line.split('|')
        if len(parts) < 6:
            continue
        row = dict(zip(('job_id', 'state', 'exit_code', 'elapsed_seconds', 'start', 'end'), parts))
        row['elapsed_seconds'] = int(row['elapsed_seconds'])
        if row['job_id'] not in entries or row['end'] >= entries[row['job_id']]['end']:
            entries[row['job_id']] = row
    assert len([k for k in entries if k.startswith('20833838_')]) == 1680
    assert len([k for k in entries if k.startswith('20833839_')]) == 24
    assert '20833840' in entries
    assert all(r['state'] == 'COMPLETED' and r['exit_code'] == '0:0' for r in entries.values())
    print(json.dumps(dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        hpc_git_commit=subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], universal_newlines=True).strip(),
        frozen_source_hashes_match=True, n_verified_trajectory_caches=len(tasks),
        n_recovery_trajectories=len(proofs), recovery_audits=proofs,
        slurm_unique_jobs=list(entries.values()), new_jobs_submitted=False), indent=2, sort_keys=True))


def verify_local(root, source):
    with source.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    for row in rows:
        relative = Path(row['path'])
        assert not relative.is_absolute() and '..' not in relative.parts
        path = root / relative
        assert path.stat().st_size == int(row['bytes']), row['path']
        assert sha(path) == row['sha256'], row['path']
        path.chmod(0o644)
    completion = json.loads((root / 'completion.json').read_text())
    assert completion['status'] == 'complete' and completion['n_tasks'] == 1680 and completion['n_seed_N_summaries'] == 24
    for name, digest in completion['artifact_sha256'].items():
        assert sha(root / name) == digest, name
    destination = root / 'completed_results_manifest.tsv'
    destination.write_bytes(source.read_bytes())
    (root / 'results_collection_receipt.json').write_text(json.dumps(dict(
        collected_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status='passed',
        n_hpc_artifacts=len(rows), total_bytes=sum(int(r['bytes']) for r in rows),
        all_hpc_sha256_matched=True, original_completion_hashes_matched=True,
        manifest_sha256=sha(destination), raw_trajectories_location='HPC only',
        analysis_complete=True, new_jobs_submitted=False), indent=2, sort_keys=True) + '\n')
    print('Verified {} collected files against HPC SHA256 and original completion hashes'.format(len(rows)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('manifest', 'audit', 'verify-local', 'gap-qc'))
    parser.add_argument('--root', type=Path)
    parser.add_argument('--repo-root', type=Path)
    parser.add_argument('--manifest', type=Path)
    args = parser.parse_args()
    if args.command == 'manifest':
        manifest(args.root)
    elif args.command == 'audit':
        audit(args.repo_root)
    elif args.command == 'gap-qc':
        gap_qc(args.repo_root)
    else:
        verify_local(args.root, args.manifest)


if __name__ == '__main__':
    main()
