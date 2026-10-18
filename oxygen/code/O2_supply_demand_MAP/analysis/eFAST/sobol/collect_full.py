#!/usr/bin/env python3
"""Collect completed N4096 Sobol deliverables using their checksum inventory."""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import tempfile


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def remote_snapshot(host, root):
    code = '''import csv,hashlib,json,pathlib,subprocess,sys
p=pathlib.Path(sys.argv[1]); rec=p/'recovery_20261018'; retry=rec/'retry5734'
s=json.loads((p/'status.json').read_text())
expected=dict(N=4096,n_neighborhoods=500,n_repeats=5,oxygen_points=201,n_pairs=91,n_response_surfaces=54)
if s.get('status')!='complete' or any(s.get(k)!=v for k,v in expected.items()):
    raise SystemExit('Full N4096 analysis is not complete')
audit=json.loads((rec/'all_cases_validated.json').read_text())
if audit.get('status')!='passed' or audit.get('n_cases')!=501:
    raise SystemExit('All 501 domains must pass the input/completion audit')
preflight=json.loads((retry/'preflight.json').read_text())
preserved=preflight['preserved_cases']
if len(preserved)!=500: raise SystemExit('Wrong preserved domain population')
for row in preserved:
    if hashlib.sha256((p/row['path']).read_bytes()).hexdigest()!=row['sha256']:
        raise SystemExit('A completed domain was replaced: '+row['path'])
with (retry/'slurm/submission_jobs.tsv').open() as handle:
    jobs=[r['job_id'] for r in csv.DictReader(handle,delimiter='\\t') if r['role']=='finish']
if len(jobs)!=1: raise SystemExit('Expected exactly one replacement final-report job')
account=subprocess.check_output(['/cm/shared/apps/slurm/current/bin/sacct','-X','-n','-P','-j',jobs[0],'-o','State,ExitCode'],universal_newlines=True).strip()
if account.splitlines()[0].split('|')[:2]!=['COMPLETED','0:0']:
    raise SystemExit('Final report job is not COMPLETED/0:0: '+account)
layout=json.loads((retry/'figure_layout_patch.json').read_text())
layout_account=subprocess.check_output(['/cm/shared/apps/slurm/current/bin/sacct','-X','-n','-P','-j',layout['job_id'],'-o','State,ExitCode'],universal_newlines=True).strip()
if layout.get('status')!='complete' or not layout.get('input_arrays_unchanged') or layout_account.splitlines()[0].split('|')[:2]!=['COMPLETED','0:0']:
    raise SystemExit('Figure layout refresh is incomplete or changed data')
print(json.dumps(dict(status=s,audit=audit,final_job_id=jobs[0],final_accounting=account,n_preserved_domains=500,
    figure_layout_job_id=layout['job_id'],figure_layout_accounting=layout_account,
    manifest_sha256=hashlib.sha256((p/'collection_manifest.tsv').read_bytes()).hexdigest())))
'''
    remote = 'python3 -c ' + shlex.quote(code) + ' ' + shlex.quote(root)
    return json.loads(subprocess.check_output(
        ['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20', host, remote],
        universal_newlines=True, timeout=120))


def validate_path(value):
    path = Path(value)
    if (not value or path.is_absolute() or '..' in path.parts
            or any(part in path.parts for part in ('runs', 'logs', 'code_snapshot'))
            or path.suffix not in ('.json', '.tsv', '.gz', '.md', '.pdf', '.png')
            or path.name in ('samples.tsv.gz', 'outputs.tsv.gz')):
        raise ValueError('Unsafe or raw-data inventory path: ' + value)
    return str(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host', default='4482173@red.moffitt.org')
    parser.add_argument('--remote-dir', default='/share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/sobol/full_N4096')
    parser.add_argument('--local-dir', default=str(Path(__file__).resolve().parents[6]/'oxygen/results/eFAST/sobol/full_N4096'))
    args = parser.parse_args()
    before = remote_snapshot(args.host, args.remote_dir)
    target = Path(args.local_dir).resolve(); target.mkdir(parents=True, exist_ok=True)
    remote = args.host + ':' + args.remote_dir.rstrip('/') + '/'
    common = ['rsync', '-azO', '--no-perms', '--no-owner', '--no-group', '-e', 'ssh']
    subprocess.run(common + [remote+'collection_manifest.tsv', str(target)+'/'], check=True)
    if sha(target/'collection_manifest.tsv') != before['manifest_sha256']:
        raise ValueError('Inventory changed before copy')
    with (target/'collection_manifest.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    paths = [validate_path(row['path']) for row in rows]
    if len(paths) != len(set(paths)):
        raise ValueError('Duplicate inventory paths')
    print(json.dumps(dict(copying_files=len(paths), bytes=sum(int(r['bytes']) for r in rows))), flush=True)
    with tempfile.NamedTemporaryFile('wb') as handle:
        handle.write(('\0'.join(paths)+'\0').encode()); handle.flush()
        subprocess.run(common + ['--from0', '--files-from='+handle.name, remote, str(target)+'/'], check=True)
    after = remote_snapshot(args.host, args.remote_dir)
    if after['manifest_sha256'] != before['manifest_sha256']:
        raise ValueError('Source inventory changed during copy')
    for row in rows:
        file = target/row['path']
        if not file.is_file() or file.stat().st_size != int(row['bytes']) or sha(file) != row['sha256']:
            raise ValueError('Copied source failed checksum validation: '+row['path'])
        file.chmod(0o644)
    receipt = dict(collected_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   source=args.host+':'+args.remote_dir,
                   source_manifest_sha256=after['manifest_sha256'],
                   final_job_id=after['final_job_id'], final_accounting=after['final_accounting'],
                   n_verified_files=len(rows), verified_bytes=sum(int(r['bytes']) for r in rows),
                   n_preexisting_domains_unchanged=after['n_preserved_domains'],
                   figure_layout_job_id=after['figure_layout_job_id'],figure_layout_accounting=after['figure_layout_accounting'],
                   status=after['status'], audit=after['audit'], collector_sha256=sha(__file__))
    (target/'artifact_collection_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__':
    main()
