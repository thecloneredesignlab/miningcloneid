#!/usr/bin/env python3
"""Capture immutable inputs and read-only scheduler evidence; never dispatch work."""
import argparse
import csv
import hashlib
import io
from pathlib import Path
import shutil
import subprocess
import tarfile
from datetime import datetime, timezone

ROOT = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'oxygen/code/O2_supply_demand_MAP').is_dir())
OLD_REL = Path('oxygen/results/structural_identifiability/20261112_invivo_continuation_141206Z')
CULTURE_REL = Path('oxygen/results/structural_identifiability/20261113_culture_symbolic_01')
SOURCE = [
    'model/model_O2_supply_demand_MAP.R',
    'model/model_O2_supply_demand_MAP.cpp',
    'util/o2_supply_demand_map_shared.R',
    'util/o2_supply_demand_map_common_semantics.R',
    'util/o2_supply_demand_map_fit_invivo_backend.R',
]
AUTHORITY = [
    'tables/resolution_closure_mode.tsv', 'tables/resolution_terminal.tsv',
    'tables/resolution_final_phase_accounting.tsv', 'tables/delivery_verified.tsv',
    'tables/baseline_policy_gates.tsv', 'tables/closure_audit.tsv',
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    args = parser.parse_args()
    run = args.run.resolve()
    expected_parent = ROOT / 'oxygen/results/structural_identifiability'
    if run.parent != expected_parent or 'invivo_structural_' not in run.name:
        raise SystemExit('Output must be an independent invivo_structural run')
    if run.exists():
        raise SystemExit('Refusing to overwrite an existing capture')
    run.mkdir()
    started = datetime.now(timezone.utc).isoformat()
    before = subprocess.check_output(['git', 'status', '--porcelain=v1', '-uall'], cwd=ROOT)
    (run / 'git_status_before.txt').write_bytes(before)
    inventory = []
    files = [(OLD_REL / 'source/numerical' / p, Path('frozen_source') / p, 'old_numerical_source') for p in SOURCE]
    files += [(OLD_REL / 'manifests' / p, Path('old_contract') / p, 'old_contract') for p in ['contract.rds', 'transformed_parameters.tsv']]
    files += [(Path('oxygen/results/analysis_documentation/frozen_inputs/in_vitro/design_and_contracts/meeting_aligned_synthetic_identifiability_master_plan_20261112.md'), Path('master_plan.md'), 'master_plan')]
    files += [(CULTURE_REL / p, Path('culture_gate') / p, 'culture_gate') for p in [
        'reports/step04_acceptance.tsv', 'reports/STOP.tsv', 'reports/reduction_summary.tsv',
        'reports/tool_access_and_compatibility.md', 'reports/parameter_reduction_map.tsv',
        'reports/structural_limitations.md', 'accounting/final_accounting.tsv',
        'accounting/jobs.tsv', 'REPORT.md', 'delivery_verification.tsv',
    ]]
    for src_rel, dest_rel, role in files:
        src = ROOT / src_rel
        dst = run / 'inputs' / dest_rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        inventory.append(dict(source=str(src), frozen_copy=str(dst.relative_to(run)), role=role,
                              sha256=sha(src), bytes=src.stat().st_size))
    for rel in SOURCE:
        current = ROOT / 'oxygen/code/O2_supply_demand_MAP' / rel
        old = ROOT / OLD_REL / 'source/numerical' / rel
        inventory.append(dict(source=str(current), frozen_copy='not_copied_not_executed',
                              role='current_checkout_comparison', sha256=sha(current), bytes=current.stat().st_size))
    tsv(run / 'inputs/source_inventory.tsv', ['source', 'frozen_copy', 'role', 'sha256', 'bytes'], inventory)
    (run / 'inputs/git_head.txt').write_text(subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, universal_newlines=True))

    # The remote program only reads files, hashes the pinned SIF and invokes sacct/squeue.
    # It writes a tar stream to stdout; no remote file, controller or counter is changed.
    remote = r'''
import csv, hashlib, io, pathlib, subprocess, sys, tarfile
from datetime import datetime, timezone
root = pathlib.Path('/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability')
old = root / 'oxygen/results/structural_identifiability/20261112_invivo_continuation_141206Z'
culture = root / 'oxygen/results/structural_identifiability/20261113_culture_symbolic_01'
image = pathlib.Path('/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif')
payload = {}
def add(name, data):
    payload[name] = data if isinstance(data, bytes) else data.encode()
def read(name, path):
    add(name, path.read_bytes())
def command(name, argv):
    p = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=45)
    if p.returncode:
        raise RuntimeError(name + ': ' + p.stderr)
    add(name, p.stdout)
    add(name + '.stderr.txt', p.stderr)
def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()
start = datetime.now(timezone.utc).isoformat()
generation = (old / 'state/CURRENT').read_text().strip()
if '/' in generation or not generation.startswith('g_'):
    raise RuntimeError('Invalid generation')
jobs_path = old / 'state' / generation / 'jobs.tsv'
jobs_data = jobs_path.read_bytes()
jobs = list(csv.DictReader(io.StringIO(jobs_data.decode()), delimiter='\t'))
ids = sorted({r['job_id'] for r in jobs if r['job_id']})
if len(ids) != 125 or any(not x.isdigit() for x in ids):
    raise RuntimeError('Owned-parent set changed; review before proceeding')
add('old_invivo/CURRENT', generation + '\n')
add('old_invivo/jobs.tsv', jobs_data)
read('old_invivo/generation_manifest.tsv', old / 'state' / generation / 'manifest.tsv')
for rel in ['tables/resolution_closure_mode.tsv','tables/resolution_terminal.tsv',
            'tables/resolution_final_phase_accounting.tsv','tables/delivery_verified.tsv',
            'tables/baseline_policy_gates.tsv','tables/closure_audit.tsv']:
    read('old_invivo/' + rel, old / rel)
slurm = '/cm/shared/apps/slurm/current/bin/'
command('old_invivo/sacct.tsv', [slurm+'sacct','-n','-X','-P','-j',','.join(ids),
        '-o','JobIDRaw,State,ExitCode,AllocCPUS,ElapsedRaw,Start,End'])
command('old_invivo/squeue.tsv', [slurm+'squeue','-h','-j',','.join(ids),'-o','%i|%T|%C|%M|%l|%j'])
for rel in ['reports/step04_acceptance.tsv','reports/STOP.tsv','accounting/jobs.tsv','accounting/final_accounting.tsv']:
    read('culture_step04/' + rel, culture / rel)
cjobs = list(csv.DictReader(io.StringIO((culture / 'accounting/jobs.tsv').read_text()), delimiter='\t'))
cidkey = 'job_id' if 'job_id' in cjobs[0] else 'job'
cids = sorted({r[cidkey] for r in cjobs if r[cidkey]})
if any(not x.isdigit() for x in cids):
    raise RuntimeError('Invalid culture job ID')
command('culture_step04/sacct.tsv', [slurm+'sacct','-n','-X','-P','-j',','.join(cids),
        '-o','JobIDRaw,State,ExitCode,AllocCPUS,ElapsedRaw,Start,End'])
command('culture_step04/squeue.tsv', [slurm+'squeue','-h','-j',','.join(cids),'-o','%i|%T|%C|%M|%l|%j'])
hashes = []
for rel in ['model/model_O2_supply_demand_MAP.R','model/model_O2_supply_demand_MAP.cpp',
            'util/o2_supply_demand_map_shared.R','util/o2_supply_demand_map_common_semantics.R',
            'util/o2_supply_demand_map_fit_invivo_backend.R']:
    p = old / 'source/numerical' / rel
    hashes.append((str(p),digest(p),p.stat().st_size))
for rel in ['manifests/contract.rds','manifests/transformed_parameters.tsv']:
    p = old / rel
    hashes.append((str(p),digest(p),p.stat().st_size))
hashes.append((str(image),digest(image),image.stat().st_size))
add('source_hashes.tsv', 'path\tsha256\tbytes\n' + ''.join('%s\t%s\t%s\n'%r for r in hashes))
after = (old / 'state/CURRENT').read_text().strip()
if after != generation or jobs_path.read_bytes() != jobs_data:
    raise RuntimeError('Incoherent old-owner snapshot')
add('capture_receipt.tsv', 'key\tvalue\nstart_utc\t%s\nend_utc\t%s\nold_generation\t%s\nold_pointer_stable\ttrue\nremote_files_written\t0\nmodel_calls\t0\ncontroller_calls\t0\nslurm_submissions\t0\n' % (start,datetime.now(timezone.utc).isoformat(),generation))
with tarfile.open(fileobj=sys.stdout.buffer, mode='w|') as tar:
    for name,data in sorted(payload.items()):
        item = tarfile.TarInfo(name)
        item.size = len(data)
        tar.addfile(item,io.BytesIO(data))
'''
    p = subprocess.run(['ssh', '-o', 'LogLevel=ERROR', '-o', 'BatchMode=yes', '-o',
                        'ConnectTimeout=15', '-o', 'ServerAliveInterval=15', '-o',
                        'ServerAliveCountMax=2', '4482173@red.moffitt.org', 'python3', '-B', '-'],
                       input=remote.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=240)
    (run / 'ssh_stderr.txt').write_bytes(p.stderr)
    if p.returncode:
        raise SystemExit('Read-only HPC capture failed: ' + p.stderr.decode())
    out = run / 'snapshots/hpc'
    out.mkdir(parents=True)
    with tarfile.open(fileobj=io.BytesIO(p.stdout), mode='r:') as archive:
        for member in archive.getmembers():
            dest = out / member.name
            if not member.isfile() or not dest.resolve().is_relative_to(out.resolve()):
                raise SystemExit('Unsafe archive member')
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(archive.extractfile(member).read())
    tsv(run / 'capture.tsv', ['key', 'value'], [
        dict(key='start_utc', value=started), dict(key='end_utc', value=datetime.now(timezone.utc).isoformat()),
        dict(key='hpc_mutation', value='false'), dict(key='model_calls', value=0),
        dict(key='symbolic_calls', value=0), dict(key='slurm_submissions', value=0)])
    print(run)


if __name__ == '__main__':
    main()
