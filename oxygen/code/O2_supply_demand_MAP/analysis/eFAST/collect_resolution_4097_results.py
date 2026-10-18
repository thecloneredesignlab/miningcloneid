#!/usr/bin/env python3
"""Audit and collect the completed six-seed N4097 pilot (Python 3.6+).

Reads existing indices/receipts only; no model evaluation or job submission.
Raw trajectories remain on HPC. Original submission/collection receipts persist.
"""
import argparse
import ast
import collections
import csv
import datetime
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import zipfile

SEEDS = ('seed25', 'seed101', 'seed138', 'seed174', 'seed66', 'seed311')
PARAMETERS = ('lam_max','p_mis_base','p_misseg','k_o_mis','buffer_smax','buffer_beta',
              'buffer_n_exp','p_wgd','alpha_o2','gamma_growth','mu_hp','gamma_mu','O2_crit','n_O')
JOB_IDS = {'smoke':'20863802','trajectory':'20864027','seed_summary':'20864028','finalize':'20864029'}


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''): digest.update(block)
    return digest.hexdigest()


def read(path):
    with path.open() as handle: return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def numeric_npy(archive, name):
    payload = archive.read(name+'.npy')
    assert payload[:6] == b'\x93NUMPY' and payload[6] in (1,2,3)
    length = 2 if payload[6] == 1 else 4
    size = int.from_bytes(payload[8:8+length], 'little'); offset = 8+length+size
    header = ast.literal_eval(payload[8+length:offset].decode('utf-8').strip())
    assert not header['fortran_order'] and header['descr'] in ('<i8','<f8')
    count = 1
    for dimension in header['shape']: count *= dimension
    code = 'q' if header['descr'] == '<i8' else 'd'
    return header['shape'], struct.unpack('<'+str(count)+code, payload[offset:])


def validate_tasks(tasks, parameters=PARAMETERS):
    assert len(tasks) == 420 and len(parameters) == 14
    for i, row in enumerate(tasks):
        seed = SEEDS[i//70]; repeat = (i%70)//14+1; p = i%14
        expected = dict(task_id=i+1,summary_id=i//70+1,fit_seed=seed,fit_seed_number=int(seed[4:]),
            N=4097,M=4,replicate=repeat,parameter_index=p,focal_parameter=parameters[p],
            phase_seed=202610080+int(seed[4:])*10+repeat-1,
            parent_sample_id_start=p*4097+1,parent_sample_id_end=(p+1)*4097)
        assert row == {k:str(v) for k,v in expected.items()}, row


def audit(repo):
    root = repo/'oxygen/results/eFAST/neighborhood10pct/pilot/N4097'
    code = repo/'oxygen/code/O2_supply_demand_MAP/analysis/eFAST'
    plan = json.loads((root/'slurm/submission_plan.json').read_text())
    assert plan['fit_seeds'] == list(SEEDS) and plan['resolutions'] == [4097]
    assert (plan['n_tasks'],plan['n_summaries'],plan['n_phases'],plan['n_parameters'],plan['n_oxygen']) == (420,6,5,14,201)
    for name,digest in plan['code_sha256'].items(): assert sha(code/name) == digest, name
    for name,digest in plan['comparison_source_sha256'].items(): assert sha(Path(name)) == digest, name
    for name,field in [('seed_selection.tsv','seed_selection_sha256'),('baseline_indices.npz','comparison_indices_sha256'),
                       ('slurm/task_manifest.tsv','task_manifest_sha256'),('slurm/summary_manifest.tsv','summary_manifest_sha256')]:
        assert sha(root/name) == plan[field], name
    tasks = read(root/'slurm/task_manifest.tsv')
    validate_tasks(tasks)
    smoke = json.loads((root/'smoke/smoke_receipt.json').read_text())
    assert smoke['status'] == 'passed' and smoke['submission_plan_sha256'] == sha(root/'slurm/submission_plan.json')
    completion = json.loads((root/'completion.json').read_text())
    assert completion['status'] == 'complete' and completion['n_tasks'] == 420 and completion['n_seed_N_summaries'] == 6
    for name,digest in completion['artifact_sha256'].items(): assert sha(root/name) == digest, name
    gaps, proofs, task_sources = {}, [], []
    for task in tasks:
        task_id = int(task['task_id'])
        marker_path = root/'slurm/task_receipts'/('task_%05d.json'%task_id)
        marker = json.loads(marker_path.read_text())
        assert marker['status'] == 'complete'
        for field,value in task.items(): assert marker[field] == value
        parent = root/'runs'/task['fit_seed']/'N4097/runs'/('N4097_R'+task['replicate'])
        part = parent/'trajectories'/('P%02d_%s'%(int(task['parameter_index'])+1,task['focal_parameter']))
        assert sha(part/'indices.npz') == marker['indices_sha256']
        assert sha(part/'outputs.tsv.gz.receipt.json') == marker['trajectory_receipt_sha256']
        cache = json.loads((part/'indices.receipt.json').read_text())
        meta = json.loads((part/'metadata.json').read_text())
        output_receipt = json.loads((part/'outputs.tsv.gz.receipt.json').read_text())
        assert cache['cache_sha256'] == marker['indices_sha256'] and cache['metadata_sha256'] == sha(part/'metadata.json')
        assert cache['analyzer_sha256'] == sha(code/'resolution_pilot.py') and cache['salib_version'] == '1.5.2'
        assert output_receipt['outputs_sha256'] == cache['outputs_sha256']
        assert output_receipt['metadata_sha256'] == cache['metadata_sha256']
        assert output_receipt['samples_sha256'] == meta['samples_sha256'] and output_receipt['n_rows'] == 4097*201
        assert (part/'outputs.tsv.gz').stat().st_size > 0
        assert sha(part/'samples.tsv.gz') == meta['samples_sha256']
        assert meta['N'] == 4097 and meta['M'] == 4 and meta['trajectory_parameter_index'] == int(task['parameter_index'])
        assert meta['parent_metadata_sha256'] == sha(parent/'metadata.json')
        task_sources.append(dict(task_id=task_id,fit_seed=task['fit_seed'],replicate=int(task['replicate']),
            parameter=task['focal_parameter'],marker_sha256=sha(marker_path),indices_sha256=marker['indices_sha256'],
            output_receipt_sha256=marker['trajectory_receipt_sha256']))
        with zipfile.ZipFile(str(part/'indices.npz')) as archive:
            shape,counts = numeric_npy(archive,'gap_counts'); assert shape == (3,4)
            shape,minimum = numeric_npy(archive,'gap_min'); assert shape == (3,)
        record = gaps.setdefault(task['fit_seed'],dict(counts=[0]*12,minimum=[float('inf')]*3,hashes=[],n=0))
        record['counts'] = [a+b for a,b in zip(record['counts'],counts)]
        record['minimum'] = [min(a,b) for a,b in zip(record['minimum'],minimum)]
        record['hashes'].append(marker['indices_sha256']);record['n'] += 1
        proof_path = part/'recovery/numerical_recovery.json'
        if proof_path.exists():
            proof = json.loads(proof_path.read_text())
            assert proof['status'] == 'passed' and proof['no_matrix_regularization'] and proof['no_sampling_change']
            for name,digest in proof['artifact_sha256'].items(): assert sha(part/'recovery'/name) == digest, name
            for name,digest in proof['code_sha256'].items(): assert sha(code/name) == digest, name
            proofs.append(dict(task_id=task_id,fit_seed=task['fit_seed'],N=4097,focal_parameter=task['focal_parameter'],
                replicate=int(task['replicate']),proof_sha256=sha(proof_path),proof=proof,
                corrections=json.loads((part/'recovery/corrections.json').read_text())))
    gap_rows = []
    for seed in SEEDS:
        record = gaps[seed];assert record['n'] == 70
        chain = hashlib.sha256(('\n'.join(record['hashes'])+'\n').encode()).hexdigest()
        for j,(label,size) in enumerate([('low_0_to_1pct',41),('middle_1_to_3pct',79),('high_3_to_5pct',81)]):
            counts = record['counts'][4*j:4*j+4];assert counts[0] == 70*4097*size
            gap_rows.append(dict(fit_seed=seed,N=4097,oxygen_band=label,n_trajectories=70,n_evaluations=counts[0],
                fraction_gap_lt_1e_6=counts[1]/counts[0],fraction_gap_lt_1e_4=counts[2]/counts[0],
                fraction_gap_lt_1e_3=counts[3]/counts[0],minimum_spectral_gap=record['minimum'][j],source_cache_hash_chain_sha256=chain))
    write(root/'spectral_gap_summary.tsv',gap_rows)
    accounting = subprocess.check_output(['sacct','-X','-S','2026-10-14','-j',','.join(JOB_IDS.values()),
        '--format=JobID%48,State%32,ExitCode,ElapsedRaw,Start,End','-n','-P'],universal_newlines=True)
    entries = {}
    for line in accounting.splitlines():
        parts = [x.strip() for x in line.split('|')]
        if len(parts) < 6: continue
        row = dict(zip(('job_id','state','exit_code','elapsed_seconds','start','end'),parts))
        row['elapsed_seconds'] = int(row['elapsed_seconds'])
        if row['job_id'] not in entries or row['end'] >= entries[row['job_id']]['end']: entries[row['job_id']] = row
    for role,count in [('trajectory',420),('seed_summary',6)]:
        selected = {int(k.split('_')[1]) for k in entries if k.startswith(JOB_IDS[role]+'_')}
        assert selected == set(range(1,count+1)), role
    assert JOB_IDS['smoke'] in entries and JOB_IDS['finalize'] in entries
    assert all(r['state']=='COMPLETED' and r['exit_code']=='0:0' for r in entries.values())
    write(root/'slurm/completed_slurm_accounting.tsv',list(entries.values()))
    write(root/'slurm/completed_trajectory_audit.tsv',task_sources)
    for seed in SEEDS:
        receipt = json.loads((root/'summaries'/(seed+'_N4097.receipt.json')).read_text())
        assert receipt['status'] == 'complete' and receipt['n_trajectories'] == 70
        assert {int(r['task_id']) for r in receipt['sources']} == {int(t['task_id']) for t in tasks if t['fit_seed']==seed}
        for row in receipt['sources']:
            source = task_sources[int(row['task_id'])-1]
            assert row['indices_sha256'] == source['indices_sha256'] and row['receipt_sha256'] == source['output_receipt_sha256']
    report = dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        hpc_git_commit=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],universal_newlines=True).strip(),
        frozen_source_hashes_match=True,n_verified_trajectory_caches=420,n_recovery_trajectories=len(proofs),
        output_validation_scope='Original output receipts cross-linked to cached indices; full raw gzip files are not rehashed except recovery source artifacts',
        recovery_audits=proofs,slurm_unique_jobs=list(entries.values()),smoke_status='passed',
        completion_sha256=sha(root/'completion.json'),n_completion_artifacts_verified=len(completion['artifact_sha256']),
        source_collector_sha256=sha(Path(__file__)),new_jobs_submitted=False)
    (root/'slurm/completed_hpc_audit.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('recovery_audits','slurm_unique_jobs')},indent=2))


def manifest(root):
    completion = json.loads((root/'completion.json').read_text())
    names = set(completion['artifact_sha256'])
    names.update(('completion.json','baseline_indices.npz','seed_selection.tsv','artifact_collection_manifest.tsv',
                  'collection_receipt.json','preflight_validation.json','spectral_gap_summary.tsv'))
    for directory,suffixes in [('summaries',('.npz','.json','.tsv')),('slurm',('.json','.tsv')),('smoke',('.json',))]:
        names.update(str(p.relative_to(root)) for p in (root/directory).iterdir() if p.is_file() and p.suffix in suffixes)
    rows = []
    for name in sorted(names):
        relative = Path(name);assert not relative.is_absolute() and '..' not in relative.parts
        path = root/relative;assert path.is_file() and not path.is_symlink()
        rows.append(dict(path=name,sha256=sha(path),bytes=path.stat().st_size))
    write(root/'completed_results_manifest.tsv',rows)
    print('Manifest: {} files, {} bytes'.format(len(rows),sum(r['bytes'] for r in rows)))


def verify(root):
    rows = read(root/'completed_results_manifest.tsv')
    for row in rows:
        relative = Path(row['path']);assert not relative.is_absolute() and '..' not in relative.parts
        path = root/relative
        assert path.is_file() and not path.is_symlink() and path.stat().st_size == int(row['bytes'])
        assert sha(path) == row['sha256'], row['path']
        # Data artifacts have no executable bits, independent of rsync defaults.
        path.chmod(0o644)
    completion = json.loads((root/'completion.json').read_text())
    assert completion['status']=='complete' and completion['n_tasks']==420 and completion['n_seed_N_summaries']==6
    for name,digest in completion['artifact_sha256'].items(): assert sha(root/name) == digest, name
    record = dict(collected_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed',
        n_hpc_artifacts=len(rows),total_bytes=sum(int(r['bytes']) for r in rows),all_hpc_sha256_matched=True,
        original_completion_hashes_matched=True,manifest_sha256=sha(root/'completed_results_manifest.tsv'),
        raw_trajectories_location='HPC only; cached indices, receipts and proofs audited on HPC',
        analysis_complete=True,new_jobs_submitted=False)
    (root/'results_collection_receipt.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(json.dumps(record,indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('audit','manifest','verify-local'))
    parser.add_argument('--repo-root',type=Path)
    parser.add_argument('--root',type=Path)
    args = parser.parse_args()
    if args.command == 'audit': audit(args.repo_root)
    elif args.command == 'manifest': manifest(args.root)
    else: verify(args.root)


if __name__ == '__main__': main()
