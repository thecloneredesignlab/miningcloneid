#!/usr/bin/env python3
"""Audit/collect completed N8193/M4 and N4097/M8 pilots; never submit jobs.

Python 3.6+ standard library. Full raw output gzip files stay on HPC and are
not rehashed except original files used by certified numerical recovery.
"""
import argparse
import collections
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

from collect_resolution_4097_results import sha, read, write, numeric_npy

PILOT = 'N8193_M4_and_N4097_M8'
CONFIGS = {'N8193_M4':(8193,4,('seed25','seed138'),'N_resolution'),
           'N4097_M8':(4097,8,('seed101','seed174'),'M_robustness')}
PARAMETERS = ('lam_max','p_mis_base','p_misseg','k_o_mis','buffer_smax','buffer_beta',
              'buffer_n_exp','p_wgd','alpha_o2','gamma_growth','mu_hp','gamma_mu','O2_crit','n_O')


def validate_tasks(tasks, config):
    n,m,seeds,_ = CONFIGS[config]
    assert len(tasks)==140
    for i,row in enumerate(tasks):
        seed=seeds[i//70];r=(i%70)//14+1;p=i%14
        expected=dict(task_id=i+1,summary_id=i//70+1,config=config,fit_seed=seed,
            fit_seed_number=int(seed[4:]),N=n,M=m,replicate=r,parameter_index=p,
            focal_parameter=PARAMETERS[p],phase_seed=202610080+int(seed[4:])*10+r-1,
            parent_sample_id_start=p*n+1,parent_sample_id_end=(p+1)*n)
        assert row=={k:str(v) for k,v in expected.items()},row


def audit(repo):
    code=repo/'oxygen/code/O2_supply_demand_MAP/analysis/eFAST'
    master=repo/'oxygen/results/eFAST/neighborhood10pct/pilot'/PILOT
    for config,(n,m,seeds,kind) in CONFIGS.items():
        root=master/config
        plan=json.loads((root/'slurm/submission_plan.json').read_text())
        assert (plan['n_tasks'],plan['n_summaries'],plan['M'],plan['resolutions'],plan['fit_seeds'])==(140,2,m,[n],list(seeds))
        assert plan['comparison_kind']==kind
        for name,digest in plan['code_sha256'].items():assert sha(code/name)==digest,name
        for name,digest in plan['comparison_source_sha256'].items():assert sha(Path(name))==digest,name
        for name,field in [('seed_selection.tsv','seed_selection_sha256'),('baseline_indices.npz','comparison_indices_sha256'),
                           ('slurm/task_manifest.tsv','task_manifest_sha256'),('slurm/summary_manifest.tsv','summary_manifest_sha256')]:
            assert sha(root/name)==plan[field],name
        for endpoint in read(root/'seed_selection.tsv'):
            for name,path in [('best_params',Path(plan['fit_root'])/endpoint['fit_seed']/'best_params.tsv'),
                              ('fit_config',Path(plan['fit_root'])/endpoint['fit_seed']/'fit_config.rds'),
                              ('bounds',Path(plan['baseline_dir'])/'bounds'/(endpoint['fit_seed']+'.tsv'))]:
                assert sha(path)==endpoint[name+'_sha256']
        smoke=json.loads((root/'smoke/smoke_receipt.json').read_text())
        assert smoke['status']=='passed' and smoke['config']==config
        assert smoke['submission_plan_sha256']==sha(root/'slurm/submission_plan.json')
        assert smoke['perron_validation']['status']=='passed'
        complete=json.loads((root/'completion.json').read_text())
        assert (complete['status'],complete['n_tasks'],complete['n_seed_summaries'],complete['comparison_kind'])==('complete',140,2,kind)
        for name,digest in complete['artifact_sha256'].items():assert sha(root/name)==digest,name
        tasks=read(root/'slurm/task_manifest.tsv');validate_tasks(tasks,config)
        gaps={};proofs=[];sources=[];parents=set()
        for task in tasks:
            seed=task['fit_seed'];p=int(task['parameter_index']);task_id=int(task['task_id'])
            marker_path=root/'slurm/task_receipts'/('task_%05d.json'%task_id)
            marker=json.loads(marker_path.read_text())
            assert marker['status']=='complete'
            for field,value in task.items():assert marker[field]==value
            parent=root/'runs'/seed/config/'runs'/('N%d_R%s'%(n,task['replicate']))
            parent_meta=json.loads((parent/'metadata.json').read_text())
            if parent not in parents:
                assert parent_meta['N']==n and parent_meta['M']==m and parent_meta['n_samples']==14*n
                assert parent_meta['n_parameters']==14 and sha(parent/'samples.tsv.gz')==parent_meta['samples_sha256']
                parents.add(parent)
            part=parent/'trajectories'/('P%02d_%s'%(p+1,PARAMETERS[p]))
            meta=json.loads((part/'metadata.json').read_text())
            cache=json.loads((part/'indices.receipt.json').read_text())
            receipt_path=part/'outputs.tsv.gz.receipt.json';receipt=json.loads(receipt_path.read_text())
            assert (meta['N'],meta['M'],meta['fit_seed'],meta['seed'],meta['replicate'],meta['trajectory_parameter_index'])==(n,m,seed,int(task['phase_seed']),int(task['replicate']),p)
            assert meta['n_samples']==n and len(meta['oxygen_pct'])==201
            assert all(abs(o-j*.025)<1e-12 for j,o in enumerate(meta['oxygen_pct']))
            assert meta['parent_metadata_sha256']==sha(parent/'metadata.json')
            assert meta['parent_samples_sha256']==parent_meta['samples_sha256']
            assert sha(part/'samples.tsv.gz')==meta['samples_sha256']
            assert sha(part/'indices.npz')==marker['indices_sha256']==cache['cache_sha256']
            assert sha(receipt_path)==marker['trajectory_receipt_sha256']
            assert cache['metadata_sha256']==receipt['metadata_sha256']==sha(part/'metadata.json')
            assert cache['analyzer_sha256']==sha(code/'resolution_nm.py') and cache['salib_version']=='1.5.2'
            assert receipt['samples_sha256']==meta['samples_sha256'] and receipt['n_rows']==n*201
            assert receipt['outputs_sha256']==cache['outputs_sha256'] and receipt['evaluator_sha256']==sha(code/'evaluate_fixed_o2.R')
            assert (part/'outputs.tsv.gz').stat().st_size>0
            sources.append(dict(task_id=task_id,config=config,fit_seed=seed,replicate=int(task['replicate']),
                parameter=task['focal_parameter'],N=n,M=m,marker_sha256=sha(marker_path),
                indices_sha256=marker['indices_sha256'],output_receipt_sha256=marker['trajectory_receipt_sha256']))
            with zipfile.ZipFile(str(part/'indices.npz')) as archive:
                shape,counts=numeric_npy(archive,'gap_counts');assert shape==(3,4)
                shape,minimum=numeric_npy(archive,'gap_min');assert shape==(3,)
            record=gaps.setdefault(seed,dict(counts=[0]*12,minimum=[float('inf')]*3,hashes=[]))
            record['counts']=[a+b for a,b in zip(record['counts'],counts)]
            record['minimum']=[min(a,b) for a,b in zip(record['minimum'],minimum)]
            record['hashes'].append(marker['indices_sha256'])
            proof_path=part/'recovery/numerical_recovery.json'
            if proof_path.exists():
                proof=json.loads(proof_path.read_text())
                assert proof['status']=='passed' and proof['no_matrix_regularization'] and proof['no_sampling_change']
                assert receipt['recovery_manifest_sha256']==sha(proof_path)
                assert proof['outputs_sha256']==receipt['outputs_sha256']
                for name,digest in proof['artifact_sha256'].items():assert sha(part/'recovery'/name)==digest,name
                for name,digest in proof['code_sha256'].items():assert sha(code/name)==digest,name
                corrections=json.loads((part/'recovery/corrections.json').read_text())
                assert corrections['status']=='passed' and corrections['n_points']==proof['n_recomputed_points']
                for row in corrections['reports']:
                    assert row['double_residual']<1e-12 and row['ploidy_precision_difference']<1e-12
                    assert row['growth_precision_difference']<1e-12 and row['vector_precision_difference']<1e-10
                proofs.append(dict(task_id=task_id,fit_seed=seed,N=n,M=m,config=config,
                    focal_parameter=task['focal_parameter'],replicate=int(task['replicate']),
                    proof_sha256=sha(proof_path),proof=proof,corrections=corrections))
        gap_rows=[]
        for seed in seeds:
            item=gaps[seed];assert len(item['hashes'])==70
            chain=hashlib.sha256(('\n'.join(item['hashes'])+'\n').encode()).hexdigest()
            for j,(label,size) in enumerate([('low_0_to_1pct',41),('middle_1_to_3pct',79),('high_3_to_5pct',81)]):
                counts=item['counts'][j*4:j*4+4];assert counts[0]==70*n*size
                gap_rows.append(dict(config=config,fit_seed=seed,N=n,M=m,oxygen_band=label,n_trajectories=70,
                    n_evaluations=counts[0],fraction_gap_lt_1e_6=counts[1]/counts[0],
                    fraction_gap_lt_1e_4=counts[2]/counts[0],fraction_gap_lt_1e_3=counts[3]/counts[0],
                    minimum_spectral_gap=item['minimum'][j],source_cache_hash_chain_sha256=chain))
        write(root/'spectral_gap_summary.tsv',gap_rows)
        jobs=read(root/'slurm/submission_jobs.tsv');ids={r['role']:r['job_id'] for r in jobs}
        raw=subprocess.check_output(['sacct','-X','-S','2026-10-18','-j',','.join(ids.values()),
            '--format=JobID%48,State%32,ExitCode,ElapsedRaw,Start,End','-n','-P'],universal_newlines=True)
        entries={}
        for line in raw.splitlines():
            fields=[f.strip() for f in line.split('|')]
            if len(fields)<6:continue
            row=dict(zip(('job_id','state','exit_code','elapsed_seconds','start','end'),fields))
            row['elapsed_seconds']=int(row['elapsed_seconds'])
            if row['job_id'] not in entries or row['end']>=entries[row['job_id']]['end']:entries[row['job_id']]=row
        for role,count in [('trajectory',140),('seed_summary',2)]:
            assert {int(k.split('_')[1]) for k in entries if k.startswith(ids[role]+'_')}==set(range(1,count+1))
        assert ids['smoke'] in entries and ids['finalize'] in entries
        assert all(r['state']=='COMPLETED' and r['exit_code']=='0:0' for r in entries.values())
        write(root/'slurm/completed_slurm_accounting.tsv',list(entries.values()))
        write(root/'slurm/completed_trajectory_audit.tsv',sources)
        for seed in seeds:
            prefix=root/'summaries'/('%s_N%d_M%d'%(seed,n,m))
            summary=json.loads(prefix.with_suffix('.receipt.json').read_text())
            assert summary['status']=='complete' and summary['n_trajectories']==70 and summary['M']==m
            assert sha(prefix.with_suffix('.npz'))==summary['indices_sha256']
            assert sha(prefix.with_name(prefix.name+'_repeat_diagnostics.tsv'))==summary['diagnostics_sha256']
            assert {int(r['task_id']) for r in summary['sources']}=={int(t['task_id']) for t in tasks if t['fit_seed']==seed}
            for row in summary['sources']:
                expected=sources[int(row['task_id'])-1]
                assert row['indices_sha256']==expected['indices_sha256'] and row['receipt_sha256']==expected['output_receipt_sha256']
        report=dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),config=config,
            hpc_git_commit=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],universal_newlines=True).strip(),
            frozen_source_hashes_match=True,n_verified_trajectory_caches=140,n_recovery_trajectories=len(proofs),
            n_recovered_points=sum(p['proof']['n_recomputed_points'] for p in proofs),
            output_validation_scope='Sample and cache SHA256 checked; output receipts cross-linked; full raw gzip rehashed only for recovery sources',
            recovery_audits=proofs,slurm_unique_jobs=list(entries.values()),smoke_status='passed',
            completion_sha256=sha(root/'completion.json'),n_completion_artifacts_verified=len(complete['artifact_sha256']),
            source_collector_sha256=sha(Path(__file__)),new_jobs_submitted=False)
        (root/'slurm/completed_hpc_audit.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
        print(json.dumps({k:v for k,v in report.items() if k not in ('recovery_audits','slurm_unique_jobs')},indent=2),flush=True)


def manifest(master):
    names=set()
    for config in CONFIGS:
        root=master/config
        complete=json.loads((root/'completion.json').read_text())
        local=set(complete['artifact_sha256']);local.update(('completion.json','spectral_gap_summary.tsv'))
        for directory,suffixes in [('summaries',('.npz','.json','.tsv')),('slurm',('.json','.tsv')),('smoke',('.json',))]:
            local.update(str(p.relative_to(root)) for p in (root/directory).iterdir() if p.is_file() and p.suffix in suffixes)
        names.update(config+'/'+p for p in local)
    rows=[]
    for name in sorted(names):
        path=master/name;assert path.is_file() and not path.is_symlink()
        rows.append(dict(path=name,sha256=sha(path),bytes=path.stat().st_size))
    write(master/'completed_results_manifest.tsv',rows)
    print('Manifest: %d files, %d bytes'%(len(rows),sum(r['bytes'] for r in rows)),flush=True)


def verify(master):
    rows=read(master/'completed_results_manifest.tsv')
    for row in rows:
        relative=Path(row['path']);assert not relative.is_absolute() and '..' not in relative.parts
        path=master/relative
        assert path.is_file() and not path.is_symlink() and path.stat().st_size==int(row['bytes'])
        assert sha(path)==row['sha256'],row['path']
        path.chmod(0o644)
    for config in CONFIGS:
        complete=json.loads((master/config/'completion.json').read_text())
        for name,digest in complete['artifact_sha256'].items():assert sha(master/config/name)==digest,name
    receipt=dict(status='passed',collected_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        n_hpc_artifacts=len(rows),total_bytes=sum(int(r['bytes']) for r in rows),all_hpc_sha256_matched=True,
        original_completion_hashes_matched=True,manifest_sha256=sha(master/'completed_results_manifest.tsv'),
        raw_trajectories_location='HPC only',new_jobs_submitted=False)
    (master/'results_collection_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('audit','manifest','verify-local'))
    parser.add_argument('--repo-root',type=Path)
    parser.add_argument('--root',type=Path)
    args=parser.parse_args()
    if args.command=='audit':audit(args.repo_root.resolve())
    elif args.command=='manifest':manifest(args.root.resolve())
    else:verify(args.root.resolve())


if __name__=='__main__':main()
