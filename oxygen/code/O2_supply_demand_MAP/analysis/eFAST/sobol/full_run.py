#!/usr/bin/env python3
"""Prepare and execute the complete fixed-O2 in-vivo second-order Sobol run."""
import argparse
import csv
import fcntl
from importlib.metadata import version
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

import numpy as np
import core as c
import full_core as f


def prepare(args):
    if version('SALib')!='1.5.2': raise ValueError('Require SALib 1.5.2')
    root=Path(args.out_dir).resolve(); root.mkdir(parents=True,exist_ok=True)
    if (root/'slurm/submission_jobs.tsv').exists(): raise ValueError('Already submitted')
    if shutil.disk_usage(root).free<1.35*1024**4: raise ValueError('Require 1.35 TiB free for compact outputs and reserve')
    old=Path(args.efast_root).resolve()/'neighborhood10pct'; fit=Path(args.fit_root).resolve()
    audit=c.load_json(old/'input_manifest.json')
    for name,field in [('seed_manifest.tsv','seed_manifest_sha256'),('neighborhood_bounds.tsv','neighborhood_bounds_sha256'),
                       ('operator_input_audit.json','operator_audit_sha256')]:
        if c.sha(old/name)!=audit[field]: raise ValueError('Existing input audit changed: '+name)
    manifest=c.read(old/'seed_manifest.tsv')
    if len(manifest)!=500 or {r['fit_seed'] for r in manifest}!={'seed%d'%i for i in range(1,501)}:
        raise ValueError('Require all 500 independent in-vivo endpoints')
    _,ranges=c.efast.get_ranges(fit/'seed25')
    inputs=root/'inputs'; inputs.mkdir(exist_ok=True); (root/'slurm').mkdir(exist_ok=True)
    oxygen,reference=c.efast.oxygen_grid(args.figure4_dir,True,'')
    if c.sha(reference)!=audit['reference_grid_sha256']: raise ValueError('Figure 4 reference changed')
    shutil.copyfile(reference,inputs/reference.name)
    c.table(inputs/'seed_manifest.tsv',manifest); c.table(inputs/'global_bounds.tsv',ranges)
    c.save_json(inputs/'source_audit.json',dict(source=audit,fit_root=str(fit),source_root=str(old)))
    common={r['parameter']:float(10**((r['encoded_lower']+r['encoded_upper'])/2) if r['transform']=='log10'
            else (r['encoded_lower']+r['encoded_upper'])/2) for r in ranges}
    units=[]; unit_sha={}
    unit_ranges=[dict(parameter=p,encoded_lower=0.,encoded_upper=1.,natural_lower=0.,natural_upper=1.,transform='identity') for p in c.ACTIVE]
    for repeat in range(1,6):
        unit=c.samples(unit_ranges,f.N,2026101300+repeat); units.append(unit)
        target=inputs/('unit_R%02d.npz'%repeat); f.atomic_npz(target,unit=unit)
        unit_sha[str(repeat)]=c.sha(target)
    rows={r['fit_seed']:r for r in manifest}; cases=[]
    definitions=[dict(scope='global',case_id='global',fit_seed='seed25')]+[
        dict(scope='neighborhood10pct',case_id='seed%d'%i,fit_seed='seed%d'%i) for i in range(1,501)]
    for number,case in enumerate(definitions,1):
        seed=case['fit_seed']; endpoint=rows[seed]
        for name,field in [('best_params.tsv','best_params_sha256'),('fit_config.rds','fit_config_sha256'),('parameter_table.csv','parameter_table_sha256')]:
            if c.sha(fit/seed/name)!=endpoint[field]: raise ValueError('Endpoint changed: '+seed+'/'+name)
        _,rr=c.efast.get_ranges(fit/seed)
        if rr!=ranges: raise ValueError('Endpoint bounds differ: '+seed)
        best={r['parameter']:float(r['value']) for r in c.read(fit/seed/'best_params.tsv')}
        bounds=c.bounds_for(ranges,None if case['scope']=='global' else best)
        target=f.case_dir(root,case); target.mkdir(parents=True,exist_ok=True)
        c.table(target/'parameter_bounds.tsv',bounds)
        c.table(target/'reference_best_params.tsv',[dict(parameter=p,value=best[p]) for p in c.ACTIVE])
        # Prove shared-unit transport matches SALib's natural design bit for bit.
        if case['scope']=='global' or seed in c.PILOT_SEEDS:
            for repeat in range(1,6):
                transformed=f.scale_unit(units[repeat-1],bounds)
                direct=c.samples(bounds,f.N,2026101300+repeat)
                if not np.array_equal(transformed,direct): raise ValueError('Transport differs from canonical SALib design')
                if not np.array_equal(direct[:2048*c.STEP],c.samples(bounds,2048,2026101300+repeat)):
                    raise ValueError('Nested diagnostic prefix differs')
        record=dict(case_index=number,**case,bounds_sha256=c.sha(target/'parameter_bounds.tsv'),
                    best_params_sha256=endpoint['best_params_sha256'],fit_config_sha256=endpoint['fit_config_sha256'],
                    parameter_table_sha256=endpoint['parameter_table_sha256'],included_in_full_analysis=True)
        c.save_json(target/'case.json',record); cases.append(record)
        if number%50==0: print('Audited/prepared cases',number,'/501',flush=True)
    c.table(root/'slurm/case_manifest.tsv',cases)
    def tasks():
        for task_id in range(1,f.N_TASKS+1):
            case_index,repeat,block=f.decode(task_id); case=cases[case_index-1]
            yield dict(task_id=task_id,case_index=case_index,scope=case['scope'],case_id=case['case_id'],
                       fit_seed=case['fit_seed'],repeat=repeat,block=c.BLOCKS[block],block_index=block)
    c.efast.write_table(root/'slurm/task_manifest.tsv',['task_id','case_index','scope','case_id','fit_seed','repeat','block','block_index'],tasks())
    sif=Path(args.sif).resolve()
    if c.sha(sif)!=c.SIF_SHA256: raise ValueError('SIF digest differs')
    template=old/'slurm/rcpp_template.tar.gz'; tpl=c.load_json(old/'slurm/rcpp_template_manifest.json')
    model=f.HERE.parents[2]/'model'
    if c.sha(template)!=tpl['archive_sha256'] or tpl['sif_sha256']!=c.SIF_SHA256 or tpl['model_cpp_sha256']!=c.sha(model/'model_O2_supply_demand_MAP.cpp'):
        raise ValueError('Runtime template differs')
    shutil.copyfile(template,root/'slurm/rcpp_template.tar.gz')
    initial=dict(status='initial',archive_sha256=c.sha(root/'slurm/rcpp_template.tar.gz'))
    c.save_json(root/'slurm/runtime_cache_manifest.json',initial)
    snapshot=Path(args.snapshot_root).resolve()
    file_hashes={str(p.relative_to(snapshot)):c.sha(p) for p in snapshot.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    plan=dict(method='Sobol_second_order',stage='full',N=f.N,resolutions=f.NS,oxygen_pct=oxygen,
        cases=cases,n_fit_seeds=500,n_cases=501,n_parameters=14,n_pairs=91,n_repeats=5,n_tasks=f.N_TASKS,
        n_operator_evaluations=f.N_TASKS*f.N*len(oxygen),n_response_tasks=54,
        fit_root=str(fit),code_snapshot_root=str(snapshot),git_commit=args.git_commit,
        code_sha256=file_hashes,engine_sha256=c.sha(f.HERE/'evaluate.R'),unit_sha256=unit_sha,
        source_audit_sha256=c.sha(old/'input_manifest.json'),reference_sha256=c.sha(inputs/reference.name),
        common_reference_params=common,sif=dict(path=str(sif),sha256=c.SIF_SHA256,bytes=sif.stat().st_size,mtime_ns=sif.stat().st_mtime_ns),
        checkpoint_oxygen_chunk=8,output_format='lossless float64 NPZ; per-point solver flags; per-oxygen eigenpair QC; exception matrices',
        convergence_policy='diagnostic_only; all 500 neighborhoods included',repeat_range_tolerance=.05,resolution_delta_tolerance=.02,
        bootstrap_resamples=200,qos='xxlarge',time='12:00:00',array_concurrency_limit=None,specified_compute_node=None,
        free_bytes_at_prepare=shutil.disk_usage(root).free,prepared_utc=c.now())
    c.save_json(root/'slurm/submission_plan.json',plan)
    c.table(root/'environment.tsv',[dict(component=k,version=version(k)) for k in ['SALib','numpy','scipy','pandas','matplotlib']]+[
        dict(component='python',version=sys.version.split()[0])])
    from structure import write_structure
    write_structure(root/'model_structure')
    c.save_json(root/'status.json',dict(status='prepared',n_cases=501,n_neighborhoods=500,n_evaluation_tasks=f.N_TASKS,
        N=f.N,repeats=5,oxygen_points=201,convergence_policy='diagnostic_only',prepared_utc=c.now()))
    print(json.dumps({k:plan[k] for k in ['n_tasks','n_operator_evaluations','N','n_cases']},indent=2),flush=True)


def evaluate_block(root,p,case,repeat,block,target,n=f.N,validate_reference=False):
    target=Path(target); target.mkdir(parents=True,exist_ok=True)
    if c.sha(f.case_dir(root,case)/'parameter_bounds.tsv')!=case['bounds_sha256']: raise ValueError('Bounds changed')
    x=f.samples(root,case,repeat,block)[:n]
    key=f.identity(p,case,repeat,block,x)
    output=target/'outputs.npz'
    with open(target/'worker.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if f.complete(output,key): return c.load_json(str(output)+'.receipt.json')
        oxygen=np.asarray(p['oxygen_pct'])
        values,flags,done,gaps,residuals,elapsed=f.checkpoint(target/'checkpoint.npz',key,n,oxygen)
        scratch=Path(os.environ['SOBOL_SCRATCH'])/'evaluation'
        scratch.mkdir(parents=True,exist_ok=True)
        c.table(scratch/'samples.tsv.gz',[dict(sample_id=i+1,**dict(zip(c.ACTIVE,map(float,row)))) for i,row in enumerate(x)])
        checks=[]; common=[]
        if (target/'checkpoint.npz').exists():
            with np.load(target/'checkpoint.npz') as z:
                checks=json.loads(str(z['reference_json'])); common=json.loads(str(z['common_json']))
        while not done.all():
            indices=np.flatnonzero(~done)[:p['checkpoint_oxygen_chunk']]
            current=oxygen[indices].tolist()
            meta=dict(method='Sobol_second_order',scope=case['scope'],case_id=case['case_id'],fit_seed=case['fit_seed'],
                n_samples=n,parameters=list(c.ACTIVE),oxygen_pct=current,samples_sha256=c.sha(scratch/'samples.tsv.gz'),
                best_params_sha256=case['best_params_sha256'],fit_config_sha256=case['fit_config_sha256'],
                validate_reference=bool(validate_reference),common_reference_params=p['common_reference_params'])
            c.save_json(scratch/'metadata.json',meta)
            subprocess.run(['Rscript','--vanilla',str(f.HERE/'evaluate.R'),str(scratch),p['fit_root'],
                            str(Path(root)/'inputs/fixed_o2_dominant_ploidy_201grid.tsv')],check=True)
            v,s,g,r=f.read_model(scratch/'outputs.tsv.gz',n,current)
            rr=c.load_json(str(scratch/'outputs.tsv.gz')+'.receipt.json')
            values[:,indices]=v; flags[:,indices]=s; gaps[indices]=g; residuals[indices]=r; done[indices]=True
            elapsed+=rr['elapsed_seconds']
            checks.extend(rr.get('reference_checks') or []); common.extend(rr.get('common_reference') or [])
            if (scratch/'recovery').exists():
                shutil.copytree(scratch/'recovery',target/'recovery',dirs_exist_ok=True)
                shutil.rmtree(scratch/'recovery')
            # All checkpoint fields are one atomically promoted, CRC-checked NPZ.
            f.atomic_npz(target/'checkpoint.npz',identity=np.array(key),oxygen=oxygen,outputs=values,solver_flags=flags,
                         done=done,minimum_gap=gaps,maximum_residual=residuals,elapsed_seconds=np.array(elapsed),
                         reference_json=np.array(json.dumps(checks)),common_json=np.array(json.dumps(common)))
            print(case['case_id'],repeat,c.BLOCKS[block],int(done.sum()),'/201 oxygen complete',flush=True)
        f.atomic_npz(output,identity=np.array(key),oxygen=oxygen,outputs=values,solver_flags=flags,
                     minimum_gap=gaps,maximum_residual=residuals)
        receipt=dict(status='complete',identity=key,sha256=c.sha(output),n_points=n*len(oxygen),
                     n_samples=n,oxygen_points=len(oxygen),n_recovered=int((flags==2).sum()),
                     minimum_spectral_gap=float(gaps.min()),maximum_residual=float(residuals.max()),
                     elapsed_seconds=elapsed,reference_checks=checks,common_reference=common,
                     completed_utc=c.now(),host=os.uname().nodename,job_id=os.environ.get('SLURM_JOB_ID',''))
        c.save_json(str(output)+'.receipt.json',receipt)
        (target/'checkpoint.npz').unlink()
        return receipt


def bootstrap(root,p):
    subprocess.run([sys.executable,str(f.HERE/'test_sobol.py')],check=True)
    subprocess.run([sys.executable,str(f.HERE/'test_full.py')],check=True)
    target=Path(root)/'slurm/runtime_validation/bootstrap'
    evaluate_block(root,p,p['cases'][0],1,0,target,n=8,validate_reference=True)
    cache=Path(os.environ['SOBOL_MODEL_CACHE'])
    archive=Path(root)/'slurm/rcpp_template.tar.gz'
    tmp=archive.with_name(archive.name+'.incomplete')
    with tarfile.open(tmp,'w:gz') as tar:
        for item in cache.iterdir(): tar.add(item,arcname=item.name)
    tmp.replace(archive)
    c.save_json(Path(root)/'slurm/runtime_cache_manifest.json',dict(status='validated_snapshot_cache',
        archive_sha256=c.sha(archive),engine_sha256=p['engine_sha256'],code_git_commit=p['git_commit'],completed_utc=c.now()))


def validate(root,p,case_id):
    case=p['cases'][case_id-1]; target=Path(root)/'slurm/runtime_validation'/case['scope']/case['case_id']
    rr=evaluate_block(root,p,case,1,0,target,n=16,validate_reference=True)
    if len(rr['common_reference'])!=201 or len(rr['reference_checks'])!=3*26:
        raise ValueError('Missing full-grid configuration/reference audit')
    # Round-trip storage must preserve exactly the canonical CSV-parsed values.
    with np.load(target/'outputs.npz') as z:
        if z['outputs'].dtype!=np.float64 or not np.isfinite(z['outputs']).all(): raise ValueError('Output precision differs')
    c.save_json(target/'validation.json',dict(status='passed',case_index=case_id,case_id=case['case_id'],
        scope=case['scope'],output_sha256=rr['sha256'],reference_checks=len(rr['reference_checks']),common_points=201,
        n_points=rr['n_points'],bytes=(target/'outputs.npz').stat().st_size,completed_utc=c.now()))


def gate(root,p):
    for name,expected in p['code_sha256'].items():
        if c.sha(f.REPO/name)!=expected: raise ValueError('Frozen source changed: '+name)
    reference=None; measured=[]
    for case in p['cases']:
        target=Path(root)/'slurm/runtime_validation'/case['scope']/case['case_id']
        v=c.load_json(target/'validation.json'); rr=c.load_json(str(target/'outputs.npz')+'.receipt.json')
        if v['status']!='passed' or not f.complete(target/'outputs.npz'): raise ValueError('Configuration validation incomplete')
        common=np.array([[r[k] for k in c.OUTPUTS] for r in rr['common_reference']])
        if reference is None: reference=common
        if not np.isfinite(common).all() or np.max(abs(common-reference))>1e-8: raise ValueError('Fitted operator configurations differ')
        measured.append(v['bytes']/v['n_points'])
    free=shutil.disk_usage(root).free
    if free<1.35*1024**4: raise ValueError('Insufficient compact-output storage reserve')
    c.save_json(Path(root)/'runtime_validation.json',dict(status='passed',cases_checked=501,oxygen_points=201,
        common_configuration_maximum_tolerance=1e-8,numerical_tests='test_sobol.py and test_full.py passed in bootstrap',
        measured_N16_bytes_per_point_p50=float(np.median(measured)),float64_outputs_and_flags_uncompressed_GiB=p['n_operator_evaluations']*17/1024**3,
        free_bytes_at_gate=free,convergence_policy='diagnostic_only',validated_utc=c.now()))
    c.save_json(Path(root)/'status.json',dict(status='validated_ready_for_full_array',n_cases=501,n_tasks=f.N_TASKS,
        all_seeds_included=True,convergence_policy='diagnostic_only',updated_utc=c.now()))


def worker(root,p,task_id):
    ci,repeat,block=f.decode(task_id); case=p['cases'][ci-1]
    target=f.case_dir(root,case)/'runs'/('R%02d'%repeat)/'blocks'/c.BLOCKS[block]
    evaluate_block(root,p,case,repeat,block,target)


def summarize(root,p,case_id):
    from full_report import summarize_case
    summarize_case(Path(root),p,p['cases'][case_id-1])


def aggregate(root,p):
    from full_report import aggregate_cases
    aggregate_cases(Path(root),p)


def response(root,p,task_id):
    tasks=c.read(Path(root)/'slurm/response_tasks.tsv'); task=tasks[task_id-1]
    from run import evaluate
    evaluate(Path(root),Path(root)/task['part'],p)


def finish(root,p):
    from full_report import finish_report
    finish_report(Path(root),p)


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('command',choices=['prepare','bootstrap','validate','gate','worker','summarize','aggregate','response','finish'])
    a.add_argument('--out-dir',required=True); a.add_argument('--task-id',type=int); a.add_argument('--case-id',type=int)
    for name in ['fit-root','efast-root','figure4-dir','sif','snapshot-root','git-commit']: a.add_argument('--'+name)
    args=a.parse_args(); root=Path(args.out_dir).resolve()
    if args.command=='prepare': prepare(args); return
    p=f.plan(root)
    try:
        if args.command in ('validate','summarize'): globals()[args.command](root,p,args.case_id)
        elif args.command in ('worker','response'): globals()[args.command](root,p,args.task_id)
        else: globals()[args.command](root,p)
    except Exception as error:
        target=root/'diagnostics/failures'/('%s_%s.json'%(args.command,args.task_id or args.case_id or 'main'))
        c.save_json(target,dict(status='failed',stage=args.command,error=str(error),recorded_utc=c.now()))
        raise


if __name__=='__main__': main()
