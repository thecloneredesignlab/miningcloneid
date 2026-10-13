#!/usr/bin/env python3
"""Prepare, evaluate and collect separate global/neighborhood Sobol pilots."""
import argparse
import fcntl
from importlib.metadata import version
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import core as c

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[5]


def file_hashes():
    files = list(HERE.glob('*.py')) + list(HERE.glob('*.R')) + list(HERE.glob('*.sh')) + list(HERE.glob('*.sbatch'))
    files += [HERE.parent / name for name in ('efast.py', 'neighborhood.py', 'perron_high_precision_safe.cpp')]
    workflow = HERE.parents[2]
    files += list((workflow / 'util').glob('*.R'))
    files += list((workflow / 'model').glob('model_O2_supply_demand_MAP.*'))
    files += [workflow / 'simulation/o2/fixed_o2/run_fixed_o2_simulation.R']
    return {str(p.relative_to(REPO)): c.sha(p) for p in sorted(set(files))}


def context(root):
    plan = c.load_json(root / 'slurm/submission_plan.json')
    for name, expected in plan['code_sha256'].items():
        if c.sha(REPO / name) != expected:
            raise ValueError('Submitted source changed: ' + name)
    for name, expected in plan['input_sha256'].items():
        if c.sha(root / name) != expected:
            raise ValueError('Submitted input changed: ' + name)
    for name in ('sif', 'template'):
        p = Path(plan[name]['path'])
        if p.stat().st_size != plan[name]['bytes'] or p.stat().st_mtime_ns != plan[name]['mtime_ns']:
            raise ValueError('Verified immutable runtime changed: ' + name)
    c.configure(plan)
    return plan


def prepare(args):
    if version('SALib') != '1.5.2':
        raise ValueError('Production workflow requires SALib 1.5.2')
    root = Path(args.out_dir).resolve()
    c.configure(dict(resolutions=args.resolutions, oxygen_pct=c.OXYGEN))
    if (root / 'slurm/submission_jobs.tsv').exists():
        raise ValueError('Submission already recorded; inspect before preparing another attempt')
    old = Path(args.efast_root).resolve() / 'neighborhood10pct'
    fit = Path(args.fit_root).resolve()
    source = c.load_json(old / 'input_manifest.json')
    for name, field in [('seed_manifest.tsv','seed_manifest_sha256'),
                        ('neighborhood_bounds.tsv','neighborhood_bounds_sha256'),
                        ('operator_input_audit.json','operator_audit_sha256')]:
        if c.sha(old / name) != source[field]:
            raise ValueError('Existing endpoint audit changed: ' + name)
    if source['n_fit_seeds'] != 500 or source['n_parameters'] != 14:
        raise ValueError('Expected audited independent in-vivo endpoints')
    manifest = {r['fit_seed']: r for r in c.read(old / 'seed_manifest.tsv')}
    _, ranges = c.efast.get_ranges(fit / 'seed25')
    root.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(root).free < 20 * 1024**3:
        raise ValueError('Require at least 20 GiB free for pilot outputs')
    inputs = root / 'inputs'; inputs.mkdir(exist_ok=True)
    selected = []
    for seed in c.PILOT_SEEDS:
        r = manifest[seed]
        for name, field in [('best_params.tsv','best_params_sha256'),('fit_config.rds','fit_config_sha256'),
                            ('parameter_table.csv','parameter_table_sha256')]:
            if c.sha(fit / seed / name) != r[field]:
                raise ValueError('Audited endpoint changed: ' + seed + '/' + name)
        _, seed_ranges = c.efast.get_ranges(fit / seed)
        if seed_ranges != ranges:
            raise ValueError('Endpoint parameter ranges differ')
        selected.append(r)
    c.table(inputs / 'representative_seeds.tsv', selected)
    c.table(inputs / 'global_bounds.tsv', ranges)
    reference = Path(args.figure4_dir).resolve() / 'fixed_o2_dominant_ploidy_201grid.tsv'
    if c.sha(reference) != source['reference_grid_sha256']:
        raise ValueError('Figure 4 reference changed')
    shutil.copyfile(reference, inputs / reference.name)
    c.save_json(inputs / 'source_audit.json', dict(source=source, source_root=str(old),
        selection='best, worst, most clipped, then farthest-point coverage in normalized encoded parameter space; first eight audited pilot endpoints',
        source_manifest_sha256=c.sha(old/'input_manifest.json')))
    sif = Path(args.sif).resolve()
    if c.sha(sif) != c.SIF_SHA256:
        raise ValueError('Production SIF hash mismatch')
    tpl = c.load_json(old / 'slurm/rcpp_template_manifest.json')
    archive = old / 'slurm/rcpp_template.tar.gz'
    model = HERE.parents[2] / 'model'
    if tpl['sif_sha256'] != c.SIF_SHA256 or tpl['model_cpp_sha256'] != c.sha(model/'model_O2_supply_demand_MAP.cpp') or tpl['model_r_sha256'] != c.sha(model/'model_O2_supply_demand_MAP.R') or tpl['archive_sha256'] != c.sha(archive):
        raise ValueError('Existing private Rcpp template is incompatible')
    (root / 'slurm').mkdir(exist_ok=True)
    shutil.copyfile(archive, root/'slurm/rcpp_template.tar.gz')
    common = {r['parameter']: float(10 ** ((r['encoded_lower']+r['encoded_upper'])/2)
              if r['transform']=='log10' else (r['encoded_lower']+r['encoded_upper'])/2) for r in ranges}
    cases = [dict(scope='global',case_id='global',fit_seed='seed25')]
    cases += [dict(scope='neighborhood10pct',case_id=s,fit_seed=s) for s in c.PILOT_SEEDS]
    reuse = None
    if args.reuse_from:
        from refine import audit_source
        reuse = audit_source(Path(args.reuse_from).resolve(), file_hashes(), cases, c.OXYGEN,
                             c.RESOLUTIONS)
        c.save_json(inputs/'reuse_source_audit.json',reuse)
    tasks = []
    for case_index, case in enumerate(cases):
        case_root = c.case_directory(root,case); case_root.mkdir(parents=True,exist_ok=True)
        best = {r['parameter']:float(r['value']) for r in c.read(fit/case['fit_seed']/'best_params.tsv')}
        bounds = c.bounds_for(ranges, None if case['scope']=='global' else best)
        c.table(case_root/'parameter_bounds.tsv',bounds)
        c.table(case_root/'reference_best_params.tsv',
                [dict(parameter=p,value=best[p]) for p in c.ACTIVE])
        endpoint = manifest[case['fit_seed']]
        for repeat in range(1,c.REPEATS+1):
            scramble_seed = 2026101300 + repeat  # paired unit-cube designs across scopes; repeats independent
            x = c.samples(bounds,max(c.RESOLUTIONS),scramble_seed)
            for n in c.RESOLUTIONS[:-1]:
                if not np.array_equal(c.prefix(x,n), c.samples(bounds,n,scramble_seed)):
                    raise ValueError('Nested Sobol prefix differs from direct sampling')
            if reuse:
                from refine import check_design_prefix
                check_design_prefix(Path(reuse['root']), case, repeat, x, bounds, scramble_seed)
            design = case_root/'runs'/('R%02d'%repeat)
            design.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(design/'ordered_samples.npz',natural=x,parameters=np.array(c.ACTIVE))
            c.save_json(design/'design.json',dict(**case,repeat=repeat,scramble_seed=scramble_seed,
                N=max(c.RESOLUTIONS),resolutions=c.RESOLUTIONS,parameters=list(c.ACTIVE),
                n_samples=len(x),oxygen_pct=c.OXYGEN,paired_scramble_across_cases=True,
                ordered_samples_sha256=c.sha(design/'ordered_samples.npz'),
                bounds_sha256=c.sha(case_root/'parameter_bounds.tsv')))
            for b, block in enumerate(c.BLOCKS):
                part = design/'blocks'/block; part.mkdir(parents=True,exist_ok=True)
                c.table(part/'samples.tsv.gz',[dict(sample_id=i+1,**dict(zip(c.ACTIVE,map(float,row))))
                      for i,row in enumerate(x[b::c.STEP])])
                meta = dict(method='Sobol_second_order',**case,
                    repeat=repeat,scramble_seed=scramble_seed,block=block,block_index=b,
                    parent_stride=c.STEP,parent_sample_id_first=b+1,N=max(c.RESOLUTIONS),
                    n_samples=max(c.RESOLUTIONS),parameters=list(c.ACTIVE),oxygen_pct=c.OXYGEN,
                    samples_sha256=c.sha(part/'samples.tsv.gz'),
                    best_params_sha256=endpoint['best_params_sha256'],fit_config_sha256=endpoint['fit_config_sha256'],
                    validate_reference=b==0,common_reference_params=common)
                if reuse:
                    from refine import source_block
                    meta['reuse'] = source_block(Path(reuse['root']), part.relative_to(root),
                                                 reuse['source_N'])
                c.save_json(part/'metadata.json',meta)
                tasks.append(dict(task_id=len(tasks)+1,case_index=case_index,scope=case['scope'],
                    case_id=case['case_id'],fit_seed=case['fit_seed'],repeat=repeat,block=block,
                    block_index=b,part=str(part.relative_to(root)),
                    metadata_sha256=c.sha(part/'metadata.json'),samples_sha256=c.sha(part/'samples.tsv.gz')))
    c.table(root/'slurm/task_manifest.tsv',tasks)
    c.table(root/'slurm/case_manifest.tsv',[dict(case_index=i+1,**case) for i,case in enumerate(cases)])
    inputs_sha = {str(p.relative_to(root)):c.sha(p) for p in inputs.iterdir() if p.is_file()}
    inputs_sha.update({str(p.relative_to(root)):c.sha(p) for p in (root/'slurm/task_manifest.tsv',root/'slurm/case_manifest.tsv')})
    for case in cases:
        for filename in ('parameter_bounds.tsv','reference_best_params.tsv'):
            p=c.case_directory(root,case)/filename
            inputs_sha[str(p.relative_to(root))]=c.sha(p)
    c.table(root/'environment.tsv',[dict(component=name,version=version(name)) for name in
        ('SALib','numpy','scipy','pandas','matplotlib')]+[dict(component='python',version=sys.version.split()[0])])
    template = root/'slurm/rcpp_template.tar.gz'
    plan = dict(method='Sobol_second_order',stage='pilot',cases=cases,fit_root=str(fit),
        n_parameters=14,n_pairs=91,n_repeats=5,resolutions=c.RESOLUTIONS,oxygen_pct=c.OXYGEN,
        n_tasks=len(tasks),max_N=max(c.RESOLUTIONS),n_operator_evaluations=len(cases)*5*max(c.RESOLUTIONS)*c.STEP*len(c.OXYGEN),
        granularity='case x scramble repeat x one A/B/AB_i/BA_i matrix block; all pilot oxygen points',
        bootstrap_resamples=200,repeat_range_tolerance=.05,resolution_delta_tolerance=.02,
        convergence_policy='diagnostic_only; retain all representative endpoints and raw negative estimates',
        source_audit_sha256=c.sha(old/'input_manifest.json'),input_sha256=inputs_sha,code_sha256=file_hashes(),
        git_commit=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip(),
        prepared_utc=c.now(),qos='xxlarge',time='12:00:00',array_concurrency_limit=None,specified_compute_node=None,
        sif=dict(path=str(sif),sha256=c.SIF_SHA256,bytes=sif.stat().st_size,mtime_ns=sif.stat().st_mtime_ns),
        template=dict(path=str(template),sha256=c.sha(template),bytes=template.stat().st_size,mtime_ns=template.stat().st_mtime_ns))
    if reuse:
        plan.update(stage='resolution_pilot',reuse_source=reuse['root'],reuse_source_N=reuse['source_N'],
                    n_reused_operator_evaluations=len(cases)*5*reuse['source_N']*c.STEP*len(c.OXYGEN))
        plan['n_new_operator_evaluations']=plan['n_operator_evaluations']-plan['n_reused_operator_evaluations']
    c.save_json(root/'slurm/submission_plan.json',plan)
    from structure import write_structure
    write_structure(root/'model_structure')
    c.save_json(root/'status.json',dict(status='prepared',n_evaluation_tasks=len(tasks),n_cases=len(cases),
        n_pairs=91,repeats=5,global_and_neighborhood_separate=True,prepared_utc=c.now()))
    print(json.dumps({k:plan[k] for k in ('n_tasks','n_operator_evaluations','cases','max_N','resolutions')},indent=2))


def evaluate(root, part, plan):
    part.mkdir(parents=True,exist_ok=True)
    with open(part/'worker.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        output=part/'outputs.tsv.gz'
        if c.complete(output):
            receipt=c.load_json(str(output)+'.receipt.json')
            if receipt['metadata_sha256']!=c.sha(part/'metadata.json') or receipt['samples_sha256']!=c.sha(part/'samples.tsv.gz'):
                raise ValueError('Completion provenance changed')
            return
        meta=c.load_json(part/'metadata.json')
        if 'reuse' in meta:
            from refine import evaluate_increment
            evaluate_increment(root,part,plan,HERE/'evaluate.R')
        else:
            subprocess.run(['Rscript','--vanilla',str(HERE/'evaluate.R'),str(part),plan['fit_root'],
                             str(root/'inputs/fixed_o2_dominant_ploidy_201grid.tsv')],check=True)
        if not c.complete(output):
            raise ValueError('Model output completion unconfirmed')


def worker(args):
    root=Path(args.out_dir).resolve(); plan=context(root)
    tasks=c.read(root/'slurm/task_manifest.tsv')
    if not 1<=args.task_id<=len(tasks): raise ValueError('Invalid task id')
    task=tasks[args.task_id-1]; part=root/task['part']
    if c.sha(part/'metadata.json')!=task['metadata_sha256'] or c.sha(part/'samples.tsv.gz')!=task['samples_sha256']:
        raise ValueError('Submitted matrix block changed')
    evaluate(root,part,plan)
    print('Completed task',args.task_id,task['scope'],task['case_id'],task['repeat'],task['block'],flush=True)


def validate(args):
    root=Path(args.out_dir).resolve(); plan=context(root)
    target=root/'slurm/runtime_validation'; target.mkdir(parents=True,exist_ok=True)
    with open(target/'numerical_tests.log','w') as log:
        subprocess.run([sys.executable,str(HERE/'test_sobol.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
    checks=[]
    for case in plan['cases']:
        source=c.case_directory(root,case)/'runs/R01/blocks/A'
        part=target/(case['scope']+'_'+case['case_id']); part.mkdir(parents=True,exist_ok=True)
        meta=c.load_json(source/'metadata.json')
        meta.pop('reuse',None)
        original_rows=c.read(source/'samples.tsv.gz')
        c.table(part/'samples.tsv.gz',original_rows[:1])
        meta.update(n_samples=1,samples_sha256=c.sha(part/'samples.tsv.gz'))
        c.save_json(part/'metadata.json',meta)
        evaluate(root,part,plan)
        _,qc=c.output_array(part/'outputs.tsv.gz',1,c.OXYGEN)
        receipt=c.load_json(str(part/'outputs.tsv.gz')+'.receipt.json')
        incremental_check=None
        if plan.get('reuse_source'):
            # Exercise the real R increment and ID restoration before the array.
            probe=target/(case['scope']+'_'+case['case_id']+'_increment')
            probe.mkdir(exist_ok=True)
            c.table(probe/'samples.tsv.gz',[original_rows[0],dict(original_rows[plan['reuse_source_N']],sample_id=2)])
            probe_meta=dict(meta,n_samples=2,samples_sha256=c.sha(probe/'samples.tsv.gz'),reuse=dict(
                part=str(part),N=1,receipt_sha256=c.sha(str(part/'outputs.tsv.gz')+'.receipt.json'),
                outputs_sha256=c.sha(part/'outputs.tsv.gz'),metadata_sha256=c.sha(part/'metadata.json'),
                samples_sha256=c.sha(part/'samples.tsv.gz')))
            c.save_json(probe/'metadata.json',probe_meta)
            evaluate(root,probe,plan)
            restored,probe_qc=c.output_array(probe/'outputs.tsv.gz',2,c.OXYGEN)
            prefix,_=c.output_array(part/'outputs.tsv.gz',1,c.OXYGEN)
            np.testing.assert_array_equal(restored[:1],prefix)
            direct=target/(case['scope']+'_'+case['case_id']+'_direct_tail')
            direct.mkdir(exist_ok=True)
            c.table(direct/'samples.tsv.gz',[dict(original_rows[plan['reuse_source_N']],sample_id=1)])
            direct_meta=dict(meta,validate_reference=False,samples_sha256=c.sha(direct/'samples.tsv.gz'))
            c.save_json(direct/'metadata.json',direct_meta)
            evaluate(root,direct,plan)
            expected,_=c.output_array(direct/'outputs.tsv.gz',1,c.OXYGEN)
            np.testing.assert_array_equal(restored[1:],expected)
            incremental_check=dict(status='passed',n_points=probe_qc['n_points'],
                                   maximum_absolute_difference=0.,prefix_sample_id=1,increment_sample_id=2)
        checks.append(dict(**case,**qc,reference_checks=receipt['reference_checks'],
                           common_reference=receipt['common_reference'],incremental_check=incremental_check))
    reference={}
    for case in checks:
        for r in case['common_reference']:
            key=float(r['O2_pct']); value=np.array([r[o] for o in c.OUTPUTS])
            if key in reference and np.max(abs(value-reference[key]))>1e-8:
                raise ValueError('Different fixed-oxygen configurations in the common-vector audit')
            reference[key]=value
    import unittest
    from test_sobol import SobolChecks
    c.save_json(root/'runtime_validation.json',dict(status='passed',validated_utc=c.now(),
        numerical_tests=unittest.defaultTestLoader.loadTestsFromTestCase(SobolChecks).countTestCases(),
        case_checks=checks,git_commit=plan['git_commit'],sif_sha256=c.SIF_SHA256,
        common_configuration_status='passed'))
    print('SIF numerical tests and nine case reference/cache checks passed',flush=True)


def summarize(args):
    root=Path(args.out_dir).resolve(); plan=context(root)
    if not 1<=args.case_id<=len(plan['cases']): raise ValueError('Invalid case id')
    case=plan['cases'][args.case_id-1]; target=c.case_directory(root,case)
    bounds=c.read(target/'parameter_bounds.tsv'); rows=[]; qc=[]; receipts=[]
    expected={r['part']:r for r in c.read(root/'slurm/task_manifest.tsv')}
    for repeat in range(1,6):
        design=target/'runs'/('R%02d'%repeat); dm=c.load_json(design/'design.json')
        if dm['bounds_sha256']!=c.sha(target/'parameter_bounds.tsv') or dm['ordered_samples_sha256']!=c.sha(design/'ordered_samples.npz'):
            raise ValueError('Ordered design/bounds changed')
        arrays=[]
        for block in c.BLOCKS:
            part=design/'blocks'/block
            task=expected[str(part.relative_to(root))]
            if c.sha(part/'metadata.json')!=task['metadata_sha256'] or c.sha(part/'samples.tsv.gz')!=task['samples_sha256']:
                raise ValueError('Submitted block provenance changed')
            if not c.complete(part/'outputs.tsv.gz'):
                c.save_json(target/'status.json',dict(status='incomplete',missing_block=str(part.relative_to(root))))
                raise ValueError('Missing valid block: '+str(part))
            values,quality=c.output_array(part/'outputs.tsv.gz',plan['max_N'],c.OXYGEN)
            arrays.append(values); qc.append(dict(repeat=repeat,block=block,**quality))
            receipts.append(dict(path=str((part/'outputs.tsv.gz').relative_to(root)),sha256=c.sha(part/'outputs.tsv.gz')))
        values=c.assemble(arrays)
        np.savez_compressed(design/'ordered_outputs.npz',outputs=values,oxygen=np.array(c.OXYGEN),output_names=np.array(c.OUTPUTS))
        for n in c.RESOLUTIONS:
            for oi,o2 in enumerate(c.OXYGEN):
                for k,output in enumerate(c.OUTPUTS):
                    si,variance=c.analyze(bounds,values[:n*c.STEP,oi,k],seed=dm['scramble_seed']+10000+n+oi*10+k)
                    common=dict(scope=case['scope'],case_id=case['case_id'],fit_seed=case['fit_seed'],
                        repeat=repeat,N=n,O2_pct=o2,output=output,output_variance=variance,
                        index_status='defined' if np.isfinite(si['S1']).all() else 'near_constant_output')
                    for index in ('S1','ST'):
                        for p,name in enumerate(c.ACTIVE):
                            value=float(si[index][p]); half=float(si[index+'_conf'][p])
                            rows.append(dict(**common,index=index,parameter_i=name,parameter_j='',value=value,
                                bootstrap_half_width=half,bootstrap_lower=value-half,bootstrap_upper=value+half))
                    for i,j in c.PAIRS:
                        value=float(si['S2'][i,j]); half=float(si['S2_conf'][i,j])
                        rows.append(dict(**common,index='S2',parameter_i=c.ACTIVE[i],parameter_j=c.ACTIVE[j],value=value,
                            bootstrap_half_width=half,bootstrap_lower=value-half,bootstrap_upper=value+half))
    c.table(target/'indices_repeats.tsv.gz',rows)
    c.table(target/'numerical_qc.tsv',qc)
    from report import summarize_rows
    summaries=summarize_rows(rows,plan)
    c.table(target/'indices_summary.tsv',summaries)
    common_checks=[]
    for repeat in range(1,6):
        receipt=c.load_json(str(target/'runs'/('R%02d'%repeat)/'blocks/A/outputs.tsv.gz')+'.receipt.json')
        common_checks.extend(receipt['common_reference'])
    c.save_json(target/'completion.json',dict(status='complete',case=case,completed_utc=c.now(),
        n_repeats=5,n_pairs=91,indices_sha256=c.sha(target/'indices_repeats.tsv.gz'),
        summary_sha256=c.sha(target/'indices_summary.tsv'),output_receipts=receipts,
        common_reference=common_checks))
    c.save_json(target/'status.json',dict(status='complete',n_pairs=91,n_repeats=5,
        n_unstable_S2=sum(r['index']=='S2' and r['N']==max(c.RESOLUTIONS) and r['convergence_status']!='passed' for r in summaries)))
    print('Summarized',case,flush=True)


def finalize(args):
    root=Path(args.out_dir).resolve(); plan=context(root)
    from report import finalize_indices
    finalize_indices(root,plan)


def response_worker(args):
    root=Path(args.out_dir).resolve(); plan=context(root)
    rp=c.load_json(root/'slurm/response_plan.json')
    if c.sha(root/'slurm/response_tasks.tsv')!=rp['tasks_sha256']:
        raise ValueError('Response plan changed')
    tasks=c.read(root/'slurm/response_tasks.tsv')
    if not 1<=args.task_id<=len(tasks): raise ValueError('Invalid response task')
    task=tasks[args.task_id-1]; part=root/task['part']
    if c.sha(part/'metadata.json')!=task['metadata_sha256']:
        raise ValueError('Response metadata changed')
    evaluate(root,part,plan)


def finish(args):
    root=Path(args.out_dir).resolve(); plan=context(root)
    from report import finish_report
    finish_report(root,plan)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','validate','worker','summarize','finalize','response-worker','finish'])
    parser.add_argument('--out-dir',required=True)
    parser.add_argument('--task-id',type=int); parser.add_argument('--case-id',type=int)
    parser.add_argument('--fit-root'); parser.add_argument('--figure4-dir'); parser.add_argument('--efast-root'); parser.add_argument('--sif')
    parser.add_argument('--resolutions',nargs='+',type=int,default=[256,512,1024])
    parser.add_argument('--reuse-from',help='Completed pilot with an identical nested design and model implementation')
    args=parser.parse_args()
    try:
        globals()[args.command.replace('-','_')](args)
    except Exception as error:
        if args.command in ('validate','finalize','finish'):
            c.save_json(Path(args.out_dir)/'status.json',dict(status='incomplete',failed_stage=args.command,
                error=str(error),recorded_utc=c.now()))
        raise


if __name__=='__main__': main()
