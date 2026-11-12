#!/usr/bin/env python3
"""Immutable global eFAST inputs, resumable model blocks and SALib estimates."""
import argparse
import importlib.metadata
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from common import (ACTIVE, OUTPUTS, METRICS, INITIAL_NS, M, CHUNK, REPEATS,
                    now, sha, load, save, read, table, locked, decode,
                    decode_analysis, model_path, model_count)


def runtime(root):
    global np, c, f
    import numpy as np
    root = Path(root)
    p = load(root / 'plan.json')
    if p['scope'] != 'global' or p['n_repeats'] != 5 or p['n_oxygen'] != 201 or p['N_hard_limit'] is not None:
        raise ValueError('Global-only five-repeat / 201-oxygen contract changed')
    if importlib.metadata.version('SALib') != '1.5.2':
        raise ValueError('Require the original SALib 1.5.2 SIF')
    for name, digest in p['controller_sha256'].items():
        if sha(root / 'slurm/code_snapshot' / name) != digest:
            raise ValueError('Frozen controller changed: ' + name)
    if sha(p['original_plan']) != p['original_plan_sha256']:
        raise ValueError('Original Sobol model provenance changed')
    original = load(p['original_plan'])
    helpers = Path(original['code_snapshot_root']) / 'oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol'
    sys.path.insert(0, str(helpers))
    import core as c
    import full_core as f
    verified = f.plan(p['source_initial_root'])
    if list(c.ACTIVE) != list(ACTIVE) or tuple(c.OUTPUTS) != OUTPUTS:
        raise ValueError('Parameter/output order differs from the certified Sobol model')
    for rel, digest in p['source_files'].items():
        if sha(root / rel) != digest:
            raise ValueError('Copied source changed: ' + rel)
    return p, verified


def design_file(root, n, repeat):
    return Path(root) / 'inputs/designs' / ('N%d_R%02d.npz' % (n, repeat))


def evaluate(root, p, original, x, oxygen, target, validate=False):
    """Use the exact canonical Sobol operator; only the accepted method label differs."""
    target.mkdir(parents=True, exist_ok=True)
    sample = target / 'samples.tsv.gz'
    c.table(sample, [dict(sample_id=i + 1, **dict(zip(ACTIVE, map(float, row)))) for i, row in enumerate(x)])
    case = original['cases'][0]
    if case['scope'] != 'global':
        raise ValueError('Only the original global reference case is permitted')
    meta = dict(method='eFAST', scope='global', case_id='global', fit_seed=case['fit_seed'],
                n_samples=len(x), parameters=list(ACTIVE), oxygen_pct=list(map(float, oxygen)),
                samples_sha256=sha(sample), best_params_sha256=case['best_params_sha256'],
                fit_config_sha256=case['fit_config_sha256'], validate_reference=validate)
    save(target / 'metadata.json', meta)
    subprocess.run(['Rscript', '--vanilla', str(Path(root) / 'slurm/code_snapshot/evaluate_efast.R'),
                    str(target), original['fit_root'], str(Path(p['source_initial_root']) / 'inputs/fixed_o2_dominant_ploidy_201grid.tsv')], check=True)
    path = target / 'outputs.tsv.gz'
    receipt = load(str(path) + '.receipt.json')
    if receipt['samples_sha256'] != meta['samples_sha256'] or receipt['metadata_sha256'] != sha(target / 'metadata.json') or receipt['outputs_sha256'] != sha(path):
        raise ValueError('Model rectangle provenance mismatch')
    return f.read_model(path, len(x), oxygen), receipt


def estimator_tests(bounds):
    from SALib.sample import fast_sampler
    from SALib.analyze import fast
    problem = dict(num_vars=3, names=['a', 'b', 'c'], bounds=[[0., 1.]] * 3)
    n = 1025
    samples = fast_sampler.sample(problem, n, M=4, seed=5826)
    y = samples @ np.array([1., 2., 3.])
    public = fast.analyze(problem, y, M=4, num_resamples=2, seed=5826)
    direct = np.array([fast.compute_orders(y[i*n:(i+1)*n], n, 4, (n-1)//8) for i in range(3)])
    error = float(np.max(np.abs(direct - np.column_stack([public['S1'], public['ST']]))))
    expected = np.array([1., 4., 9.]) / 14
    analytic = float(np.max(np.abs(direct - expected[:, None])))
    if error > 1e-14 or analytic > .02:
        raise ValueError('FAST public-estimator parity or additive reference failed')
    # Verify unit transport against SALib's encoded-bound sampling before exponentiation.
    unit_problem = dict(num_vars=14, names=list(ACTIVE), bounds=[[0., 1.]] * 14)
    enc_problem = dict(num_vars=14, names=list(ACTIVE), bounds=[[float(r['encoded_lower']), float(r['encoded_upper'])] for r in bounds])
    unit = fast_sampler.sample(unit_problem, 2049, M=4, seed=20260923)
    enc = fast_sampler.sample(enc_problem, 2049, M=4, seed=20260923)
    lo = np.array([float(r['encoded_lower']) for r in bounds])
    hi = np.array([float(r['encoded_upper']) for r in bounds])
    if not np.array_equal(unit * (hi-lo) + lo, enc):
        raise ValueError('Transport does not match SALib transformed bounds')
    return dict(public_compute_orders_max_error=error, additive_reference_max_error=analytic,
                transformed_bound_transport_bitwise_equal=True, bootstrap_confidence_not_used=True)


def prepare(root):
    from SALib.sample import fast_sampler
    p, original = runtime(root)
    root = Path(root)
    if (root / 'preparation.json').exists():
        if load(root / 'preparation.json')['plan_sha256'] != sha(root / 'plan.json'):
            raise ValueError('Preparation identity changed')
        return
    subprocess.run([sys.executable,str(root/'slurm/code_snapshot/test_workflow.py')],check=True)
    if sha(original['sif']['path']) != original['sif']['sha256']:
        raise ValueError('SIF digest differs')
    for name, expected in original['code_sha256'].items():
        if sha(Path(original['code_snapshot_root']) / name) != expected:
            raise ValueError('Frozen canonical source differs: ' + name)
    expected_engine=(f.HERE/'evaluate.R').read_text().replace(
        '# Canonical fixed-O2 matrix, checked dominant eigenpair, unchanged-matrix recovery.',
        '# eFAST inputs; identical frozen Sobol fixed-O2 matrix and eigenpair recovery.')
    expected_engine=expected_engine.replace('here <- dirname(normalizePath(script))',
        "here <- normalizePath(Sys.getenv('EFAST_CANONICAL_HELPERS', unset='__missing_canonical_helpers__'))")
    expected_engine=expected_engine.replace("stopifnot(meta$method %in% c('Sobol_second_order', 'Sobol_response_surface'))",
        "stopifnot(identical(meta$method, 'eFAST'))")
    if (root/'slurm/code_snapshot/evaluate_efast.R').read_text()!=expected_engine:
        raise ValueError('The eFAST model body differs from the frozen Sobol evaluator')
    case = original['cases'][0]
    for name, field in [('best_params.tsv', 'best_params_sha256'), ('fit_config.rds', 'fit_config_sha256'), ('parameter_table.csv', 'parameter_table_sha256')]:
        if sha(Path(original['fit_root']) / case['fit_seed'] / name) != case[field]:
            raise ValueError('Original fit input differs')
    bounds = read(root / 'inputs/parameter_bounds.tsv')
    legacy_bounds = read(root / 'inputs/legacy_parameter_ranges.tsv')
    if bounds != legacy_bounds or [r['parameter'] for r in bounds] != list(ACTIVE):
        raise ValueError('Global eFAST and Sobol bounds differ')
    results = estimator_tests(bounds)
    oxygen = np.asarray(p['oxygen_pct'])
    if not np.allclose(oxygen, np.linspace(0, 5, 201), rtol=0, atol=1e-12):
        raise ValueError('Oxygen design changed')
    problem = dict(num_vars=14, names=list(ACTIVE), bounds=[[0., 1.]] * 14)
    designs = {}
    for n in p['resolutions']:
        for repeat, seed in enumerate(p['phase_seeds'], 1):
            target = design_file(root, n, repeat)
            unit = fast_sampler.sample(problem, n, M=p['M'], seed=seed)
            x = f.scale_unit(unit, bounds)
            for j, row in enumerate(bounds):
                if np.min(x[:, j]) < float(row['natural_lower']) - 1e-12 or np.max(x[:, j]) > float(row['natural_upper']) + 1e-12:
                    raise ValueError('Natural sample out of bounds')
            if target.exists():
                with np.load(target) as z:
                    if not np.array_equal(z['natural'], x):
                        raise ValueError('Existing eFAST design differs')
            else:
                f.atomic_npz(target, natural=x, N=np.array(n), M=np.array(p['M']), phase_seed=np.array(seed))
            designs[str(target.relative_to(root))] = sha(target)
    # Method-alignment checks reuse cached Sobol outputs, not another Sobol run.
    directory = f.case_dir(p['source_initial_root'], case)
    x = f.samples(p['source_initial_root'], case, 1, 0)
    with np.load(directory / 'runs/R01/blocks/A/outputs.npz') as z:
        old_values = z['outputs']; old_flags = z['solver_flags']
        selected = sorted(set([0, 1, 2] + list(np.flatnonzero((old_flags == 2).any(axis=1))[:3])))
        expected = old_values[selected]
    (actual, flags, gaps, residuals), receipt = evaluate(root, p, original, x[selected], oxygen,
        Path(os.environ['SOBOL_SCRATCH']) / 'alignment', validate=True)
    error = np.max(np.abs(actual - expected), axis=(0, 1))
    if np.any(error > 1e-8):
        raise ValueError('New eFAST evaluator differs from cached Sobol outputs')
    results.update(cached_sobol_samples=selected, checked_model_points=len(selected)*201,
                   output_max_absolute_error=dict(zip(OUTPUTS, map(float, error))),
                   maximum_eigen_residual=float(np.max(residuals)), minimum_spectral_gap=float(np.min(gaps)))
    save(root / 'validation/method_alignment.json', dict(status='passed', completed_utc=now(), **results))
    save(root / 'preparation.json', dict(status='complete', plan_sha256=sha(root / 'plan.json'),
        design_sha256=designs, model_alignment_validation_points=len(selected)*201,
        model_alignment_elapsed_seconds=receipt['elapsed_seconds'], neighborhood_files_read=0,
        completed_utc=now()))


def model(root, n, task, smoke=False):
    p, original = runtime(root)
    if n not in p['resolutions']:
        raise ValueError('N is not in this submitted batch; append another batch explicitly')
    root = Path(root)
    prep = load(root / 'preparation.json')
    repeat, pi, begin, end = decode(n, task)
    design = design_file(root, n, repeat)
    if sha(design) != prep['design_sha256'][str(design.relative_to(root))]:
        raise ValueError('eFAST design changed')
    with np.load(design) as z:
        x = z['natural'][pi*n+begin:pi*n+end]
    key = __import__('hashlib').sha256(('%s/%d/%d/%d/%d/%d' % (sha(root/'plan.json'), n, repeat, pi, begin, end)).encode()).hexdigest()
    path = model_path(root, n, repeat, pi, begin)
    with locked(path.parent / 'worker.lock'):
        if f.complete(path, key):
            return
        oxygen = np.asarray(p['oxygen_pct'])
        values, flags, done, gaps, residuals, elapsed = f.checkpoint(path.parent/'checkpoint.npz', key, len(x), oxygen)
        scratch = Path(os.environ['SOBOL_SCRATCH']) / 'evaluation'
        while not done.all():
            positions = np.flatnonzero(~done)[:8]
            if shutil.disk_usage(root).free < 32 * 1024**3:
                raise ValueError('Storage reserve below 32 GiB; retain checkpoint')
            (v, s, g, residual), receipt = evaluate(root, p, original, x, oxygen[positions], scratch)
            values[:, positions], flags[:, positions] = v, s
            gaps[positions], residuals[positions], done[positions] = g, residual, True
            elapsed += receipt['elapsed_seconds']
            if (scratch/'recovery').exists():
                shutil.copytree(scratch/'recovery', path.parent/'recovery', dirs_exist_ok=True)
                shutil.rmtree(scratch/'recovery')
            f.atomic_npz(path.parent/'checkpoint.npz', identity=np.array(key), oxygen=oxygen, outputs=values,
                         solver_flags=flags, done=done, minimum_gap=gaps, maximum_residual=residuals,
                         elapsed_seconds=np.array(elapsed))
            print('N', n, 'R', repeat, ACTIVE[pi], begin, 'oxygen', int(done.sum()), '/201', flush=True)
        f.atomic_npz(path, identity=np.array(key), oxygen=oxygen, outputs=values, solver_flags=flags,
                     minimum_gap=gaps, maximum_residual=residuals)
        save(str(path)+'.receipt.json', dict(status='complete', identity=key, sha256=sha(path), N=n,
             repeat=repeat, parameter=ACTIVE[pi], sample_begin=begin, sample_end=end,
             n_points=len(x)*201, n_recovered=int((flags==2).sum()), elapsed_seconds=elapsed,
             completed_utc=now(), job_id=os.environ.get('SLURM_JOB_ID')))
        (path.parent/'checkpoint.npz').unlink()
        if smoke:
            save(root/'smoke.json', dict(status='complete', path=str(path), sha256=sha(path), n_points=len(x)*201,
                 full_oxygen_grid=True, production_reuses_completed_block=True, elapsed_seconds=elapsed, completed_utc=now()))


def analyze(root, n, task):
    from SALib.analyze.fast import compute_orders
    p, original = runtime(root)
    root = Path(root)
    repeat, pi = decode_analysis(task)
    target = root/'efast'/('N%d'%n)/'analysis'/('R%02d_P%02d.npz'%(repeat,pi))
    key = sha(root/'plan.json')+'/%d/%d/%d'%(n,repeat,pi)
    with locked(target.parent/(target.stem+'.lock')):
        if f.complete(target,key):
            return
        arrays, receipts = [], []
        for begin in range(0,n,CHUNK):
            path=model_path(root,n,repeat,pi,begin)
            if not f.complete(path):
                raise ValueError('Incomplete eFAST trajectory; never discard rows')
            rec=load(str(path)+'.receipt.json')
            if (rec['N'],rec['repeat'],rec['parameter'],rec['sample_begin'],rec['sample_end']) != (n,repeat,ACTIVE[pi],begin,min(begin+CHUNK,n)):
                raise ValueError('Trajectory block position differs')
            with np.load(path) as z:
                if not np.array_equal(z['oxygen'],p['oxygen_pct']):
                    raise ValueError('Trajectory oxygen differs')
                arrays.append(z['outputs'])
            receipts.append(dict(path=str(path),sha256=rec['sha256']))
        y=np.concatenate(arrays,axis=0)
        if y.shape!=(n,201,2) or not np.isfinite(y).all():
            raise ValueError('Incomplete/nonfinite ordered eFAST trajectory')
        result=np.full((201,2,2),np.nan)
        variance=np.var(y,axis=0)
        for oi in range(201):
            for k in range(2):
                if variance[oi,k] > (1e-12*max(1.,float(np.max(np.abs(y[:,oi,k])))))**2:
                    result[oi,k]=compute_orders(y[:,oi,k],n,p['M'],(n-1)//(2*p['M']))
        f.atomic_npz(target,identity=np.array(key),values=result,variance=variance,output_mean=np.mean(y,axis=0),oxygen=np.asarray(p['oxygen_pct']))
        save(str(target)+'.receipt.json',dict(status='complete',identity=key,sha256=sha(target),
             source_receipts=receipts,estimator='SALib.analyze.fast.compute_orders',
             public_analyze_parity_verified=True,completed_utc=now()))


def summarize(root,n):
    p, original=runtime(root)
    root=Path(root); dest=root/'efast'/('N%d'%n)
    values=np.full((5,201,2,28),np.nan); variances=np.full((5,201,2,14),np.nan)
    output_means=np.full_like(variances,np.nan)
    for task in range(1,71):
        repeat,pi=decode_analysis(task)
        path=dest/'analysis'/('R%02d_P%02d.npz'%(repeat,pi))
        if not f.complete(path,sha(root/'plan.json')+'/%d/%d/%d'%(n,repeat,pi)):
            raise ValueError('Missing phase/trajectory estimate')
        with np.load(path) as z:
            values[repeat-1,:,:,pi]=z['values'][:,:,0]
            values[repeat-1,:,:,pi+14]=z['values'][:,:,1]
            variances[repeat-1,:,:,pi]=z['variance']; output_means[repeat-1,:,:,pi]=z['output_mean']
    mean=values.mean(axis=0); spread=np.ptp(values,axis=0)
    low=p['resolutions'][p['resolutions'].index(n)-1] if n!=p['resolutions'][0] else None
    delta=np.full_like(mean,np.nan)
    if low is not None:
        with np.load(root/'efast'/('N%d'%low)/'indices_summary.npz') as z:
            delta=np.abs(mean-z['mean'])
    rows=[]
    for oi,o2 in enumerate(p['oxygen_pct']):
        for k,out in enumerate(OUTPUTS):
            for j,(idx,pi,pj) in enumerate(METRICS):
                valid=bool(np.isfinite(values[:,oi,k,j]).all())
                status='undefined' if not valid else ('no_adjacent_new_design' if low is None else ('passed' if spread[oi,k,j]<=.05 and delta[oi,k,j]<=.02 else 'not_passed'))
                rows.append(dict(method='eFAST',scope='global',N=n,M=p['M'],O2_pct=o2,output=out,index=idx,
                    parameter_i=pi,parameter_j='',value_mean=mean[oi,k,j],repeat_min=np.min(values[:,oi,k,j]),
                    repeat_max=np.max(values[:,oi,k,j]),repeat_range=spread[oi,k,j],resolution_low=low or '',
                    resolution_delta=delta[oi,k,j],convergence_status=status,included=True))
    table(dest/'indices_summary.tsv.gz',rows)
    f.atomic_npz(dest/'indices_repeats.npz',values=values,variance=variances,output_mean=output_means,oxygen=np.asarray(p['oxygen_pct']))
    f.atomic_npz(dest/'indices_summary.npz',mean=mean,repeat_range=spread,resolution_delta=delta,oxygen=np.asarray(p['oxygen_pct']))
    table(dest/'not_passed.tsv.gz',[r for r in rows if r['convergence_status']!='passed'],list(rows[0]))
    save(dest/'completion.json',dict(status='complete',method='eFAST',N=n,M=p['M'],n_repeats=5,n_oxygen=201,n_indices=len(rows),
        summary_sha256=sha(dest/'indices_summary.tsv.gz'),convergence_policy='diagnostic_only',completed_utc=now()))
    from report import efast_report
    efast_report(root,n)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['prepare','smoke','model','analyze','report','sobol_report','final_report'])
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--n',type=int,default=2049)
    parser.add_argument('--task',type=int,default=1)
    args=parser.parse_args()
    try:
        if args.mode=='prepare':prepare(args.root)
        elif args.mode in ['smoke','model']:model(args.root,args.n,args.task,smoke=args.mode=='smoke')
        elif args.mode=='analyze':analyze(args.root,args.n,args.task)
        elif args.mode=='report':summarize(args.root,args.n)
        else:
            runtime(args.root)
            from report import sobol_report,final_report
            (sobol_report if args.mode=='sobol_report' else final_report)(args.root)
    except Exception as exc:
        save(args.root/'failures'/('%s_N%d_T%d.json'%(args.mode,args.n,args.task)),
             dict(error=str(exc),utc=now(),job_id=os.environ.get('SLURM_JOB_ID')))
        raise


if __name__=='__main__':main()
