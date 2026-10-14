#!/usr/bin/env python3
"""Independent full-bound Sobol post-processing; the model run is read-only."""
import argparse
import csv
import fcntl
import gzip
from importlib.metadata import version
import json
import os
from pathlib import Path
import shutil
import sys
import time

import numpy as np
import core as c
import full_core as f
import full_report as report


def inventory(snapshot):
    return {str(p.relative_to(snapshot)): c.sha(p)
            for p in sorted(Path(snapshot).rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts}


def source_manifest(root, check_outputs=False):
    root = Path(root).resolve()
    p = c.load_json(root/'slurm/submission_plan.json')
    if (p['N'], p['n_repeats'], p['n_parameters'], p['n_pairs']) != (4096, 5, 14, 91):
        raise ValueError('Require the submitted N4096/five-repeat/14-factor design')
    if p['resolutions'] != [2048, 4096] or len(p['oxygen_pct']) != 201:
        raise ValueError('Resolution or oxygen axes changed')
    np.testing.assert_allclose(p['oxygen_pct'], np.arange(201)*.025, atol=1e-12, rtol=0)
    c.configure(p)
    case = p['cases'][0]
    if case['scope'] != 'global' or case['case_id'] != 'global':
        raise ValueError('Require the full-bound global domain')
    source = f.case_dir(root, case)
    if c.load_json(source/'case.json') != case:
        raise ValueError('Global domain metadata changed')
    if c.sha(source/'parameter_bounds.tsv') != case['bounds_sha256']:
        raise ValueError('Global bounds changed')
    c.problem(c.read(source/'parameter_bounds.tsv'))
    original = Path(p['code_snapshot_root'])
    if c.sha(original/'oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/evaluate.R') != p['engine_sha256']:
        raise ValueError('Original model engine changed')
    for path, digest in p['code_sha256'].items():
        if c.sha(original/path) != digest:
            raise ValueError('Original frozen code changed: '+path)
    for repeat, digest in p['unit_sha256'].items():
        if c.sha(root/'inputs'/('unit_R%02d.npz'%int(repeat))) != digest:
            raise ValueError('Original unit design changed')
    if c.sha(root/'inputs/fixed_o2_dominant_ploidy_201grid.tsv') != p['reference_sha256']:
        raise ValueError('Original Figure 4 oxygen reference changed')
    for filename, field in [('best_params.tsv', 'best_params_sha256'),
                            ('fit_config.rds', 'fit_config_sha256'),
                            ('parameter_table.csv', 'parameter_table_sha256')]:
        if c.sha(Path(p['fit_root'])/case['fit_seed']/filename) != case[field]:
            raise ValueError('Original in-vivo fit configuration changed')
    sif = Path(p['sif']['path'])
    if sif.stat().st_size != p['sif']['bytes'] or sif.stat().st_mtime_ns != p['sif']['mtime_ns']:
        raise ValueError('Verified SIF changed')
    actual = {r['component']: (sys.version.split()[0] if r['component']=='python'
                              else version(r['component'])) for r in c.read(root/'environment.tsv')}
    expected = {r['component']: r['version'] for r in c.read(root/'environment.tsv')}
    if actual != expected or actual['SALib'] != '1.5.2':
        raise ValueError('Use the original verified SIF Python environment')
    blocks = []
    for repeat in range(1, 6):
        for block in c.BLOCKS:
            path = source/'runs'/('R%02d'%repeat)/'blocks'/block/'outputs.npz'
            receipt = Path(str(path)+'.receipt.json')
            r = c.load_json(receipt)
            if r['status'] != 'complete' or r['n_samples'] != 4096 or r['oxygen_points'] != 201:
                raise ValueError('Global model block incomplete: '+str(path))
            if not path.is_file() or (check_outputs and c.sha(path) != r['sha256']):
                raise ValueError('Global model output changed: '+str(path))
            blocks.append(dict(repeat=repeat, block=block, path=str(path.relative_to(root)),
                               bytes=path.stat().st_size, sha256=r['sha256'],
                               receipt_sha256=c.sha(receipt)))
    return p, dict(source_root=str(root), source_plan_sha256=c.sha(root/'slurm/submission_plan.json'),
                   source_code_git_commit=p['git_commit'], case=case, blocks=blocks,
                   engine_sha256=p['engine_sha256'], unit_sha256=p['unit_sha256'],
                   sif_sha256=p['sif']['sha256'], environment=actual)


def prepare(root, target, snapshot, revision):
    target = Path(target).resolve(); snapshot = Path(snapshot).resolve()
    if target != Path(root).resolve()/'global/standalone_N4096':
        raise ValueError('Use the isolated global/standalone_N4096 result directory')
    if f.REPO != snapshot:
        raise ValueError('Prepare from the independent code snapshot')
    p, sources = source_manifest(root)
    target.mkdir(parents=True, exist_ok=True)
    c.save_json(target/'inputs/source_manifest.json', sources)
    original = f.case_dir(root, p['cases'][0])
    for name in ['parameter_bounds.tsv', 'reference_best_params.tsv', 'case.json']:
        shutil.copyfile(original/name, target/'inputs'/name)
    shutil.copyfile(Path(root)/'environment.tsv', target/'environment.tsv')
    plan = dict(method='Sobol_second_order', scope='global', N=4096, repeats=5,
                oxygen_points=201, parameters=14, parameter_pairs=91,
                source_root=str(Path(root).resolve()), result_root=str(target),
                code_snapshot_root=str(snapshot), postprocess_git_commit=revision,
                code_sha256=inventory(snapshot),
                source_manifest_sha256=c.sha(target/'inputs/source_manifest.json'),
                source_code_git_commit=p['git_commit'], sif=p['sif'],
                source_environment=sources['environment'], bootstrap_resamples=200,
                confidence_level=.95, resolutions=[2048,4096],
                convergence_policy='diagnostic_only', qos='xxlarge', time='12:00:00',
                cpus=1, memory='16G', specified_node=None,
                dependency_on_neighborhood=None, prepared_utc=c.now())
    c.save_json(target/'slurm/postprocess_plan.json', plan)
    c.save_json(target/'chart_contract.json', dict(
        question='Global parameter main effects, total effects and named pair interactions across fixed oxygen',
        renderer='Matplotlib PDF/PNG; quantitative heatmaps',
        population='one global full-bound domain; five independently scrambled repeats',
        oxygen_axis='0 to 5 percent by 0.025 percent; 201 points',
        matrices={'S1':[14,201], 'ST':[14,201], 'S2':[91,201]},
        outputs=list(c.OUTPUTS), independent_output_color_scales=True,
        palette={'S1':'white to deep purple', 'ST':'white to orange',
                 'S2':'white to deep purple', 'negative':'blue', 'undefined':'gray'},
        source='N4096 five-repeat means; full signed estimates; all diagnostic states retained',
        export='three PDF/PNG pairs in figures/', qa='inspect PDFs and PNGs after collection'))
    c.save_json(target/'status.json', dict(status='prepared', scope='global', analysis_complete=False))
    print('Prepared independent global analysis; 150 completed model blocks', flush=True)


def repeat_rows(values, widths, variances, oxygen, resolutions):
    for ni, n in enumerate(resolutions):
        for ri in range(5):
            for oi, o2 in enumerate(oxygen):
                for k, output in enumerate(c.OUTPUTS):
                    for mi, (index, pi, pj) in enumerate(f.METRICS):
                        v = values[ni,ri,oi,k,mi]; w = widths[ni,ri,oi,k,mi]
                        yield dict(scope='global', N=int(n), repeat=ri+1, O2_pct=float(o2),
                                   output=output, index=index, parameter_i=pi, parameter_j=pj,
                                   group_i=c.efast.GROUP[pi], group_j=c.efast.GROUP[pj] if pj else '',
                                   value=v, bootstrap_half_width=w,
                                   bootstrap_confidence_low=v-w, bootstrap_confidence_high=v+w,
                                   confidence_level=.95, output_variance=variances[ni,ri,oi,k],
                                   included='TRUE')


def write_stream(path, rows):
    from itertools import chain
    rows = iter(rows); first = next(rows)
    c.efast.write_table(path, list(first), chain([first], rows))


def band_tables(target, oxygen, mean, status, spread, delta):
    rows = []
    for band, lo, hi in [('low_0_to_1pct',0,1), ('high_3_to_5pct',3,5)]:
        weights = c.efast.normalized_trapezoid_weights(np.asarray(oxygen),lo,hi)
        selected = weights > 0
        for k, output in enumerate(c.OUTPUTS):
            for mi, (index, pi, pj) in enumerate(f.METRICS):
                curve = mean[selected,k,mi]
                score = float(np.dot(curve, weights[selected])) if np.isfinite(curve).all() else float('nan')
                rows.append(dict(scope='global', output=output, band=band, index=index,
                                 parameter_i=pi, parameter_j=pj, group_i=c.efast.GROUP[pi],
                                 group_j=c.efast.GROUP[pj] if pj else '', value_mean=score,
                                 n_oxygen_points=int(selected.sum()),
                                 diagnostic_pass_fraction=float(np.mean(status[selected,k,mi]==1)),
                                 repeat_range_max=float(np.max(spread[selected,k,mi])),
                                 resolution_delta_max=float(np.max(delta[selected,k,mi])),
                                 included='TRUE'))
    c.table(target/'parameter_band_summary.tsv',[r for r in rows if r['index']!='S2'])
    pairs = [r for r in rows if r['index']=='S2']
    for output in c.OUTPUTS:
        for band in ['low_0_to_1pct','high_3_to_5pct']:
            selected = sorted([r for r in pairs if r['output']==output and r['band']==band],
                              key=lambda r:r['value_mean'] if np.isfinite(r['value_mean']) else -np.inf,
                              reverse=True)
            for rank, row in enumerate(selected,1): row['rank_within_output_band']=rank
    c.table(target/'pair_band_summary.tsv',pairs)
    groups = []
    for output in c.OUTPUTS:
        for band in ['low_0_to_1pct','high_3_to_5pct']:
            for index in ['S1','ST']:
                for group in sorted(set(c.efast.GROUP.values())):
                    rr = [r for r in rows if (r['output'],r['band'],r['index'],r['group_i'])==(output,band,index,group)]
                    groups.append(dict(output=output,band=band,index=index,process_group=group,
                                       n_parameters=len(rr),
                                       mean_parameter_index=float(np.mean([r['value_mean'] for r in rr])),
                                       diagnostic_pass_fraction=float(np.mean([r['diagnostic_pass_fraction'] for r in rr])),
                                       aggregation='mean of individual parameter indices; not group variance'))
    c.table(target/'process_band_summary.tsv',groups)
    return rows, pairs, groups


def finalize(target, plan):
    with np.load(target/'indices_repeats.npz') as z:
        values=z['values']; widths=z['bootstrap_half_width']; variances=z['variance']
        oxygen=z['oxygen']; resolutions=z['resolutions']
    write_stream(target/'indices_repeats.tsv.gz',repeat_rows(values,widths,variances,oxygen,resolutions))
    with np.load(target/'indices_summary.npz') as z:
        mean=z['mean']; status=z['status']; spread=z['repeat_range']; delta=z['resolution_delta']
    diagnostic=[]
    for k, output in enumerate(c.OUTPUTS):
        for index in ['S1','ST','S2']:
            ii=[mi for mi,m in enumerate(f.METRICS) if m[0]==index]
            states=status[:,k,ii]
            diagnostic.append(dict(output=output,index=index,n_estimates=int(states.size),
                                   passed=int(np.sum(states==1)),not_passed=int(np.sum(states==2)),
                                   undefined=int(np.sum(states==0)),
                                   negative_estimates=int(np.sum(mean[:,k,ii]<0)),
                                   convergence_policy='diagnostic_only'))
    c.table(target/'diagnostic_summary.tsv',diagnostic)
    rows,pairs,groups=band_tables(target,oxygen,mean,status,spread,delta)
    paragraphs=[]
    for output in c.OUTPUTS:
        paragraphs.append('## '+output)
        for index in ['S1','ST']:
            for band in ['low_0_to_1pct','high_3_to_5pct']:
                ranked=sorted([r for r in groups if (r['output'],r['index'],r['band'])==(output,index,band)],
                              key=lambda r:r['mean_parameter_index'] if np.isfinite(r['mean_parameter_index']) else -np.inf,reverse=True)
                if ranked and np.isfinite(ranked[0]['mean_parameter_index']):
                    r=ranked[0]
                    paragraphs.append('- '+index+' / '+band+': largest process mean is '+r['process_group']+
                                      ' ('+format(r['mean_parameter_index'],'.4g')+
                                      '; diagnostic pass fraction '+format(r['diagnostic_pass_fraction'],'.1%')+').')
        for band in ['low_0_to_1pct','high_3_to_5pct']:
            ranked=[r for r in pairs if r['output']==output and r['band']==band and r['rank_within_output_band']==1]
            if ranked and np.isfinite(ranked[0]['value_mean']):
                r=ranked[0]
                paragraphs.append('- Largest estimated S2 / '+band+': '+r['parameter_i']+' × '+r['parameter_j']+
                                  ' ('+format(r['value_mean'],'.4g')+
                                  '; diagnostic pass fraction '+format(r['diagnostic_pass_fraction'],'.1%')+').')
    (target/'interpretation.md').write_text(
        '# Global full-bound Sobol sensitivity\n\n'
        'N4096, five independent scrambles, 201 fixed oxygen conditions, 14 parameters and all 91 pairs. '
        'Both model outputs use the completed canonical in-vivo operator evaluations. '
        'S1 is the independent main effect; ST includes all interactions involving that parameter; '
        'S2 identifies a two-parameter variance contribution and does not include all higher-order interactions.\n\n'+
        '\n\n'.join(paragraphs)+'\n\n'
        'These are global full-bound estimates. Neighborhood robustness across 500 fitted endpoints is pending. '
        'Process means average parameter indices and are not formal group variance fractions; ST must not be summed as a variance fraction. '
        'Ranks retain all estimates regardless of diagnostics. A high estimate with poor repeat/resolution agreement is an exploratory signal, '
        'not a stable mechanistic ranking. Negative estimates reflect sampling error and are retained. '
        'Gray cells are undefined. Neither S1, ST nor S2 gives an effect direction. '
        'Per-repeat bootstrap intervals describe numerical estimation uncertainty, not biological or fitted-endpoint uncertainty; '
        'the median half-width in the summary is not a confidence interval for the five-repeat mean.\n')
    summary=c.read(target/'indices_summary.tsv.gz')
    expected=len(oxygen)*2*len(f.METRICS)
    if len(summary)!=expected or any(r['included_in_full_analysis']!='TRUE' for r in summary):
        raise ValueError('Summary export incomplete or filtered')
    with gzip.open(target/'indices_repeats.tsv.gz','rt',newline='') as handle:
        n_repeated=sum(1 for _ in csv.DictReader(handle,delimiter='\t'))
    if n_repeated!=expected*10:
        raise ValueError('Repeated export incomplete')
    for name in ['s1_parameter_by_oxygen','st_parameter_by_oxygen','s2_pair_by_oxygen']:
        for ext in ['pdf','png']:
            if not (target/'figures'/(name+'.'+ext)).is_file():
                raise ValueError('Figure missing')
    qc=dict(status='passed',summary_rows=expected,repeat_rows=expected*10,
            raw_model_blocks=150,oxygen_points=len(oxygen),N=4096,repeats=5,n_parameters=14,n_pairs=91,
            all_estimates_retained=True,source_outputs_verified_by_sha256=True,
            mathematical_validation='additive/product model; interrupted resume; idempotence tests',
            figures_generated=6,figure_visual_review='pending collection',validated_utc=c.now())
    c.save_json(target/'validation.json',qc)
    (target/'README.md').write_text(
        '# Independent global N4096 Sobol results\n\n'
        'This directory contains the global full-bound analysis only. The existing 500-neighborhood pipeline continues separately.\n\n'
        '- figures/: S1 and ST parameter-by-oxygen heatmaps; S2 pair-by-oxygen heatmaps. Each has separate ploidy/growth color scales.\n'
        '- indices_summary.tsv.gz: N4096 five-repeat means, repeat ranges, N2048→4096 changes and diagnostics.\n'
        '- indices_repeats.tsv.gz: every estimate and per-repeat 95% bootstrap interval at N2048 and N4096.\n'
        '- diagnostic_summary.tsv and numerical_qc.tsv: sensitivity and original solver diagnostics.\n'
        '- parameter_band_summary.tsv, pair_band_summary.tsv, process_band_summary.tsv: low 0–1% and high 3–5% O2 summaries.\n'
        '- inputs/source_manifest.json: original code, input designs and 150 raw model-block hashes. Raw blocks remain at the source HPC root.\n'
        '- interpretation.md, validation.json and collection_manifest.tsv: interpretation, export checks and file inventory.\n\n'
        'All convergence statuses and negative estimates are retained. Five scrambled sampling repeats are distinct from 200 bootstrap resamples. '
        'N2048 is an existing design prefix; no model samples were rerun.\n')
    c.save_json(target/'status.json',dict(status='complete',scope='global',analysis_complete=True,N=4096,
               repeats=5,oxygen_points=201,n_parameters=14,n_pairs=91,convergence_policy='diagnostic_only',
               source_code_git_commit=plan['source_code_git_commit'],
               postprocess_git_commit=plan['postprocess_git_commit'],
               job_id=os.environ.get('SLURM_JOB_ID',''),completed_utc=c.now()))
    excluded={'collection_manifest.tsv','postprocess.lock','indices_checkpoint.npz'}
    files=[]
    for path in sorted(target.rglob('*')):
        rel=path.relative_to(target)
        if (not path.is_file() or path.name in excluded or
            'code_snapshot' in rel.parts or 'logs' in rel.parts or
            path.suffix not in ['.json','.tsv','.gz','.md','.pdf','.png']):
            continue
        files.append(dict(path=str(rel),bytes=path.stat().st_size,sha256=c.sha(path)))
    c.table(target/'collection_manifest.tsv',files)
    print('Global figures and source tables complete; files:',len(files),flush=True)


def run(target):
    target=Path(target).resolve(); plan=c.load_json(target/'slurm/postprocess_plan.json')
    if f.REPO!=Path(plan['code_snapshot_root']) or inventory(f.REPO)!=plan['code_sha256']:
        raise ValueError('Independent frozen post-processing code changed')
    with open(target/'postprocess.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        start=time.monotonic()
        p,sources=source_manifest(plan['source_root'],check_outputs=True)
        saved=c.load_json(target/'inputs/source_manifest.json')
        if sources!=saved or c.sha(target/'inputs/source_manifest.json')!=plan['source_manifest_sha256']:
            raise ValueError('Read-only source inventory changed')
        metadata=dict(postprocess_git_commit=plan['postprocess_git_commit'],
                      source_manifest_sha256=plan['source_manifest_sha256'],code_sha256=plan['code_sha256'])
        already_complete=c.load_json(target/'status.json').get('status')=='complete'
        if already_complete:
            report.summarize_case(Path(plan['source_root']),p,p['cases'][0],
                                  target_dir=target,analysis_metadata=metadata)
            for row in c.read(target/'collection_manifest.tsv'):
                artifact=target/row['path']
                if artifact.stat().st_size!=int(row['bytes']) or c.sha(artifact)!=row['sha256']:
                    raise ValueError('Completed artifact changed: '+row['path'])
            print('Independent global analysis already complete and verified',flush=True)
            return
        c.save_json(target/'status.json',dict(status='computing_global_indices',scope='global',
                   analysis_complete=False,started_utc=c.now(),job_id=os.environ.get('SLURM_JOB_ID','')))
        report.summarize_case(Path(plan['source_root']),p,p['cases'][0],
                              target_dir=target,analysis_metadata=metadata)
        finalize(target,plan)
        print('Elapsed global post-processing seconds:',time.monotonic()-start,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','run'])
    parser.add_argument('--out-dir',required=True)
    parser.add_argument('--source-root')
    parser.add_argument('--snapshot-root')
    parser.add_argument('--git-commit')
    args=parser.parse_args()
    if args.command=='prepare':
        prepare(args.source_root,args.out_dir,args.snapshot_root,args.git_commit)
    else:
        run(args.out_dir)


if __name__=='__main__':
    main()
