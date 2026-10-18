#!/usr/bin/env python3
"""Independently reconcile and interpret six-seed N4097 results.

No model execution, FAST resampling or Slurm submission. Statistics are
recomputed from phase arrays, independently of the frozen calculation core.
"""
import argparse
import csv
import gzip
import json
import platform
from pathlib import Path
import numpy as np
import efast
from report_resolution_pilot import (KEYS, WINDOWS, INDICES, Verifier, sha,
                                    read, write, stats, band, diagnostic, mechanism)
from collect_resolution_4097_results import SEEDS, PARAMETERS, validate_tasks
NS = (513,1025,2049,4097)
PAIRS = tuple(zip(NS[:-1],NS[1:]))


def load(root):
    selection = read(root/'seed_selection.tsv')
    seeds = [r['fit_seed'] for r in selection]
    assert seeds == list(SEEDS) and tuple(efast.ACTIVE) == PARAMETERS
    plan = json.loads((root/'slurm/submission_plan.json').read_text())
    assert plan['fit_seeds'] == seeds and plan['primary_comparison'] == [2049,4097]
    assert sha(root/'baseline_indices.npz') == plan['comparison_indices_sha256']
    for name,digest in plan['code_sha256'].items(): assert sha(Path(__file__).parent/name) == digest
    previous = root.parent/'N1025_N2049'
    prefix = plan['previous_pilot_dir']+'/'
    for name,digest in plan['comparison_source_sha256'].items():
        assert name.startswith(prefix) and sha(previous/name[len(prefix):]) == digest, name
    tasks = read(root/'slurm/task_manifest.tsv');validate_tasks(tasks)
    assert len(read(root/'slurm/summary_manifest.tsv')) == 6
    values,variance = {},{}
    with np.load(root/'baseline_indices.npz',allow_pickle=False) as data:
        assert list(data['fit_seeds']) == seeds and list(data['resolutions']) == list(NS[:-1])
        assert list(data['parameters']) == list(efast.ACTIVE) and list(data['outputs']) == list(efast.OUTPUTS)
        assert list(data['index_names']) == list(INDICES) and data['indices'].shape == (6,3,5,2,14,201,2)
        oxygen = data['oxygen'].copy()
        for i,seed in enumerate(seeds):
            for j,n in enumerate(NS[:-1]): values[(seed,n)] = data['indices'][i,j].copy()
    np.testing.assert_allclose(oxygen,np.arange(201)*.025,rtol=0,atol=1e-12)
    with np.load(previous/'baseline_N513.npz',allow_pickle=False) as data:
        for seed in seeds:
            position = list(data['fit_seeds']).index(seed)
            np.testing.assert_array_equal(values[(seed,513)],data['indices'][position])
    for seed in seeds:
        for n in NS[1:]:
            directory = root if n == 4097 else previous
            stem = directory/'summaries'/('%s_N%d'%(seed,n))
            receipt = json.loads(stem.with_suffix('.receipt.json').read_text())
            assert receipt['status']=='complete' and receipt['n_trajectories']==70
            assert receipt['fit_seed']==seed and receipt['N']==n and len(receipt['sources'])==70
            assert len({int(r['task_id']) for r in receipt['sources']})==70
            assert sha(stem.with_suffix('.npz'))==receipt['indices_sha256']
            assert sha(stem.with_name(stem.name+'_repeat_diagnostics.tsv'))==receipt['diagnostics_sha256']
            with np.load(stem.with_suffix('.npz'),allow_pickle=False) as data:
                assert str(data['fit_seed'].item())==seed and int(data['N'])==n
                assert list(data['parameters'])==list(efast.ACTIVE) and list(data['outputs'])==list(efast.OUTPUTS)
                assert list(data['index_names'])==list(INDICES)
                np.testing.assert_allclose(data['oxygen'],oxygen,rtol=0,atol=1e-12)
                if n!=4097: np.testing.assert_array_equal(values[(seed,n)],data['indices'])
                values[(seed,n)]=data['indices'].copy();variance[(seed,n)]=data['variance'].copy()
    assert len(values)==24 and all(x.shape==(5,2,14,201,2) for x in values.values())
    assert len(variance)==18 and all(x.shape==(5,14,201,2) and np.isfinite(x).all() and (x>=0).all() for x in variance.values())
    return selection,values,variance,oxygen


def validate(root, values, oxygen):
    v=Verifier()
    cached={key:stats(x) for key,x in values.items()}
    repeats=read(root/'repeat_diagnostics.tsv')
    resolutions=read(root/'resolution_diagnostics.tsv')
    assert len(repeats)==4032 and len(resolutions)==3024
    seen=set()
    for row in repeats:
        seed,n=row['fit_seed'],int(row['N'])
        t,p,k=INDICES.index(row['index']),efast.ACTIVE.index(row['parameter']),efast.OUTPUTS.index(row['output'])
        _,lo,hi=next(w for w in WINDOWS if w[0]==row['oxygen_window'])
        key=(seed,n,t,p,k,lo,hi)
        assert key not in seen;seen.add(key)
        mean,sd,span,count=cached[(seed,n)]
        q,m,at,valid,passed,status=diagnostic(span[t,p,:,k],oxygen,lo,hi,.10)
        mask=(oxygen>=lo)&(oxygen<=hi)
        curves=values[(seed,n)][:,t,p,:,k]
        bands=band(curves,oxygen,lo,hi)
        peaks=np.max(curves[:,mask],axis=1)
        expected=dict(phase_range_p90=q,phase_range_max=m,phase_range_max_O2_pct=at,
            n_valid_oxygen=valid,n_expected_oxygen=int(mask.sum()),
            n_cells_with_five_defined_phases=int((count[t,p,mask,k]==5).sum()),
            phase_sd_p90=diagnostic(sd[t,p,:,k],oxygen,lo,hi,.10)[0],
            phase_band_mean=np.mean(bands),phase_band_range=np.ptp(bands),phase_band_sd=np.std(bands),
            phase_peak_range=np.ptp(peaks),threshold_absolute_index_units=.10,p90_pass=passed,
            diagnostic_status=status,key_parameter=row['parameter'] in KEYS)
        for field,value in expected.items():v.check(row[field],value,(key,field))
    seen=set()
    for row in resolutions:
        seed,a,b=row['fit_seed'],int(row['N_low']),int(row['N_high'])
        assert (a,b) in PAIRS
        t,p,k=INDICES.index(row['index']),efast.ACTIVE.index(row['parameter']),efast.OUTPUTS.index(row['output'])
        _,lo,hi=next(w for w in WINDOWS if w[0]==row['oxygen_window'])
        key=(seed,a,b,t,p,k,lo,hi);assert key not in seen;seen.add(key)
        x,_,_,cx=cached[(seed,a)];y,_,_,cy=cached[(seed,b)]
        diff=np.where((cx[t,p,:,k]==5)&(cy[t,p,:,k]==5),np.abs(y[t,p,:,k]-x[t,p,:,k]),np.nan)
        q,m,at,valid,passed,status=diagnostic(diff,oxygen,lo,hi,.05)
        mask=(oxygen>=lo)&(oxygen<=hi);complete=valid==int(mask.sum())
        expected=dict(absolute_delta_p90=q,absolute_delta_max=m,absolute_delta_max_O2_pct=at,
            n_valid_oxygen=valid,n_expected_oxygen=int(mask.sum()),p90_pass=passed,diagnostic_status=status,
            threshold_absolute_index_units=.05,primary_comparison=(a,b)==(2049,4097),key_parameter=row['parameter'] in KEYS,
            band_delta=band(y[t,p,:,k],oxygen,lo,hi)-band(x[t,p,:,k],oxygen,lo,hi) if complete else None,
            peak_delta=np.max(y[t,p,mask,k])-np.max(x[t,p,mask,k]) if complete else None,
            peak_O2_low=oxygen[mask][np.argmax(x[t,p,mask,k])] if complete else None,
            peak_O2_high=oxygen[mask][np.argmax(y[t,p,mask,k])] if complete else None)
        for field,value in expected.items():v.check(row[field],value,(key,field))
    curve_count=0
    curve_seen={key:np.zeros((2,14,201,2),dtype=bool) for key in values}
    with gzip.open(root/'oxygen_curve_summary.tsv.gz','rt',newline='') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=(row['fit_seed'],int(row['N']))
            t,p,k=INDICES.index(row['index']),efast.ACTIVE.index(row['parameter']),efast.OUTPUTS.index(row['output'])
            o=float(row['O2_pct']);j=int(round(o/.025));assert abs(oxygen[j]-o)<1e-12
            assert not curve_seen[key][t,p,j,k], 'Duplicate oxygen-curve source row'
            curve_seen[key][t,p,j,k]=True
            mean,sd,span,count=cached[key]
            for field,x in (('phase_mean',mean),('phase_sd',sd),('phase_range',span),('n_defined_phases',count)):
                v.check(row[field],x[t,p,j,k],(key,t,p,j,k,field))
            curve_count+=1
    assert curve_count==270144
    assert all(x.all() for x in curve_seen.values())
    # Verify each endpoint's band scores and submitted deterministic parameter ranks.
    bands_table=read(root/'band_rank_mechanism_comparison.tsv');assert len(bands_table)==1344
    lookup={};band_lookup={}
    for row in bands_table:
        key=(row['fit_seed'],int(row['N']));t,k=INDICES.index(row['index']),efast.OUTPUTS.index(row['output'])
        p=efast.ACTIVE.index(row['parameter'])
        band_key=(key,t,k,p);assert band_key not in band_lookup
        band_lookup[band_key]=row
        if (key,t,k) not in lookup:lookup[(key,t,k)]=mechanism(values[key],oxygen,t,k)
        summary=lookup[(key,t,k)]
        expected=dict(low_band_mean=summary['low'][p],high_band_mean=summary['high'][p],
            high_minus_low=summary['high'][p]-summary['low'][p],
            low_band_rank=summary['low_ranks'][p],high_band_rank=summary['high_ranks'][p],
            n_oxygen_with_five_phases=int((cached[key][3][t,p,:,k]==5).sum()),
            low_death_largest_group=summary['conditions'][0],high_missegregation_increase=summary['conditions'][1],
            high_buffering_increase=summary['conditions'][2],noemi_all_three_supported=summary['supported'])
        for field,x in expected.items():v.check(row[field],x,(key,t,k,p,field))
    stability=read(root/'band_rank_mechanism_stability.tsv');assert len(stability)==1008
    stability_seen=set()
    for row in stability:
        seed,a,b=row['fit_seed'],int(row['N_low']),int(row['N_high'])
        t,p,k=INDICES.index(row['index']),efast.ACTIVE.index(row['parameter']),efast.OUTPUTS.index(row['output'])
        key=(seed,a,b,t,p,k);assert key not in stability_seen;stability_seen.add(key)
        x,y=band_lookup[((seed,a),t,k,p)],band_lookup[((seed,b),t,k,p)]
        expected=dict(high_minus_low_at_lower_N=float(x['high_minus_low']),high_minus_low_at_higher_N=float(y['high_minus_low']),
            contrast_direction_unchanged=bool(np.sign(float(x['high_minus_low']))==np.sign(float(y['high_minus_low']))),
            low_band_rank_delta=int(y['low_band_rank'])-int(x['low_band_rank']),
            high_band_rank_delta=int(y['high_band_rank'])-int(x['high_band_rank']),
            noemi_support_at_lower_N=x['noemi_all_three_supported'],noemi_support_at_higher_N=y['noemi_all_three_supported'],
            noemi_support_unchanged=x['noemi_all_three_supported']==y['noemi_all_three_supported'])
        for field,value in expected.items():v.check(row[field],value,(key,field))
    # Focus subsets must be exact subsets, not independently rounded copies.
    for name,parent in [('key_ploidy_ST_repeat_diagnostics.tsv',repeats),('key_ploidy_ST_resolution_diagnostics.tsv',resolutions)]:
        expected=[r for r in parent if r['parameter'] in KEYS and r['index']=='ST' and r['output']=='dominant_mean_ploidy']
        assert read(root/name)==expected
    return dict(status='passed',independent_numeric_or_categorical_checks=v.n,
        maximum_absolute_numeric_difference=v.max_error,repeat_rows=len(repeats),resolution_rows=len(resolutions),
        curve_rows=curve_count,band_mechanism_rows=len(bands_table),stability_rows=len(stability),n_source_arrays=len(values),
        undefined_index_cells=sum(int((~np.isfinite(x)).sum()) for x in values.values()),
        min_index=min(float(np.nanmin(x)) for x in values.values()),
        max_index=max(float(np.nanmax(x)) for x in values.values()),
        scope='NPZ-to-table reconciliation; no new model evaluations or FAST resampling'),cached


def boolean(value):
    assert value in ('True','False',True,False), value
    return value is True or value == 'True'


def detail_rows(repeats, resolutions):
    rmap = {(r['fit_seed'],int(r['N']),r['parameter'],r['oxygen_window']):r for r in repeats
            if r['index']=='ST' and r['output']==efast.OUTPUTS[0]}
    dmap = {(r['fit_seed'],int(r['N_low']),int(r['N_high']),r['parameter'],r['oxygen_window']):r for r in resolutions
            if r['index']=='ST' and r['output']==efast.OUTPUTS[0]}
    records = []
    for seed in SEEDS:
        for p in efast.ACTIVE:
            for window,_,_ in WINDOWS:
                row = dict(fit_seed=seed,parameter=p,group=efast.GROUP[p],key_parameter=p in KEYS,oxygen_window=window)
                for n in NS:
                    r = rmap[(seed,n,p,window)]
                    for field in ('phase_range_p90','phase_range_max','phase_range_max_O2_pct','p90_pass','diagnostic_status'):
                        row['N%d_%s'%(n,field)] = r[field]
                for a,b in PAIRS:
                    r = dmap[(seed,a,b,p,window)]
                    for field in ('absolute_delta_p90','absolute_delta_max','absolute_delta_max_O2_pct','p90_pass','diagnostic_status'):
                        row['%d_to_%d_%s'%(a,b,field)] = r[field]
                    row['N%d_both_p90_pass'%b] = boolean(row['N%d_p90_pass'%b]) and boolean(row['%d_to_%d_p90_pass'%(a,b)])
                    row['N%d_both_max_pass'%b] = row['N%d_diagnostic_status'%b]=='passed_by_threshold' and row['%d_to_%d_diagnostic_status'%(a,b)]=='passed_by_threshold'
                row['repeat_p90_improved_2049_to_4097'] = float(row['N4097_phase_range_p90']) < float(row['N2049_phase_range_p90'])
                row['delta_p90_improved_vs_previous_step'] = float(row['2049_to_4097_absolute_delta_p90']) < float(row['1025_to_2049_absolute_delta_p90'])
                records.append(row)
    return records


def make_tables(root, selection, values, oxygen):
    repeats = read(root/'repeat_diagnostics.tsv'); resolutions = read(root/'resolution_diagnostics.tsv')
    details = detail_rows(repeats,resolutions)
    focus = [r for r in details if r['key_parameter']]
    write(root/'ploidy_ST_parameter_details.tsv',details)
    write(root/'ploidy_ST_key_parameter_details.tsv',focus)
    summary = []
    for seed in SEEDS:
        f = [r for r in focus if r['fit_seed']==seed and r['oxygen_window']==WINDOWS[0][0]]
        assert len(f)==5
        row = dict(fit_seed=seed,n_key_parameters=5)
        for n in NS:
            row['N%d_repeat_p90_pass_of_5'%n] = sum(boolean(r['N%d_p90_pass'%n]) for r in f)
            row['N%d_repeat_max_pass_of_5'%n] = sum(r['N%d_diagnostic_status'%n]=='passed_by_threshold' for r in f)
            row['N%d_worst_repeat_p90'%n] = max(float(r['N%d_phase_range_p90'%n]) for r in f)
        for a,b in PAIRS:
            row['%d_to_%d_delta_p90_pass_of_5'%(a,b)] = sum(boolean(r['%d_to_%d_p90_pass'%(a,b)]) for r in f)
            row['%d_to_%d_delta_max_pass_of_5'%(a,b)] = sum(r['%d_to_%d_diagnostic_status'%(a,b)]=='passed_by_threshold' for r in f)
            row['N%d_both_p90_pass_of_5'%b] = sum(r['N%d_both_p90_pass'%b] for r in f)
            row['N%d_both_max_pass_of_5'%b] = sum(r['N%d_both_max_pass'%b] for r in f)
        worst = max(f,key=lambda r:float(r['2049_to_4097_absolute_delta_p90']))
        row.update(worst_resolution_parameter=worst['parameter'],worst_resolution_p90=float(worst['2049_to_4097_absolute_delta_p90']),
            improved_repeat_of_5=sum(r['repeat_p90_improved_2049_to_4097'] for r in f),
            improved_delta_of_5=sum(r['delta_p90_improved_vs_previous_step'] for r in f))
        summary.append(row)
    write(root/'seed_convergence_summary.tsv',summary)
    overall = []
    for kind,table,nfields in [('repeat',repeats,['N']),('resolution',resolutions,['N_low','N_high'])]:
        grouped = {}
        for r in table:
            key = tuple(r[f] for f in [*nfields,'index','output','oxygen_window'])
            grouped.setdefault(key,[]).append(r)
        for key,items in grouped.items():
            for scope in ('all14','key5'):
                rows = [r for r in items if scope=='all14' or r['parameter'] in KEYS]
                row = dict(diagnostic=kind,parameter_scope=scope,N_low=key[0],N_high=key[1] if kind=='resolution' else '',
                    index=key[-3],output=key[-2],oxygen_window=key[-1],n_parameter_seed_cells=len(rows),
                    p90_pass=sum(boolean(r['p90_pass']) for r in rows))
                for label in ('passed_by_threshold','localized_review','not_passed','undefined_or_incomplete'):
                    row[label] = sum(r['diagnostic_status']==label for r in rows)
                overall.append(row)
    write(root/'convergence_pass_summary.tsv',overall)
    mechanisms,phase_mechanisms,profiles = [],[],[]
    for seed in SEEDS:
        for n in NS:
            x = values[(seed,n)]
            for t,index in enumerate(INDICES):
                for k,output in enumerate(efast.OUTPUTS):
                    m = mechanism(x,oxygen,t,k)
                    mechanisms.append(dict(fit_seed=seed,N=n,index=index,output=output,
                        low_death_largest=m['conditions'][0],high_missegregation_increase=m['conditions'][1],
                        high_buffering_increase=m['conditions'][2],all_three_supported=m['supported'],
                        low_dominant_parameter=efast.ACTIVE[int(np.nanargmax(m['low']))],
                        high_dominant_parameter=efast.ACTIVE[int(np.nanargmax(m['high']))],
                        low_top_group=max(m['groups']['low'],key=m['groups']['low'].get),high_top_group=max(m['groups']['high'],key=m['groups']['high'].get),
                        low_death_mean=m['groups']['low']['death'],high_death_mean=m['groups']['high']['death'],
                        missegregation_high_minus_low=m['groups']['high']['missegregation']-m['groups']['low']['missegregation'],
                        buffering_high_minus_low=m['groups']['high']['buffering']-m['groups']['low']['buffering']))
                    for r in range(5):
                        phase = mechanism(np.repeat(x[r:r+1],5,axis=0),oxygen,t,k)
                        phase_mechanisms.append(dict(fit_seed=seed,N=n,replicate=r+1,index=index,output=output,
                            low_death_largest=phase['conditions'][0],high_missegregation_increase=phase['conditions'][1],
                            high_buffering_increase=phase['conditions'][2],all_three_supported=phase['supported']))
            for parameter in KEYS:
                for r in range(5):
                    y = x[r,1,efast.ACTIVE.index(parameter),:,0]
                    profiles.append(dict(fit_seed=seed,N=n,parameter=parameter,replicate=r+1,
                        whole_band_mean=float(band(y,oxygen,0,5)),low_band_mean=float(band(y,oxygen,0,1)),
                        high_band_mean=float(band(y,oxygen,3,5)),maximum_ST=float(np.nanmax(y)),
                        peak_O2_pct=float(oxygen[np.nanargmax(y)])))
    write(root/'mechanism_endpoint_summary.tsv',mechanisms)
    write(root/'mechanism_phase_summary.tsv',phase_mechanisms)
    write(root/'key_parameter_phase_profiles.tsv',profiles)
    audit = json.loads((root/'slurm/completed_hpc_audit.json').read_text())
    numerical = []
    for entry in audit['recovery_audits']:
        for row in entry['corrections']['reports']:
            numerical.append(dict(fit_seed=entry['fit_seed'],N=4097,focal_parameter=entry['focal_parameter'],replicate=entry['replicate'],
                sample_id=row['sample_id'],O2_pct=row['O2_pct'],spectral_gap=float(row['original']['spectral_gap']),
                ploidy_change=row['ploidy_change'],growth_change=row['growth_change'],double_residual=row['double_residual'],
                leading_root_bound_width=row['leading_root_bound_width'],ploidy_precision_difference=row['ploidy_precision_difference'],
                growth_precision_difference=row['growth_precision_difference']))
    write(root/'numerical_recovery_summary.tsv',numerical)
    # Compare the same six endpoints in the historical and current accounting.
    previous = root.parent/'N1025_N2049'
    old_audit = json.loads((previous/'slurm/completed_hpc_audit.json').read_text())
    old_tasks = {r['task_id']:r for r in read(previous/'slurm/task_manifest.tsv')}
    observations = []
    for a,native_job in [(old_audit,'20833838'),(audit,'20864027')]:
        for row in a['slurm_unique_jobs']:
            if not row['job_id'].startswith(native_job+'_'): continue
            task = old_tasks[row['job_id'].split('_')[1]] if native_job=='20833838' else read_task(root,int(row['job_id'].split('_')[1]))
            if task['fit_seed'] in SEEDS:
                observations.append(dict(fit_seed=task['fit_seed'],N=int(task['N']),task_id=int(task['task_id']),
                    elapsed_hours=row['elapsed_seconds']/3600,job_id=row['job_id']))
    write(root/'runtime_observations.tsv',observations)
    per_seed,costs = [],[]
    for n in NS[1:]:
        all_times = [r['elapsed_hours'] for r in observations if r['N']==n]
        assert len(all_times)==420
        for seed in SEEDS:
            times = [r['elapsed_hours'] for r in observations if r['N']==n and r['fit_seed']==seed]
            assert len(times)==70
            per_seed.append(dict(fit_seed=seed,N=n,n_tasks=70,median_hours=float(np.median(times)),
                p95_hours=float(np.quantile(times,.95)),maximum_hours=float(max(times))))
        costs.append(dict(N=n,type='observed_same_six_seeds',n_tasks=420,median_hours=float(np.median(all_times)),
            p95_hours=float(np.quantile(all_times,.95)),maximum_hours=float(max(all_times)),
            n_seeds=6,n_phases=5,n_parameters=14,n_oxygen=201,operator_evaluation_rows=420*n*201,new_jobs_submitted=False))
    base = costs[-1]; factor = 8193/4097
    costs.append(dict(N=8193,type='linear_projection_not_measured',n_tasks=420,
        median_hours=base['median_hours']*factor,p95_hours=base['p95_hours']*factor,maximum_hours=base['maximum_hours']*factor,
        n_seeds=6,n_phases=5,n_parameters=14,n_oxygen=201,operator_evaluation_rows=420*8193*201,new_jobs_submitted=False))
    write(root/'runtime_by_seed.tsv',per_seed);write(root/'runtime_and_N_projection.tsv',costs)
    comparison = []
    overall_focus = [r for r in focus if r['oxygen_window']==WINDOWS[0][0]]
    for a,b in ((1025,2049),(2049,4097)):
        counts = [sum(boolean(r['N%d_p90_pass'%b]) for r in overall_focus),
                  sum(boolean(r['%d_to_%d_p90_pass'%(a,b)]) for r in overall_focus),
                  sum(r['N%d_both_p90_pass'%b] for r in overall_focus),sum(r['N%d_both_max_pass'%b] for r in overall_focus)]
        for criterion,count in zip(('repeat_p90','resolution_p90','both_p90','both_maximum'),counts):
            comparison.append(dict(N=b,comparison_N=a,criterion=criterion,n_passed=count,n_parameter_seed_pairs=30,
                fit_seeds=','.join(SEEDS),parameter_scope='key5',oxygen_window='overall_0_to_5pct',repeat_threshold=.10,resolution_threshold=.05))
    write(root/'ploidy_ST_pass_comparison.tsv',comparison)
    return summary,focus,overall,mechanisms


_TASK_CACHE = {}
def read_task(root, task_id):
    if root not in _TASK_CACHE: _TASK_CACHE[root] = read(root/'slurm/task_manifest.tsv')
    return _TASK_CACHE[root][task_id-1]


def plot_pass_counts(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':10})
    records = read(root/'ploidy_ST_pass_comparison.tsv')
    fig,ax = plt.subplots(figsize=(10,5.4))
    for j,(n,color,hatch,label) in enumerate(((2049,'#aaaaaa','//','N2049 (vs N1025)'),(4097,'#087e8b','','N4097 (vs N2049)'))):
        rows = [r for r in records if int(r['N'])==n]
        y = [int(r['n_passed']) for r in rows];x = np.arange(4)+(j-.5)*.34
        ax.bar(x,y,width=.31,color=color,edgecolor='#303030',linewidth=.7,hatch=hatch,label=label)
        for p,value in zip(x,y):ax.text(p,value+.6,str(value),ha='center',fontsize=10)
    ax.set(xticks=range(4),xticklabels=['Repeat p90','Resolution p90','Both p90','Both maxima'],ylim=(0,33),
        ylabel='Passing parameter / seed pairs (of 30)',title='Ploidy ST convergence checks in the same six fitted neighborhoods')
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    ax.legend(frameon=False,loc='upper right',fontsize=9)
    fig.text(.12,.03,'Five key parameters x six endpoints; 201 oxygen points (0-5%). Repeat threshold 0.10; resolution threshold 0.05.\n"Both maxima" requires the full-curve maxima of both diagnostics to pass. These are practical checks, not confidence levels.',fontsize=8)
    fig.tight_layout(rect=(0,.13,1,1))
    for ext in ('pdf','png'):fig.savefig(root/'figures'/('ploidy_ST_pass_comparison.'+ext),dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    args = parser.parse_args();root = args.root.resolve()
    selection,values,variance,oxygen = load(root)
    validation,cached = validate(root,values,oxygen)
    gaps = read(root/'spectral_gap_summary.tsv')
    assert len(gaps)==18 and sum(int(r['n_evaluations']) for r in gaps)==345868740
    assert len({(r['fit_seed'],r['oxygen_band']) for r in gaps})==18
    for r in gaps:
        assert r['fit_seed'] in SEEDS and int(r['N'])==4097 and int(r['n_trajectories'])==70
        size={'low_0_to_1pct':41,'middle_1_to_3pct':79,'high_3_to_5pct':81}[r['oxygen_band']]
        assert int(r['n_evaluations'])==70*4097*size
        assert 0<=float(r['fraction_gap_lt_1e_6'])<=float(r['fraction_gap_lt_1e_4'])<=float(r['fraction_gap_lt_1e_3'])<=1
        assert np.isfinite(float(r['minimum_spectral_gap'])) and float(r['minimum_spectral_gap'])>=0
    summary,focus,overall,mechanisms = make_tables(root,selection,values,oxygen)
    plot_pass_counts(root)
    audit = json.loads((root/'slurm/completed_hpc_audit.json').read_text())
    assert audit['source_collector_sha256']==sha(Path(__file__).with_name('collect_resolution_4097_results.py'))
    validation.update(python=platform.python_version(),numpy=np.__version__,script_sha256=sha(Path(__file__)),
        independent_helpers_sha256=sha(Path(__file__).with_name('report_resolution_pilot.py')),
        collector_sha256=sha(Path(__file__).with_name('collect_resolution_4097_results.py')),
        completed_manifest_sha256=sha(root/'completed_results_manifest.tsv'),hpc_audit_sha256=sha(root/'slurm/completed_hpc_audit.json'),
        n_variance_cells=sum(x.size for x in variance.values()),min_output_variance=min(float(x.min()) for x in variance.values()),
        n_seed_parameter_key_pairs=30,n_all_parameter_seed_pairs=84,spectral_gap_rows=18,new_jobs_submitted=False)
    (root/'analysis_validation.json').write_text(json.dumps(validation,indent=2,sort_keys=True)+'\n')
    print(json.dumps(validation,indent=2))
    print('SEED_SUMMARY '+json.dumps(summary,indent=2))


if __name__ == '__main__': main()
