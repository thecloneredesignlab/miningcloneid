#!/usr/bin/env python3
"""Independently audit completed pilot tables against five-phase NPZ sources.

Uses NumPy and Matplotlib; never evaluates the model or submits a job.
The frozen calculation module resolution_pilot.py is not imported or changed.
"""
import argparse
import csv
import gzip
import hashlib
import json
import platform
from pathlib import Path
import warnings

import numpy as np
import efast

KEYS = ('buffer_beta', 'buffer_smax', 'p_mis_base', 'mu_hp', 'O2_crit')
NS = (513, 1025, 2049)
WINDOWS = (('overall_0_to_5pct', 0., 5.), ('low_0_to_1pct', 0., 1.), ('high_3_to_5pct', 3., 5.))
INDICES = ('S1', 'ST')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows):
    efast.write_table(path, list(rows[0]), rows)


def stats(x):
    count = np.isfinite(x).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        mean, sd = np.nanmean(x, axis=0), np.nanstd(x, axis=0, ddof=0)
        span = np.nanmax(x, axis=0) - np.nanmin(x, axis=0)
    return mean, np.where(count == 5, sd, np.nan), np.where(count == 5, span, np.nan), count


def band(x, oxygen, lo, hi):
    selected = (oxygen >= lo) & (oxygen <= hi)
    # Independent integration, not the submitted weighted-sum helper.
    return np.trapezoid(x[..., selected], x=oxygen[selected], axis=-1) / (hi-lo)


def diagnostic(curve, oxygen, lo, hi, threshold):
    mask = (oxygen >= lo) & (oxygen <= hi)
    valid = mask & np.isfinite(curve)
    ids = np.flatnonzero(valid)
    if not len(ids):
        return None, None, None, 0, False, 'undefined_or_incomplete'
    q = float(np.quantile(curve[valid], .9))
    at = ids[np.argmax(curve[valid])]
    maximum = float(curve[at])
    complete = int(valid.sum()) == int(mask.sum())
    passed = complete and q <= threshold
    status = 'undefined_or_incomplete' if not complete else (
        'not_passed' if not passed else 'localized_review' if maximum > threshold else 'passed_by_threshold')
    return q, maximum, float(oxygen[at]), int(valid.sum()), passed, status


class Verifier:
    def __init__(self):
        self.n = 0
        self.max_error = 0.

    def check(self, actual, expected, context):
        self.n += 1
        if isinstance(expected, (str, bool)):
            if str(actual) != str(expected):
                raise AssertionError((context, actual, expected))
        elif expected is None or not np.isfinite(expected):
            if actual not in ('', 'None', 'nan', 'NaN', 'NA'):
                raise AssertionError((context, actual, expected))
        else:
            value = float(actual)
            if not np.isfinite(value) or not np.isclose(value, expected, rtol=1e-12, atol=1e-12):
                raise AssertionError((context, actual, expected))
            self.max_error = max(self.max_error, abs(value-float(expected)))


def load(root):
    selection = read(root/'seed_selection.tsv')
    seeds = [r['fit_seed'] for r in selection]
    assert len(seeds) == len(set(seeds)) == 12
    values, variance = {}, {}
    plan = json.loads((root/'slurm/submission_plan.json').read_text())
    assert sha(root/'baseline_N513.npz') == plan['baseline_selected_indices_sha256']
    with np.load(root/'baseline_N513.npz') as data:
        assert list(data['fit_seeds']) == seeds
        oxygen = data['oxygen'].copy()
        assert list(data['parameters']) == list(efast.ACTIVE)
        assert list(data['outputs']) == list(efast.OUTPUTS)
        assert list(data['index_names']) == list(INDICES)
        for seed, x in zip(seeds, data['indices']):
            values[(seed,513)] = x.copy()
    for seed in seeds:
        for n in NS[1:]:
            prefix=root/'summaries'/('%s_N%d'%(seed,n))
            receipt=json.loads(prefix.with_suffix('.receipt.json').read_text())
            assert receipt['status']=='complete' and receipt['n_trajectories']==70
            assert receipt['fit_seed']==seed and receipt['N']==n
            assert len(receipt['sources'])==70 and len({r['task_id'] for r in receipt['sources']})==70
            assert sha(prefix.with_suffix('.npz'))==receipt['indices_sha256']
            assert sha(prefix.with_name(prefix.name+'_repeat_diagnostics.tsv'))==receipt['diagnostics_sha256']
            with np.load(prefix.with_suffix('.npz')) as data:
                assert data['fit_seed'].item()==seed and data['N'].item()==n
                assert list(data['parameters'])==list(efast.ACTIVE)
                assert list(data['outputs'])==list(efast.OUTPUTS)
                assert list(data['index_names'])==list(INDICES)
                np.testing.assert_allclose(data['oxygen'],oxygen,atol=1e-12,rtol=0)
                values[(seed,n)]=data['indices'].copy()
                variance[(seed,n)]=data['variance'].copy()
    np.testing.assert_allclose(oxygen,np.arange(201)*.025,atol=1e-12,rtol=0)
    assert all(x.shape==(5,2,14,201,2) for x in values.values())
    assert all(x.shape==(5,14,201,2) and np.isfinite(x).all() and (x>=0).all() for x in variance.values())
    return selection, values, variance, oxygen


def validate(root, values, oxygen):
    v=Verifier()
    cached={key:stats(x) for key,x in values.items()}
    repeats=read(root/'repeat_diagnostics.tsv')
    resolutions=read(root/'resolution_diagnostics.tsv')
    assert len(repeats)==6048 and len(resolutions)==4032
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
        assert (a,b) in ((513,1025),(1025,2049))
        t,p,k=INDICES.index(row['index']),efast.ACTIVE.index(row['parameter']),efast.OUTPUTS.index(row['output'])
        _,lo,hi=next(w for w in WINDOWS if w[0]==row['oxygen_window'])
        key=(seed,a,b,t,p,k,lo,hi);assert key not in seen;seen.add(key)
        x,_,_,cx=cached[(seed,a)];y,_,_,cy=cached[(seed,b)]
        diff=np.where((cx[t,p,:,k]==5)&(cy[t,p,:,k]==5),np.abs(y[t,p,:,k]-x[t,p,:,k]),np.nan)
        q,m,at,valid,passed,status=diagnostic(diff,oxygen,lo,hi,.05)
        mask=(oxygen>=lo)&(oxygen<=hi);complete=valid==int(mask.sum())
        expected=dict(absolute_delta_p90=q,absolute_delta_max=m,absolute_delta_max_O2_pct=at,
            n_valid_oxygen=valid,n_expected_oxygen=int(mask.sum()),p90_pass=passed,diagnostic_status=status,
            threshold_absolute_index_units=.05,primary_comparison=a==1025,key_parameter=row['parameter'] in KEYS,
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
    assert curve_count==405216
    assert all(x.all() for x in curve_seen.values())
    # Verify each endpoint's band scores and submitted deterministic parameter ranks.
    bands_table=read(root/'band_rank_mechanism_comparison.tsv');assert len(bands_table)==2016
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
    stability=read(root/'band_rank_mechanism_stability.tsv');assert len(stability)==1344
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


def mechanism(x,oxygen,t,k):
    mean,_,_,count=stats(x)
    mean=np.where(count==5,mean,np.nan)
    low=band(mean[t,:,:,k],oxygen,0,1);high=band(mean[t,:,:,k],oxygen,3,5)
    groups={}
    for label,curve in [('low',low),('high',high)]:
        groups[label]={g:float(np.mean([curve[p] for p,name in enumerate(efast.ACTIVE) if efast.GROUP[name]==g]))
                       for g in sorted(set(efast.GROUP.values()))}
    a,b=groups['low'],groups['high'];defined=np.isfinite(list(a.values())+list(b.values())).all()
    conditions=[bool(a['death']>=max(a.values())),bool(b['missegregation']>a['missegregation']),bool(b['buffering']>a['buffering'])] if defined else [None]*3
    # Stable sort is the documented tie policy of the published tables.
    rank=lambda c:np.argsort(np.argsort(-np.where(np.isfinite(c),c,-np.inf),kind='stable'),kind='stable')+1
    return dict(low=low,high=high,groups=groups,conditions=conditions,supported=all(conditions) if defined else None,
                low_ranks=rank(low),high_ranks=rank(high))


def make_tables(root,selection,values,oxygen,cached):
    repeats=read(root/'repeat_diagnostics.tsv');resolutions=read(root/'resolution_diagnostics.tsv')
    rlookup={(r['fit_seed'],int(r['N']),r['parameter'],r['oxygen_window'],r['index'],r['output']):r for r in repeats}
    dlookup={(r['fit_seed'],int(r['N_low']),r['parameter'],r['oxygen_window'],r['index'],r['output']):r for r in resolutions}
    focus=[];summary=[];mechanisms=[];mechanism_phases=[];overall=[];numerical=[]
    for endpoint in selection:
        seed=endpoint['fit_seed']
        for name in KEYS:
            row=dict(fit_seed=seed,parameter=name,group=efast.GROUP[name],pilot_order=int(endpoint['pilot_order']))
            for n in NS:
                r=rlookup[(seed,n,name,WINDOWS[0][0],'ST',efast.OUTPUTS[0])]
                for field in ('phase_range_p90','phase_range_max','phase_range_max_O2_pct','p90_pass','diagnostic_status'):
                    row['N%d_%s'%(n,field)]=r[field]
            for n in NS[:-1]:
                r=dlookup[(seed,n,name,WINDOWS[0][0],'ST',efast.OUTPUTS[0])]
                for field in ('absolute_delta_p90','absolute_delta_max','absolute_delta_max_O2_pct','p90_pass','diagnostic_status'):
                    row['%d_to_%d_%s'%(n,1025 if n==513 else 2049,field)]=r[field]
            row['both_repeat_and_resolution_p90_pass']=row['N2049_p90_pass']=='True' and row['1025_to_2049_p90_pass']=='True'
            row['both_repeat_and_resolution_max_pass']=row['N2049_diagnostic_status']=='passed_by_threshold' and row['1025_to_2049_diagnostic_status']=='passed_by_threshold'
            row['phase_p90_improved_1025_to_2049']=float(row['N2049_phase_range_p90'])<float(row['N1025_phase_range_p90'])
            focus.append(row)
        f=[r for r in focus if r['fit_seed']==seed]
        s=dict(fit_seed=seed,pilot_order=int(endpoint['pilot_order']))
        for n in NS:
            s['N%d_repeat_p90_pass_of_5'%n]=sum(r['N%d_p90_pass'%n]=='True' for r in f)
            s['N%d_repeat_max_pass_of_5'%n]=sum(r['N%d_diagnostic_status'%n]=='passed_by_threshold' for r in f)
            s['N%d_worst_repeat_p90'%n]=max(float(r['N%d_phase_range_p90'%n]) for r in f)
        s['resolution_p90_pass_of_5']=sum(r['1025_to_2049_p90_pass']=='True' for r in f)
        s['resolution_max_pass_of_5']=sum(r['1025_to_2049_diagnostic_status']=='passed_by_threshold' for r in f)
        worst=max(f,key=lambda r:float(r['1025_to_2049_absolute_delta_p90']))
        s['worst_resolution_parameter']=worst['parameter'];s['worst_resolution_p90']=float(worst['1025_to_2049_absolute_delta_p90'])
        s['both_p90_pass_of_5']=sum(r['both_repeat_and_resolution_p90_pass'] for r in f)
        s['both_max_pass_of_5']=sum(r['both_repeat_and_resolution_max_pass'] for r in f)
        s['improved_repeat_of_5']=sum(r['phase_p90_improved_1025_to_2049'] for r in f)
        summary.append(s)
        for n in NS:
            for t,index in enumerate(INDICES):
                for k,output in enumerate(efast.OUTPUTS):
                    result=mechanism(values[(seed,n)],oxygen,t,k)
                    mechanisms.append(dict(fit_seed=seed,N=n,index=index,output=output,
                        low_death_largest=result['conditions'][0],high_missegregation_increase=result['conditions'][1],
                        high_buffering_increase=result['conditions'][2],all_three_supported=result['supported'],
                        low_dominant_parameter=efast.ACTIVE[int(np.nanargmax(result['low']))],
                        high_dominant_parameter=efast.ACTIVE[int(np.nanargmax(result['high']))],
                        low_top_group=max(result['groups']['low'],key=result['groups']['low'].get),
                        high_top_group=max(result['groups']['high'],key=result['groups']['high'].get),
                        low_death_mean=result['groups']['low']['death'],high_death_mean=result['groups']['high']['death'],
                        missegregation_high_minus_low=result['groups']['high']['missegregation']-result['groups']['low']['missegregation'],
                        buffering_high_minus_low=result['groups']['high']['buffering']-result['groups']['low']['buffering']))
                    for r in range(5):
                        # Repeat a phase along the phase axis so the same integration
                        # and process definition can be used for phase-level robustness.
                        phase=mechanism(np.repeat(values[(seed,n)][r:r+1],5,axis=0),oxygen,t,k)
                        mechanism_phases.append(dict(fit_seed=seed,N=n,replicate=r+1,index=index,output=output,
                            low_death_largest=phase['conditions'][0],high_missegregation_increase=phase['conditions'][1],
                            high_buffering_increase=phase['conditions'][2],all_three_supported=phase['supported']))
    for kind,table,nfields in [('repeat',repeats,['N']),('resolution',resolutions,['N_low','N_high'])]:
        grouped={}
        for row in table:
            key=tuple(row[f] for f in [*nfields,'index','output','oxygen_window'])
            grouped.setdefault(key,[]).append(row)
        for key,items in grouped.items():
            for scope in ('all14','key5'):
                rows=[r for r in items if scope=='all14' or r['parameter'] in KEYS]
                record=dict(diagnostic=kind,parameter_scope=scope,N_low=key[0],N_high=key[1] if kind=='resolution' else '',
                    index=key[-3],output=key[-2],oxygen_window=key[-1],n_parameter_seed_cells=len(rows),
                    p90_pass=sum(r['p90_pass']=='True' for r in rows))
                for status in ('passed_by_threshold','localized_review','not_passed','undefined_or_incomplete'):
                    record[status]=sum(r['diagnostic_status']==status for r in rows)
                overall.append(record)
    for name,rows in [('ploidy_ST_key_parameter_details.tsv',focus),('seed_convergence_summary.tsv',summary),
                      ('mechanism_endpoint_summary.tsv',mechanisms),('mechanism_phase_summary.tsv',mechanism_phases),
                      ('convergence_pass_summary.tsv',overall)]:write(root/name,rows)
    audit=json.loads((root/'slurm/completed_hpc_audit.json').read_text())
    for entry in audit['recovery_audits']:
        for row in entry['corrections']['reports']:
            numerical.append(dict(fit_seed=entry['fit_seed'],N=entry['N'],focal_parameter=entry['focal_parameter'],replicate=entry['replicate'],
                sample_id=row['sample_id'],O2_pct=row['O2_pct'],spectral_gap=float(row['original']['spectral_gap']),
                ploidy_change=row['ploidy_change'],growth_change=row['growth_change'],double_residual=row['double_residual'],
                leading_root_bound_width=row['leading_root_bound_width'],ploidy_precision_difference=row['ploidy_precision_difference'],
                growth_precision_difference=row['growth_precision_difference']))
    write(root/'numerical_recovery_summary.tsv',numerical)
    taskmap={str(r['task_id']):int(r['N']) for r in read(root/'slurm/task_manifest.tsv')}
    times={n:[r['elapsed_seconds']/3600 for r in audit['slurm_unique_jobs']
              if r['job_id'].startswith('20833838_') and taskmap[r['job_id'].split('_')[-1]]==n] for n in NS[1:]}
    costs=[]
    for n in (1025,2049,4097,8193):
        base=times[n] if n in times else times[2049]
        factor=1. if n in times else n/2049
        costs.append(dict(N=n,type='observed' if n in times else 'linear_projection_not_measured',
            n_observed_trajectory_jobs=len(base),trajectory_median_hours=float(np.median(base)*factor),
            trajectory_p95_hours=float(np.quantile(base,.95)*factor),trajectory_max_hours=float(np.max(base)*factor),
            n_seeds_proposed=6 if n>2049 else 12,n_phases=5,n_sampled_parameters=14,n_oxygen=201,
            projected_operator_evaluations=(6 if n>2049 else 12)*5*14*n*201,
            new_jobs_submitted=False))
    write(root/'runtime_and_N_projection.tsv',costs)
    return focus,summary,overall,mechanisms


def plots(root,selection,focus):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'font.family':'DejaVu Sans'})
    seeds=[r['fit_seed'] for r in selection]
    lookup={(r['fit_seed'],r['parameter']):r for r in focus}
    fig,axes=plt.subplots(1,3,figsize=(13.5,6.7),sharey=True)
    cmap=ListedColormap(['#EEEEEE','#B8A2D3','#E8A45F','#FFF3D0'])
    panels=[('N1025_phase_range_p90','N1025_diagnostic_status','A  Phase repeat range, N=1025',.10),
            ('N2049_phase_range_p90','N2049_diagnostic_status','B  Phase repeat range, N=2049',.10),
            ('1025_to_2049_absolute_delta_p90','1025_to_2049_diagnostic_status','C  Mean-curve change, 1025 to 2049',.05)]
    statuses={'passed_by_threshold':0,'localized_review':1,'not_passed':2,'undefined_or_incomplete':3}
    for ax,(field,status,title,threshold) in zip(axes,panels):
        z=np.array([[statuses[lookup[(seed,p)][status]] for p in KEYS] for seed in seeds])
        ax.imshow(z,aspect='auto',cmap=cmap,vmin=-.5,vmax=3.5,interpolation='nearest')
        for i,seed in enumerate(seeds):
            for j,p in enumerate(KEYS):
                r=lookup[(seed,p)];value=float(r[field]);star='*' if r[status]=='localized_review' else ''
                ax.text(j,i,'%.3f%s'%(value,star),ha='center',va='center',fontsize=8)
        ax.set_xticks(range(5),KEYS,rotation=40,ha='right')
        ax.set_yticks(range(12),seeds)
        ax.set_title(title,fontsize=10,pad=12)
        ax.set_xticks(np.arange(-.5,5),minor=True);ax.set_yticks(np.arange(-.5,12),minor=True)
        ax.grid(which='minor',color='white',lw=1.5);ax.tick_params(which='minor',bottom=False,left=False)
        for spine in ax.spines.values():spine.set_visible(False)
        ax.set_xlabel('p90 threshold = %.2f index units'%threshold,labelpad=10)
    fig.suptitle('Ploidy ST convergence diagnostics: 12 selected endpoints, five phases, oxygen 0-5%',fontsize=13,y=.98)
    fig.text(.06,.035,'Numbers: p90 across 201 O2 points. Grey: p90 and maximum pass; purple *: p90 passes, local maximum fails; orange: p90 fails.\nPhase range = max - min of five repeats; mean-curve change = absolute difference of five-phase means. These are practical thresholds, not confidence intervals.',fontsize=8,va='bottom')
    fig.subplots_adjust(left=.075,right=.99,top=.87,bottom=.24,wspace=.13)
    for ext in ('pdf','png'):fig.savefig(root/'figures'/('ploidy_ST_convergence_overview.'+ext),dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args();root=args.root.resolve()
    selection,values,variance,oxygen=load(root)
    validation,cached=validate(root,values,oxygen)
    gap_rows=read(root/'spectral_gap_summary.tsv')
    assert len(gap_rows)==72 and len({(r['fit_seed'],r['N'],r['oxygen_band']) for r in gap_rows})==72
    assert sum(int(r['n_evaluations']) for r in gap_rows)==519014160
    for row in gap_rows:
        assert (row['fit_seed'],int(row['N'])) in values and int(row['N']) in NS[1:]
        size={'low_0_to_1pct':41,'middle_1_to_3pct':79,'high_3_to_5pct':81}[row['oxygen_band']]
        assert int(row['n_trajectories'])==70 and int(row['n_evaluations'])==70*int(row['N'])*size
        assert 0<=float(row['fraction_gap_lt_1e_6'])<=float(row['fraction_gap_lt_1e_4'])<=float(row['fraction_gap_lt_1e_3'])<=1
        assert np.isfinite(float(row['minimum_spectral_gap'])) and float(row['minimum_spectral_gap'])>=0
    focus,summary,overall,mechanisms=make_tables(root,selection,values,oxygen,cached)
    plots(root,selection,focus)
    validation.update(python=platform.python_version(),numpy=np.__version__,
        script_sha256=sha(Path(__file__)),new_jobs_submitted=False,
        hpc_collection_manifest_sha256=sha(root/'completed_results_manifest.tsv'),
        spectral_gap_summary_sha256=sha(root/'spectral_gap_summary.tsv'),spectral_gap_rows=len(gap_rows),
        hpc_audit_sha256=sha(root/'slurm/completed_hpc_audit.json'),
        collector_sha256=sha(Path(__file__).with_name('collect_resolution_results.py')),
        n_high_N_variance_cells=sum(x.size for x in variance.values()),
        min_high_N_output_variance=min(float(x.min()) for x in variance.values()))
    (root/'analysis_validation.json').write_text(json.dumps(validation,indent=2,sort_keys=True)+'\n')
    print(json.dumps(validation,indent=2))
    print('SEED_SUMMARY',json.dumps(summary,indent=2))


if __name__=='__main__':main()
