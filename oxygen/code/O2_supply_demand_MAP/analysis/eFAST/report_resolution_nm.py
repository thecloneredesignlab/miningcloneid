#!/usr/bin/env python3
"""Independent reconciliation and interpretation of the four-endpoint N/M pilots.

Only phase arrays, diagnostic tables and audit receipts are read. No model
execution, resampling or job submission. N and M comparisons stay separate.
"""
import argparse
import collections
import csv
import gzip
import json
import platform
from pathlib import Path

import numpy as np
import efast
from report_resolution_pilot import KEYS, WINDOWS, INDICES, Verifier, sha, read, write, stats, band, diagnostic, mechanism
from collect_resolution_nm_results import CONFIGS, PILOT, validate_tasks


def load(master,config):
    n,m,seeds,kind=CONFIGS[config]
    root=master/config
    plan=json.loads((root/'slurm/submission_plan.json').read_text())
    assert (plan['fit_seeds'],plan['resolutions'],plan['M'],plan['comparison_kind'])==(list(seeds),[n],m,kind)
    code=Path(__file__).parent
    for name,digest in plan['code_sha256'].items():assert sha(code/name)==digest,name
    previous=master.parent/'N4097'
    prefix=plan['previous_pilot_dir']+'/'
    for name,digest in plan['comparison_source_sha256'].items():
        assert name.startswith(prefix) and sha(previous/name[len(prefix):])==digest,name
    assert sha(root/'baseline_indices.npz')==plan['comparison_indices_sha256']
    validate_tasks(read(root/'slurm/task_manifest.tsv'),config)
    values,variance={},{}
    with np.load(root/'baseline_indices.npz',allow_pickle=False) as data:
        assert list(data['fit_seeds'])==list(seeds) and list(data['resolutions'])==[2049,4097] and int(data['M'])==4
        assert data['indices'].shape==(2,2,5,2,14,201,2)
        assert list(data['parameters'])==list(efast.ACTIVE) and list(data['outputs'])==list(efast.OUTPUTS)
        assert list(data['index_names'])==list(INDICES)
        oxygen=data['oxygen'].copy()
        for i,seed in enumerate(seeds):
            for j,level in enumerate((2049,4097)):values[(seed,level,4)]=data['indices'][i,j].copy()
    np.testing.assert_allclose(oxygen,np.arange(201)*.025,rtol=0,atol=1e-12)
    with np.load(previous/'baseline_indices.npz',allow_pickle=False) as historical:
        for seed in seeds:
            i=list(historical['fit_seeds']).index(seed);j=list(historical['resolutions']).index(2049)
            np.testing.assert_array_equal(values[(seed,2049,4)],historical['indices'][i,j])
    for seed in seeds:
        for directory,level,harmonics in ((previous,4097,4),(root,n,m)):
            stem=seed+'_N%d'%level+('_M%d'%harmonics if directory==root else '')
            prefix_path=directory/'summaries'/stem
            receipt=json.loads(prefix_path.with_suffix('.receipt.json').read_text())
            assert receipt['status']=='complete' and receipt['n_trajectories']==70 and len(receipt['sources'])==70
            assert sha(prefix_path.with_suffix('.npz'))==receipt['indices_sha256']
            assert sha(prefix_path.with_name(stem+'_repeat_diagnostics.tsv'))==receipt['diagnostics_sha256']
            with np.load(prefix_path.with_suffix('.npz'),allow_pickle=False) as data:
                assert str(data['fit_seed'].item())==seed and int(data['N'])==level
                if directory==root:assert int(data['M'])==harmonics
                assert list(data['parameters'])==list(efast.ACTIVE) and list(data['outputs'])==list(efast.OUTPUTS)
                assert list(data['index_names'])==list(INDICES)
                np.testing.assert_allclose(data['oxygen'],oxygen,rtol=0,atol=1e-12)
                if directory==previous:np.testing.assert_array_equal(values[(seed,4097,4)],data['indices'])
                values[(seed,level,harmonics)]=data['indices'].copy()
                variance[(seed,level,harmonics)]=data['variance'].copy()
    assert len(values)==6 and all(v.shape==(5,2,14,201,2) for v in values.values())
    assert all(v.shape==(5,14,201,2) and np.isfinite(v).all() and (v>=0).all() for v in variance.values())
    return values,variance,oxygen


def validate(root, values, oxygen, config):
    new_n,new_m,seeds,kind=CONFIGS[config]
    v=Verifier()
    cached={key:stats(x) for key,x in values.items()}
    repeats=read(root/'repeat_diagnostics.tsv')
    resolutions=read(root/'configuration_comparison.tsv')
    assert len(repeats)==1008 and len(resolutions)==336
    seen=set()
    for row in repeats:
        seed,n,M=row['fit_seed'],int(row['N']),int(row['M'])
        t,p,k=INDICES.index(row['index']),efast.ACTIVE.index(row['parameter']),efast.OUTPUTS.index(row['output'])
        _,lo,hi=next(w for w in WINDOWS if w[0]==row['oxygen_window'])
        key=(seed,n,M,t,p,k,lo,hi)
        assert key not in seen;seen.add(key)
        mean,sd,span,count=cached[(seed,n,M)]
        q,m,at,valid,passed,status=diagnostic(span[t,p,:,k],oxygen,lo,hi,.10)
        mask=(oxygen>=lo)&(oxygen<=hi)
        curves=values[(seed,n,M)][:,t,p,:,k]
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
        seed,a,b,M=row['fit_seed'],int(row['N_reference']),int(row['N_test']),int(row['M_test'])
        assert (a,b,M)==(4097,new_n,new_m) and int(row['M_reference'])==4
        assert row['comparison_kind']==kind
        t,p,k=INDICES.index(row['index']),efast.ACTIVE.index(row['parameter']),efast.OUTPUTS.index(row['output'])
        _,lo,hi=next(w for w in WINDOWS if w[0]==row['oxygen_window'])
        key=(seed,a,b,M,t,p,k,lo,hi);assert key not in seen;seen.add(key)
        x,_,_,cx=cached[(seed,a,4)];y,_,_,cy=cached[(seed,b,M)]
        diff=np.where((cx[t,p,:,k]==5)&(cy[t,p,:,k]==5),np.abs(y[t,p,:,k]-x[t,p,:,k]),np.nan)
        q,m,at,valid,passed,status=diagnostic(diff,oxygen,lo,hi,.05)
        mask=(oxygen>=lo)&(oxygen<=hi);complete=valid==int(mask.sum())
        expected=dict(absolute_delta_p90=q,absolute_delta_max=m,absolute_delta_max_O2_pct=at,
            n_valid_oxygen=valid,n_expected_oxygen=int(mask.sum()),p90_pass=passed,diagnostic_status=status,
            threshold_absolute_index_units=.05,primary_comparison=True,key_parameter=row['parameter'] in KEYS,
            band_delta=band(y[t,p,:,k],oxygen,lo,hi)-band(x[t,p,:,k],oxygen,lo,hi) if complete else None,
            peak_delta=np.max(y[t,p,mask,k])-np.max(x[t,p,mask,k]) if complete else None,
            peak_O2_low=oxygen[mask][np.argmax(x[t,p,mask,k])] if complete else None,
            peak_O2_high=oxygen[mask][np.argmax(y[t,p,mask,k])] if complete else None)
        for field,value in expected.items():v.check(row[field],value,(key,field))
    curve_count=0
    curve_seen={key:np.zeros((2,14,201,2),dtype=bool) for key in values}
    with gzip.open(root/'oxygen_curve_summary.tsv.gz','rt',newline='') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            key=(row['fit_seed'],int(row['N']),int(row['M']))
            t,p,k=INDICES.index(row['index']),efast.ACTIVE.index(row['parameter']),efast.OUTPUTS.index(row['output'])
            o=float(row['O2_pct']);j=int(round(o/.025));assert abs(oxygen[j]-o)<1e-12
            assert not curve_seen[key][t,p,j,k], 'Duplicate oxygen-curve source row'
            curve_seen[key][t,p,j,k]=True
            mean,sd,span,count=cached[key]
            for field,x in (('phase_mean',mean),('phase_sd',sd),('phase_range',span),('n_defined_phases',count)):
                v.check(row[field],x[t,p,j,k],(key,t,p,j,k,field))
            curve_count+=1
    assert curve_count==67536
    assert all(x.all() for x in curve_seen.values())
    # Verify each endpoint's band scores and submitted deterministic parameter ranks.
    bands_table=read(root/'band_rank_mechanism_comparison.tsv');assert len(bands_table)==336
    lookup={};band_lookup={}
    for row in bands_table:
        key=(row['fit_seed'],int(row['N']),int(row['M']));t,k=INDICES.index(row['index']),efast.OUTPUTS.index(row['output'])
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
    # Focus subsets must be exact subsets, not independently rounded copies.
    for name,parent in [('key_ploidy_ST_repeat_diagnostics.tsv',repeats),('key_ploidy_ST_configuration_comparison.tsv',resolutions)]:
        expected=[r for r in parent if r['parameter'] in KEYS and r['index']=='ST' and r['output']=='dominant_mean_ploidy']
        assert read(root/name)==expected
    return dict(status='passed',independent_numeric_or_categorical_checks=v.n,
        maximum_absolute_numeric_difference=v.max_error,repeat_rows=len(repeats),resolution_rows=len(resolutions),
        curve_rows=curve_count,band_mechanism_rows=len(bands_table),n_source_arrays=len(values),
        undefined_index_cells=sum(int((~np.isfinite(x)).sum()) for x in values.values()),
        min_index=min(float(np.nanmin(x)) for x in values.values()),
        max_index=max(float(np.nanmax(x)) for x in values.values()),
        scope='NPZ-to-table reconciliation; no new model evaluations or FAST resampling'),cached


def metric_fields(prefix,curve,oxygen,lo,hi,threshold):
    q,maximum,at,count,passed,status=diagnostic(curve,oxygen,lo,hi,threshold)
    return {prefix+'90':q,prefix+'max':maximum,prefix+'max_O2_pct':at,
            prefix+'90_pass':passed,prefix+'status':status}


def detail_rows(config,values,oxygen):
    n,m,seeds,kind=CONFIGS[config]
    records=[]
    for seed in seeds:
        a,_,ra,ca=stats(values[(seed,4097,4)])
        b,_,rb,cb=stats(values[(seed,n,m)])
        older,_,_,co=stats(values[(seed,2049,4)])
        delta=np.where((ca==5)&(cb==5),np.abs(b-a),np.nan)
        history=np.where((co==5)&(ca==5),np.abs(a-older),np.nan)
        for t,index in enumerate(INDICES):
            for k,output in enumerate(efast.OUTPUTS):
                for p,parameter in enumerate(efast.ACTIVE):
                    for window,lo,hi in WINDOWS:
                        row=dict(config=config,comparison_kind=kind,fit_seed=seed,N_reference=4097,M_reference=4,
                            N_test=n,M_test=m,index=index,output=output,parameter=parameter,
                            group=efast.GROUP[parameter],key_parameter=parameter in KEYS,oxygen_window=window)
                        row.update(metric_fields('R_reference_',ra[t,p,:,k],oxygen,lo,hi,.10))
                        row.update(metric_fields('R_test_',rb[t,p,:,k],oxygen,lo,hi,.10))
                        row.update(metric_fields('D_',delta[t,p,:,k],oxygen,lo,hi,.05))
                        row.update(metric_fields('D_historical_N2049_to_N4097_',history[t,p,:,k],oxygen,lo,hi,.05))
                        row['test_repeat_and_difference_p90_pass']=row['R_test_90_pass'] and row['D_90_pass']
                        row['test_repeat_and_difference_max_pass']=row['R_test_status']=='passed_by_threshold' and row['D_status']=='passed_by_threshold'
                        row['both_configurations_R90_and_D90_pass']=row['test_repeat_and_difference_p90_pass'] and row['R_reference_90_pass']
                        row['repeat_p90_improved']=row['R_test_90']<row['R_reference_90'] if row['R_test_90'] is not None and row['R_reference_90'] is not None else None
                        row['delta_p90_improved_vs_previous_N_step']=(row['D_90']<row['D_historical_N2049_to_N4097_90']) if kind=='N_resolution' and row['D_90'] is not None else None
                        records.append(row)
    return records


def make_tables(master,config,values,variance,oxygen):
    root=master/config
    n,m,seeds,kind=CONFIGS[config]
    details=detail_rows(config,values,oxygen)
    ploidy=[r for r in details if r['index']=='ST' and r['output']==efast.OUTPUTS[0]]
    focus=[r for r in ploidy if r['key_parameter']]
    write(root/'all_output_parameter_diagnostics.tsv',details)
    write(root/'ploidy_ST_parameter_details.tsv',ploidy)
    write(root/'ploidy_ST_key_parameter_details.tsv',focus)
    summaries=[]
    for seed in seeds:
        for scope in ('key5','all14'):
            items=[r for r in ploidy if r['fit_seed']==seed and r['oxygen_window']==WINDOWS[0][0]
                   and (scope=='all14' or r['key_parameter'])]
            summaries.append(summarize(items,config,seed,scope,WINDOWS[0][0],'ST',efast.OUTPUTS[0]))
    write(root/'seed_convergence_summary.tsv',summaries)
    aggregate=[]
    for scope in ('key5','all14'):
        for window,_,_ in WINDOWS:
            for index in INDICES:
                for output in efast.OUTPUTS:
                    items=[r for r in details if r['oxygen_window']==window and r['index']==index and r['output']==output
                           and (scope=='all14' or r['key_parameter'])]
                    aggregate.append(summarize(items,config,'pooled_selected_two',scope,window,index,output))
    write(root/'convergence_pass_summary.tsv',aggregate)
    mechanisms,phases,profiles,stability=[],[],[],[]
    for seed in seeds:
        for level,harmonics in ((2049,4),(4097,4),(n,m)):
            x=values[(seed,level,harmonics)]
            for t,index in enumerate(INDICES):
                for k,output in enumerate(efast.OUTPUTS):
                    mechanism_result=mechanism(x,oxygen,t,k)
                    low,high=mechanism_result['groups']['low'],mechanism_result['groups']['high']
                    mechanisms.append(dict(config=config,fit_seed=seed,N=level,M=harmonics,index=index,output=output,
                        low_death_largest=mechanism_result['conditions'][0],high_missegregation_increase=mechanism_result['conditions'][1],
                        high_buffering_increase=mechanism_result['conditions'][2],all_three_supported=mechanism_result['supported'],
                        low_dominant_parameter=efast.ACTIVE[int(np.nanargmax(mechanism_result['low']))],
                        high_dominant_parameter=efast.ACTIVE[int(np.nanargmax(mechanism_result['high']))],
                        low_top_group=max(low,key=low.get),high_top_group=max(high,key=high.get),
                        low_death_mean=low['death'],high_death_mean=high['death'],
                        missegregation_high_minus_low=high['missegregation']-low['missegregation'],
                        buffering_high_minus_low=high['buffering']-low['buffering']))
                    for r in range(5):
                        phase=mechanism(np.repeat(x[r:r+1],5,axis=0),oxygen,t,k)
                        phases.append(dict(config=config,fit_seed=seed,N=level,M=harmonics,replicate=r+1,index=index,output=output,
                            low_death_largest=phase['conditions'][0],high_missegregation_increase=phase['conditions'][1],
                            high_buffering_increase=phase['conditions'][2],all_three_supported=phase['supported']))
            for parameter in KEYS:
                for r in range(5):
                    y=x[r,1,efast.ACTIVE.index(parameter),:,0]
                    profiles.append(dict(config=config,fit_seed=seed,N=level,M=harmonics,parameter=parameter,replicate=r+1,
                        whole_band_mean=float(band(y,oxygen,0,5)),low_band_mean=float(band(y,oxygen,0,1)),
                        high_band_mean=float(band(y,oxygen,3,5)),maximum_ST=float(np.nanmax(y)),peak_O2_pct=float(oxygen[np.nanargmax(y)])))
        for t,index in enumerate(INDICES):
            for k,output in enumerate(efast.OUTPUTS):
                a=mechanism(values[(seed,4097,4)],oxygen,t,k);b=mechanism(values[(seed,n,m)],oxygen,t,k)
                for p,parameter in enumerate(efast.ACTIVE):
                    da=a['high'][p]-a['low'][p];db=b['high'][p]-b['low'][p]
                    stability.append(dict(config=config,comparison_kind=kind,fit_seed=seed,index=index,output=output,parameter=parameter,
                        contrast_reference=float(da),contrast_test=float(db),contrast_sign_unchanged=bool(np.sign(da)==np.sign(db)),
                        low_rank_delta=int(b['low_ranks'][p]-a['low_ranks'][p]),high_rank_delta=int(b['high_ranks'][p]-a['high_ranks'][p])))
    write(root/'mechanism_endpoint_summary.tsv',mechanisms)
    write(root/'mechanism_phase_summary.tsv',phases)
    write(root/'mechanism_setting_stability.tsv',stability)
    write(root/'key_parameter_phase_profiles.tsv',profiles)
    audit=json.loads((root/'slurm/completed_hpc_audit.json').read_text())
    assert audit['source_collector_sha256']==sha(Path(__file__).with_name('collect_resolution_nm_results.py'))
    numerical=[]
    for entry in audit['recovery_audits']:
        for row in entry['corrections']['reports']:
            numerical.append(dict(config=config,fit_seed=entry['fit_seed'],N=n,M=m,
                focal_parameter=entry['focal_parameter'],replicate=entry['replicate'],sample_id=row['sample_id'],O2_pct=row['O2_pct'],
                spectral_gap=float(row['original']['spectral_gap']),ploidy_change=row['ploidy_change'],growth_change=row['growth_change'],
                double_residual=row['double_residual'],leading_root_bound_width=row['leading_root_bound_width'],
                ploidy_precision_difference=row['ploidy_precision_difference'],growth_precision_difference=row['growth_precision_difference']))
    if numerical:write(root/'numerical_recovery_summary.tsv',numerical)
    observations=[]
    tasks={r['task_id']:r for r in read(root/'slurm/task_manifest.tsv')}
    trajectory=next(r['job_id'] for r in read(root/'slurm/submission_jobs.tsv') if r['role']=='trajectory')
    for job in audit['slurm_unique_jobs']:
        if not job['job_id'].startswith(trajectory+'_'):continue
        task=tasks[job['job_id'].split('_')[1]]
        observations.append(dict(config=config,fit_seed=task['fit_seed'],N=n,M=m,task_id=int(task['task_id']),
            elapsed_hours=job['elapsed_seconds']/3600,job_id=job['job_id']))
    assert len(observations)==140
    write(root/'runtime_observations.tsv',observations)
    runtime=[]
    baseline_runtime=read(master.parent/'N4097/runtime_by_seed.tsv')
    for seed in seeds:
        hours=[r['elapsed_hours'] for r in observations if r['fit_seed']==seed]
        previous=next(r for r in baseline_runtime if r['fit_seed']==seed and r['N']=='4097')
        runtime.append(dict(config=config,fit_seed=seed,N=n,M=m,n_tasks=70,
            median_hours=float(np.median(hours)),p95_hours=float(np.quantile(hours,.95)),maximum_hours=max(hours),
            baseline_N4097_M4_median_hours=float(previous['median_hours']),
            observed_median_ratio_to_N4097_M4=float(np.median(hours))/float(previous['median_hours']),
            projected_next_N=2*n-1,projected_next_M=m,
            projected_next_median_hours=float(np.median(hours))*(2*n-1)/n,
            projected_next_p95_hours=float(np.quantile(hours,.95))*(2*n-1)/n,
            projected_next_maximum_hours=max(hours)*(2*n-1)/n,
            projection_is_not_measured=True,new_jobs_submitted=False))
    write(root/'runtime_by_seed.tsv',runtime)
    variance_rows=[]
    for (seed,level,harmonics),v in variance.items():
        for k,output in enumerate(efast.OUTPUTS):
            variance_rows.append(dict(config=config,fit_seed=seed,N=level,M=harmonics,output=output,
                minimum_variance=float(v[:,:,:,k].min()),maximum_variance=float(v[:,:,:,k].max()),
                nonfinite_cells=int((~np.isfinite(v[:,:,:,k])).sum()),zero_variance_cells=int((v[:,:,:,k]==0).sum())))
    write(root/'output_variance_summary.tsv',variance_rows)
    return details,summaries,aggregate


def summarize(items,config,seed,scope,window,index,output):
    return dict(config=config,comparison_kind=CONFIGS[config][3],fit_seed=seed,parameter_scope=scope,
        oxygen_window=window,index=index,output=output,n_parameter_seed_cells=len(items),
        reference_repeat_p90_pass=sum(r['R_reference_90_pass'] for r in items),
        test_repeat_p90_pass=sum(r['R_test_90_pass'] for r in items),
        difference_p90_pass=sum(r['D_90_pass'] for r in items),
        test_repeat_and_difference_p90_pass=sum(r['test_repeat_and_difference_p90_pass'] for r in items),
        test_repeat_and_difference_max_pass=sum(r['test_repeat_and_difference_max_pass'] for r in items),
        both_configurations_R90_and_D90_pass=sum(r['both_configurations_R90_and_D90_pass'] for r in items),
        repeat_p90_improved=sum(r['repeat_p90_improved'] is True for r in items),
        reference_repeat_max_pass=sum(r['R_reference_status']=='passed_by_threshold' for r in items),
        test_repeat_max_pass=sum(r['R_test_status']=='passed_by_threshold' for r in items),
        difference_max_pass=sum(r['D_status']=='passed_by_threshold' for r in items))


def plot_overview(master,config,details):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    n,m,seeds,kind=CONFIGS[config]
    records={(r['fit_seed'],r['parameter']):r for r in details if r['index']=='ST'
             and r['output']==efast.OUTPUTS[0] and r['oxygen_window']==WINDOWS[0][0]}
    fig,axes=plt.subplots(1,3,figsize=(14.8,4.6),sharey=True)
    cmap=ListedColormap(['#eeeeee','#d5e3ed','#f3bd8e','#777777'])
    states={'passed_by_threshold':0,'localized_review':1,'not_passed':2,'undefined_or_incomplete':3}
    titles=('A  N4097/M4 phase range','B  N%d/M%d phase range'%(n,m),
            'C  Absolute mean difference')
    for ax,prefix,title in zip(axes,('R_reference_','R_test_','D_'),titles):
        cells=np.array([[states[records[(s,p)][prefix+'status']] for p in KEYS] for s in seeds])
        ax.imshow(cells,cmap=cmap,vmin=0,vmax=3,aspect='auto')
        for i,seed in enumerate(seeds):
            for j,p in enumerate(KEYS):
                row=records[(seed,p)];value=row[prefix+'90']
                ax.text(j,i,'NA' if value is None else '%.3f%s'%(value,'*' if cells[i,j]==1 else ''),ha='center',va='center',fontsize=10)
        ax.set(xticks=range(5),yticks=range(2),xticklabels=KEYS,yticklabels=seeds,title=title)
        plt.setp(ax.get_xticklabels(),rotation=45,ha='right')
    label='N resolution at fixed M=4' if kind=='N_resolution' else 'M setting robustness at fixed N=4097'
    fig.suptitle('Ploidy ST numerical diagnostics: '+label+'\nFive key parameters; five phases; 201 O2 points (0-5%)',fontsize=12)
    fig.text(.05,.02,'Cells show p90 of the diagnostic, not ST. Grey: p90 and maximum pass; blue *: p90 passes, local maximum fails; orange: p90 fails.\nPhase-range threshold 0.10; mean-difference threshold 0.05. Full phase ranges are not confidence intervals.',fontsize=9)
    fig.tight_layout(rect=(0,.13,1,.88))
    for ext in ('pdf','png'):fig.savefig(master/config/'figures'/('ploidy_ST_diagnostic_overview.'+ext),dpi=170)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args();master=args.root.resolve()
    validations=[];all_details=[];all_summaries=[]
    for config in CONFIGS:
        values,variance,oxygen=load(master,config)
        validation,_=validate(master/config,values,oxygen,config)
        details,summaries,aggregate=make_tables(master,config,values,variance,oxygen)
        plot_overview(master,config,details)
        gaps=read(master/config/'spectral_gap_summary.tsv')
        assert len(gaps)==6 and sum(int(r['n_evaluations']) for r in gaps)==140*CONFIGS[config][0]*201
        validation.update(config=config,n_variance_cells=sum(x.size for x in variance.values()),
            minimum_variance=min(float(v.min()) for v in variance.values()),
            python=platform.python_version(),numpy=np.__version__,script_sha256=sha(Path(__file__)),
            independent_helpers_sha256=sha(Path(__file__).with_name('report_resolution_pilot.py')),
            collector_sha256=sha(Path(__file__).with_name('collect_resolution_nm_results.py')),
            hpc_audit_sha256=sha(master/config/'slurm/completed_hpc_audit.json'),new_jobs_submitted=False)
        (master/config/'analysis_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
        validations.append(validation);all_details.extend(details);all_summaries.extend(summaries)
        print('VALIDATION '+json.dumps(validation),flush=True)
        print('SEED_SUMMARY '+json.dumps(summaries),flush=True)
    write(master/'ploidy_ST_parameter_details.tsv',[r for r in all_details if r['index']=='ST' and r['output']==efast.OUTPUTS[0]])
    write(master/'ploidy_ST_key_parameter_details.tsv',[r for r in all_details if r['index']=='ST' and r['output']==efast.OUTPUTS[0] and r['key_parameter']])
    write(master/'seed_convergence_summary.tsv',all_summaries)
    (master/'analysis_validation.json').write_text(json.dumps(dict(status='passed',configurations=validations,
        independent_checks=sum(v['independent_numeric_or_categorical_checks'] for v in validations),
        max_numeric_difference=max(v['maximum_absolute_numeric_difference'] for v in validations),new_jobs_submitted=False),indent=2)+'\n')


if __name__=='__main__':main()
