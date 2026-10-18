#!/usr/bin/env python3
"""Read existing outputs only: phase, response and Fourier diagnostics.

No model invocation, no new samples, no cluster submission. Formal indices
use unmodified periodic response sequences and the original design's N/M.
Auxiliary cutoff sweeps are not sensitivity estimates for a different M.
"""
import argparse
import csv
import gzip
import itertools
import json
import platform
from pathlib import Path

import numpy as np
from collect_phase_spectrum import FOCUS, KEYS
from collect_resolution_nm_results import PARAMETERS, sha, read

OXYGEN=np.arange(201)*.025
WINDOWS={'seed25':(0,5),'seed138':(1.5,2.8),'seed101':(3,5),'seed174':(0,1)}
CASE_O2=np.array([0,.2,.425,.55,.75,1,1.575,1.625,1.7,1.775,1.925,1.975,2.1,2.25,2.475,2.575,2.8,3.35,4.025,4.35,4.75,4.95,5])
CASE_J=np.rint(CASE_O2/.025).astype(int)
OUTPUTS=('ploidy','growth')
MULTIPLIERS=(.5,.75,.8,.9,1.,1.1,1.2,1.25,1.5,2.)
METRICS=('mean','variance','low_power','high_power','S1','ST','near_cutoff_power_fraction',
         'nyquist_tail_fraction','top1pct_variance_fraction','q01','q10','q50','q90','q99','minimum','maximum',
         'periodic_jump_q99','min_gap','gap_below_1e6_fraction','gap_below_1e4_fraction','top1pct_lowgap_fraction',
         'spectral_entropy','k90_over_omega','cutoff_sweep_span_20pct')


def table(path,rows):
    rows=list(rows)
    assert rows,path
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


def frequencies(n,m,p,d=14):
    assert n%2==1 and n>4*m*m and 0<=p<d
    omega=(n-1)//(2*m);budget=omega//(2*m)
    other=np.floor(np.linspace(1,budget,d-1)).astype(int) if budget>=d-1 else np.arange(d-1)%budget+1
    result=np.zeros(d,dtype=int);result[p]=omega;result[np.arange(d)!=p]=other
    return result


def infer_phase(u,n,m,p):
    """Infer an existing phase from stored coordinates; never draw new points."""
    omega=frequencies(n,m,p)
    a=np.pi*(u[0,p]-.5)
    candidates=(a%(2*np.pi),(np.pi-a)%(2*np.pi))
    s=2*np.pi*np.arange(n)/n
    errors=[]
    for phi in candidates:
        restored=.5+np.arcsin(np.sin(s[:,None]*omega+phi))/np.pi
        errors.append(float(np.max(np.abs(restored-u))))
    best=int(np.argmin(errors))
    assert errors[best]<2e-7,errors
    return candidates[best],errors[best],omega


def spectrum(y,m):
    n=y.shape[0]
    assert n%2==1 and n>4*m*m and np.isfinite(y).all()
    omega=(n-1)//(2*m);cut=omega//2
    p=2*np.abs(np.fft.rfft(y,axis=0)[1:]/n)**2
    v=p.sum(axis=0)
    if np.any(v<=1e-20*np.maximum(1,np.mean(y,axis=0)**2)):
        raise ValueError('Undefined index: zero or negligible variance')
    low=p[:cut].sum(axis=0)
    first=p[np.arange(1,m+1)*omega-1].sum(axis=0)/v
    return p,v,low,first,1-low/v,omega,cut


def transformed_inputs(samples,ranges):
    u=np.empty_like(samples)
    for j,row in enumerate(ranges):
        x=np.log10(samples[:,j]) if row['transform']=='log10' else samples[:,j]
        u[:,j]=(x-float(row['encoded_lower']))/(float(row['encoded_upper'])-float(row['encoded_lower']))
    assert np.isfinite(u).all() and u.min()>-1e-8 and u.max()<1+1e-8
    return u


def process(root):
    receipt=json.loads((root/'raw_collection_receipt.json').read_text());assert receipt['status']=='passed'
    assert receipt['raw_manifest_sha256']==sha(root/'raw_file_manifest.tsv')
    cache=root/'processed_cache';cache.mkdir(exist_ok=True)
    records=read(root/'trajectory_manifest.tsv');metrics=[];design=[];validation=[];peaks=[];sweeps=[]
    signature=sha(Path(__file__))
    for number,row in enumerate(records,1):
        identity={k:row[k] for k in ('trajectory_id','fit_seed','N','M','replicate','parameter')}
        n,m,p=int(row['N']),int(row['M']),int(row['parameter_index'])
        path=root/'raw_cache'/row['part_path'];target=cache/('trajectory_%03d.npz'%int(row['trajectory_id']))
        # The full collection is SHA checked; each parse also checks its raw content.
        assert sha(path/'outputs.tsv.gz')==row['output_sha256']
        print('%d/%d %s N%d M%d R%s %s'%(number,len(records),row['fit_seed'],n,m,row['replicate'],row['parameter']),flush=True)
        with np.load(path/'indices.npz',allow_pickle=False) as z:
            expected=z['indices'].copy();expected_v=z['variance'].copy()
            np.testing.assert_allclose(z['oxygen'],OXYGEN,rtol=0,atol=1e-12)
        # Canonical order is validated before reshape: sample_id x O2, never sort by output.
        with gzip.open(path/'outputs.tsv.gz','rt') as f:
            data=np.loadtxt(f,skiprows=1,delimiter='\t',usecols=(0,1,2,3,4))
        assert data.shape==(n*201,5) and np.isfinite(data).all()
        data=data.reshape(n,201,5)
        np.testing.assert_array_equal(data[:,:,0],np.broadcast_to(np.arange(1,n+1)[:,None],(n,201)))
        np.testing.assert_allclose(data[:,:,1],np.broadcast_to(OXYGEN,(n,201)),rtol=0,atol=1e-12)
        y=data[:,:,2:4];gap=data[:,:,4];assert (gap>=0).all()
        with gzip.open(path/'samples.tsv.gz','rt') as f:samples=np.loadtxt(f,skiprows=1,delimiter='\t')
        assert samples.shape==(n,15)
        np.testing.assert_array_equal(samples[:,0],np.arange(1,n+1))
        u=transformed_inputs(samples[:,1:],read(root/'raw_cache'/row['ranges_path']))
        phi,error,omega_all=infer_phase(u,n,m,p)
        corr=np.corrcoef(u,rowvar=False);np.fill_diagonal(corr,0)
        other=np.arange(14)!=p;other_corr=corr[np.ix_(other,other)]
        standardized=(u-u.mean(axis=0))/u.std(axis=0)
        triples=[(abs(float(np.mean(standardized[:,a]*standardized[:,b]*standardized[:,c]))),a,b,c)
                 for a,b,c in itertools.combinations(np.where(other)[0],3)]
        triple=max(triples)
        design.append(dict(identity,phase_radians=phi,coordinate_reconstruction_max_error=error,
            omega_main=int(omega_all[p]),formal_cutoff=int(omega_all[p]//2),
            complementary_frequencies=','.join(map(str,omega_all[other])),
            max_abs_pair_correlation=float(np.abs(corr).max()),max_abs_nonfocal_pair_correlation=float(np.abs(other_corr).max()),
            max_abs_nonfocal_third_moment=triple[0],third_moment_parameters=','.join(PARAMETERS[j] for j in triple[1:]),
            normalized_marginal_mean_max_error=float(np.max(np.abs(u.mean(axis=0)-.5))),
            normalized_marginal_variance_max_error=float(np.max(np.abs(u.var(axis=0)-1/12)))))
        power,v,low,s1,st,omega,cut=spectrum(y,m)
        parseval=float(np.max(np.abs(v-np.var(y,axis=0))/v))
        cache_index=float(np.max(np.abs(np.stack((s1,st))-expected)))
        variance_error=float(np.max(np.abs(v-expected_v)/expected_v))
        assert parseval<1e-11 and variance_error<1e-11 and cache_index<1e-10
        validation.append(dict(identity,n_rows=n*201,index_max_absolute_error=cache_index,
            variance_max_relative_error=variance_error,parseval_max_relative_error=parseval,
            coordinate_reconstruction_max_error=error,output_sha256=row['output_sha256'],
            formal_index_cells=804,source_flags_scope='Original certified analyzer checked every status/nonnegative flag; complete raw SHA matched'))
        normalized=power/v;cum=np.cumsum(normalized,axis=0)
        near=power[int(np.floor(.8*cut)):int(np.ceil(1.2*cut))].sum(axis=0)/v
        tail=power[int(.9*len(power)):].sum(axis=0)/v
        sq=(y-y.mean(axis=0))**2;ntop=int(np.ceil(.01*n))
        top=np.argpartition(sq,-ntop,axis=0)[-ntop:]
        leverage=np.take_along_axis(sq,top,axis=0).sum(axis=0)/sq.sum(axis=0)
        lowgap=np.mean(np.take_along_axis(np.broadcast_to(gap[:,:,None],y.shape),top,axis=0)<1e-6,axis=0)
        q=np.quantile(y,[.01,.1,.5,.9,.99],axis=0)
        entropy=-np.sum(normalized*np.log(np.maximum(normalized,1e-300)),axis=0)/np.log(len(power))
        k90=(np.argmax(cum>=.9,axis=0)+1)/omega
        aux=np.array([1-power[:int(np.floor(cut*factor))].sum(axis=0)/v for factor in MULTIPLIERS])
        spread20=aux[MULTIPLIERS.index(.8)]-aux[MULTIPLIERS.index(1.2)]
        arrays=[y.mean(axis=0),v,low,v-low,s1,st,near,tail,leverage,*q,y.min(axis=0),y.max(axis=0),
                np.quantile(np.abs(np.roll(y,-1,axis=0)-y),.99,axis=0),
                np.broadcast_to(gap.min(axis=0)[:,None],v.shape),
                np.broadcast_to((gap<1e-6).mean(axis=0)[:,None],v.shape),
                np.broadcast_to((gap<1e-4).mean(axis=0)[:,None],v.shape),lowgap,entropy,k90,spread20]
        assert len(arrays)==len(METRICS)
        matrix=np.stack(arrays,axis=-1)
        for j,o in enumerate(OXYGEN):
            for k,out in enumerate(OUTPUTS):
                metrics.append(dict(identity,O2_pct=float(o),output=out,**{name:float(matrix[j,k,t]) for t,name in enumerate(METRICS)}))
        for ci,j in enumerate(CASE_J):
            for k,out in enumerate(OUTPUTS):
                for rank,frequency in enumerate(np.argsort(normalized[:,j,k])[-8:][::-1],1):
                    peaks.append(dict(identity,O2_pct=float(OXYGEN[j]),output=out,rank=rank,frequency=int(frequency+1),
                        frequency_over_omega=float((frequency+1)/omega),variance_fraction=float(normalized[frequency,j,k])))
                for a,factor in enumerate(MULTIPLIERS):
                    sweeps.append(dict(identity,O2_pct=float(OXYGEN[j]),output=out,cutoff_multiplier=factor,
                        diagnostic_high_frequency_fraction=float(aux[a,j,k]),formal_ST=bool(factor==1),
                        interpretation='auxiliary partition only; not a new M design'))
        np.savez_compressed(target,code_sha256=signature,metrics=matrix,case_oxygen=OXYGEN[CASE_J],
            case_y=y[:,CASE_J],case_gap=gap[:,CASE_J],case_power=normalized[:,CASE_J],
            u=u,phase=phi,omega=omega,cutoff=cut)
        del data,y,gap,power,normalized,sq,top,cum
    table(root/'phase_metrics.tsv.gz',metrics);table(root/'input_design_diagnostics.tsv',design)
    table(root/'trajectory_validation.tsv',validation);table(root/'top_spectral_peaks.tsv.gz',peaks)
    table(root/'auxiliary_cutoff_sweep.tsv.gz',sweeps)
    record=dict(status='passed',n_trajectories=len(records),n_raw_response_rows=sum(int(r['N'])*201 for r in records),
        n_formal_index_cells=804*len(records),n_parameter_coordinates=sum(int(r['N'])*14 for r in records),
        max_cached_index_absolute_error=max(r['index_max_absolute_error'] for r in validation),
        max_variance_relative_error=max(r['variance_max_relative_error'] for r in validation),
        max_parseval_relative_error=max(r['parseval_max_relative_error'] for r in validation),
        max_design_coordinate_absolute_error=max(r['coordinate_reconstruction_max_error'] for r in validation),
        code_sha256=signature,python=platform.python_version(),numpy=np.__version__,new_model_evaluations=0,new_jobs_submitted=0,
        original_certified_status_flags_rechecked=False,
        original_status_flags_scope='Original analyzer certification plus complete raw SHA256; local numeric parse uses first five columns')
    (root/'numerical_validation.json').write_text(json.dumps(record,indent=2)+'\n')
    summarize(root)


def summarize(root):
    with gzip.open(root/'phase_metrics.tsv.gz','rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    groups={}
    for r in rows:
        key=tuple(r[k] for k in ('fit_seed','N','M','parameter','output','O2_pct'))
        groups.setdefault(key,[]).append(r)
    aggregate=[]
    for key,rs in groups.items():
        assert len(rs)==5
        identity=dict(zip(('fit_seed','N','M','parameter','output','O2_pct'),key))
        data={name:np.array([float(r[name]) for r in rs]) for name in METRICS}
        result=dict(identity,ST_mean=float(data['ST'].mean()),ST_range=float(np.ptp(data['ST'])),S1_range=float(np.ptp(data['S1'])),
            variance_max_over_min=float(data['variance'].max()/data['variance'].min()),
            low_power_max_over_min=float(data['low_power'].max()/data['low_power'].min()),
            high_power_max_over_min=float(data['high_power'].max()/data['high_power'].min()))
        for name in ('near_cutoff_power_fraction','top1pct_variance_fraction','min_gap','gap_below_1e6_fraction','cutoff_sweep_span_20pct'):
            result[name+'_mean']=float(data[name].mean());result[name+'_max']=float(data[name].max())
        result['min_gap_min']=float(data['min_gap'].min())
        # Descriptive correlations across five phases, not a significance test.
        for name in ('variance','low_power','high_power','q99','top1pct_variance_fraction'):
            a=data[name];b=data['ST']
            result['five_phase_corr_ST_'+name]=float(np.corrcoef(a,b)[0,1]) if np.std(a)>0 and np.std(b)>0 else float('nan')
        aggregate.append(result)
    table(root/'phase_spread_by_oxygen.tsv.gz',aggregate)
    cases=[];bands=[]
    for seed,parameters in FOCUS.items():
        configs=sorted({(int(r['N']),int(r['M'])) for r in aggregate if r['fit_seed']==seed})
        for parameter in parameters:
            for n,m in configs:
                for output in OUTPUTS:
                    candidates=[r for r in aggregate if (r['fit_seed'],r['N'],r['M'],r['parameter'],r['output'])==(seed,str(n),str(m),parameter,output)]
                    for band,lo,hi in [('focus',*WINDOWS[seed]),('low',0,1),('high',3,5),('outside_focus',0,5)]:
                        if band=='outside_focus':
                            chosen=[r for r in candidates if float(r['O2_pct'])<WINDOWS[seed][0]-1e-10 or float(r['O2_pct'])>WINDOWS[seed][1]+1e-10]
                        else:chosen=[r for r in candidates if lo-1e-10<=float(r['O2_pct'])<=hi+1e-10]
                        if not chosen:continue
                        maximum=max(chosen,key=lambda r:r['ST_range'])
                        bands.append(dict(fit_seed=seed,N=n,M=m,parameter=parameter,output=output,oxygen_window=band,n_oxygen=len(chosen),
                            ST_range_p90=float(np.quantile([r['ST_range'] for r in chosen],.9)),ST_range_max=maximum['ST_range'],
                            ST_range_max_O2_pct=float(maximum['O2_pct']),S1_range_p90=float(np.quantile([r['S1_range'] for r in chosen],.9)),
                            variance_ratio_p90=float(np.quantile([r['variance_max_over_min'] for r in chosen],.9)),
                            cutoff_near_fraction_p90=float(np.quantile([r['near_cutoff_power_fraction_mean'] for r in chosen],.9)),
                            cutoff_sweep20_max= max(r['cutoff_sweep_span_20pct_max'] for r in chosen),
                            tail1pct_variance_fraction_max=max(r['top1pct_variance_fraction_max'] for r in chosen),
                            min_gap=min(float(r['min_gap_min']) for r in chosen)))
            # One cached representative O2 per parameter: maximum test-setting spread in focus among disclosed candidate grid.
            n,m=configs[-1]
            chosen=[r for r in aggregate if (r['fit_seed'],r['N'],r['M'],r['parameter'],r['output'])==(seed,str(n),str(m),parameter,'ploidy')
                    and WINDOWS[seed][0]-1e-9<=float(r['O2_pct'])<=WINDOWS[seed][1]+1e-9 and np.min(np.abs(CASE_O2-float(r['O2_pct'])))<1e-9]
            winner=max(chosen,key=lambda r:r['ST_range'])
            cases.append(dict(fit_seed=seed,parameter=parameter,O2_pct=float(winner['O2_pct']),N_test=n,M_test=m,
                              test_ST_range=winner['ST_range'],selection='maximum test phase range within focus among fixed cached candidate O2 values'))
    table(root/'representative_cases.tsv',cases);table(root/'band_diagnostics.tsv',bands)


def figures(root):
    """Chart contract: x is O2, normalized input coordinate, or k/omega, never mixed.
    Five fixed phase colours; reference dashed, test solid. Spectra are raw
    unwindowed power/variance. Binned response means are display summaries
    conditional on the existing joint path, not one-at-a-time perturbations.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    palette=['#0072B2','#D55E00','#E69F00','#009E73','#CC79A7']
    manifest=read(root/'trajectory_manifest.tsv');cases=read(root/'representative_cases.tsv')
    by_id={int(r['trajectory_id']):r for r in manifest};source_spectra=[];source_response=[];source_hist=[];anomalies=[]
    figdir=root/'figures';figdir.mkdir(exist_ok=True)
    for seed in FOCUS:
        selected=[r for r in cases if r['fit_seed']==seed]
        fig,axes=plt.subplots(len(selected),4,figsize=(17,2.7*len(selected)),squeeze=False)
        for i,case in enumerate(selected):
            parameter=case['parameter'];o=float(case['O2_pct']);ci=int(np.argmin(np.abs(CASE_O2-o)))
            trajectories=[r for r in manifest if r['fit_seed']==seed and r['parameter']==parameter]
            configs=sorted({(int(r['N']),int(r['M'])) for r in trajectories})
            ax0,ax1,ax2,ax3=axes[i];bin_edges=np.linspace(0,1,25)
            for row in trajectories:
                n,m,rep=map(int,(row['N'],row['M'],row['replicate']));test=(n,m)==configs[-1]
                color=palette[rep-1];style='-' if test else '--';alpha=.9 if test else .45
                with np.load(root/'processed_cache'/('trajectory_%03d.npz'%int(row['trajectory_id'])),allow_pickle=False) as z:
                    y=z['case_y'][:,ci,0];gap=z['case_gap'][:,ci];u=z['u'][:,int(row['parameter_index'])]
                    power=z['case_power'][:,ci,0];omega=int(z['omega']);cut=int(z['cutoff']);met=z['metrics'][:,0,:]
                    assert abs(float(z['case_oxygen'][ci])-o)<1e-10
                    digit=np.clip(np.digitize(u,bin_edges)-1,0,23)
                    means=np.array([np.mean(y[digit==j]) for j in range(24)])
                    ax0.plot((bin_edges[1:]+bin_edges[:-1])/2,means,style,color=color,alpha=alpha,lw=1)
                    # Exact CDF rather than a smoothed spectrum makes cutoff partition visible.
                    k=np.arange(1,len(power)+1)/omega
                    ax1.plot(k,np.cumsum(power),style,color=color,alpha=alpha,lw=1)
                    ax2.plot(OXYGEN,met[:,METRICS.index('ST')],style,color=color,alpha=alpha,lw=.85)
                    # Same power, independent growth response serves as control.
                    ax3.plot(k,power,style,color=color,alpha=alpha,lw=.65)
                    identity={key:row[key] for key in ('trajectory_id','fit_seed','N','M','replicate','parameter')}
                    for j,val in enumerate(power):source_spectra.append(dict(identity,O2_pct=o,frequency=j+1,frequency_over_omega=float(k[j]),ploidy_variance_fraction=float(val)))
                    for j,mean in enumerate(means):source_response.append(dict(identity,O2_pct=o,input_bin_center=float((bin_edges[j]+bin_edges[j+1])/2),n_in_bin=int((digit==j).sum()),conditional_ploidy_mean=float(mean)))
                    hist,edges=np.histogram(y,bins=48,density=False)
                    for j,count in enumerate(hist):source_hist.append(dict(identity,O2_pct=o,lower=float(edges[j]),upper=float(edges[j+1]),count=int(count),fraction=float(count/n)))
                    top=np.argsort((y-y.mean())**2)[-3:][::-1]
                    # Actual input vector and spectral gap for the most influential response points.
                    for rank,j in enumerate(top,1):
                        rec=dict(identity,O2_pct=o,rank=rank,sample_id=int(j+1),ploidy=float(y[j]),spectral_gap=float(gap[j]),
                                 squared_deviation_fraction=float((y[j]-y.mean())**2/np.sum((y-y.mean())**2)))
                        rec.update({'normalized_'+name:float(z['u'][j,t]) for t,name in enumerate(PARAMETERS)});anomalies.append(rec)
            ax0.set(ylabel=parameter+'\nPloidy (joint-path bin mean)',xlabel='Focal input quantile coordinate (0–1)')
            ax0.set_title('%s | O₂ %.3g%%'%(parameter,o),loc='left',fontsize=10)
            ax1.axvline(.5,color='0.3',ls=':',lw=1);ax1.set(xlabel='Frequency / focal frequency (k/ω)',ylabel='Cumulative variance fraction',xlim=(0,2),ylim=(0,1.02))
            ax2.axvspan(*WINDOWS[seed],color='0.85',alpha=.4);ax2.axvline(o,color='0.3',ls=':',lw=.7)
            ax2.set(xlabel='O₂ (%)',ylabel='Formal ploidy ST',xlim=(0,5),ylim=(0,1))
            ax3.set_yscale('log');ax3.set(xlabel='Frequency / focal frequency (k/ω)',ylabel='Variance fraction / Fourier bin',xlim=(0,2),ylim=(1e-8,1))
            ax3.axvline(.5,color='0.3',ls=':',lw=1);ax3.axvline(1,color='0.3',ls=':',lw=.6)
        fig.suptitle(seed+' | Existing phase responses and raw Fourier power',fontsize=15,y=.995)
        handles=[Line2D([0],[0],color=c,label='Phase %d'%(j+1)) for j,c in enumerate(palette)]
        handles.extend([Line2D([0],[0],color='0.2',ls='--',label='Reference N%d/M%d'%configs[0]),Line2D([0],[0],color='0.2',ls='-',label='Test N%d/M%d'%configs[-1])])
        fig.legend(handles=handles,loc='lower center',ncol=7,bbox_to_anchor=(.5,.005),frameon=False)
        fig.tight_layout(rect=(0,.035,1,.975));fig.savefig(figdir/(seed+'_phase_spectrum.pdf'));fig.savefig(figdir/(seed+'_phase_spectrum.png'),dpi=150);plt.close(fig)
    table(root/'representative_spectra.tsv.gz',source_spectra);table(root/'representative_response_bins.tsv',source_response)
    table(root/'representative_histograms.tsv.gz',source_hist);table(root/'influential_response_points.tsv',anomalies)
    overview(root,plt,palette)


def overview(root,plt,palette):
    manifest=read(root/'trajectory_manifest.tsv');cases=read(root/'representative_cases.tsv')
    focus={'seed25':'mu_hp','seed138':'p_mis_base','seed101':'p_mis_base','seed174':'O2_crit'}
    fig,axes=plt.subplots(4,3,figsize=(13,12))
    for i,(seed,parameter) in enumerate(focus.items()):
        case=next(r for r in cases if r['fit_seed']==seed and r['parameter']==parameter);o=float(case['O2_pct']);ci=int(np.argmin(np.abs(CASE_O2-o)))
        trajectories=[r for r in manifest if r['fit_seed']==seed and r['parameter']==parameter and int(r['N'])==int(case['N_test']) and int(r['M'])==int(case['M_test'])]
        for row in trajectories:
            rep=int(row['replicate']);color=palette[rep-1]
            with np.load(root/'processed_cache'/('trajectory_%03d.npz'%int(row['trajectory_id']))) as z:
                y=z['case_y'][:,ci,0];metrics=z['metrics'];j=int(round(o/.025));omega=int(z['omega'])
                hist,edges=np.histogram(y,bins=48,density=True)
                axes[i,0].stairs(hist,edges,color=color,lw=.9,label='Phase %d'%rep)
                v=metrics[j,0,METRICS.index('variance')];l=metrics[j,0,METRICS.index('low_power')];h=v-l
                axes[i,1].bar(rep,l,color=color,alpha=.35,width=.65)
                axes[i,1].bar(rep,h,bottom=l,color=color,width=.65)
                axes[i,1].text(rep,v,'%.2f'%metrics[j,0,METRICS.index('ST')],ha='center',va='bottom',fontsize=8)
                axes[i,2].plot(OXYGEN,metrics[:,0,METRICS.index('ST')],color=color,lw=1)
                axes[i,2].plot(OXYGEN,metrics[:,1,METRICS.index('ST')],color=color,ls=':',lw=.9)
        axes[i,0].set(title='%s: %s at %.3g%% O₂'%(seed,parameter,o),xlabel='Ploidy',ylabel='Response density')
        axes[i,1].set(title='V = L + H; labels are ST = H/V',xlabel='Existing phase repeat',ylabel='Ploidy variance',xticks=range(1,6));axes[i,1].margins(y=.25)
        axes[i,2].set(title='Solid: ploidy; dotted: growth',xlabel='O₂ (%)',ylabel='Formal ST',xlim=(0,5),ylim=(0,1));axes[i,2].axvspan(*WINDOWS[seed],color='0.85',alpha=.4)
    fig.suptitle('Phase-dependent output variance and ST: selected diagnostic examples',fontsize=14)
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=5,frameon=False)
    fig.tight_layout(rect=(0,.035,1,.97));fig.savefig(root/'figures/diagnostic_overview.pdf');fig.savefig(root/'figures/diagnostic_overview.png',dpi=160);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('process','summarize','figures'));p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();globals()[a.command](a.root.resolve())


if __name__=='__main__':main()
