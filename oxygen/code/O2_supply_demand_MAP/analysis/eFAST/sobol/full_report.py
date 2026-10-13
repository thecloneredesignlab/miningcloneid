"""Seed-resolved full Sobol summaries, diagnostic-only ensemble reports and plots."""
import csv
from collections import defaultdict
from pathlib import Path
import warnings
import hashlib
import json

import numpy as np
import core as c
import full_core as f


def summary_rows(case,oxygen,mean,spread,delta,status,variance,width):
    for oi,o2 in enumerate(oxygen):
        for k,output in enumerate(c.OUTPUTS):
            for m,(index,pi,pj) in enumerate(f.METRICS):
                yield dict(scope=case['scope'],case_id=case['case_id'],fit_seed=case['fit_seed'],N=f.N,
                    O2_pct=o2,output=output,index=index,parameter_i=pi,parameter_j=pj,
                    group_i=c.efast.GROUP[pi],group_j=c.efast.GROUP[pj] if pj else '',
                    value_mean=mean[oi,k,m],repeat_range=spread[oi,k,m],resolution_delta=delta[oi,k,m],
                    resolution_low=2048,resolution_high=4096,bootstrap_half_width_median=width[oi,k,m],
                    output_variance_mean=variance[oi,k],
                    convergence_status=['undefined','passed','not_passed'][status[oi,k,m]],
                    included_in_full_analysis='TRUE')


def summarize_case(root,p,case):
    target=f.case_dir(root,case); target.mkdir(parents=True,exist_ok=True)
    key=hashlib.sha256(json.dumps(dict(case=case,unit_sha256=p['unit_sha256'],engine=p['engine_sha256'],
        git_commit=p['git_commit'],oxygen=c.OXYGEN,resolutions=f.NS),sort_keys=True).encode()).hexdigest()
    if (target/'completion.json').exists():
        old=c.load_json(target/'completion.json')
        if old['analysis_identity']!=key or old['summary_sha256']!=c.sha(target/'indices_summary.npz'): raise ValueError('Summary changed')
        return
    bounds=c.read(target/'parameter_bounds.tsv'); shape=(2,5,len(c.OXYGEN),2,len(f.METRICS))
    checkpoint=target/'indices_checkpoint.npz'
    if checkpoint.exists():
        with np.load(checkpoint) as z:
            if str(z['identity'])!=key: raise ValueError('Analysis checkpoint identity changed')
            values,widths,variance,done=z['values'],z['widths'],z['variance'],z['done']
        if values.shape!=shape: raise ValueError('Analysis checkpoint axes differ')
    else:
        values=np.full(shape,np.nan); widths=np.full(shape,np.nan)
        variance=np.full(shape[:-1],np.nan); done=np.zeros((5,len(c.OXYGEN),2),bool)
    receipts=[]; numerical=[]
    for repeat in range(1,6):
        arrays=[]
        for bi,block in enumerate(c.BLOCKS):
            path=target/'runs'/('R%02d'%repeat)/'blocks'/block/'outputs.npz'
            x=f.samples(root,case,repeat,bi); key=f.identity(p,case,repeat,bi,x)
            if not f.complete(path,key): raise ValueError('Missing/invalid block: '+str(path))
            rr=c.load_json(str(path)+'.receipt.json'); receipts.append(dict(path=str(path.relative_to(root)),sha256=rr['sha256']))
            numerical.append(dict(repeat=repeat,block=block,n_points=rr['n_points'],n_recovered=rr['n_recovered'],
                minimum_spectral_gap=rr['minimum_spectral_gap'],maximum_residual=rr['maximum_residual']))
            if not done[repeat-1].all():
                with np.load(path) as z:
                    if not np.array_equal(z['oxygen'],c.OXYGEN) or z['outputs'].shape!=(f.N,len(c.OXYGEN),2):
                        raise ValueError('Output axes differ')
                    arrays.append(z['outputs'])
        if done[repeat-1].all(): continue
        ordered=c.assemble(arrays); del arrays
        for oi,o2 in enumerate(c.OXYGEN):
            for k,output in enumerate(c.OUTPUTS):
                if done[repeat-1,oi,k]: continue
                for ni,n in enumerate(f.NS):
                    si,var=c.analyze(bounds,ordered[:n*c.STEP,oi,k],seed=2026101300+repeat+10000+n+oi*10+k)
                    vv,ww=f.metric_arrays(si); values[ni,repeat-1,oi,k]=vv; widths[ni,repeat-1,oi,k]=ww
                    variance[ni,repeat-1,oi,k]=var
                done[repeat-1,oi,k]=True
            if (oi+1)%8==0 or oi==len(c.OXYGEN)-1:
                f.atomic_npz(checkpoint,identity=np.array(key),values=values,widths=widths,variance=variance,done=done)
                print('Indices',case['case_id'],'repeat',repeat,'oxygen',oi+1,'/201',flush=True)
        del ordered
    if not done.all(): raise ValueError('Incomplete five-repeat indices')
    f.atomic_npz(target/'indices_repeats.npz',values=values,bootstrap_half_width=widths,variance=variance,
        resolutions=np.array(f.NS),oxygen=np.asarray(c.OXYGEN),metrics=np.array(f.METRICS),outputs=np.array(c.OUTPUTS))
    mean,spread,delta,status=f.diagnostics(values)
    outvar=np.mean(variance[-1],axis=0); width=np.median(widths[-1],axis=0)
    f.atomic_npz(target/'indices_summary.npz',mean=mean,repeat_range=spread,resolution_delta=delta,
                 status=status,variance=outvar,bootstrap_half_width_median=width,oxygen=np.asarray(c.OXYGEN))
    rows=summary_rows(case,c.OXYGEN,mean,spread,delta,status,outvar,width)
    first=next(rows); fields=list(first)
    from itertools import chain
    c.efast.write_table(target/'indices_summary.tsv.gz',fields,chain([first],rows))
    c.table(target/'numerical_qc.tsv',numerical)
    plot_indices(mean,c.OXYGEN,target/'figures',case['scope']+'/'+case['case_id']+'; mean of five repeats')
    c.save_json(target/'completion.json',dict(status='complete',analysis_identity=key,case=case,N=f.N,repeats=5,oxygen_points=201,
        summary_sha256=c.sha(target/'indices_summary.npz'),summary_table_sha256=c.sha(target/'indices_summary.tsv.gz'),
        repeat_indices_sha256=c.sha(target/'indices_repeats.npz'),output_receipts=receipts,
        n_undefined_indices=int((status==0).sum()),n_not_passed_indices=int((status==2).sum()),
        n_recovered=sum(r['n_recovered'] for r in numerical),completed_utc=c.now(),
        included_in_full_analysis=True,convergence_policy='diagnostic_only'))
    checkpoint.unlink()


def plot_indices(values,oxygen,directory,title,fraction=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    for index in ('S1','ST','S2'):
        selected=[i for i,m in enumerate(f.METRICS) if m[0]==index]
        labels=[f.METRICS[i][1]+(' × '+f.METRICS[i][2] if index=='S2' else ' ['+c.efast.GROUP[f.METRICS[i][1]]+']') for i in selected]
        height=23 if index=='S2' else 6
        fig,axes=plt.subplots(1,2,figsize=(16,height),sharey=True,layout='constrained')
        for k,ax in enumerate(axes):
            data=values[:,k,selected].T; finite=data[np.isfinite(data)]
            lo=min(0,float(finite.min())) if len(finite) else 0.; hi=max(.001,float(finite.max())) if len(finite) else .001
            high='#e66101' if index=='ST' else '#54278f'
            cmap=LinearSegmentedColormap.from_list('sensitivity',['#2166ac','white',high] if lo<0 else ['white',high]); cmap.set_bad('#d9d9d9')
            if lo<0: image=ax.imshow(np.ma.masked_invalid(data),aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=lo,vcenter=0,vmax=hi))
            else: image=ax.imshow(np.ma.masked_invalid(data),aspect='auto',cmap=cmap,vmin=0,vmax=hi)
            ticks=np.arange(0,len(oxygen),20)
            ax.set_xticks(ticks,[format(oxygen[i],'.2g') for i in ticks]); ax.set_xlabel('Oxygen (%)')
            ax.set_yticks(range(len(labels)),labels,fontsize=7)
            ax.set_title('Dominant mean ploidy' if k==0 else 'Net live-cell growth rate')
            fig.colorbar(image,ax=ax,orientation='horizontal',pad=.025,fraction=.025,label=index+' (raw estimates; blue = negative sampling estimates)')
        fig.suptitle(title+'; N4096\nAll endpoints retained; gray = undefined; diagnostics in source tables',fontsize=11)
        name='s2_pair_by_oxygen' if index=='S2' else index.lower()+'_parameter_by_oxygen'
        fig.savefig(directory/(name+'.pdf')); fig.savefig(directory/(name+'.png'),dpi=120); plt.close(fig)
        if index=='S2' and fraction is not None:
            fig,axes=plt.subplots(1,2,figsize=(16,23),sharey=True,layout='constrained')
            for k,ax in enumerate(axes):
                image=ax.imshow(fraction[:,k,selected].T,aspect='auto',cmap='viridis',vmin=0,vmax=1)
                ax.set_xticks(np.arange(0,len(oxygen),20),[format(oxygen[i],'.2g') for i in range(0,len(oxygen),20)])
                ax.set_yticks(range(len(labels)),labels,fontsize=7); ax.set_xlabel('Oxygen (%)')
                ax.set_title(c.OUTPUTS[k]); fig.colorbar(image,ax=ax,orientation='horizontal',label='Fraction of all 500 endpoints passing both diagnostics')
            fig.suptitle('500-neighborhood diagnostic coverage; no endpoint exclusion')
            fig.savefig(directory/'s2_diagnostic_fraction.pdf'); fig.savefig(directory/'s2_diagnostic_fraction.png',dpi=120); plt.close(fig)


def aggregate_cases(root,p):
    cube=[]; statuses=[]; variances=[]; numerical=[]
    for case in p['cases']:
        target=f.case_dir(root,case); rr=c.load_json(target/'completion.json')
        if rr['status']!='complete' or rr['summary_sha256']!=c.sha(target/'indices_summary.npz') or rr['repeat_indices_sha256']!=c.sha(target/'indices_repeats.npz'):
            raise ValueError('Missing/changed case summary: '+case['case_id'])
        with np.load(target/'indices_summary.npz') as z:
            cube.append(z['mean']); statuses.append(z['status']); variances.append(z['variance'])
        numerical.append(dict(case_id=case['case_id'],scope=case['scope'],n_recovered=rr['n_recovered'],
                              n_undefined=rr['n_undefined_indices'],n_not_passed=rr['n_not_passed_indices'],included=True))
    cube=np.stack(cube); statuses=np.stack(statuses); variances=np.stack(variances)
    c.table(root/'diagnostics/case_diagnostics.tsv',numerical)
    scope_values={}; scope_records={}
    for scope in ['global','neighborhood10pct']:
        xx=cube[:1] if scope=='global' else cube[1:]; ss=statuses[:1] if scope=='global' else statuses[1:]
        count=np.isfinite(xx).sum(axis=0); n=len(xx)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            median=np.nanmedian(xx,axis=0); q10=np.nanquantile(xx,.1,axis=0); q90=np.nanquantile(xx,.9,axis=0)
        passed=(ss==1).sum(axis=0); fraction=passed/n
        scope_values[scope]=median; target=root/scope/'summaries'; target.mkdir(parents=True,exist_ok=True)
        f.atomic_npz(target/'scope_summary.npz',median=median,endpoint_q10=q10,endpoint_q90=q90,n_defined=count,
                     n_passed=passed,pass_fraction_all_endpoints=fraction,oxygen=np.asarray(c.OXYGEN))
        records=[]
        for oi,o in enumerate(c.OXYGEN):
            for k,output in enumerate(c.OUTPUTS):
                for m,(index,pi,pj) in enumerate(f.METRICS):
                    records.append(dict(scope=scope,N=f.N,O2_pct=o,output=output,index=index,parameter_i=pi,parameter_j=pj,
                        value_mean=median[oi,k,m],endpoint_q10=q10[oi,k,m],endpoint_q90=q90[oi,k,m],
                        n_endpoints=n,n_defined=int(count[oi,k,m]),n_passed=int(passed[oi,k,m]),
                        pass_fraction_all_endpoints=fraction[oi,k,m],
                        convergence_status='passed' if passed[oi,k,m]==n else 'mixed_or_not_passed',
                        summary_population='all endpoints; median over defined values; undefined count explicit',included='TRUE'))
        c.table(target/'indices_scope_summary.tsv.gz',records); scope_records[scope]=records
        plot_indices(median,c.OXYGEN,root/scope/'figures',
            'Global full-bound sensitivity' if scope=='global' else 'Median of 500 endpoint five-repeat means',
            None if scope=='global' else fraction)
    mechanism_report(root,cube[1:],statuses[1:],p)
    prepare_responses(root,p,scope_records)
    c.save_json(root/'status.json',dict(status='all_indices_complete_responses_prepared',N=f.N,n_neighborhoods=500,
        repeats=5,oxygen_points=201,n_pairs=91,convergence_policy='diagnostic_only',updated_utc=c.now()))


def mechanism_report(root,cube,statuses,p):
    low=c.efast.normalized_trapezoid_weights(np.asarray(c.OXYGEN),0,1)
    high=c.efast.normalized_trapezoid_weights(np.asarray(c.OXYGEN),3,5)
    groups=sorted(set(c.efast.GROUP.values())); rows=[]; comparisons=[]
    for index,offset in [('S1',0),('ST',14)]:
        values=cube[:,:,:,offset:offset+14]
        bands={}
        for band,w in [('low_0_to_1pct',low),('high_3_to_5pct',high)]:
            selected=w>0
            scores=np.einsum('sokp,o->skp',values[:,selected],w[selected])
            bands[band]=np.stack([np.mean(scores[:,:,[i for i,param in enumerate(c.ACTIVE) if c.efast.GROUP[param]==g]],axis=2) for g in groups],axis=2)
            for si in range(500):
                for k,output in enumerate(c.OUTPUTS):
                    for gi,g in enumerate(groups):
                        rows.append(dict(fit_seed='seed%d'%(si+1),index=index,output=output,oxygen_band=band,group=g,
                            group_parameter_mean_index=bands[band][si,k,gi],
                            interpretation='group mean of parameter indices; not a group variance decomposition'))
        a,b=bands['low_0_to_1pct'],bands['high_3_to_5pct']
        for k,output in enumerate(c.OUTPUTS):
            valid=np.isfinite(a[:,k]).all(axis=1)&np.isfinite(b[:,k]).all(axis=1)
            death=groups.index('death'); mis=groups.index('missegregation'); buf=groups.index('buffering')
            criterion=(a[:,k,death]>np.max(np.delete(a[:,k],death,axis=1),axis=1))&(b[:,k,mis]>a[:,k,mis])&(b[:,k,buf]>a[:,k,buf])&valid
            comparisons.append(dict(index=index,output=output,n_total=500,n_evaluable=int(valid.sum()),
                n_hypothesis_pattern=int(criterion.sum()),fraction_of_all_500=float(criterion.sum()/500),
                convergence_exclusion='none',interpretation='optimizer-endpoint robustness; not biological confidence or formal group Sobol indices'))
    c.table(root/'neighborhood10pct/summaries/mechanism_band_seed_summary.tsv.gz',rows)
    c.table(root/'neighborhood10pct/summaries/hypothesis_pattern_fraction.tsv',comparisons)


def prepare_responses(root,p,scope_records):
    candidates=[]; selected=[]; used=set()
    for scope,records in scope_records.items():
        for output in c.OUTPUTS:
            rr=[r for r in records if r['index']=='S2' and r['output']==output and np.isfinite(r['value_mean'])]
            rr.sort(key=lambda r:r['value_mean'],reverse=True); candidates.extend(rr)
            for r in rr:
                key=(r['parameter_i'],r['parameter_j'])
                if key not in used: selected.append(r); used.add(key); break
    for r in sorted(candidates,key=lambda r:r['value_mean'],reverse=True):
        key=(r['parameter_i'],r['parameter_j'])
        if key not in used: selected.append(r); used.add(key)
        if len(selected)==6: break
    if len(selected)!=6: raise ValueError('Cannot select six response pairs')
    keys=[dict(pair_id=i+1,parameter_i=r['parameter_i'],parameter_j=r['parameter_j'],
               selection_origin=r['scope']+'/'+r['output'],peak_oxygen=r['O2_pct'],peak_estimate=r['value_mean'],
               peak_pass_fraction=r['pass_fraction_all_endpoints'],interpretation='candidate; peak stability shown explicitly') for i,r in enumerate(selected)]
    c.table(root/'key_pairs.tsv',keys)
    from report import anchored_grid
    representative=[p['cases'][0]]+[case for case in p['cases'][1:] if case['fit_seed'] in c.PILOT_SEEDS]
    tasks=[]
    for case in representative:
        target=f.case_dir(root,case)
        bounds={r['parameter']:r for r in c.read(target/'parameter_bounds.tsv')}
        best={r['parameter']:float(r['value']) for r in c.read(target/'reference_best_params.tsv')}
        for pair in keys:
            pi,pj=pair['parameter_i'],pair['parameter_j']; bi,bj=bounds[pi],bounds[pj]
            gi=anchored_grid(bi['natural_lower'],bi['natural_upper'],best[pi],bi['transform']=='log10')
            gj=anchored_grid(bj['natural_lower'],bj['natural_upper'],best[pj],bj['transform']=='log10')
            part=target/'responses'/(pi+'__'+pj); part.mkdir(parents=True,exist_ok=True)
            samples=[dict(sample_id=i*len(gj)+j+1,**dict(best,**{pi:float(vi),pj:float(vj)})) for i,vi in enumerate(gi) for j,vj in enumerate(gj)]
            c.table(part/'samples.tsv.gz',samples,fields=['sample_id']+list(c.ACTIVE))
            c.save_json(part/'metadata.json',dict(method='Sobol_response_surface',scope=case['scope'],case_id=case['case_id'],fit_seed=case['fit_seed'],
                parameter_i=pi,parameter_j=pj,grid_i=gi.tolist(),grid_j=gj.tolist(),anchor_i=best[pi],anchor_j=best[pj],
                transform_i=bi['transform'],transform_j=bj['transform'],n_samples=len(samples),parameters=list(c.ACTIVE),
                oxygen_pct=[0.,.5,2.5,5.],samples_sha256=c.sha(part/'samples.tsv.gz'),
                best_params_sha256=case['best_params_sha256'],fit_config_sha256=case['fit_config_sha256'],validate_reference=False))
            tasks.append(dict(task_id=len(tasks)+1,part=str(part.relative_to(root)),scope=case['scope'],case_id=case['case_id'],fit_seed=case['fit_seed'],
                              parameter_i=pi,parameter_j=pj,metadata_sha256=c.sha(part/'metadata.json')))
    c.table(root/'slurm/response_tasks.tsv',tasks)


def finish_report(root,p):
    from report import plot_surface
    c.OXYGEN=[0.,.5,2.5,5.]
    tasks=c.read(root/'slurm/response_tasks.tsv'); surface_rows=[]
    for task in tasks:
        part=root/task['part']; meta=c.load_json(part/'metadata.json')
        if c.sha(part/'metadata.json')!=task['metadata_sha256'] or not c.complete(part/'outputs.tsv.gz'): raise ValueError('Response incomplete/changed')
        values,qc=c.output_array(part/'outputs.tsv.gz',meta['n_samples'],c.OXYGEN)
        directory=part.parents[1]/'figures/responses'
        data,contrasts=plot_surface(meta,values,directory/(meta['parameter_i']+'__'+meta['parameter_j']))
        rows=[]
        for i,vi in enumerate(meta['grid_i']):
            for j,vj in enumerate(meta['grid_j']):
                for oi,o2 in enumerate(c.OXYGEN):
                    rows.append(dict(value_i=vi,value_j=vj,O2_pct=o2,dominant_mean_ploidy=data[i,j,oi,0],
                        dominant_growth_rate=data[i,j,oi,1],ploidy_nonadditive_contrast=contrasts[oi][i,j],growth_nonadditive_contrast=contrasts[4+oi][i,j]))
        c.table(part/'response_surface_source.tsv.gz',rows); surface_rows.append(dict(**task,**qc))
    c.table(root/'response_surface_inventory.tsv',surface_rows)
    (root/'analysis_notes.md').write_text('''# Complete independent in-vivo second-order Sobol analysis

One full-bound domain and all 500 clipped natural-span 10% neighborhoods are retained.
N4096, five independently scrambled repeats and the complete Figure 4 201-point grid.
Main indices use N4096; N2048 nested prefixes provide resolution diagnostics without extra model evaluations.
Repeat range <=0.05 and resolution delta <=0.02 are operational diagnostics, not proof of convergence.
No neighborhood is excluded for failing a diagnostic. Undefined indices remain undefined; endpoint medians
use defined values and report their counts explicitly. Negative estimates are sampling error, not antagonism.
S2 excludes main effects and interactions of order three or higher. S1/ST complement mechanistic interpretation;
ST cannot be summed into a process variance decomposition. Oxygen-band group summaries are group means
of parameter indices, not formal group Sobol indices. Endpoint fractions describe optimizer-solution robustness,
not biological confidence. Conditional response surfaces fix other parameters at the indicated optimum;
their nonadditive contrasts are not marginalized S2. Separate output/scope color scales must be read independently.

Source tables and diagnostic coverage: global/summaries and neighborhood10pct/summaries.
Individual 500-endpoint sources: neighborhood10pct/full/seed*/indices_summary.tsv.gz.
Main figures: global/figures and neighborhood10pct/figures. Per-endpoint figures remain separately available.
Raw float64 outputs, full repeated index arrays, checkpoints and exception matrices remain on HPC.
''')
    c.save_json(root/'status.json',dict(status='complete',method='Sobol_second_order',N=f.N,n_repeats=5,
        oxygen_points=201,n_global_designs=1,n_neighborhoods=500,n_pairs=91,n_response_surfaces=54,
        all_500_neighborhoods_run=True,convergence_policy='diagnostic_only',completed_utc=c.now()))
    def inventory():
        for path in sorted(root.rglob('*')):
            relative=path.relative_to(root)
            if not path.is_file() or 'runs' in relative.parts or 'recovery' in relative.parts or 'logs' in relative.parts or 'code_snapshot' in relative.parts or 'runtime_validation' in relative.parts: continue
            if path.name in ['collection_manifest.tsv','samples.tsv.gz','outputs.tsv.gz'] or path.suffix not in ['.json','.tsv','.gz','.md','.pdf','.png']: continue
            if path.name.startswith('outputs.tsv.gz') or path.name.endswith('.tar.gz'): continue
            yield dict(path=str(relative),bytes=path.stat().st_size,sha256=c.sha(path))
    c.efast.write_table(root/'collection_manifest.tsv',['path','bytes','sha256'],inventory())
