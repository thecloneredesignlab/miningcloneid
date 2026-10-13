"""Per-case diagnostics, separate scope reports, and conditional response surfaces."""
from collections import defaultdict
import json
from pathlib import Path
import time

import numpy as np
import core as c


def summarize_rows(rows, plan):
    resolutions=plan.get('resolutions',c.RESOLUTIONS)
    lower,highest_n=resolutions[-2:]
    groups=defaultdict(list)
    fields=('scope','case_id','fit_seed','N','O2_pct','output','index','parameter_i','parameter_j')
    for row in rows: groups[tuple(row[k] for k in fields)].append(row)
    means={key:float(np.mean([r['value'] for r in rr])) for key,rr in groups.items()}
    result=[]
    for key, rr in groups.items():
        if len(rr)!=5 or len({r['repeat'] for r in rr})!=5:
            raise ValueError('Require five distinct scrambles for every index')
        values=np.array([r['value'] for r in rr]); widths=np.array([r['bootstrap_half_width'] for r in rr])
        valid=bool(np.isfinite(values).all()); repeat_range=float(np.ptp(values)) if valid else np.nan
        k=list(key); k[3]=lower; mid=means.get(tuple(k),np.nan)
        k[3]=highest_n; highest=means.get(tuple(k),np.nan)
        delta=abs(highest-mid)
        status='undefined' if not valid or not np.isfinite(delta) else (
            'passed' if repeat_range<=plan['repeat_range_tolerance'] and delta<=plan['resolution_delta_tolerance'] else 'not_passed')
        result.append(dict(zip(fields,key),value_mean=means[key],value_median=float(np.median(values)),
            repeat_min=float(np.min(values)),repeat_max=float(np.max(values)),repeat_range=repeat_range,
            n_valid_repeats=int(np.isfinite(values).sum()),bootstrap_half_width_median=float(np.median(widths)),
            n_bootstrap_intervals_above_zero=sum(r['bootstrap_lower']>0 for r in rr),
            resolution_low=lower,resolution_high=highest_n,resolution_delta=delta,
            output_variance_mean=float(np.mean([r['output_variance'] for r in rr])),
            convergence_status=status,retained='TRUE',**{'delta_N%d_N%d'%(lower,highest_n):delta}))
    return result


def scope_summary(case_summaries, scope):
    selected=[r for rr in case_summaries for r in rr if r['scope']==scope and int(r['N'])==max(c.RESOLUTIONS)]
    if scope=='global': return selected
    groups=defaultdict(list)
    for r in selected: groups[(r['O2_pct'],r['output'],r['index'],r['parameter_i'],r['parameter_j'])].append(r)
    result=[]
    for key,rr in groups.items():
        if len(rr)!=8: raise ValueError('Require all eight representative neighborhoods')
        values=np.array([float(r['value_mean']) for r in rr]); finite=np.isfinite(values)
        # Undefined endpoints are retained; an incomplete eight-endpoint median is undefined.
        median=float(np.median(values)) if finite.all() else np.nan
        result.append(dict(scope=scope,case_id='representative8_median',fit_seed='representative8',N=max(c.RESOLUTIONS),
            O2_pct=key[0],output=key[1],index=key[2],parameter_i=key[3],parameter_j=key[4],
            value_mean=median,value_median=median,n_endpoints=8,n_defined_endpoints=int(finite.sum()),
            endpoint_q25=float(np.quantile(values,.25)) if finite.all() else np.nan,
            endpoint_q75=float(np.quantile(values,.75)) if finite.all() else np.nan,
            n_convergence_passed=sum(r['convergence_status']=='passed' for r in rr),
            convergence_status='passed' if all(r['convergence_status']=='passed' for r in rr) else 'mixed_or_not_passed',
            output_variance_mean=float(np.median([float(r['output_variance_mean']) for r in rr])),retained='TRUE'))
    return result


def plot_pair_heatmap(rows, path, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    lookup={(r['output'],r['parameter_i'],r['parameter_j'],float(r['O2_pct'])):r for r in rows if r['index']=='S2' and int(r['N'])==max(c.RESOLUTIONS)}
    fig,axes=plt.subplots(1,2,figsize=(10,23),sharey=True,layout='constrained')
    for ax,output in zip(axes,c.OUTPUTS):
        data=np.array([[float(lookup[(output,c.ACTIVE[i],c.ACTIVE[j],o)]['value_mean']) for o in c.OXYGEN] for i,j in c.PAIRS])
        finite=data[np.isfinite(data)]; vmin=min(0.,float(finite.min())) if len(finite) else 0.
        vmax=max(.001,float(finite.max())) if len(finite) else .001
        if vmin<0:
            cmap=LinearSegmentedColormap.from_list('signed_s2',['#2166ac','white','#54278f'])
            cmap.set_bad('#d9d9d9')
            im=ax.imshow(np.ma.masked_invalid(data),aspect='auto',cmap=cmap,norm=TwoSlopeNorm(vmin=vmin,vcenter=0,vmax=vmax))
        else:
            cmap=LinearSegmentedColormap.from_list('positive_s2',['white','#54278f'])
            cmap.set_bad('#d9d9d9')
            im=ax.imshow(np.ma.masked_invalid(data),aspect='auto',cmap=cmap,vmin=0,vmax=vmax)
        ax.set_xticks(range(len(c.OXYGEN)),[str(o) for o in c.OXYGEN]); ax.set_xlabel('Oxygen (%)')
        ax.set_title('Dominant mean ploidy' if output==c.OUTPUTS[0] else 'Net live-cell growth rate')
        ax.set_yticks(range(len(c.PAIRS)),[c.ACTIVE[i]+' × '+c.ACTIVE[j] for i,j in c.PAIRS],fontsize=7)
        unstable=[]
        for p,(i,j) in enumerate(c.PAIRS):
            for oi,o in enumerate(c.OXYGEN):
                if lookup[(output,c.ACTIVE[i],c.ACTIVE[j],o)]['convergence_status']!='passed': unstable.append((oi,p))
        if unstable:
            xx,yy=zip(*unstable); ax.scatter(xx,yy,s=6,c='#555555',marker='x',linewidths=.5)
        bar=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.025,pad=.025,label='S2 (raw estimates; negative values retained)')
        if vmin<0: bar.set_ticks([vmin,0,vmax/2,vmax])
    fig.suptitle(title+'; N='+str(max(c.RESOLUTIONS))+'\n× = diagnostic not passed/undefined; gray cells = undefined',fontsize=11)
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path.with_suffix('.pdf')); fig.savefig(path.with_suffix('.png'),dpi=120); plt.close(fig)


def choose_pairs(scope_rows):
    groups=[]
    for scope,rows in scope_rows.items():
        for output in c.OUTPUTS:
            candidates=[]
            for i,j in c.PAIRS:
                rr=[r for r in rows if r['index']=='S2' and r['output']==output and r['parameter_i']==c.ACTIVE[i] and r['parameter_j']==c.ACTIVE[j]]
                values=[float(r['value_mean']) for r in rr if np.isfinite(float(r['value_mean']))]
                strength=max(values) if values else -np.inf
                stable=any(r['convergence_status']=='passed' and float(r['value_mean'])>0 for r in rr)
                candidates.append(dict(parameter_i=c.ACTIVE[i],parameter_j=c.ACTIVE[j],strength=strength,
                    stable=stable,origin=scope+'/'+output))
            candidates.sort(key=lambda r:(r['stable'],r['strength']),reverse=True)
            groups.append(candidates)
    selected=[]; used=set()
    # Cover both outputs and both scopes before filling by estimated strength.
    for group in groups:
        for r in group:
            key=(r['parameter_i'],r['parameter_j'])
            if key not in used:
                selected.append(r); used.add(key); break
    for r in sorted([r for g in groups for r in g],key=lambda r:(r['stable'],r['strength']),reverse=True):
        key=(r['parameter_i'],r['parameter_j'])
        if key not in used:
            selected.append(r); used.add(key)
        if len(selected)==6: break
    return [dict(pair_id=i+1,parameter_i=r['parameter_i'],parameter_j=r['parameter_j'],
                 selection_origin=r['origin'],selection_strength=r['strength'],
                 stable_somewhere=str(r['stable']).upper(),
                 interpretation='candidate response surface; stability is reported separately') for i,r in enumerate(selected)]


def anchored_grid(lo,hi,best,logarithmic,size=41):
    lo,hi,best=map(float,(lo,hi,best))
    tolerance=64*np.finfo(float).eps*max(1.,abs(lo),abs(hi),abs(best))
    if best<lo-tolerance or best>hi+tolerance: raise ValueError('Response anchor outside bounds')
    # Preserve the actual optimum when decimal/log bound conversions differ only
    # by roundoff. This applies to the deterministic surface grid, not Sobol bounds.
    lo,hi=min(lo,best),max(hi,best)
    values=np.geomspace(lo,hi,size) if logarithmic else np.linspace(lo,hi,size)
    close=np.flatnonzero(np.isclose(values,best,rtol=1e-12,atol=0))
    index=int(close[0]) if len(close) else int(np.argmin(np.abs(values[1:-1]-best)))+1
    values[index]=best
    return np.sort(values)


def interaction_contrast(z,ai,aj):
    return z-z[:,aj][:,None]-z[ai,:][None,:]+z[ai,aj]


def finalize_indices(root,plan):
    all_summaries=[]; common=[]
    for case in plan['cases']:
        target=c.case_directory(root,case)
        receipt=c.load_json(target/'completion.json')
        if receipt['status']!='complete' or receipt['summary_sha256']!=c.sha(target/'indices_summary.tsv') or receipt['indices_sha256']!=c.sha(target/'indices_repeats.tsv.gz'):
            raise ValueError('Incomplete/changed case summary')
        all_summaries.append(c.read(target/'indices_summary.tsv'))
        common.extend(receipt['common_reference'])
    reference={}
    for r in common:
        key=float(r['O2_pct']); value=np.array([r[o] for o in c.OUTPUTS])
        if key in reference and np.max(np.abs(value-reference[key]))>1e-8:
            raise ValueError('Common-vector fixed-O2 output differs among fitted configurations')
        reference[key]=value
    c.save_json(root/'inputs/common_configuration_check.json',dict(status='passed',
        n_points=len(common),tolerance=1e-8,oxygen_pct=c.OXYGEN,
        description='Identical midpoint parameter vector produces matching outputs in all eight endpoint configurations and the global reference'))
    scope_rows={}
    for scope in ('global','neighborhood10pct'):
        rows=scope_summary(all_summaries,scope); scope_rows[scope]=rows
        target=root/scope/'pilot'
        c.table(target/'indices_scope_summary.tsv',rows)
        title='Global full-bound sensitivity: mean of five scrambles' if scope=='global' else 'Neighborhood sensitivity: median of eight endpoint five-scramble means'
        plot_pair_heatmap(rows,target/'figures/s2_pair_by_oxygen',title)
    for case,rows in zip(plan['cases'],all_summaries):
        if case['scope']=='neighborhood10pct':
            plot_pair_heatmap(rows,c.case_directory(root,case)/'figures/s2_pair_by_oxygen',case['case_id']+': neighborhood, mean of five scrambles')
    comparison=[]
    fields=('scope','case_id','fit_seed','O2_pct','output','index','parameter_i','parameter_j')
    groups=defaultdict(dict)
    for rows in all_summaries:
        for row in rows: groups[tuple(row[f] for f in fields)][int(row['N'])]=row
    for key,by_n in groups.items():
        if set(by_n)!=set(c.RESOLUTIONS): raise ValueError('Incomplete resolution comparison')
        record=dict(zip(fields,key))
        for n in c.RESOLUTIONS:
            record.update({'value_N%d'%n:by_n[n]['value_mean'],
                           'repeat_range_N%d'%n:by_n[n]['repeat_range'],
                           'diagnostic_N%d'%n:by_n[n]['convergence_status']})
        record.update(resolution_low=c.RESOLUTIONS[-2],resolution_high=c.RESOLUTIONS[-1],
                      resolution_delta=by_n[c.RESOLUTIONS[-1]]['resolution_delta'])
        comparison.append(record)
    c.table(root/'resolution_comparison.tsv',comparison)
    if plan.get('reuse_source'):
        # Same N1024 vectors and outputs must reproduce the old indices exactly.
        source=Path(plan['reuse_source']); checks=[]
        audit=c.load_json(root/'inputs/reuse_source_audit.json')
        for name,expected in audit['summary_sha256'].items():
            if c.sha(name)!=expected: raise ValueError('Reuse summary changed after submission: '+name)
        for case,rows in zip(plan['cases'],all_summaries):
            old=c.read(c.case_directory(source,case)/'indices_summary.tsv')
            lookup={tuple(r[f] for f in fields):r for r in old if int(r['N'])==plan['reuse_source_N']}
            maximum=0.; count=0
            for row in rows:
                if int(row['N'])!=plan['reuse_source_N']: continue
                prior=lookup[tuple(row[f] for f in fields)]
                a,b=float(row['value_mean']),float(prior['value_mean'])
                if np.isnan(a) and np.isnan(b): count+=1; continue
                if not np.isfinite([a,b]).all() or a!=b:
                    raise ValueError('Reused N indices differ: '+case['case_id'])
                maximum=max(maximum,abs(a-b)); count+=1
            if count!=len(lookup): raise ValueError('Missing source-N index comparison')
            checks.append(dict(**case,n_indices=count,maximum_absolute_difference=maximum,status='passed'))
        c.save_json(root/'reused_indices_validation.json',dict(status='passed',N=plan['reuse_source_N'],cases=checks))
    keys=choose_pairs(scope_rows)
    c.table(root/'key_pairs.tsv',keys)
    tasks=[]
    source_manifest={r['fit_seed']:r for r in c.read(root/'inputs/representative_seeds.tsv')}
    for case in plan['cases']:
        case_root=c.case_directory(root,case)
        bounds={r['parameter']:r for r in c.read(case_root/'parameter_bounds.tsv')}
        best={r['parameter']:float(r['value']) for r in c.read(case_root/'reference_best_params.tsv')}
        endpoint=source_manifest[case['fit_seed']]
        for pair in keys:
            pi,pj=pair['parameter_i'],pair['parameter_j']; bi,bj=bounds[pi],bounds[pj]
            gi=anchored_grid(bi['natural_lower'],bi['natural_upper'],best[pi],bi['transform']=='log10')
            gj=anchored_grid(bj['natural_lower'],bj['natural_upper'],best[pj],bj['transform']=='log10')
            part=case_root/'responses'/(pi+'__'+pj); part.mkdir(parents=True,exist_ok=True)
            samples=[]
            for vi in gi:
                for vj in gj:
                    values=dict(best); values[pi]=float(vi); values[pj]=float(vj)
                    samples.append(dict(sample_id=len(samples)+1,**values))
            c.table(part/'samples.tsv.gz',samples,fields=['sample_id']+list(c.ACTIVE))
            c.save_json(part/'metadata.json',dict(method='Sobol_response_surface',**case,
                parameter_i=pi,parameter_j=pj,grid_i=gi.tolist(),grid_j=gj.tolist(),
                anchor_i=best[pi],anchor_j=best[pj],transform_i=bi['transform'],transform_j=bj['transform'],
                original_bounds_i=[float(bi['natural_lower']),float(bi['natural_upper'])],
                original_bounds_j=[float(bj['natural_lower']),float(bj['natural_upper'])],
                anchor_boundary_policy='exact optimum retained within 64 machine epsilons; Sobol sampling bounds unchanged',
                n_samples=len(samples),parameters=list(c.ACTIVE),
                oxygen_pct=c.OXYGEN,samples_sha256=c.sha(part/'samples.tsv.gz'),
                best_params_sha256=endpoint['best_params_sha256'],fit_config_sha256=endpoint['fit_config_sha256'],
                validate_reference=False,other_parameters='held at this fitted optimum; conditional surface, not marginalized S2'))
            tasks.append(dict(task_id=len(tasks)+1,**case,parameter_i=pi,parameter_j=pj,
                              part=str(part.relative_to(root)),metadata_sha256=c.sha(part/'metadata.json')))
    c.table(root/'slurm/response_tasks.tsv',tasks)
    c.save_json(root/'slurm/response_plan.json',dict(n_tasks=len(tasks),key_pairs_sha256=c.sha(root/'key_pairs.tsv'),
        tasks_sha256=c.sha(root/'slurm/response_tasks.tsv'),nominal_grid_size=41,anchor_inserted_exactly=True,
        selection='strongest stable candidates, cover both outputs and scopes first; retain unstable candidates with explicit labels'))
    c.save_json(root/'status.json',dict(status='indices_complete_responses_prepared',n_cases=9,n_pairs=91,
        n_response_tasks=len(tasks),completed_utc=c.now()))
    print('Indices complete; prepared',len(tasks),'conditional response tasks',flush=True)


def plot_surface(meta,values,path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    gi,gj=np.array(meta['grid_i']),np.array(meta['grid_j'])
    ai=int(np.argmin(abs(gi-meta['anchor_i']))); aj=int(np.argmin(abs(gj-meta['anchor_j'])))
    data=values.reshape(len(gi),len(gj),len(c.OXYGEN),2)
    contrasts=[interaction_contrast(data[:,:,oi,k],ai,aj) for k in range(2) for oi in range(len(c.OXYGEN))]
    fig,axes=plt.subplots(4,4,figsize=(15,13),layout='constrained')
    for k,output in enumerate(c.OUTPUTS):
        finite=data[:,:,:,k]; lo,hi=float(finite.min()),float(finite.max())
        if hi==lo: hi=lo+max(1.,abs(lo))*1e-12
        limit=max(1e-12,max(float(np.max(abs(contrasts[k*4+oi]))) for oi in range(4)))
        for oi,o in enumerate(c.OXYGEN):
            ax=axes[k*2,oi]; ax2=axes[k*2+1,oi]
            image=ax.pcolormesh(gi,gj,data[:,:,oi,k].T,shading='auto',cmap='viridis',vmin=lo,vmax=hi,rasterized=True)
            contrast_image=ax2.pcolormesh(gi,gj,contrasts[k*4+oi].T,shading='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True)
            for a in (ax,ax2):
                a.axvline(meta['anchor_i'],color='black',ls='--',lw=.6)
                a.axhline(meta['anchor_j'],color='black',ls='--',lw=.6)
                a.set_xlabel(meta['parameter_i']); a.set_ylabel(meta['parameter_j'])
            if oi==0:
                ax.set_title(('Mean ploidy' if k==0 else 'Growth (day⁻¹)')+'; O₂='+str(o)+'%')
                ax2.set_title('Nonadditive contrast; O₂='+str(o)+'%')
            else:
                ax.set_title('O₂='+str(o)+'%'); ax2.set_title('O₂='+str(o)+'%')
            if meta['transform_i']=='log10':
                ax.set_xscale('log'); ax2.set_xscale('log')
            if meta['transform_j']=='log10':
                ax.set_yscale('log'); ax2.set_yscale('log')
        fig.colorbar(image,ax=axes[k*2,:],label='Ploidy' if k==0 else 'day⁻¹',shrink=.8)
        fig.colorbar(contrast_image,ax=axes[k*2+1,:],label='Contrast (ploidy)' if k==0 else 'Contrast (day⁻¹)',shrink=.8)
    fig.suptitle(meta['scope']+' / '+meta['fit_seed']+': conditional response surfaces\nOther parameters held at optimum; dashed lines = reference values',fontsize=12)
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path.with_suffix('.pdf')); fig.savefig(path.with_suffix('.png'),dpi=120); plt.close(fig)
    return data,contrasts


def finish_report(root,plan):
    tasks=c.read(root/'slurm/response_tasks.tsv'); inventory=[]
    for task in tasks:
        part=root/task['part']; meta=c.load_json(part/'metadata.json')
        if not c.complete(part/'outputs.tsv.gz'): raise ValueError('Incomplete response surface: '+task['part'])
        values,qc=c.output_array(part/'outputs.tsv.gz',meta['n_samples'],c.OXYGEN)
        case=dict(scope=task['scope'],case_id=task['case_id'],fit_seed=task['fit_seed'])
        case_root=c.case_directory(root,case)
        data,contrasts=plot_surface(meta,values,case_root/'figures/responses'/(meta['parameter_i']+'__'+meta['parameter_j']))
        rows=[]
        for i,vi in enumerate(meta['grid_i']):
            for j,vj in enumerate(meta['grid_j']):
                for oi,o in enumerate(c.OXYGEN):
                    rows.append(dict(scope=task['scope'],fit_seed=task['fit_seed'],parameter_i=meta['parameter_i'],parameter_j=meta['parameter_j'],
                        value_i=vi,value_j=vj,O2_pct=o,dominant_mean_ploidy=data[i,j,oi,0],dominant_growth_rate=data[i,j,oi,1],
                        ploidy_nonadditive_contrast=contrasts[oi][i,j],growth_nonadditive_contrast=contrasts[4+oi][i,j]))
        c.table(part/'response_surface_source.tsv.gz',rows)
        inventory.append(dict(**task,**qc,outputs_sha256=c.sha(part/'outputs.tsv.gz')))
    c.table(root/'response_surface_inventory.tsv',inventory)
    for scope in ('global','neighborhood10pct'):
        target=root/scope/'pilot'; rr=c.read(target/'indices_scope_summary.tsv')
        pairs=[r for r in rr if r['index']=='S2']
        passed=sum(r['convergence_status']=='passed' for r in pairs)
        top=sorted([r for r in pairs if np.isfinite(float(r['value_mean']))],key=lambda r:float(r['value_mean']),reverse=True)[:12]
        c.table(target/'top_pair_estimates.tsv',top)
        lines=['# '+scope+' second-order Sobol pilot','',
            'Scope: '+('full independent in-vivo bounds; five randomized Sobol scrambles.' if scope=='global' else
            'eight representative fitted endpoints; natural-span ±10% clipped neighborhoods; five scrambles per endpoint. This is not a 500-endpoint result.'),'',
            'Oxygen: '+', '.join(map(str,c.OXYGEN))+'%; 14 parameters, all 91 pairs; N='+','.join(map(str,c.RESOLUTIONS))+'.',
            'Each neighborhood is averaged over five scrambles before taking the median across eight endpoints.',
            'Input distributions retain the original independent uniform/log-uniform choices.',
            'Indices describe relative variance within this scope; output variances are retained for comparison.',
            'Raw negative estimates and undefined/unstable cells are retained. A negative estimate is sampling error, not antagonism.',
            'SALib 95% bootstrap intervals are estimator diagnostics; randomized-QMC repeat/resolution checks are reported separately.',
            'Diagnostic tolerances: five-repeat range <=0.05 and mean change N'+str(c.RESOLUTIONS[-2])+' to N'+str(c.RESOLUTIONS[-1])+' <=0.02. These are operational checks, not a proof of convergence.',
            'Conditional surfaces hold the other 12 parameters at the indicated optimum. Their nonadditive contrast is not the marginalized Sobol S2.',
            'Pairwise S2 excludes third- and higher-order interactions. S2/ST estimates do not give biological effect direction.',
            '',f'S2 cells passing all applicable checks: {passed}/{len(pairs)}.',
            'For the neighborhood median, this count requires all eight endpoint cells to pass.',
            '', '## Strongest estimated pair/output/oxygen cells','',
            '| Pair | Output | O2 (%) | S2 | Diagnostic |','|---|---|---:|---:|---|']
        for r in top:
            lines.append('| '+r['parameter_i']+' × '+r['parameter_j']+' | '+r['output']+' | '+r['O2_pct']+' | '+format(float(r['value_mean']),'.4g')+' | '+r['convergence_status']+' |')
        lines+=['','These trial oxygen points do not estimate a continuous low/high-band AUC. Claims about the full oxygen range require the 201-point extension.']
        (target/'interpretation.md').write_text('\n'.join(lines)+'\n')
    maximum=plan['max_N']
    # Counts exclude optional recomputations and scheduler overhead; measured time kept per task.
    elapsed=[]; output_bytes=[]
    for task in c.read(root/'slurm/task_manifest.tsv'):
        p=root/task['part']/'outputs.tsv.gz'; receipt=c.load_json(str(p)+'.receipt.json')
        elapsed.append(float(receipt['elapsed_seconds'])); output_bytes.append(p.stat().st_size)
    p50,p90=np.quantile(elapsed,[.5,.9]); bytes_per_evaluation=sum(output_bytes)/plan['n_operator_evaluations']
    full_evaluations=(500+1)*5*maximum*c.STEP*201
    projection=dict(N=maximum,n_full_neighborhoods=500,n_global_designs=1,oxygen_points=201,
        full_operator_evaluations=full_evaluations,
        measured_task_seconds_p50=float(p50),measured_task_seconds_p90=float(p90),
        estimated_serial_core_hours_p50=float(p50*full_evaluations/(maximum*4)/3600),
        estimated_serial_core_hours_p90=float(p90*full_evaluations/(maximum*4)/3600),
        estimated_compressed_output_GiB=bytes_per_evaluation*full_evaluations/1024**3,
        limitation='Linear pilot extrapolation; no queue wait, compilation, analysis or storage overhead; target N needs review of convergence diagnostics.')
    c.save_json(root/'full_scale_projection.json',projection)
    c.save_json(root/'status.json',dict(status='complete',stage=plan['stage'],resolutions=c.RESOLUTIONS,
        n_global_designs=1,n_neighborhoods=8,
        n_pairs=91,n_repeats=5,n_response_surfaces=len(tasks),completed_utc=c.now(),
        convergence_policy='diagnostic_only',all_500_neighborhoods_run=False))
    paths=[]
    for p in sorted(root.rglob('*')):
        if p.is_file() and p.suffix not in ('.lock','.incomplete') and '/logs/' not in str(p) and p.name not in ('collection_manifest.tsv',):
            paths.append(dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,sha256=c.sha(p)))
    c.table(root/'collection_manifest.tsv',paths)
    print('Completed separate global and eight-neighborhood pilots',flush=True)
