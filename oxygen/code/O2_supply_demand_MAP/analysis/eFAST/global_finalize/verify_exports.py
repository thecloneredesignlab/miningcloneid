#!/usr/bin/env python3
"""Independently validate collected global exports; never evaluate the model."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import numpy as np

from common import ACTIVE,OUTPUTS,METRICS,read,table,load,save,sha,now


def validate(root):
    root=Path(root);p=load(root/'plan.json')
    if p['scope']!='global' or p['n_repeats']!=5 or p['n_oxygen']!=201 or p['N_hard_limit'] is not None:
        raise ValueError('Scope differs')
    files=load(root/'collection_manifest.json')['files']
    for row in files:
        path=root/row['path']
        if path.stat().st_size!=row['size'] or sha(path)!=row['sha256']:
            raise ValueError('Collected file differs: '+str(path))
    check=load(root/'validation/method_alignment.json')
    if check['status']!='passed' or max(check['output_max_absolute_error'].values())>1e-8:
        raise ValueError('Model alignment not verified')
    with np.load(root/'inputs/sobol_N131072/indices_repeats.npz') as z:
        sv=z['values'].copy();oxygen=z['oxygen'].copy();svar=z['variance'][1].copy()
    rows=read(root/'inputs/sobol_N131072/indices_summary.tsv.gz')
    metrics=[(r['index'],r['parameter_i'],r['parameter_j']) for r in rows if float(r['O2_pct'])==0 and r['output']==OUTPUTS[0]]
    if len(rows)!=47838 or len(set(metrics))!=119 or sv.shape!=(2,5,201,2,119):
        raise ValueError('Sobol completeness differs')
    passed=np.isfinite(sv).all(axis=(0,1))&(np.ptp(sv[1],axis=0)<=.05)&(np.abs(sv[1].mean(axis=0)-sv[0].mean(axis=0))<=.02)
    if not passed.all():raise ValueError('Certified Sobol convergence count differs')
    refs=sv[1].mean(axis=0)[:,:,[metrics.index(m) for m in METRICS]]
    repeats={};checks={};moment_rows=[];class_rows=[]
    for n in p['resolutions']:
        directory=root/'efast'/('N%d'%n)
        done=load(directory/'completion.json')
        if done['N']!=n or done['n_repeats']!=5 or done['n_oxygen']!=201 or done['n_indices']!=11256:
            raise ValueError('eFAST completion differs')
        if sha(directory/'indices_summary.tsv.gz')!=done['summary_sha256']:
            raise ValueError('eFAST summary receipt differs')
        with np.load(directory/'indices_repeats.npz') as z:
            v=z['values'].copy();var=z['variance'].copy()
            if not np.array_equal(z['oxygen'],oxygen):raise ValueError('Oxygen axes differ')
        if v.shape!=(5,201,2,28) or var.shape!=(5,201,2,14):raise ValueError('Repeat/trajectory dimensions differ')
        mean=np.add.reduce(v,axis=0)/5
        spread=np.max(v,axis=0)-np.min(v,axis=0)
        i=p['resolutions'].index(n);previous=p['resolutions'][i-1] if i else None
        delta=np.abs(mean-repeats[previous].mean(axis=0)) if previous else np.full_like(mean,np.nan)
        data=read(directory/'indices_summary.tsv.gz')
        keys={(r['output'],round(float(r['O2_pct']),12),r['index'],r['parameter_i']) for r in data}
        if len(data)!=11256 or len(keys)!=11256:raise ValueError('Missing or duplicate eFAST rows')
        statuses={}
        for row in data:
            k=OUTPUTS.index(row['output']);oi=int(round(float(row['O2_pct'])*40));j=METRICS.index((row['index'],row['parameter_i'],''))
            for field,expected in [('value_mean',mean[oi,k,j]),('repeat_range',spread[oi,k,j]),('resolution_delta',delta[oi,k,j])]:
                if not np.isclose(float(row[field]),expected,rtol=0,atol=1e-12,equal_nan=True):raise ValueError('Independent '+field+' differs')
            finite=np.isfinite(v[:,oi,k,j]).all()
            status='undefined' if not finite else ('no_adjacent_new_design' if previous is None else ('passed' if spread[oi,k,j]<=.05 and delta[oi,k,j]<=.02 else 'not_passed'))
            if row['convergence_status']!=status or row['included']!='True':raise ValueError('Diagnostic inclusion differs')
            statuses[status]=statuses.get(status,0)+1
        repeats[n]=v
        checks[str(n)]=dict(rows=len(data),repeats=5,oxygen_points=201,parameters=14,status_counts=statuses,
                             independent_mean_range_delta_verified=True)
        classes=read(directory/'classification.tsv')
        if len(classes)!=28 or any(r['peak_threshold']!='none' or 'phase' not in r['bootstrap_unit'] for r in classes):
            raise ValueError('eFAST classification semantics differ')
        for oi,o2 in enumerate(oxygen):
            for k,out in enumerate(OUTPUTS):
                moment_rows.append(dict(efast_N=n,O2_pct=o2,output=out,sobol_variance_mean=float(svar[:,oi,k].mean()),
                    efast_trajectory_variance_mean=float(var[:,oi,k,:].mean()),
                    efast_trajectory_variance_min=float(var[:,oi,k,:].min()),efast_trajectory_variance_max=float(var[:,oi,k,:].max()),
                    interpretation='finite-design output variances; phase/trajectory spread is not a confidence interval'))
    # Cross-method rows are matched by labels, independently of report array order.
    comparative=read(root/'comparison/method_pointwise_comparison.tsv.gz')
    expected_rows=(1+len(p['resolutions']))*11256
    if len(comparative)!=expected_rows:raise ValueError('Comparison row count differs')
    selected=[r for r in comparative if int(r['efast_N']) in repeats]
    keys=set()
    for row in selected:
        n=int(row['efast_N']);oi=int(round(float(row['O2_pct'])*40));k=OUTPUTS.index(row['output']);j=METRICS.index((row['index'],row['parameter'],''))
        key=(n,oi,k,j)
        if key in keys:raise ValueError('Duplicate comparison row')
        keys.add(key)
        e=float(repeats[n][:,oi,k,j].mean());s=float(refs[oi,k,j])
        if not all(np.isclose(float(row[field]),value,rtol=0,atol=1e-12,equal_nan=True) for field,value in [('efast_mean',e),('sobol_mean',s),('difference',e-s)]):
            raise ValueError('Comparison values differ')
    if len(keys)!=len(p['resolutions'])*11256:raise ValueError('Comparison coverage differs')
    rankings=read(root/'sobol/all_parameter_and_pair_rankings.tsv')
    if sum(r['index']=='S2' for r in rankings)!=91*2*4:raise ValueError('91-pair ranking coverage differs')
    if len(read(root/'sobol/focal_36_pair_rankings.tsv'))!=36*2*4:raise ValueError('36 focal-pair ranking coverage differs')
    limits=read(root/'comparison/comparison_color_limits.tsv')
    for out in OUTPUTS:
        for idx in ['S1','ST']:
            pair=[float(r['vmax']) for r in limits if r['output']==out and r['index']==idx]
            if len(pair)!=2 or pair[0]!=pair[1]:raise ValueError('Cross-method color scale mismatch')
    table(root/'comparison/output_variance_comparison.tsv',moment_rows)
    diagnostic_rows=[]
    bands=[('low',0.,1.),('middle',1.,3.),('high',3.,5.)]
    masks=[('all',np.ones(201,dtype=bool)),('low_0_to_1_inclusive',oxygen<=1),
           ('middle_above1_below3',(oxygen>1)&(oxygen<3)),('high_3_to_5_inclusive',oxygen>=3)]
    configurations=[('Sobol',131072,sv[1],sv[0],metrics,root/'sobol')]+[
        ('eFAST',n,repeats[n],None if i==0 else repeats[p['resolutions'][i-1]],METRICS,root/'efast'/('N%d'%n))
        for i,n in enumerate(p['resolutions'])]
    for method,n,v,previous,mm,directory in configurations:
        # Independent trapezoidal integration, without calling the plot helper.
        bands_source=read(directory/'oxygen_band_indices.tsv.gz')
        expected={}
        for name,lo,hi in bands:
            selected=(oxygen>=lo)&(oxygen<=hi)
            by_repeat=np.trapezoid(v[:,selected,:,:],x=oxygen[selected],axis=1)/(hi-lo)
            expected[name]=by_repeat
        if len(bands_source)!=3*2*len(mm):raise ValueError('Band summary coverage differs')
        for row in bands_source:
            j=mm.index((row['index'],row['parameter_i'],row['parameter_j']));k=OUTPUTS.index(row['output'])
            b=expected[row['oxygen_band']][:,k,j]
            for field,value in [('band_mean',b.mean()),('repeat_min',b.min()),('repeat_max',b.max())]:
                if not np.isclose(float(row[field]),value,rtol=0,atol=1e-12,equal_nan=True):raise ValueError('Independent band integration differs')
        spread=np.ptp(v,axis=0)
        delta=np.abs(v.mean(axis=0)-previous.mean(axis=0)) if previous is not None else np.full_like(spread,np.nan)
        for k,out in enumerate(OUTPUTS):
            for idx in ['S1','ST']+(['S2'] if method=='Sobol' else []):
                jj=[j for j,m in enumerate(mm) if m[0]==idx]
                for name,mask in masks:
                    s=spread[:,k,:][mask][:,jj];d=delta[:,k,:][mask][:,jj]
                    finite=np.isfinite(v[:,:,k,:][:,:,jj]).all(axis=0)[mask]
                    passed=finite&(s<=.05)&(d<=.02)
                    diagnostic_rows.append(dict(method=method,N=n,output=out,index=idx,oxygen_band=name,n_total=int(s.size),
                        n_passed=int(passed.sum()) if previous is not None else '',
                        n_not_passed=int((finite&~passed).sum()) if previous is not None else '',n_undefined=int((~finite).sum()),
                        pass_fraction=float(passed.mean()) if previous is not None else '',
                        n_repeat_range_failed=int((s>.05).sum()),n_resolution_delta_failed=int((d>.02).sum()) if previous is not None else '',
                        repeat_range_p95=float(np.quantile(s,.95)),repeat_range_max=float(np.max(s)),
                        resolution_delta_max=float(np.max(d)) if previous is not None else '',
                        diagnostic_status='assessed' if previous is not None else 'no_adjacent_new_design',
                        bands_for_counts='disjoint; physical-O2 integration bands are separately inclusive'))
    table(root/'comparison/diagnostic_summary_disjoint_bands.tsv',diagnostic_rows)
    # At O2=0, h=C**n/(C**n+0**n)=1 for the allowed C>0,n>0.
    # Both parameters enter this fixed-O2 model only through h. This is a
    # structural null, independent of the empirical convergence thresholds.
    null_rows=[]
    for method,n,v,mm in [('Sobol',131072,sv[1],metrics)]+[
            ('eFAST',n,v,METRICS) for n,v in repeats.items()]:
        for k,out in enumerate(OUTPUTS):
            for param in ['O2_crit','n_O']:
                for idx in ['S1','ST']:
                    j=mm.index((idx,param,''))
                    for ri,value in enumerate(v[:,0,k,j],1):
                        null_rows.append(dict(method=method,N=n,repeat=ri,O2_pct=0.,output=out,parameter=param,index=idx,
                            estimate=float(value),expected_structural_value=0.,absolute_error=float(abs(value)),
                            interpretation='structurally inactive input at O2=0; not a significance test'))
    table(root/'comparison/structural_null_O2_zero.tsv',null_rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for k,out in enumerate(OUTPUTS):
        for col,param in enumerate(['O2_crit','n_O']):
            ax=axes[k,col];ns=p['resolutions']
            for ri in range(5):
                j=METRICS.index(('ST',param,''))
                ax.plot(ns,[repeats[n][ri,0,k,j] for n in ns],'.-',color='#762a83',alpha=.5,lw=.7)
            ax.axhline(0,color='#222222',lw=.8,label='Structural null / Sobol ST=0')
            ax.set_xscale('log',base=2);ax.set_xticks(ns,[str(n) for n in ns]);ax.minorticks_off()
            ax.set_title(('Ploidy' if k==0 else 'Net live-cell growth')+'; '+param)
            ax.set_xlabel('eFAST N');ax.set_ylabel('Estimated ST at O2=0');ax.legend(fontsize=8)
    fig.suptitle('Structural zero-effect check; five independent eFAST phase curves')
    fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(root/'comparison/figures/structural_null_O2_zero.pdf',bbox_inches='tight')
    fig.savefig(root/'comparison/figures/structural_null_O2_zero.png',dpi=140,bbox_inches='tight');plt.close(fig)
    # Explicit methods comparison receipt; old N513 remains an uncertified historical reference.
    alignment=[dict(item='parameter_order',Sobol=';'.join(ACTIVE),eFAST=';'.join(ACTIVE),status='matched'),
        dict(item='bounds_distributions',Sobol=sha(root/'inputs/parameter_bounds.tsv'),eFAST=sha(root/'inputs/legacy_parameter_ranges.tsv'),status='matched'),
        dict(item='oxygen',Sobol='201 points; 0:0.025:5 percent',eFAST='201 points; 0:0.025:5 percent',status='matched'),
        dict(item='outputs',Sobol=';'.join(OUTPUTS),eFAST=';'.join(OUTPUTS),status='matched'),
        dict(item='model_operator',Sobol=p['original_plan_sha256'],eFAST='same frozen source; wrapper body identity and cached-output parity checked',status='matched'),
        dict(item='runtime',Sobol=p['sif']['sha256'],eFAST=p['sif']['sha256'],status='matched'),
        dict(item='repetitions',Sobol='5 independent scrambles',eFAST='5 independent phase repeats',status='different randomization; same repeat count'),
        dict(item='estimator',Sobol='SALib.sample.sobol; analyze.sobol; calc_second_order=True',eFAST='SALib.sample.fast_sampler; analyze.fast.compute_orders; M4; public-interface parity verified',status='intentional method difference'),
        dict(item='N',Sobol='131072',eFAST=';'.join(map(str,p['resolutions'])),status='different design definitions; not matched by numeric N'),
        dict(item='legacy_N513',Sobol='not applicable',eFAST='old best/config identities incomplete',status='historical reference only')]
    if read(root/'inputs/parameter_bounds.tsv')!=read(root/'inputs/legacy_parameter_ranges.tsv'):
        raise ValueError('Input bounds/distributions differ')
    table(root/'method_alignment.tsv',alignment)
    source_plan=root/'inputs/source_sobol_submission_plan.json'
    if sha(source_plan)!=p['original_plan_sha256']:
        raise ValueError('Copied original Sobol provenance differs')
    original=load(source_plan);case=original['cases'][0]
    for name,field in [('best_params.tsv','best_params_sha256'),('fit_config.rds','fit_config_sha256'),('parameter_table.csv','parameter_table_sha256')]:
        if sha(root/'inputs/reference_seed25'/name)!=case[field]:raise ValueError('Global reference input differs')
    if sha(root/'inputs/source_sobol_evaluate.R')!=original['engine_sha256']:
        raise ValueError('Copied canonical evaluator differs')
    process={r['parameter']:r['parameter_group'] for r in read(root/'inputs/figure4_parameter_groups_source.tsv')}
    equation={r['parameter']:r for r in read(root/'inputs/model_structure/parameter_equation_map.tsv')}
    table(root/'process_group_membership.tsv',[dict(parameter=param,figure4_display_group=process[param],
        primary_group_for_descriptive_summary=equation[param]['primary_group'],shared_pathways=equation[param]['pathways']) for param in ACTIVE])
    result=dict(status='passed',verified_utc=now(),source_files_verified=len(files),sobol_entries_verified=47838,
                sobol_all_operational_diagnostics_passed=True,efast=checks,method_comparison_rows_verified=len(keys),
                all_pairs=91,focal_pairs=36,matched_color_limits_verified=True,global_reference_identity_verified=True,
                structural_null_rows=len(null_rows),band_integrals_independently_verified=True,
                diagnostic_rows_disjoint_bands=len(diagnostic_rows),model_calls=0,estimator_calls=0,
                neighborhood_files_read=0,PDF_visual_validation='recorded separately')
    save(root/'validation/independent_exports.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    validate(parser.parse_args().root)
