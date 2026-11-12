"""Source-backed global sensitivity tables, figures and scientific interpretation.

This module reads indices only; it never calls a sensitivity estimator or model.
"""
from pathlib import Path
import shutil
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch

from common import ACTIVE, OUTPUTS, METRICS, now, sha, load, save, read, table
import plot_helpers as gp

FOCAL = {'p_mis_base', 'mu_hp', 'gamma_mu'}
LABELS = {'dominant_mean_ploidy':'Ploidy', 'dominant_growth_rate':'Net live-cell growth rate'}
BANDS = [('low',0.,1.),('middle',1.,3.),('high',3.,5.)]
NOTE = ('P: biological process; C: Figure 4B correlation O2 group; S: ploidy sensitivity O2 preference.\n'
        'S is also ploidy-based on growth panels; exploratory grouping is separate from convergence.\n'
        'O2 axis: linear 0-0.025%, logarithmic above; dashed boundaries: 1% and 3%.')


def sobol_data(root):
    source=Path(root)/'inputs/sobol_N131072'
    done=load(source/'completion.json')
    if sha(source/'indices_summary.tsv.gz')!=done['summary_sha256']:
        raise ValueError('Sobol source checksum differs')
    rows=read(source/'indices_summary.tsv.gz')
    metrics=[(r['index'],r['parameter_i'],r['parameter_j']) for r in rows if float(r['O2_pct'])==0 and r['output']==OUTPUTS[0]]
    if len(rows)!=47838 or len(metrics)!=119 or len(set(metrics))!=119:
        raise ValueError('Incomplete/duplicate Sobol sources')
    with np.load(source/'indices_repeats.npz') as z:
        v=z['values']; oxygen=z['oxygen']; resolutions=z['resolutions']
    if v.shape!=(2,5,201,2,119) or not np.array_equal(resolutions,[65536,131072]):
        raise ValueError('Sobol resolution/repeat axes differ')
    mean=np.mean(v[1],axis=0)
    exported=np.array([float(r['value_mean']) for r in rows]).reshape(201,2,119)
    if not np.isfinite(v).all() or not np.allclose(mean,exported,rtol=0,atol=1e-12):
        raise ValueError('Sobol repeats disagree with exported indices')
    return v,oxygen,metrics


def annotations(root):
    root=Path(root)
    process={r['parameter']:r['parameter_group'] for r in read(root/'inputs/figure4_parameter_groups_source.tsv')}
    correlation={r['parameter']:r['o2_association_group'] for r in read(root/'inputs/figure4b_o2_classification_source.tsv')}
    palette={r['parameter_group']:r['color'] for r in read(root/'inputs/figure4_parameter_group_palette_source.tsv')}
    return process,correlation,palette


def heatmaps(root,dest,values,metrics,oxygen,title,phase=False):
    dest=Path(dest);figures=dest/'figures';figures.mkdir(parents=True,exist_ok=True)
    mean=np.mean(values,axis=0)
    classes=gp.classify(values,metrics,oxygen)
    for row in classes:
        if phase:row['bootstrap_unit']='whole independent phase-repeat curve'
    table(dest/'classification.tsv',classes)
    process,correlation,palette=annotations(root)
    selected={kind:gp.row_order([i for i,m in enumerate(metrics) if m[0]==kind],classes) for kind in ['S1','ST']}
    if any(m[0]=='S2' for m in metrics):
        selected['S2']=gp.row_order([i for i,m in enumerate(metrics) if m[0]=='S2'],classes)
        selected['S2_focal']=gp.row_order([i for i,m in enumerate(metrics) if m[0]=='S2' and bool(set(m[1:])&FOCAL)],classes)
    table(dest/'row_order.tsv',[dict(panel=kind,rank=j+1,**classes[i]) for kind,ids in selected.items() for j,i in enumerate(ids)])
    for kind,ids in selected.items():
        fig=plt.figure(figsize=(20,max(6,.27*len(ids)+2)))
        grid=fig.add_gridspec(1,2,wspace=.65)
        for k,out in enumerate(OUTPUTS):
            gp.panel(fig,grid[0,k],mean[:,k],ids,kind.split('_')[0],metrics,classes,oxygen,
                     process,correlation,palette,LABELS[out])
        fig.suptitle(title+'; '+kind)
        fig.text(.5,.015,NOTE,ha='center',fontsize=8)
        fig.subplots_adjust(top=.85,bottom=.19 if len(ids)<20 else .09)
        gp.savefig(fig,figures/(kind.lower()+'_parameter_by_oxygen'))
    kinds=['S1','ST']+(['S2_focal'] if 'S2_focal' in selected else [])
    fig=plt.figure(figsize=(20,7*len(kinds)+2))
    grid=fig.add_gridspec(len(kinds),2,wspace=.65,hspace=.4,
                         height_ratios=[max(14,len(selected[k])) for k in kinds])
    for row,kind in enumerate(kinds):
        for k,out in enumerate(OUTPUTS):
            gp.panel(fig,grid[row,k],mean[:,k],selected[kind],kind.split('_')[0],metrics,classes,
                     oxygen,process,correlation,palette,kind+': '+LABELS[out])
    fig.suptitle(title+'; five independent repeats')
    fig.text(.5,.012,NOTE,ha='center',fontsize=9)
    fig.subplots_adjust(top=.94,bottom=.09)
    gp.savefig(fig,figures/'sensitivity_summary')
    palette_rows=read(Path(root)/'inputs/figure4_parameter_group_palette_source.tsv')
    handles=[Patch(facecolor=r['color'],label=r['display_label']) for r in palette_rows
             if r['parameter_group'] in {process[p] for p in ACTIVE}]
    handles += [Patch(facecolor=gp.COLORS[label],label=label) for label in ['High O2','Low O2','No resolved O2 preference']]
    fig=plt.figure(figsize=(15,3.5))
    fig.legend(handles=handles,loc='center',ncol=3,frameon=False)
    fig.text(.5,.025,NOTE+'\nS2 rows have P/C only; each split cell represents the two pair members.',ha='center',fontsize=9)
    gp.savefig(fig,figures/'annotation_legend')
    return classes,selected


def scientific_tables(root,dest,values,oxygen,metrics,method,n):
    """Each band is integrated in physical O2; ST process means are descriptive."""
    dest=Path(dest);dest.mkdir(parents=True,exist_ok=True)
    mean=np.mean(values,axis=0)
    groups={r['parameter']:r['primary_group'] for r in read(Path(root)/'inputs/model_structure/parameter_equation_map.tsv')}
    bandrows=[];group_rows=[];hypotheses=[]
    for k,out in enumerate(OUTPUTS):
        for band,lo,hi in BANDS:
            w=gp.weights(oxygen,lo,hi)
            each=np.einsum('rom,o->rm',values[:,:,k],w)
            for j,(idx,pi,pj) in enumerate(metrics):
                bandrows.append(dict(method=method,N=n,output=out,oxygen_band=band,index=idx,parameter_i=pi,parameter_j=pj,
                    focal=bool({pi,pj}&FOCAL),band_mean=float(each[:,j].mean()),repeat_min=float(each[:,j].min()),
                    repeat_max=float(each[:,j].max()),positive_repeats=int((each[:,j]>0).sum())))
            for idx in ['S1','ST']:
                for group in sorted(set(groups[p] for p in ACTIVE)):
                    ids=[j for j,m in enumerate(metrics) if m[0]==idx and not m[2] and groups[m[1]]==group]
                    v=each[:,ids].mean(axis=1)
                    group_rows.append(dict(method=method,N=n,output=out,index=idx,oxygen_band=band,group=group,n_parameters=len(ids),
                        descriptive_parameter_mean=float(v.mean()),repeat_min=float(v.min()),repeat_max=float(v.max()),
                        interpretation='mean of parameter indices; not a group variance fraction'))
        for idx in ['S1','ST']:
            grouped={}
            for band,lo,hi in [BANDS[0],BANDS[-1]]:
                w=gp.weights(oxygen,lo,hi)
                for group in sorted(set(groups[p] for p in ACTIVE)):
                    ids=[j for j,m in enumerate(metrics) if m[0]==idx and not m[2] and groups[m[1]]==group]
                    grouped[band,group]=np.einsum('rom,o->rm',values[:,:,k,:][:,:,ids],w).mean(axis=1)
            death=grouped['low','death'];others=np.stack([v for (band,g),v in grouped.items() if band=='low' and g!='death'])
            for name,test,estimate in [
                ('low_O2_death_largest_descriptive_group_mean',death>=others.max(axis=0),death.mean()-others.mean(axis=1).max()),
                ('missegregation_high_above_low',grouped['high','missegregation']>grouped['low','missegregation'],
                 (grouped['high','missegregation']-grouped['low','missegregation']).mean()),
                ('buffering_high_above_low',grouped['high','buffering']>grouped['low','buffering'],
                 (grouped['high','buffering']-grouped['low','buffering']).mean())]:
                hypotheses.append(dict(method=method,N=n,output=out,index=idx,condition=name,contrast=float(estimate),
                    supporting_repeats=int(test.sum()),n_repeats=5,
                    interpretation='descriptive within-model check; repeats are numerical replicates, not biological confidence'))
    table(dest/'oxygen_band_indices.tsv.gz',bandrows)
    table(dest/'process_descriptive_summary.tsv',group_rows)
    table(dest/'noemi_hypothesis_checks.tsv',hypotheses)
    ranking=[]
    for k,out in enumerate(OUTPUTS):
        for idx in ['S1','ST','S2']:
            ids=[j for j,m in enumerate(metrics) if m[0]==idx]
            for band,lo,hi in BANDS+[('all',0.,5.)]:
                score=np.einsum('om,o->m',mean[:,k],gp.weights(oxygen,lo,hi))
                for rank,j in enumerate(sorted(ids,key=lambda j:(-score[j],j)),1):
                    series=mean[:,k,j];peak=int(np.argmax(series))
                    ranking.append(dict(method=method,N=n,output=out,index=idx,oxygen_band=band,rank=rank,
                        parameter_i=metrics[j][1],parameter_j=metrics[j][2],focal=bool(set(metrics[j][1:])&FOCAL),
                        band_mean=float(score[j]),peak=float(series[peak]),peak_O2_pct=float(oxygen[peak])))
    table(dest/'all_parameter_and_pair_rankings.tsv',ranking)
    table(dest/'focal_36_pair_rankings.tsv',[r for r in ranking if r['index']=='S2' and r['focal']],list(ranking[0]))
    gaps=[]
    for k,out in enumerate(OUTPUTS):
        for pi,param in enumerate(ACTIVE):
            i=metrics.index(('S1',param,''));j=metrics.index(('ST',param,''))
            delta=values[:,:,k,j]-values[:,:,k,i]
            for oi,o2 in enumerate(oxygen):
                gaps.append(dict(output=out,parameter=param,O2_pct=o2,ST_minus_S1=float(delta[:,oi].mean()),
                    repeat_min=float(delta[:,oi].min()),repeat_max=float(delta[:,oi].max()),
                    interpretation='all interaction orders involving this parameter; not a named pair'))
    table(dest/'ST_minus_S1.tsv.gz',gaps)
    return ranking,group_rows,hypotheses


def key_pair_curves(dest,values,oxygen,metrics):
    mean=values.mean(axis=0);dest=Path(dest)
    fig,axes=plt.subplots(6,2,figsize=(14,17))
    selected=[]
    for k,out in enumerate(OUTPUTS):
        ids=[j for j,m in enumerate(metrics) if m[0]=='S2']
        ids=sorted(ids,key=lambda j:-float(np.max(mean[:,k,j])))[:6]
        for row,j in enumerate(ids):
            ax=axes[row,k];ax.plot(oxygen,mean[:,k,j],color='#54278f',lw=1.1)
            ax.fill_between(oxygen,values[:,:,k,j].min(axis=0),values[:,:,k,j].max(axis=0),color='#54278f',alpha=.17)
            ax.axhline(0,color='#555555',lw=.5);gp.oxygen_axis(ax)
            ax.set_title(' x '.join(metrics[j][1:]),fontsize=9);ax.set_ylabel('S2')
            selected.append(dict(output=out,rank_by_peak=row+1,parameter_i=metrics[j][1],parameter_j=metrics[j][2],
                                 selection='six largest peaks among all 91 pairs; not just focal pairs'))
    fig.text(.25,.985,'Ploidy',ha='center');fig.text(.75,.985,'Net live-cell growth rate',ha='center')
    fig.suptitle('Sobol N131072: strongest pairs by peak S2',y=1.005)
    fig.text(.5,.008,'Shading: min-max of five scramble repeats, not a 95% confidence interval. Negative estimates are retained.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.035,1,.97]);gp.savefig(fig,dest/'figures/strongest_91_pair_curves')
    table(dest/'strongest_pair_curve_selection.tsv',selected)


def sobol_report(root):
    root=Path(root);dest=root/'sobol';dest.mkdir(exist_ok=True)
    v,oxygen,metrics=sobol_data(root)
    heatmaps(root,dest,v[1],metrics,oxygen,'Global Sobol N131072')
    rankings,groups,hypotheses=scientific_tables(root,dest,v[1],oxygen,metrics,'Sobol',131072)
    key_pair_curves(dest,v[1],oxygen,metrics)
    stability=[]
    for oi,o2 in enumerate(oxygen):
        for k,out in enumerate(OUTPUTS):
            for idx in ['S1','ST','S2']:
                ids=[j for j,m in enumerate(metrics) if m[0]==idx]
                a,b=v[0,:,oi,k,:].mean(axis=0)[ids],v[1,:,oi,k,:].mean(axis=0)[ids]
                stability.append(dict(O2_pct=o2,output=out,index=idx,rank_spearman=float(spearmanr(a,b).statistic),
                    max_absolute_change=float(np.max(np.abs(a-b))),N_low=65536,N_high=131072))
    table(dest/'ranking_stability_N65536_N131072.tsv',stability)
    lines=['# 全域 Sobol N131072 科学解读','',
        '使用既有5次独立scramble、201个固定氧点、14个独立输入。没有重算模型或Sobol估计。',
        '全部47,838个单元通过既定操作性收敛诊断；通过不是交互显著性或机制方向证明。','',
        '## Noemi假设的分项检查','',
        '以下是参数索引的过程组平均，不能解释为过程方差份额；ST包含重叠交互。mu_hp、gamma_mu同时影响死亡和死亡耦合的误分离，因此过程并不彼此独立。','']
    for out in OUTPUTS:
        lines += ['### '+LABELS[out],'']
        for idx in ['S1','ST']:
            low=sorted([r for r in groups if r['output']==out and r['index']==idx and r['oxygen_band']=='low'],key=lambda r:-r['descriptive_parameter_mean'])
            lines.append('- %s：低氧描述性过程均值最高为%s（%.5f）。'%(idx,low[0]['group'],low[0]['descriptive_parameter_mean']))
            for group in ['death','missegregation','buffering']:
                scores={r['oxygen_band']:r['descriptive_parameter_mean'] for r in groups if r['output']==out and r['index']==idx and r['group']==group}
                lines.append('  - %s：低氧%.5f，高氧%.5f，高−低%+.5f。'%(group,scores['low'],scores['high'],scores['high']-scores['low']))
        lines += ['', '按全氧区物理氧浓度积分的S2前六对：','']
        for r in [r for r in rankings if r['output']==out and r['index']=='S2' and r['oxygen_band']=='all'][:6]:
            lines.append('- %s × %s：氧区平均S2 %.5f；峰值%.5f，O2=%.3f%%。'%(r['parameter_i'],r['parameter_j'],r['band_mean'],r['peak'],r['peak_O2_pct']))
        lines.append('')
    bands=read(dest/'oxygen_band_indices.tsv.gz')
    lines += ['## 对原假设的综合判断','',
        '倍性结果支持部分氧依赖转变：死亡参数在低氧更敏感，基线误分离与缓冲在高氧增强。'
        '但低氧死亡组与误分离组ST均值非常接近，不能写成死亡独占主导；整个误分离参数组也没有一致随高氧增强。'
        '增长率必须单独解释，其缓冲组贡献没有显示相同的高氧增强。','',
        '以下逐参数ST带平均值说明为什么不能用一个过程组概括所有成员：','']
    for out in OUTPUTS:
        for param in ['p_mis_base','p_misseg','mu_hp','gamma_mu','buffer_smax','buffer_beta','buffer_n_exp','lam_max']:
            scores={r['oxygen_band']:float(r['band_mean']) for r in bands if r['output']==out and r['index']=='ST' and r['parameter_i']==param}
            lines.append('- %s / %s：低氧%.5f，高氧%.5f。'%(LABELS[out],param,scores['low'],scores['high']))
    lines += ['## 解释边界','',
        '- 高ST但较低S1说明该参数主要通过非加性作用影响输出，不能单凭ST−S1确定搭档。',
        '- S2量化纯两两交互方差份额，不提供协同、拮抗或补偿方向，也不包括三阶以上作用。',
        '- 关键响应面沿用既有full_N4096/global/full/global/responses；每个面固定其他参数于seed25参考值，不能当成对其他参数边缘化的S2。',
        '- p_wgd范围为1e-5至1e-3，未检验精确零；接近拟合下界不能证明WGD生物学缺失或唯一可辨识。',
        '- 当前eFAST比较另行生成；不再用会议中未定位参数/氧点的0.6对0.8作为最终证据。',
        '- 邻域继续封存；本报告不恢复或重新汇总邻域。','',
        '完整91对排名、焦点36对、各氧区和重复间差异见同目录源表。']
    (dest/'interpretation_zh.md').write_text('\n'.join(lines)+'\n')
    save(dest/'completion.json',dict(status='complete',N=131072,source_summary_sha256=sha(root/'inputs/sobol_N131072/indices_summary.tsv.gz'),
        model_calls=0,sobol_estimator_calls=0,all_pairs=91,focal_pairs=36,completed_utc=now()))


def efast_report(root,n):
    root=Path(root);dest=root/'efast'/('N%d'%n)
    with np.load(dest/'indices_repeats.npz') as z:values=z['values'];oxygen=z['oxygen']
    heatmaps(root,dest,values,METRICS,oxygen,'Global eFAST N%d M4'%n,phase=True)
    scientific_tables(root,dest,values,oxygen,METRICS,'eFAST',n)


def legacy_values(root):
    rows=read(Path(root)/'inputs/legacy_efast_indices.tsv')
    result=np.full((5,201,2,28),np.nan)
    for r in rows:
        if int(r['N'])!=513:continue
        ri=int(r['replicate'])-1;oi=int(round(float(r['O2_pct'])*40));k=OUTPUTS.index(r['output']);pi=ACTIVE.index(r['parameter'])
        if np.isfinite(result[ri,oi,k,pi]):raise ValueError('Duplicate legacy estimate')
        result[ri,oi,k,pi]=float(r['S1']);result[ri,oi,k,pi+14]=float(r['ST'])
    if not np.isfinite(result).all():raise ValueError('Legacy five-phase data incomplete')
    return result


def comparison_figures(root,dest,series,sobol,oxygen):
    dest=Path(dest);dest.joinpath('figures').mkdir(exist_ok=True)
    colors={n:color for n,color in zip(sorted(series),['#9b9b9b','#c5b0d5','#9970ab','#762a83','#40004b','#2166ac'])}
    with PdfPages(dest/'figures/parameter_comparison_curves.pdf') as pdf:
        for start in range(0,14,7):
            fig,axes=plt.subplots(7,2,figsize=(15,19))
            for row,pi in enumerate(range(start,min(start+7,14))):
                for k,out in enumerate(OUTPUTS):
                    ax=axes[row,k]
                    for idx,offset,style in [('S1',0,'--'),('ST',14,'-')]:
                        j=pi+offset
                        ref=sobol[:,:,k,j];ax.plot(oxygen,ref.mean(axis=0),style,color='#222222',lw=1.25,label='Sobol '+idx)
                        for n in sorted(series)[-2:]:
                            v=series[n][:,:,k,j]
                            ax.plot(oxygen,v.mean(axis=0),style,color=colors[n],lw=1,label='eFAST N%d %s'%(n,idx))
                            if idx=='ST':ax.fill_between(oxygen,v.min(axis=0),v.max(axis=0),color=colors[n],alpha=.1)
                    gp.oxygen_axis(ax);ax.set_title(ACTIVE[pi]+'; '+LABELS[out],fontsize=9)
                    ax.axhline(0,color='#777777',lw=.4)
                    if row==0:ax.legend(fontsize=7,ncol=2)
            fig.suptitle('S1 dashed / ST solid; highest two new eFAST designs versus Sobol N131072')
            fig.text(.5,.008,'ST shading: five independent phase repeats min-max; not a confidence interval.',ha='center')
            fig.tight_layout(rect=[0,.03,1,.97]);pdf.savefig(fig);plt.close(fig)
    last=max(series);smean=sobol.mean(axis=0);emean=series[last].mean(axis=0)
    classes=gp.classify(sobol,METRICS,oxygen);process,correlation,palette=annotations(root)
    limits=[]
    for kind in ['S1','ST']:
        ids=gp.row_order([j for j,m in enumerate(METRICS) if m[0]==kind],classes)
        fig=plt.figure(figsize=(20,13));grid=fig.add_gridspec(2,2,wspace=.65,hspace=.5)
        for row,(method,data,rep) in enumerate([('Sobol N131072',smean,sobol),('eFAST N%d'%last,emean,series[last])]):
            own=gp.classify(rep,METRICS,oxygen)
            for k,out in enumerate(OUTPUTS):
                vmax=float(max(np.abs(smean[:,k,ids]).max(),np.abs(emean[:,k,ids]).max()))
                gp.panel(fig,grid[row,k],data[:,k],ids,kind,METRICS,own,oxygen,process,correlation,palette,
                         method+': '+LABELS[out],vmax_override=vmax)
                limits.append(dict(index=kind,output=out,method=method,vmax=vmax,row_order_source='Sobol N131072 ploidy'))
        fig.suptitle(kind+': matched row order and shared positive color limits within each output')
        fig.text(.5,.015,NOTE,ha='center',fontsize=9);fig.subplots_adjust(top=.92,bottom=.13)
        gp.savefig(fig,dest/'figures'/('method_comparison_'+kind.lower()))
    table(dest/'comparison_color_limits.tsv',limits)
    fig,axes=plt.subplots(2,2,figsize=(13,9))
    for k,out in enumerate(OUTPUTS):
        for col,idx in enumerate(['S1','ST']):
            ax=axes[k,col];ids=list(range(0,14)) if idx=='S1' else list(range(14,28))
            ns=sorted(series)
            spread=[float(np.quantile(np.ptp(series[n][:,:,k,ids],axis=0),.95)) for n in ns]
            difference=[float(np.quantile(np.abs(series[n][:,:,k,:].mean(axis=0)[:,ids]-smean[:,k,ids]),.95)) for n in ns]
            ax.plot(ns,spread,'o-',color='#762a83',label='p95 phase range')
            ax.plot(ns,difference,'s--',color='#2166ac',label='p95 absolute difference from Sobol')
            ax.axhline(.05,color='#555555',ls=':',lw=.7,label='phase range reference .05')
            ax.set_xscale('log',base=2);ax.set_xticks(ns,[str(n) for n in ns]);ax.minorticks_off()
            ax.set_title(LABELS[out]+'; '+idx);ax.set_xlabel('eFAST N');ax.legend(fontsize=8)
    fig.tight_layout();gp.savefig(fig,dest/'figures/method_resolution_summary')


def final_report(root):
    root=Path(root);p=load(root/'plan.json');dest=root/'comparison';dest.mkdir(exist_ok=True)
    sv,oxygen,sm=sobol_data(root)
    ids=[sm.index(m) for m in METRICS];sobol=sv[1][:,:,:,ids]
    series={513:legacy_values(root)}
    for n in p['resolutions']:
        directory=root/'efast'/('N%d'%n)
        if load(directory/'completion.json')['status']!='complete':raise ValueError('Missing eFAST resolution')
        with np.load(directory/'indices_repeats.npz') as z:
            if z['values'].shape!=(5,201,2,28) or not np.array_equal(z['oxygen'],oxygen):raise ValueError('eFAST axes differ')
            series[n]=z['values']
    rows=[];scores=[];diagnostics=[];ranks=[]
    reference=sobol.mean(axis=0)
    for n,values in sorted(series.items()):
        mean=values.mean(axis=0);spread=np.ptp(values,axis=0)
        previous=None if n in [513,p['resolutions'][0]] else p['resolutions'][p['resolutions'].index(n)-1]
        delta=np.abs(mean-series[previous].mean(axis=0)) if previous else np.full_like(mean,np.nan)
        for k,out in enumerate(OUTPUTS):
            for idx in ['S1','ST']:
                mids=[j for j,m in enumerate(METRICS) if m[0]==idx]
                diff=mean[:,k,mids]-reference[:,k,mids]
                scores.append(dict(efast_N=n,output=out,index=idx,n_cells=2814,mean_absolute_difference=float(np.abs(diff).mean()),
                    rmse=float(np.sqrt((diff**2).mean())),p95_absolute_difference=float(np.quantile(np.abs(diff),.95)),
                    max_absolute_difference=float(np.abs(diff).max()),legacy_reference_only=n==513))
                for band,lo,hi in [('all',0.,5.)]+BANDS:
                    oxy=np.flatnonzero((oxygen>=lo)&(oxygen<=hi))
                    s=spread[:,k,:][np.ix_(oxy,mids)];d=delta[:,k,:][np.ix_(oxy,mids)]
                    finite=np.isfinite(values[:,:,k,:][:,:,mids]).all(axis=0)[oxy]
                    passed=finite&(s<=.05)&(d<=.02)
                    diagnostics.append(dict(efast_N=n,output=out,index=idx,oxygen_band=band,n_cells=int(s.size),
                        repeat_range_p95=float(np.nanquantile(s,.95)),repeat_range_max=float(np.nanmax(s)),
                        resolution_low=previous or '',resolution_delta_max=float(np.nanmax(d)) if previous else '',
                        n_passed=int(passed.sum()) if previous else '',n_not_passed=int((finite&~passed).sum()) if previous else '',
                        n_undefined=int((~finite).sum()),diagnostic_policy='all results included',legacy_reference_only=n==513))
            for oi,o2 in enumerate(oxygen):
                for idx in ['S1','ST']:
                    mids=[j for j,m in enumerate(METRICS) if m[0]==idx]
                    ranks.append(dict(efast_N=n,O2_pct=o2,output=out,index=idx,rank_spearman=float(spearmanr(mean[oi,k,mids],reference[oi,k,mids]).statistic)))
                for j,(idx,pi,pj) in enumerate(METRICS):
                    valid=np.isfinite(values[:,oi,k,j]).all()
                    status='undefined' if not valid else ('historical_reference_only' if n==513 else ('no_adjacent_new_design' if not previous else ('passed' if spread[oi,k,j]<=.05 and delta[oi,k,j]<=.02 else 'not_passed')))
                    rows.append(dict(efast_N=n,sobol_N=131072,output=out,O2_pct=o2,index=idx,parameter=pi,
                        efast_mean=float(mean[oi,k,j]),sobol_mean=float(reference[oi,k,j]),difference=float(mean[oi,k,j]-reference[oi,k,j]),
                        efast_repeat_min=float(values[:,oi,k,j].min()),efast_repeat_max=float(values[:,oi,k,j].max()),
                        sobol_repeat_min=float(sobol[:,oi,k,j].min()),sobol_repeat_max=float(sobol[:,oi,k,j].max()),
                        efast_repeat_range=float(spread[oi,k,j]),resolution_delta=float(delta[oi,k,j]),diagnostic_status=status,included=True))
    table(dest/'method_pointwise_comparison.tsv.gz',rows)
    table(dest/'method_difference_summary.tsv',scores)
    table(dest/'efast_diagnostics.tsv',diagnostics)
    table(dest/'method_rank_comparison.tsv.gz',ranks)
    reference_hyp={tuple(r[k] for k in ['output','index','condition']):r for r in read(root/'sobol/noemi_hypothesis_checks.tsv')}
    hypothesis_comparison=[]
    for n in p['resolutions']:
        for row in read(root/'efast'/('N%d'%n)/'noemi_hypothesis_checks.tsv'):
            ref=reference_hyp[tuple(row[k] for k in ['output','index','condition'])]
            hypothesis_comparison.append(dict(efast_N=n,output=row['output'],index=row['index'],condition=row['condition'],
                sobol_contrast=float(ref['contrast']),efast_contrast=float(row['contrast']),
                sobol_supporting_repeats=int(ref['supporting_repeats']),efast_supporting_repeats=int(row['supporting_repeats']),
                same_contrast_direction=bool(np.sign(float(ref['contrast']))==np.sign(float(row['contrast']))),
                interpretation='descriptive process means, not formal group indices or biological confidence'))
    table(dest/'hypothesis_method_comparison.tsv',hypothesis_comparison)
    comparison_figures(root,dest,series,sobol,oxygen)
    last=max(series)
    interpretation=['# 全域eFAST与Sobol最终比较','',
        'Sobol固定为N131072；新eFAST分辨率为%s，每档5次独立相位、201个固定氧点。没有预设最大N。'%p['resolutions'],
        '新eFAST与Sobol使用相同参数范围/分布、冻结固定氧矩阵和异常恢复；方法一致性检查见validation/method_alignment.json。',
        '旧N513仅作为历史参考；旧元数据未完整保存best_params与fit_config哈希，不能赋予与新运行同等的身份认证。','',
        '## S1与ST差异','']
    for row in [r for r in scores if r['efast_N']==last]:
        interpretation.append('- %s / %s：平均绝对差%.5f，95分位绝对差%.5f，最大绝对差%.5f。'%(LABELS[row['output']],row['index'],row['mean_absolute_difference'],row['p95_absolute_difference'],row['max_absolute_difference']))
    interpretation+=['','## Noemi假设的方法一致程度','']
    for out in OUTPUTS:
        selected=[r for r in hypothesis_comparison if r['efast_N']==last and r['output']==out and r['index']=='ST']
        interpretation.append('- %s：最高新eFAST分辨率与Sobol在三项ST描述性检验中，%d/3项对比方向一致。'%(LABELS[out],sum(r['same_contrast_direction'] for r in selected)))
        for row in selected:
            interpretation.append('  - %s：Sobol对比值%+.5f，eFAST对比值%+.5f；各自支持重复数%d/5、%d/5。'%(row['condition'],row['sobol_contrast'],row['efast_contrast'],row['sobol_supporting_repeats'],row['efast_supporting_repeats']))
    interpretation += ['', '## 科学解释与使用边界','',
        'Noemi假设的Sobol逐项解释见../sobol/interpretation_zh.md；各eFAST分辨率的同口径过程描述和假设检查分别保存在efast/N*/源表中。',
        'S1/ST比较评价方法稳健性；S2仍只来自Sobol。不得把eFAST的ST−S1标为S2。',
        '收敛仅作诊断，所有参数、氧点和重复保留；未通过不能被隐藏为缺失。通过也不等于交互显著或生物学机制被证实。',
        '五次重复范围不是95%置信区间。微弱S2的细微排名不作为机制证据。',
        '过程均值不等于正式组Sobol指数；死亡参数也参与误分离通路。响应面固定背景，不是全域交互方向。',
        '本轮只完成批准的四档；更高N可作为后续独立不可变批次追加，无N=8193硬上限，也未无限自动提交未指定的计算。',
        '邻域保持封存。','', '全部逐点差值、重复范围及未通过项见本目录源表。']
    (dest/'interpretation_zh.md').write_text('\n'.join(interpretation)+'\n')
    english=(
        '# Global variance-based sensitivity analysis\n\n'
        'We analyzed the separate in-vivo model at 201 externally fixed oxygen levels (0-5%, 0.025% increments). '
        'Fourteen active inputs were varied independently over the documented fitting bounds, uniformly in encoded space '
        '(log-uniform on the natural scale for log10-encoded inputs). Outputs were mean ploidy from the normalized dominant '
        'eigenvector and asymptotic net live-cell growth from the leading eigenvalue real part.\n\n'
        'The completed Sobol/Saltelli design used N=131072 and five independent scrambles, with second-order estimation enabled. '
        'It supplied S1, ST and all 91 unordered S2 pairs. Extended FAST used SALib 1.5.2 fast_sampler and its compute_orders '
        'estimator, verified against the public fast.analyze interface, at N='+', '.join(map(str,p['resolutions']))+
        ', M=4 and five independent phase repeats. Identical model and numerical-recovery rules were used. '
        'FAST bootstrap intervals were not used as convergence evidence.\n\n'
        'Operational diagnostics were a five-repeat range <=0.05 and an adjacent-resolution mean change <=0.02. '
        'They were reported per output, parameter/index and oxygen level without excluding entries. '
        'The first new FAST resolution has no adjacent new-design diagnostic. Numerical convergence, nonzero interaction '
        'evidence, conditional response-surface direction and biological interpretation were treated separately. '
        'Endpoint neighborhoods remained archived.\n\n'
        '## Results\n\n'
        'All 47,838 Sobol index-by-output-by-oxygen entries met the operational diagnostics. '
        'Ploidy sensitivities showed stronger low-oxygen effects of death parameters and stronger high-oxygen '
        'effects of baseline missegregation and buffering. This supported parts of the proposed shift, '
        'while the complete missegregation group did not uniformly increase at high oxygen and death '
        'was not uniquely dominant over other processes. Net growth showed a different buffering pattern. '
        'Among all 91 pairs, lam_max with p_misseg had the largest peak S2 for both outputs '
        '(ploidy 0.10495 at 0% oxygen; growth 0.23646 at 0.075% oxygen). '
        'Rankings based on integrated oxygen bands were reported separately from peak rankings.\n\n')
    for row in [r for r in scores if r['efast_N']==last]:
        english+='For %s %s, eFAST N%d versus Sobol N131072 had mean absolute difference %.5f and maximum absolute difference %.5f.\n\n'%(LABELS[row['output']],row['index'],last,row['mean_absolute_difference'],row['max_absolute_difference'])
    (root/'methods_results_en.md').write_text(english)
    table(root/'meeting_update.tsv',[
        dict(item='final_Sobol_N',old='unclear / approximately 4000',updated='131072; 5 scrambles; 201 oxygen points'),
        dict(item='pair_count',old='92',updated='91 unique unordered pairs; 36 focal pairs'),
        dict(item='convergence',old='approximately 92 percent',updated='47838/47838 Sobol entries passed; operational, not confidence level'),
        dict(item='completion_date',old='just finished Oct 29',updated='HPC report Oct 23; collection and independent validation Oct 29'),
        dict(item='method_assignment',old='check whether N or M / source of S2',updated='Sobol S1/ST/S2; eFAST S1/ST, M4 and increasing N'),
        dict(item='neighborhood',old='complementary near-fit analysis',updated='archived; not resumed; old q10/q90 and diagnostics remain available')])
    files=[]
    for sub in ['sobol','comparison']:
        for path in sorted((root/sub).rglob('*')):
            if path.is_file():files.append(dict(path=str(path.relative_to(root)),size=path.stat().st_size,sha256=sha(path)))
    for n in p['resolutions']:
        directory=root/'efast'/('N%d'%n)
        paths=list(directory.iterdir())+list((directory/'figures').rglob('*'))
        for path in sorted(paths):
            if path.is_file():files.append(dict(path=str(path.relative_to(root)),size=path.stat().st_size,sha256=sha(path)))
    save(root/'export_manifest.json',dict(files=files,created_utc=now(),raw_model_blocks_excluded=True))
    save(root/'completion.json',dict(status='complete',scope='global',sobol_N=131072,efast_resolutions=p['resolutions'],
        n_repeats=5,n_oxygen=201,N_hard_limit=None,convergence_policy='diagnostic_only',
        neighborhood_jobs_submitted=0,completed_utc=now()))
