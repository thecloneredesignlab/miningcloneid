#!/usr/bin/env python3
"""Assemble final figures without overwriting checksum-certified HPC exports."""
import argparse
from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from common import load,save,sha,now,read,table,METRICS,OUTPUTS,ACTIVE
import plot_helpers as gp
from report import sobol_data,LABELS


def strongest_pairs(root,dest):
    values,oxygen,metrics=sobol_data(root);values=values[1];mean=values.mean(axis=0)
    fig,axes=plt.subplots(6,2,figsize=(14,19))
    for k,out in enumerate(OUTPUTS):
        ids=sorted([j for j,m in enumerate(metrics) if m[0]=='S2'],key=lambda j:-float(np.max(mean[:,k,j])))[:6]
        for row,j in enumerate(ids):
            ax=axes[row,k];curve=values[:,:,k,j]
            ax.plot(oxygen,curve.mean(axis=0),color='#54278f',lw=1.1)
            ax.fill_between(oxygen,curve.min(axis=0),curve.max(axis=0),color='#54278f',alpha=.17)
            ax.axhline(0,color='#555555',lw=.5);gp.oxygen_axis(ax)
            ax.set_title(' x '.join(metrics[j][1:]),fontsize=10,pad=22);ax.set_ylabel('S2')
    fig.text(.25,.985,'Ploidy',ha='center');fig.text(.75,.985,'Net live-cell growth rate',ha='center')
    fig.suptitle('Sobol N131072: strongest pairs by peak S2',y=1.005)
    fig.text(.5,.006,'Shading: min-max of five scramble repeats, not a 95% confidence interval. Negative estimates retained.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.035,1,.96]);gp.savefig(fig,dest/'sobol/strongest_91_pair_curves')


def comparison_curves(root,dest):
    p=load(root/'plan.json');ns=p['resolutions'][-2:]
    v,oxygen,metrics=sobol_data(root);ids=[metrics.index(m) for m in METRICS];sv=v[1][:,:,:,ids]
    series={}
    for n in ns:
        with np.load(root/'efast'/('N%d'%n)/'indices_repeats.npz') as z:series[n]=z['values'].copy()
    colors=dict(zip(ns,['#9970ab','#54278f']))
    with PdfPages(dest/'comparison/parameter_comparison_curves.pdf') as pdf:
        for start in range(0,14,7):
            fig,axes=plt.subplots(7,2,figsize=(15,21))
            for row,pi in enumerate(range(start,start+7)):
                for k,out in enumerate(OUTPUTS):
                    ax=axes[row,k]
                    for idx,off,style in [('S1',0,'--'),('ST',14,'-')]:
                        j=pi+off;ref=sv[:,:,k,j]
                        ax.plot(oxygen,ref.mean(axis=0),style,color='#222222',lw=1.25,label='Sobol '+idx)
                        for n in ns:
                            curve=series[n][:,:,k,j]
                            ax.plot(oxygen,curve.mean(axis=0),style,color=colors[n],lw=1,label='eFAST N%d %s'%(n,idx))
                            if idx=='ST':ax.fill_between(oxygen,curve.min(axis=0),curve.max(axis=0),color=colors[n],alpha=.1)
                    gp.oxygen_axis(ax);ax.set_title(ACTIVE[pi]+'; '+LABELS[out],fontsize=10,pad=22)
                    ax.axhline(0,color='#777777',lw=.4)
            handles,labels=axes[0,0].get_legend_handles_labels()
            fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.978),ncol=3,frameon=False,fontsize=9)
            fig.suptitle('S1 dashed / ST solid; highest two new eFAST designs versus Sobol N131072',y=.997)
            fig.text(.5,.008,'ST shading: min-max of five independent phase repeats; not a confidence interval.',ha='center')
            fig.tight_layout(rect=[0,.03,1,.935]);pdf.savefig(fig,bbox_inches='tight');plt.close(fig)


def publish(root,sobol_only=False):
    root=Path(root);dest=root/'publication'
    for name in ['sobol','comparison','efast']:(dest/name).mkdir(parents=True,exist_ok=True)
    sources=[]
    def copy(src,target):
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(str(src),str(target));sources.append(dict(source=str(src.relative_to(root)),source_sha256=sha(src),
            path=str(target.relative_to(root)),sha256=sha(target),action='unchanged copy'))
    for source in (root/'sobol/figures').iterdir():
        if source.suffix in ['.pdf','.png'] and source.stem!='strongest_91_pair_curves':copy(source,dest/'sobol'/source.name)
    strongest_pairs(root,dest)
    for source in (root/'inputs').glob('figure4b_spearman_direction.*'):copy(source,dest/source.name)
    if not sobol_only:
        if load(root/'completion.json')['status']!='complete':raise ValueError('Final comparison not complete')
        for source in (root/'comparison/figures').iterdir():
            if source.suffix in ['.pdf','.png'] and source.stem!='parameter_comparison_curves':copy(source,dest/'comparison'/source.name)
        comparison_curves(root,dest)
        for n in load(root/'plan.json')['resolutions']:
            for source in (root/'efast'/('N%d'%n)/'figures').iterdir():
                if source.suffix in ['.pdf','.png']:copy(source,dest/'efast'/('N%d'%n)/source.name)
    manifest=[]
    for path in sorted(dest.rglob('*')):
        if path.is_file() and path.name!='manifest.json':manifest.append(dict(path=str(path.relative_to(root)),sha256=sha(path),size=path.stat().st_size))
    save(dest/'manifest.json',dict(status='Sobol_only' if sobol_only else 'complete',created_utc=now(),files=manifest,
        unchanged_source_copies=sources,redrawn_figures=['sobol/strongest_91_pair_curves','comparison/parameter_comparison_curves'] if not sobol_only else ['sobol/strongest_91_pair_curves'],
        correction='separate panel titles from Low/Mid/High oxygen labels; separate curve legend from data',
        model_calls=0,estimator_calls=0,original_HPC_exports_overwritten=0))
    (dest/'README.md').write_text(
        '# Figures for the final sensitivity report\n\n'
        'Use this directory for the final figure package. Original HPC exports remain unchanged for SHA256 provenance. '
        'The strongest-pair and method-comparison curve figures were redrawn from the same saved repeated indices '
        'to separate parameter titles from oxygen-band labels and place the comparison legend outside the data panels. '
        'No model or sensitivity estimator was rerun. All other figures are unchanged copies.\n\n'
        'Five-repeat shading means min-max, not a calibrated confidence interval. S1/ST/S2 have no effect direction. '
        'The retained Figure 4B correlation panel supplies direction; C annotations alone do not. '
        'P display groups and primary groups in descriptive process summaries are mapped in ../process_group_membership.tsv. '
        'Existing conditional response surfaces and their sources are in ../response_surfaces.\n')
    print('Publication figures assembled; Sobol-only=',sobol_only)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--sobol-only',action='store_true')
    a=parser.parse_args();publish(a.root,a.sobol_only)
