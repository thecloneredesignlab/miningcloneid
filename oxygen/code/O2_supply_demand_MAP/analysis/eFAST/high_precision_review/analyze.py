#!/usr/bin/env python3
"""Independent decimal matrix residuals, source tables and review figures."""
import argparse
import collections
from decimal import Decimal,localcontext
import gzip
import json
from pathlib import Path
import numpy as np
from review import sha,read,table,dump,audit_point


def verify(root):
    points=read(root/'point_manifest.tsv');summaries=[];checks=[];vector_rows=[]
    selection=json.loads((root/'selection.json').read_text())
    for row in points:
        path=root/'points'/('point_%04d'%int(row['point_id']));s=audit_point(path)
        report=json.loads((path/'result.json').read_text())
        assert report['identity']==row
        assert report['runtime']['selection_sha256']==sha(root/'selection.json')
        assert report['runtime']['code_sha256']==selection['code_sha256']
        # Recreate the exact binary-double values, not truncated decimal rates.
        matrix=np.loadtxt(gzip.open(path/'matrix.tsv.gz','rt'),delimiter='\t')
        vectors=np.genfromtxt(gzip.open(path/'vectors_double.tsv.gz','rt'),delimiter='\t',names=True)
        states=vectors['N'];np.testing.assert_array_equal(states,np.arange(22,155))
        assert abs(np.max(np.abs(matrix).sum(axis=1))-report['matrix']['norm_inf'])<1e-12
        check=dict(point_id=row['point_id'])
        with localcontext() as ctx:
            ctx.prec=160
            entries=[(i,j,Decimal.from_float(float(matrix[i,j]))) for i,j in zip(*np.nonzero(matrix))]
            for digits in (50,100):
                p=json.loads((path/('precision_%d.json'%digits)).read_text())
                for side in ('right','left'):
                    z=p[side]
                    if z.get('error') or not z.get('converged',False):continue
                    x=[Decimal(a) for a in z['vector_decimal']];lam=Decimal(z['lambda']);mx=[Decimal(0)]*len(x)
                    rowsum=[Decimal(0)]*len(x)
                    for i,j,a in entries:
                        if side=='left':i,j=j,i
                        mx[i]+=a*x[j];rowsum[i]+=abs(a)
                    residue=max(abs(a-lam*b) for a,b in zip(mx,x))
                    scaled=residue/((max(rowsum)+abs(lam))*max(abs(a) for a in x))
                    assert scaled<Decimal(z['relative_tolerance'])*10,(row['point_id'],digits,side,str(scaled))
                    check['independent_%s_residual_%d'%(side,digits)]=str(scaled)
                    upper=max(a/b for a,b in zip(mx,x));bound_gap=abs(upper-lam)/max(rowsum)
                    assert bound_gap<Decimal(z['relative_tolerance'])*10,(row['point_id'],digits,side,'upper bound',str(bound_gap))
                    check['independent_%s_upper_gap_%d'%(side,digits)]=str(bound_gap)
                    if side=='right':
                        ploidy=sum(Decimal(int(n))*a/22 for n,a in zip(states,x))
                        assert abs(ploidy-Decimal(z['ploidy']))<Decimal('1e-45' if digits==50 else '1e-95')
            z=json.loads((path/'precision_100.json').read_text())['right']
            if not z.get('error') and z.get('converged',False):
                x=[Decimal(a) for a in z['vector_decimal']]
                xd=np.array([float(a) for a in x]);old=vectors['normalized']
                s.update(vector_L1_double_vs_100=float(sum(abs(Decimal.from_float(float(a))-b) for a,b in zip(old,x))),
                    upper_boundary_mass_double=float(old[-1]),upper_boundary_mass_100=float(x[-1]),
                    lower_boundary_mass_double=float(old[0]),lower_boundary_mass_100=float(x[0]),
                    max_component_log10_100=float(max(x).log10()),min_component_log10_100=float(min(x).log10()))
                for j,n in enumerate(states):vector_rows.append(dict(point_id=row['point_id'],fit_seed=row['fit_seed'],N=int(n),
                    double_raw_real=float(vectors['raw_real'][j]),double_signed=float(vectors['signed_real'][j]),
                    double_normalized=float(old[j]),precision100_decimal=str(x[j]),precision100_log10=float(x[j].log10())))
        summaries.append(s);checks.append(check)
        if len(summaries)%25==0:print('Independently verified',len(summaries),'/',len(points),flush=True)
    fields=list(dict.fromkeys(k for r in summaries for k in r));table(root/'point_evidence.tsv',[{k:r.get(k,'') for k in fields} for r in summaries])
    fields=list(dict.fromkeys(k for r in checks for k in r));table(root/'independent_residuals.tsv',[{k:r.get(k,'') for k in fields} for r in checks])
    table(root/'full_vector_comparison.tsv',vector_rows)
    groups=[]
    for seed in ('seed138','seed101','seed174'):
        for role in ('abnormal','ordinary_control','oxygen_neighbor'):
            rs=[r for r in summaries if r['fit_seed']==seed and role in r['roles'].split(',')]
            good=[r for r in rs if r['precision_verified']]
            groups.append(dict(fit_seed=seed,role=role,n_points=len(rs),n_precision_verified=len(good),
                n_double_response_changed=sum(r['classification']=='double_response_changed' for r in rs),
                n_reducible=sum(not r['strongly_connected'] for r in rs),
                max_absolute_ploidy_change=max([abs(r['ploidy_change']) for r in good] or [0]),
                max_absolute_growth_change=max([abs(r['growth_change']) for r in good] or [0]),
                n_negative_growth_100=sum(Decimal(r['growth_100'])<0 for r in good),
                n_upper_boundary_mass_ge99pct=sum(r.get('upper_boundary_mass_100',0)>=.99 for r in good),
                max_vector_L1_double_vs_100=max([r.get('vector_L1_double_vs_100',0) for r in good] or [0]),
                min_component_log10_100=min([r.get('min_component_log10_100',0) for r in good] or [0]),
                max_double_normalized_residual=max(r['normalized_relative_residual'] for r in rs)))
    table(root/'group_summary.tsv',groups)
    tails(root)
    dump(root/'independent_validation.json',dict(status='passed',n_points=len(points),n_matrix_residual_checks=sum('_residual_' in key for c in checks for key in c),
        n_independent_leading_root_bound_checks=sum('_upper_gap_' in key for c in checks for key in c),
        matrix_source='17-digit round-trip text parsed to binary double, then Decimal.from_float at 160 decimal digits',
        source_hashes_verified=True,n_states_per_vector=133,n_vector_rows=len(vector_rows),new_FAST_samples=0))


def tails(root):
    points=read(root/'point_manifest.tsv')
    with gzip.open(root.parent/'raw_representative_responses.tsv.gz','rt') as f:
        import csv
        original=list(csv.DictReader(f,delimiter='\t'))
    verified={(r['fit_seed'],r['N'],r['M'],r['replicate'],r['parameter'],r['sample_id'],round(float(r['O2_pct']),10)) for r in points}
    tail=[]
    for seed in ('seed138','seed101','seed174'):
        for phase in range(1,6):
            subset=[r for r in original if r['fit_seed']==seed and r['replicate']==str(phase)]
            extreme=[r for r in subset if (float(r['ploidy'])>=6.5 if seed!='seed174' else float(r['ploidy'])>1.01)]
            count=sum((r['fit_seed'],r['N'],r['M'],r['replicate'],r['parameter'],r['sample_id'],round(float(r['O2_pct']),10)) in verified for r in extreme)
            first=subset[0]
            tail.append(dict(fit_seed=seed,N=first['N'],M=first['M'],replicate=phase,parameter=first['parameter'],O2_pct=first['O2_pct'],
                tail_definition='ploidy>=6.5' if seed!='seed174' else 'ploidy>1.01',n_samples=len(subset),n_tail=len(extreme),
                tail_fraction=len(extreme)/len(subset),n_tail_points_in_this_precision_review=count,
                scope='Descriptive existing sequence counts; only stated covered points were precision reevaluated'))
    table(root/'representative_tail_counts.tsv',tail)


def figures(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    rs=read(root/'point_evidence.tsv');figdir=root/'figures';figdir.mkdir(exist_ok=True)
    colors={'abnormal':'#D55E00','ordinary_control':'#0072B2','oxygen_neighbor':'#999999'}
    fig,axes=plt.subplots(1,3,figsize=(13,4.4))
    for ax,seed in zip(axes,('seed138','seed101','seed174')):
        for role in ('oxygen_neighbor','ordinary_control','abnormal'):
            chosen=[r for r in rs if r['fit_seed']==seed and role in r['roles'].split(',') and r['precision_verified']=='True']
            ax.scatter([float(r['original_ploidy']) for r in chosen],[float(r['ploidy_100']) for r in chosen],s=15,alpha=.7,color=colors[role],label=role.replace('_',' '))
        ax.plot([.8,7.2],[.8,7.2],':',color='0.2',lw=1);ax.set(title=seed,xlabel='Original double ploidy',ylabel='100-digit ploidy',xlim=(.8,7.2),ylim=(.8,7.2));ax.set_aspect('equal')
    h,l=axes[-1].get_legend_handles_labels();fig.legend(h,l,loc='lower center',ncol=3,frameon=False)
    fig.suptitle('Same parameter vectors: original versus high-precision ploidy')
    fig.tight_layout(rect=(0,.08,1,.94));fig.savefig(figdir/'ploidy_before_after.pdf');fig.savefig(figdir/'ploidy_before_after.png',dpi=160);plt.close(fig)
    chosen=[]
    for seed in ('seed138','seed101','seed174'):
        for role in ('abnormal','ordinary_control'):
            candidates=[r for r in rs if r['fit_seed']==seed and role in r['roles'].split(',') and r['precision_verified']=='True']
            if candidates:chosen.append(max(candidates,key=lambda r:abs(float(r['ploidy_change']))))
    fig,axes=plt.subplots(len(chosen),2,figsize=(13,2.6*len(chosen)),squeeze=False)
    vector=read(root/'full_vector_comparison.tsv')
    for i,row in enumerate(chosen):
        selected=[v for v in vector if v['point_id']==row['point_id']];n=np.array([int(v['N']) for v in selected])
        raw=np.array([float(v['double_signed']) for v in selected]);old=np.array([float(v['double_normalized']) for v in selected]);new=np.array([float(v['precision100_decimal']) for v in selected])
        ax,logax=axes[i];ax.plot(n,raw,color='#D55E00',ls=':',lw=1,label='Signed double (raw L2 scale)')
        ax.plot(n,old,color='0.4',lw=1,label='Clipped double (L1 normalized)');ax.plot(n,new,color='#6A3D9A',lw=1,label='100 digits (L1 normalized)')
        ax.set(title='%s | point %s | %s | O₂ %.4g%%'%(row['fit_seed'],row['point_id'],row['roles'],float(row['O2_pct'])),xlabel='Chromosome state N',ylabel='Component')
        logax.plot(n,np.log10(np.maximum(old,1e-300)),color='0.4',label='Double (zeros at -300)')
        logax.plot(n,[float(v['precision100_log10']) for v in selected],color='#6A3D9A',label='100-digit decimal log10')
        logax.set(title='Ploidy %.6g → %.6g; ||Δv||₁=%.3g'%(float(row['original_ploidy']),float(row['ploidy_100']),float(row['vector_L1_double_vs_100'])),xlabel='Chromosome state N',ylabel='log₁₀ component')
    for ax in axes[0]:ax.legend(fontsize=7,frameon=False)
    fig.suptitle('Complete right eigenvectors: selected anomaly and ordinary examples',fontsize=14)
    fig.tight_layout(rect=(0,0,1,.98),h_pad=1.4);fig.savefig(figdir/'complete_vectors.pdf');fig.savefig(figdir/'complete_vectors.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4));relations=read(root/'relationships.tsv');manifest=read(root/'point_manifest.tsv');keyid={r['point_key']:r['point_id'] for r in manifest};byid={r['point_id']:r for r in rs}
    for ax,seed in zip(axes,('seed138','seed101','seed174')):
        row=max((r for r in rs if r['fit_seed']==seed and 'abnormal' in r['roles'] and r['precision_verified']=='True'),key=lambda r:abs(float(r['ploidy_change'])))
        original=next(r for r in manifest if r['point_id']==row['point_id']);keys=[original['point_key']]+[r['related_key'] for r in relations if r['abnormal_key']==original['point_key'] and r['relationship']=='oxygen_neighbor']
        neighbor=sorted([byid[keyid[k]] for k in keys],key=lambda r:float(r['O2_pct']))
        ax.plot([float(r['O2_pct']) for r in neighbor],[float(r['original_ploidy']) for r in neighbor],'o--',color='0.4',label='Original double')
        ax.plot([float(r['O2_pct']) for r in neighbor],[float(r['ploidy_100']) for r in neighbor],'o-',color='#6A3D9A',label='100 digits')
        ax.set(title=seed+' | same parameter vector',xlabel='O₂ (%)',ylabel='Ploidy',ylim=(.8,7.2));ax.legend(frameon=False,fontsize=8)
        ax.ticklabel_format(axis='y',style='plain',useOffset=False)
    fig.suptitle('Neighboring existing O₂ points: common ploidy scale; N and M unchanged')
    fig.tight_layout(rect=(0,0,1,.94));fig.savefig(figdir/'oxygen_neighbors.pdf');fig.savefig(figdir/'oxygen_neighbors.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,8))
    independent={r['point_id']:r for r in read(root/'independent_residuals.tsv')}
    for i,seed in enumerate(('seed138','seed101','seed174')):
        for role in ('oxygen_neighbor','ordinary_control','abnormal'):
            selected=[r for r in rs if r['fit_seed']==seed and role in r['roles'].split(',') and r['precision_verified']=='True']
            axes[0,i].scatter([np.log10(max(float(r['normalized_relative_residual']),1e-300)) for r in selected],
                [float(Decimal(independent[r['point_id']]['independent_right_residual_100']).log10()) if Decimal(independent[r['point_id']]['independent_right_residual_100'])>0 else -300 for r in selected],
                color=colors[role],alpha=.6,s=12)
            axes[1,i].scatter([np.log10(max(abs(float(r['normalization_ploidy_change'])),1e-18)) for r in selected],
                [np.log10(max(abs(float(r['ploidy_change'])),1e-18)) for r in selected],color=colors[role],alpha=.6,s=12)
        axes[0,i].set(title=seed,xlabel='log₁₀ double normalized relative residual',ylabel='log₁₀ saved-vector relative residual')
        axes[1,i].set(xlabel='log₁₀ |clipping-induced ploidy change|',ylabel='log₁₀ |high-precision ploidy change|')
    fig.suptitle('Residuals recomputed from saved 100-digit vectors; normalization diagnostics\nResidual zeros display at -300; ploidy differences below 10⁻¹⁸ display at -18',fontsize=13)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([0],[0],marker='o',ls='',color=colors[r],label=r.replace('_',' ')) for r in colors],loc='lower center',ncol=3,frameon=False)
    fig.tight_layout(rect=(0,.04,1,.94));fig.savefig(figdir/'residuals_normalization.pdf');fig.savefig(figdir/'residuals_normalization.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('verify','figures','tails'));p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();globals()[a.command](a.root.resolve())
