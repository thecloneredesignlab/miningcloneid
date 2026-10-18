#!/usr/bin/env python3
"""Reconcile diagnostic tables and summarize evidence; no model evaluations."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import numpy as np
import phase_spectrum as diag
from collect_resolution_nm_results import PILOT, PARAMETERS, sha, read


def gzread(path):
    with gzip.open(path,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))


def run(root):
    parent=root.parent;repo=Path(__file__).resolve().parents[5]
    checked=[];finder_metadata=[]
    for name in ('submission_collection_manifest.tsv','completed_results_manifest.tsv','review_artifact_manifest.tsv'):
        for r in read(parent/name):
            p=repo/r['path'] if name=='review_artifact_manifest.tsv' else parent/r['path']
            if p.name=='.DS_Store':
                finder_metadata.append(dict(path=str(p),matches_historical_hash=p.is_file() and sha(p)==r['sha256'],
                                            scope='Finder folder metadata; excluded from scientific preservation gate'))
                continue
            assert p.is_file() and sha(p)==r['sha256'],str(p)
            checked.append(str(p))
    raw=read(root/'trajectory_manifest.tsv');cases=read(root/'representative_cases.tsv')
    design={r['trajectory_id']:r for r in read(root/'input_design_diagnostics.tsv')}
    metrics=gzread(root/'phase_metrics.tsv.gz')
    lookup={(r['trajectory_id'],round(float(r['O2_pct']),10),r['output']):r for r in metrics}
    result=[];comparisons=[];max_identity_error=0;phase_alignment_error=0
    for case in cases:
        seed,parameter=case['fit_seed'],case['parameter'];o=float(case['O2_pct']);j=int(round(o/.025))
        subset=[r for r in raw if r['fit_seed']==seed and r['parameter']==parameter]
        config=sorted({(int(r['N']),int(r['M'])) for r in subset})
        vs={}
        for n,m in config:
            phase_rows=[]
            for r in sorted([r for r in subset if int(r['N'])==n and int(r['M'])==m],key=lambda r:int(r['replicate'])):
                met=lookup[(r['trajectory_id'],round(o,10),'ploidy')]
                v,l=float(met['variance']),float(met['low_power']);h=v-l
                with np.load(root/'processed_cache'/('trajectory_%03d.npz'%int(r['trajectory_id']))) as z:
                    assert str(z['code_sha256'].item())==sha(Path(diag.__file__))
                    ci=int(np.argmin(np.abs(z['case_oxygen']-o)));assert abs(float(z['case_oxygen'][ci])-o)<1e-10
                    y=z['case_y'][:,ci,0];g=z['case_gap'][:,ci]
                    assert np.max(np.abs(z['metrics'][j,0]-np.array([float(met[k]) for k in diag.METRICS])))<1e-12
                    record={k:r[k] for k in ('trajectory_id','fit_seed','N','M','replicate','parameter')}
                    record.update(O2_pct=o,**{k:float(met[k]) for k in diag.METRICS})
                    record.update(phase_radians=float(z['phase']),fraction_ploidy_ge2=float(np.mean(y>=2)),fraction_ploidy_ge4=float(np.mean(y>=4)),
                        fraction_ploidy_ge8=float(np.mean(y>=8)),fraction_ploidy_ge16=float(np.mean(y>=16)),
                        median_gap=float(np.median(g)),
                        growth_ST=float(lookup[(r['trajectory_id'],round(o,10),'growth')]['ST']),
                        max_abs_nonfocal_pair_correlation=float(design[r['trajectory_id']]['max_abs_nonfocal_pair_correlation']),
                        max_abs_nonfocal_third_moment=float(design[r['trajectory_id']]['max_abs_nonfocal_third_moment']))
                    result.append(record);phase_rows.append(record)
            vs[n,m]=phase_rows
            # Exact ratio decomposition about phase-mean V and L. Two terms can cancel.
            vm=np.mean([r['variance'] for r in phase_rows]);lm=np.mean([r['low_power'] for r in phase_rows]);reference=1-lm/vm
            for r in phase_rows:
                a=-(r['low_power']-lm)/r['variance'];b=lm*(r['variance']-vm)/(r['variance']*vm)
                error=abs((r['ST']-reference)-(a+b));max_identity_error=max(max_identity_error,error)
                r.update(ST_at_ratio_of_phase_mean_powers=float(reference),exact_low_power_term=float(a),exact_variance_term=float(b))
        for a,b in zip(vs[config[0]],vs[config[1]]):
            assert a['replicate']==b['replicate']
            phase_alignment_error=max(phase_alignment_error,abs(a['phase_radians']-b['phase_radians']))
            comparisons.append(dict(fit_seed=seed,parameter=parameter,O2_pct=o,replicate=a['replicate'],
                N_reference=config[0][0],M_reference=config[0][1],N_test=config[1][0],M_test=config[1][1],
                reference_ST=a['ST'],test_ST=b['ST'],ST_difference=b['ST']-a['ST'],
                variance_test_over_reference=b['variance']/a['variance'],
                tail1pct_variance_reference=a['top1pct_variance_fraction'],tail1pct_variance_test=b['top1pct_variance_fraction'],
                phase_radians_difference=b['phase_radians']-a['phase_radians'],
                caution='same phase, different joint input paths; sample IDs are not matched parameter vectors'))
    assert max_identity_error<1e-12 and phase_alignment_error<1e-8
    diag.table(root/'representative_phase_evidence.tsv',result);diag.table(root/'paired_phase_configuration_changes.tsv',comparisons)
    # Screen all 14 cached parameters, including those outside the requested raw-trajectory focus.
    screen=[]
    for seed in diag.FOCUS:
        folder=parent/('N8193_M4' if seed in ('seed25','seed138') else 'N4097_M8')
        for r in read(folder/'repeat_diagnostics.tsv'):
            if r['fit_seed']!=seed or r['output']!='dominant_mean_ploidy' or r['index']!='ST' or int(r['N'])==2049:continue
            if r['oxygen_window']!='overall_0_to_5pct':continue
            screen.append({**r,'raw_spectrum_collected':r['parameter'] in diag.FOCUS[seed]})
    assert len(screen)==112
    diag.table(root/'all_parameter_cached_ST_screen.tsv',screen)
    record=dict(status='passed',n_preserved_historical_manifest_entries=len(checked),
        excluded_finder_metadata=finder_metadata,
        n_all_parameter_screen_rows=len(screen),n_representative_phase_rows=len(result),
        n_paired_configuration_phase_rows=len(comparisons),exact_ST_decomposition_max_error=max_identity_error,
        matched_phase_radians_max_difference=phase_alignment_error,new_model_evaluations=0,new_jobs_submitted=0,
        phase_code_sha256=sha(Path(diag.__file__)),report_code_sha256=sha(Path(__file__)))
    (root/'evidence_validation.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


def response_figures(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    palette=['#0072B2','#D55E00','#E69F00','#009E73','#CC79A7']
    focus={'seed25':'mu_hp','seed138':'p_mis_base','seed101':'p_mis_base','seed174':'O2_crit'}
    manifest=read(root/'trajectory_manifest.tsv');cases=read(root/'representative_cases.tsv');source=[]
    fig,axes=plt.subplots(4,2,figsize=(13,12))
    for i,(seed,parameter) in enumerate(focus.items()):
        case=next(r for r in cases if r['fit_seed']==seed and r['parameter']==parameter);o=float(case['O2_pct'])
        rows=[r for r in manifest if r['fit_seed']==seed and r['parameter']==parameter and r['N']==case['N_test'] and r['M']==case['M_test']]
        for row in rows:
            n,rep=int(row['N']),int(row['replicate']);color=palette[rep-1]
            with np.load(root/'processed_cache'/('trajectory_%03d.npz'%int(row['trajectory_id']))) as z:
                ci=int(np.argmin(np.abs(z['case_oxygen']-o)));y=z['case_y'][:,ci,0];gap=z['case_gap'][:,ci]
                s=np.arange(n)/n
                axes[i,0].plot(s,y,color=color,alpha=.55,lw=.45,rasterized=True)
                axes[i,1].scatter(np.maximum(gap,1e-20),y,color=color,s=2,alpha=.4,rasterized=True,label='Phase %d'%rep)
                for j in range(n):source.append(dict(fit_seed=seed,parameter=parameter,N=n,M=int(row['M']),replicate=rep,
                    O2_pct=o,sample_id=j+1,FAST_sequence_coordinate=float(s[j]),ploidy=float(y[j]),spectral_gap=float(gap[j])))
        axes[i,0].set(title='%s: %s, O₂ %.3g%%'%(seed,parameter,o),xlabel='FAST sequence coordinate s/(2π)',ylabel='Raw ploidy response',xlim=(0,1))
        axes[i,1].set_xscale('log');axes[i,1].axvline(1e-6,color='0.2',ls=':',lw=.8)
        axes[i,1].set(title='Same stored response points; no reevaluation',xlabel='Spectral gap (log scale; zeros displayed at 10⁻²⁰)',ylabel='Raw ploidy response')
    handles,labels=axes[0,1].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=5,frameon=False)
    fig.suptitle('Original FAST sequence responses and spectral gaps',fontsize=14)
    fig.tight_layout(rect=(0,.035,1,.97));fig.savefig(root/'figures/phase_response_sequences.pdf',dpi=180)
    fig.savefig(root/'figures/phase_response_sequences.png',dpi=160);plt.close(fig)
    diag.table(root/'raw_representative_responses.tsv.gz',source)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    root=p.parse_args().root.resolve();run(root);response_figures(root)
