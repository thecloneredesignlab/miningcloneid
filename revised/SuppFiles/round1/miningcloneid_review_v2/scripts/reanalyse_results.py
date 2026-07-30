#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, spearmanr


def locate_root(candidate: Path, marker: str, label: str) -> Path:
    """Resolve an extracted archive root from a direct or one-level-nested path."""
    candidate = candidate.expanduser().resolve()
    if (candidate / marker).exists():
        return candidate
    hits = sorted(candidate.glob(f"*/{marker}")) if candidate.is_dir() else []
    if len(hits) == 1:
        return hits[0].parent.resolve()
    raise SystemExit(
        f"Could not locate {label} root under {candidate}. "
        f"Expected marker: {marker}"
    )


def parse_cli() -> argparse.Namespace:
    package_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=(
            "Recompute the manuscript audit tables from extracted top10.zip and "
            "joint_pre.zip directories."
        )
    )
    parser.add_argument(
        "--top10-root",
        type=Path,
        default=Path(os.environ["TOP10_ROOT"]) if os.environ.get("TOP10_ROOT") else None,
        help="Extracted top10 archive root containing top10_index.tsv.",
    )
    parser.add_argument(
        "--joint-pre-root",
        type=Path,
        default=Path(os.environ["JOINT_PRE_ROOT"]) if os.environ.get("JOINT_PRE_ROOT") else None,
        help="Extracted joint_pre archive root containing cross_validation/ and landscape_subcluster/.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=package_root / "analysis",
        help="Output directory; defaults to the package analysis/ directory.",
    )
    args = parser.parse_args()
    if args.top10_root is None or args.joint_pre_root is None:
        parser.error(
            "Provide --top10-root and --joint-pre-root, or set TOP10_ROOT and JOINT_PRE_ROOT."
        )
    return args


ARGS = parse_cli()
TOP = locate_root(ARGS.top10_root, "top10_index.tsv", "top10")
JP = locate_root(ARGS.joint_pre_root, "cross_validation", "joint_pre")
if not (JP / "landscape_subcluster").exists():
    raise SystemExit(f"joint_pre root lacks landscape_subcluster/: {JP}")
OUT = ARGS.out_dir.expanduser().resolve()
OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep='\t')

def metric_map(p: Path) -> dict[str, str]:
    d=read_tsv(p)
    return dict(zip(d.iloc[:,0].astype(str), d.iloc[:,1].astype(str)))

def fnum(x, default=np.nan):
    try:
        v=float(x)
        return v if np.isfinite(v) else default
    except Exception:
        return default

def best_map(seed_dir: Path):
    d=read_tsv(seed_dir/'best_params.tsv')
    return dict(zip(d['parameter'].astype(str), pd.to_numeric(d['value'], errors='coerce')))

def qstats(s):
    s=pd.to_numeric(pd.Series(s),errors='coerce').dropna()
    if not len(s): return {}
    return {'min':float(s.min()),'median':float(s.median()),'max':float(s.max()),'n':int(len(s))}

idx=read_tsv(TOP/'top10_index.tsv')
# normalize extracted paths

def extracted_dir(row):
    group=row.result_group
    if group.startswith('fit_joint'):
        return TOP/group/row.pair/row.seed
    return TOP/group/row.seed

# ---------------- in vitro ----------------
vit_rows=[]
for row in idx[idx.result_group=='fit_invitro_O2_buffering_500seed'].itertuples(index=False):
    sd=extracted_dir(row)
    ls=read_tsv(sd/'invitro_lineage_summary.tsv')
    dc=read_tsv(sd/'invitro_daily_counts.tsv')
    u4=ls[ls.cohort.eq('4N')].drop_duplicates('segment_id').copy()
    start=float(u4.loc[u4.passage_index.idxmin(),'predicted_mean_kary_N'])
    terminal=float(u4.loc[u4.passage_index.idxmax(),'predicted_mean_kary_N'])
    d4=dc[dc.cohort.eq('4N')].copy()
    d4['hyp_frac']=d4['burden_dead_hypoxia']/d4['burden_total']
    d4['cin_frac']=d4['burden_dead_buffer']/d4['burden_total']
    # late 2N karyotype comparison
    k2=ls[(ls.cohort=='2N') & (pd.to_numeric(ls.observed_n_kary, errors='coerce')>0)].copy()
    late=k2[k2.passage_index==k2.passage_index.max()]
    vit_rows.append({
        'rank':int(row.rank),'seed':row.seed,'objective':float(row.objective),
        'fourN_start_pred_mean_N':start,'fourN_terminal_pred_mean_N':terminal,
        'fourN_decline_N':start-terminal,
        'fourN_max_hypoxia_dead_fraction':float(d4.hyp_frac.max()),
        'fourN_max_CIN_dead_fraction':float(d4.cin_frac.max()),
        'late_2N_predicted_mean_N':float(late.predicted_mean_kary_N.iloc[0]),
        'late_2N_observed_mean_min':float(late.observed_mean_kary_N.min()),
        'late_2N_observed_mean_max':float(late.observed_mean_kary_N.max()),
    })
vit=pd.DataFrame(vit_rows).sort_values('rank')
vit.to_csv(OUT/'invitro_top10_robustness.csv',index=False)

# ---------------- in vivo ----------------
def logit(x, eps=1e-4):
    x=np.clip(np.asarray(x,dtype=float),eps,1-eps)
    return np.log(x/(1-x))

vivo_rows=[]
param_rows=[]
bound_rows=[]
prov_rows=[]
for row in idx[idx.result_group=='fit_invivo_O2_buffering_500seed'].itertuples(index=False):
    sd=extracted_dir(row)
    fm=metric_map(sd/'fit_summary.tsv')
    b=read_tsv(sd/'burden_fit.tsv')
    use=b[(b.day>0)&(b.obs_burden>0)].copy()
    res=use.pred_log_burden-use.obs_log_burden
    rmse=float(np.sqrt(np.mean(res**2)))
    per_h=use.assign(sq=res**2).groupby('harvest').sq.mean().pow(0.5)
    tumor_bal=float(np.sqrt(use.assign(sq=res**2).groupby('harvest').sq.mean().mean()))
    tp=read_tsv(sd/'terminal_ploidy_fit.tsv')
    hm=[]
    for h,g in tp.groupby('harvest'):
        pred=g.pred_fraction.to_numpy(float); pred=pred/pred.sum()
        obs=g.obs_count.to_numpy(float); obs=obs/obs.sum()
        N=g.N.to_numpy(float)
        pm=float((N*pred).sum()); om=float((N*obs).sum())
        tv=0.5*float(np.abs(pred-obs).sum())
        w1=float(np.abs(np.cumsum(pred)-np.cumsum(obs)).sum())
        hm.append((h,pm,om,abs(pm-om),tv,w1))
    hmd=pd.DataFrame(hm,columns=['harvest','pred_mean','obs_mean','mean_abs_error','tv','w1'])
    nf=read_tsv(sd/'necrosis_fit.tsv')
    nrows=[]
    for r in nf.itertuples(index=False):
        z=b[(b.harvest==r.harvest)&(b.day==r.day)]
        if z.empty: continue
        z=z.iloc[0]
        pred=float(z.pred_burden_dead_total_volume_mm3/z.pred_burden_volume_mm3)
        resid=float((logit(pred)-logit(float(r.obs_necrosis_fraction)))/0.75)
        nrows.append((r.harvest,pred,float(r.obs_necrosis_fraction),abs(pred-float(r.obs_necrosis_fraction)),resid**2,0.5*resid**2))
    nd=pd.DataFrame(nrows,columns=['harvest','pred','obs','abs_error','objective_no_half','objective_half'])
    vivo_rows.append({
        'rank':int(row.rank),'seed':row.seed,'objective':float(row.objective),
        'objective_data':fnum(fm.get('objective_data')),
        'burden_log_rmse':rmse,'burden_tumor_balanced_log_rmse':tumor_bal,
        'terminal_mean_N_MAE':float(hmd.mean_abs_error.mean()),
        'terminal_distribution_TV_mean':float(hmd.tv.mean()),
        'terminal_distribution_W1_mean':float(hmd.w1.mean()),
        'necrosis_fraction_MAE':float(nd.abs_error.mean()),
        'necrosis_objective_reconstructed_no_half':float(nd.objective_no_half.mean()),
        'necrosis_objective_reconstructed_half':float(nd.objective_half.mean()),
        'necrosis_objective_reported':fnum(fm.get('objective_necrosis')),
        'necrosis_export_predicted_all_missing':bool(read_tsv(sd/'necrosis_fit.tsv').pred_necrosis_fraction.isna().all()),
        'local_attempted':str(fm.get('optimizer_local_attempted','')).upper()=='TRUE',
        'local_accepted':str(fm.get('optimizer_local_accepted','')).upper()=='TRUE',
    })
    bp=best_map(sd)
    for p,v in bp.items(): param_rows.append({'fit':'invivo','seed':row.seed,'rank':int(row.rank),'parameter':p,'value':v})
    ptab=pd.read_csv(sd/'parameter_table_input.csv')
    for pr in ptab.itertuples(index=False):
        if pr.param_symbol in bp and bool(pr.estimate):
            v=float(bp[pr.param_symbol]); lo=float(pr.lower_bound); hi=float(pr.upper_bound)
            tol=max(1e-10,1e-7*max(abs(lo),abs(hi),1))
            status='lower' if abs(v-lo)<=tol else ('upper' if abs(v-hi)<=tol else 'interior')
            bound_rows.append({'fit':'invivo','seed':row.seed,'rank':int(row.rank),'parameter':pr.param_symbol,'value':v,'lower':lo,'upper':hi,'status':status})
    rp=sd/'run_provenance.tsv'
    if rp.exists():
        pdv=read_tsv(rp)
        for key in ['branch','commit','dirty_status']:
            q=pdv[(pdv.section=='git')&(pdv.key==key)]
            prov_rows.append({'fit':'invivo','seed':row.seed,'key':key,'value':q.value.iloc[-1] if len(q) else np.nan})

vivo=pd.DataFrame(vivo_rows).sort_values('rank')
vivo.to_csv(OUT/'invivo_top10_fit_quality.csv',index=False)
pd.DataFrame(param_rows).to_csv(OUT/'independent_parameter_values_long.csv',index=False)
pd.DataFrame(bound_rows).to_csv(OUT/'independent_parameter_boundary_hits.csv',index=False)
pd.DataFrame(prov_rows).to_csv(OUT/'provenance_audit.csv',index=False)

# in vitro bounds / provenance presence
for row in idx[idx.result_group=='fit_invitro_O2_buffering_500seed'].itertuples(index=False):
    sd=extracted_dir(row); bp=best_map(sd); ptab=pd.read_csv(sd/'parameter_table_input.csv')
    for pr in ptab.itertuples(index=False):
        if pr.param_symbol in bp and bool(pr.estimate) and str(getattr(pr,'use_invitro_fit',True)).upper()!='FALSE':
            v=float(bp[pr.param_symbol]); lo=float(pr.lower_bound); hi=float(pr.upper_bound)
            tol=max(1e-10,1e-7*max(abs(lo),abs(hi),1))
            status='lower' if abs(v-lo)<=tol else ('upper' if abs(v-hi)<=tol else 'interior')
            bound_rows.append({'fit':'invitro','seed':row.seed,'rank':int(row.rank),'parameter':pr.param_symbol,'value':v,'lower':lo,'upper':hi,'status':status})
    prov_rows.append({'fit':'invitro','seed':row.seed,'key':'run_provenance_present','value':(sd/'run_provenance.tsv').exists()})
pd.DataFrame(bound_rows).to_csv(OUT/'independent_parameter_boundary_hits.csv',index=False)
pd.DataFrame(prov_rows).to_csv(OUT/'provenance_audit.csv',index=False)

# ---------------- joint ----------------
joint_idx=idx[idx.result_group.str.startswith('fit_joint')].copy()
joint_rows=[]; soft_rows=[]; surv_rows=[]; mis_rows=[]; move_rows=[]; jprov=[]
O2_refs=[0,0.1,0.5,1,2,5]
for row in joint_idx.itertuples(index=False):
    sd=extracted_dir(row); fm=metric_map(sd/'fit_summary.tsv'); comp=read_tsv(sd/'joint_components.tsv')
    cm=dict(zip(comp.component.astype(str),pd.to_numeric(comp.value,errors='coerce')))
    sc=read_tsv(sd/'joint_soft_coupling.tsv')
    sc.insert(0,'rank',int(row.rank)); sc.insert(0,'seed',row.seed); sc.insert(0,'pair',row.pair)
    soft_rows.append(sc)
    pen=float(cm.get('objective_soft_coupling',np.nan))
    joint_rows.append({
        'pair':row.pair,'rank':int(row.rank),'seed':row.seed,'objective':float(row.objective),
        'objective_invivo':float(cm.get('objective_invivo',np.nan)),
        'objective_invitro':float(cm.get('objective_invitro',np.nan)),
        'objective_soft_coupling':pen,'soft_penalty_fraction_of_cap':pen/(14*0.08),
        'local_attempted':str(fm.get('optimizer_local_attempted','')).upper()=='TRUE',
        'local_accepted':str(fm.get('optimizer_local_accepted','')).upper()=='TRUE',
        'deoptim_iterations':fnum(fm.get('optimizer_iter_completed')),
        'deoptim_stop_reason':fm.get('deoptim_stop_reason'),
    })
    vals={r.parameter:{'vivo':float(r.vivo_natural),'vitro':float(r.vitro_natural)} for r in sc.itertuples(index=False)}
    for N in [44,88]:
        for ctx in ['vivo','vitro']:
            smax=vals['buffer_smax'][ctx]; beta=vals['buffer_beta'][ctx]; nexp=vals['buffer_n_exp'][ctx]
            s=smax*math.exp(-beta*(44/N)**nexp)
            surv_rows.append({'pair':row.pair,'rank':int(row.rank),'seed':row.seed,'context':ctx,'N':N,'s_per_copy':s})
    for O2 in O2_refs:
        for N in [44,88]:
            for ctx in ['vivo','vitro']:
                o2c=vals['O2_crit'][ctx]; no=vals['n_O'][ctx]
                h=(o2c**no)/(o2c**no+O2**no) if O2>0 else 1.0
                mu=vals['mu_hp'][ctx]*h*(N/44)**vals['gamma_mu'][ctx]
                p=vals['p_mis_base'][ctx]+vals['p_misseg'][ctx]*(mu/(mu+vals['k_o_mis'][ctx]))
                p=min(max(p,0),1)
                mis_rows.append({'pair':row.pair,'rank':int(row.rank),'seed':row.seed,'context':ctx,'O2':O2,'N':N,'p_mis':p})
    initp=sd/'joint_soft_coupling_initial_values.tsv'
    if initp.exists():
        ini=read_tsv(initp)
        final={r.center_name:float(r.center_transformed) for r in sc.itertuples(index=False)}
        final.update({r.delta_name:float(r.delta_transformed) for r in sc.itertuples(index=False)})
        for r in ini.itertuples(index=False):
            if r.optimizer_name in final:
                move_rows.append({'pair':row.pair,'rank':int(row.rank),'seed':row.seed,'optimizer_name':r.optimizer_name,'initial':float(r.optimizer_value),'final':final[r.optimizer_name],'abs_move':abs(final[r.optimizer_name]-float(r.optimizer_value))})
    rp=sd/'run_provenance.tsv'
    if rp.exists():
        d=read_tsv(rp)
        for key in ['branch','commit','dirty_status']:
            q=d[(d.section=='git')&(d.key==key)]
            jprov.append({'fit':'joint','pair':row.pair,'seed':row.seed,'key':key,'value':q.value.iloc[-1] if len(q) else np.nan})

joint=pd.DataFrame(joint_rows).sort_values(['pair','rank'])
soft=pd.concat(soft_rows,ignore_index=True)
surv=pd.DataFrame(surv_rows)
mis=pd.DataFrame(mis_rows)
move=pd.DataFrame(move_rows)
joint.to_csv(OUT/'joint_selected60_fit_summary.csv',index=False)
soft.to_csv(OUT/'joint_soft_coupling_all60.csv',index=False)
surv.to_csv(OUT/'joint_survival_function_all60.csv',index=False)
mis.to_csv(OUT/'joint_missegregation_function_all60.csv',index=False)
move.to_csv(OUT/'joint_warm_start_movement_all60.csv',index=False)
pd.concat([pd.DataFrame(prov_rows),pd.DataFrame(jprov)],ignore_index=True).to_csv(OUT/'provenance_audit.csv',index=False)

# joint pair summary
pair_summary=joint.groupby('pair').agg(n=('seed','size'),objective_min=('objective','min'),objective_median=('objective','median'),objective_max=('objective','max'),soft_penalty_fraction_median=('soft_penalty_fraction_of_cap','median'),local_accepted_n=('local_accepted','sum')).reset_index()
pair_summary.to_csv(OUT/'joint_pair_summary.csv',index=False)
# soft parameter summary
ss=[]
for p,g in soft.groupby('parameter'):
    ss.append({'parameter':p,'n':len(g),'vivo_median':g.vivo_natural.median(),'vitro_median':g.vitro_natural.median(),'ratio_median':g.ratio_vivo_to_vitro.median(),'ratio_min':g.ratio_vivo_to_vitro.min(),'ratio_max':g.ratio_vivo_to_vitro.max(),'n_vivo_gt_vitro':int((g.vivo_natural>g.vitro_natural).sum()),'n_vivo_lt_vitro':int((g.vivo_natural<g.vitro_natural).sum()),'n_saturation_ge_0p95':int((g.welsch_saturation_fraction>=.95).sum()),'n_vivo_at_union_bound':int((np.isclose(g.vivo_transformed,g.joint_union_lower_transformed,rtol=0,atol=1e-8)|np.isclose(g.vivo_transformed,g.joint_union_upper_transformed,rtol=0,atol=1e-8)).sum()),'n_vitro_at_union_bound':int((np.isclose(g.vitro_transformed,g.joint_union_lower_transformed,rtol=0,atol=1e-8)|np.isclose(g.vitro_transformed,g.joint_union_upper_transformed,rtol=0,atol=1e-8)).sum())})
pd.DataFrame(ss).to_csv(OUT/'joint_soft_coupling_parameter_summary.csv',index=False)
# survival wide derived
sw0=surv.pivot_table(index=['pair','rank','seed'],columns=['context','N'],values='s_per_copy')
sw=pd.DataFrame({
    'vivo_s44':sw0[('vivo',44)], 'vivo_s88':sw0[('vivo',88)],
    'vitro_s44':sw0[('vitro',44)], 'vitro_s88':sw0[('vitro',88)]
}).reset_index()
sw['vivo_ratio_s88_s44']=sw.vivo_s88/sw.vivo_s44; sw['vitro_ratio_s88_s44']=sw.vitro_s88/sw.vitro_s44
sw['vivo_abs_gradient']=sw.vivo_s88-sw.vivo_s44; sw['vitro_abs_gradient']=sw.vitro_s88-sw.vitro_s44
sw.to_csv(OUT/'joint_survival_function_summary_all60.csv',index=False)
# misseg wide ratios
mw=mis.pivot_table(index=['pair','rank','seed','O2','N'],columns='context',values='p_mis').reset_index(); mw['ratio_vivo_vitro']=mw.vivo/mw.vitro
mw.to_csv(OUT/'joint_missegregation_function_comparison_all60.csv',index=False)
# movement by seed
if len(move):
    mv=move.groupby(['pair','rank','seed']).agg(n_dims=('optimizer_name','size'),n_exact=('abs_move',lambda x:int((x<1e-12).sum())),max_abs_move=('abs_move','max'),median_abs_move=('abs_move','median')).reset_index()
    mv['all_exact']=mv.n_exact==mv.n_dims
    mv.to_csv(OUT/'joint_warm_start_movement_by_seed.csv',index=False)

# ---------------- joint_pre / full-500 ----------------
F=JP/'cross_validation/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/monotonicity_classification/dense-grid_monotonicity_classification/tables'
by=read_tsv(F/'fixed_o2_ploidy_monotonicity_by_seed.tsv')
counts=read_tsv(F/'fixed_o2_ploidy_monotonicity_class_counts.tsv')
runargs=read_tsv(F/'fixed_o2_ploidy_monotonicity_run_arguments.tsv')
by.to_csv(OUT/'fixed_o2_by_seed_500.csv',index=False)
counts.to_csv(OUT/'fixed_o2_class_counts.csv',index=False)
runargs.to_csv(OUT/'fixed_o2_run_arguments.csv',index=False)
rel=by.monotonicity_reliability.value_counts(dropna=False).rename_axis('monotonicity_reliability').reset_index(name='n_seed')
rel['fraction']=rel.n_seed/len(by); rel.to_csv(OUT/'fixed_o2_reliability_counts.csv',index=False)
finalc=by.final_interpretation_class.value_counts(dropna=False).rename_axis('final_interpretation_class').reset_index(name='n_seed'); finalc['fraction']=finalc.n_seed/len(by); finalc.to_csv(OUT/'fixed_o2_final_interpretation_counts.csv',index=False)
# stable class by reliability
pd.crosstab(by.curve_class,by.monotonicity_reliability).reset_index().to_csv(OUT/'fixed_o2_class_by_reliability.csv',index=False)
# objective source audit and correlations with actual separate objective table
invpar=pd.read_csv(JP/'landscape_subcluster/SeedParameterTables/invivo_best_params_by_seed.csv')
obj=by[['seed_number','seed_id','objective','objective_source','objective_total','objective_data','objective_burden','objective_ploidy']].merge(invpar[['seed','objective']].rename(columns={'objective':'separate_MAP_objective'}),left_on='seed_number',right_on='seed',how='left')
obj['raw_minus_MAP']=obj.objective-obj.separate_MAP_objective
obj.to_csv(OUT/'fixed_o2_objective_definition_audit.csv',index=False)
corr=[]
for a,b in [('objective','separate_MAP_objective'),('objective','objective_data'),('objective_data','separate_MAP_objective')]:
    x=obj[[a,b]].dropna(); rho,p=spearmanr(x[a],x[b]); corr.append({'variable_a':a,'variable_b':b,'spearman_rho':rho,'p_value':p,'n':len(x)})
pd.DataFrame(corr).to_csv(OUT/'fixed_o2_objective_spearman_correlations.csv',index=False)
# landscape and subclusters
S=JP/'landscape_subcluster/pooled_invivo_invitro/full_data_in_vivo_clustring/Tables'
cluster=pd.read_csv(S/'pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_cluster_summary.csv')
sil=pd.read_csv(S/'pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_silhouette.csv')
subsum=pd.read_csv(S/'pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_subcluster_summary.csv')
subsil=pd.read_csv(S/'pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_subcluster_silhouette.csv')
suball=pd.read_csv(S/'pooled_invivo_invitro_initial_vs_best_tsne_best_clusters_best_subclusters.csv')
cluster.to_csv(OUT/'landscape_primary_cluster_summary.csv',index=False); sil.to_csv(OUT/'landscape_primary_cluster_silhouette.csv',index=False); subsum.to_csv(OUT/'landscape_subcluster_summary.csv',index=False); subsil.to_csv(OUT/'landscape_subcluster_silhouette.csv',index=False)
# invivo cluster/class association
ig=suball[suball.dataset=='invivo'][['seed','primary_cluster_id','subcluster_id','objective']].copy()
ig=ig.merge(by[['seed_number','curve_class','final_interpretation_class','monotonicity_reliability']],left_on='seed',right_on='seed_number',how='left')
ig.to_csv(OUT/'fixed_o2_landscape_cluster_crosswalk.csv',index=False)
ct=pd.crosstab(ig.primary_cluster_id,ig.curve_class)
chi2,p,dof,exp=chi2_contingency(ct)
n=ct.to_numpy().sum(); r,k=ct.shape; cramer=math.sqrt(chi2/(n*min(k-1,r-1)))
pd.DataFrame([{'chi2':chi2,'p_value':p,'dof':dof,'n':n,'cramers_v':cramer,'n_primary_clusters':r,'n_curve_classes':k}]).to_csv(OUT/'fixed_o2_cluster_class_association.csv',index=False)
# warm start reliability crosswalk
warm=[]
for pair in sorted(joint_idx.pair.unique()):
    m=re.search(r'vi_seed(\d+)_([A-Za-z0-9]+)_vt_seed(\d+)',pair)
    if m:
        vi=int(m.group(1)); vitro=int(m.group(3)); q=by[by.seed_number==vi].iloc[0]
        sr=ig[ig.seed==vi].iloc[0]
        warm.append({'pair':pair,'invivo_seed':vi,'invitro_seed':vitro,'primary_cluster_id':sr.primary_cluster_id,'subcluster_id':sr.subcluster_id,'curve_class':q.curve_class,'final_interpretation_class':q.final_interpretation_class,'monotonicity_reliability':q.monotonicity_reliability,'min_spectral_gap':q.min_spectral_gap,'fraction_gap_below_0p005':q.fraction_o2_gap_below_0p005,'fraction_gap_below_0p01':q.fraction_o2_gap_below_0p01,'separate_MAP_objective':float(invpar.loc[invpar.seed==vi,'objective'].iloc[0])})
pd.DataFrame(warm).to_csv(OUT/'joint_warm_start_fixed_o2_crosswalk.csv',index=False)
# AUC inventory
auc_files=[str(p.relative_to(JP)) for p in JP.rglob('*') if p.is_file() and 'auc' in p.name.lower()]
(Path(OUT/'joint_pre_AUC_file_inventory.txt')).write_text('\n'.join(auc_files) if auc_files else 'No file whose basename contains AUC was found in joint_pre.zip.\n')
(Path(OUT/'joint_pre_file_inventory.txt')).write_text('\n'.join(sorted(str(p.relative_to(JP)) for p in JP.rglob('*') if p.is_file()))+'\n')

summary={
 'invitro':{
   'objective':qstats(vit.objective),'fourN_decline_N':qstats(vit.fourN_decline_N),
   'fourN_terminal_pred_mean_N':qstats(vit.fourN_terminal_pred_mean_N),
   'fourN_max_hypoxia_dead_fraction':qstats(vit.fourN_max_hypoxia_dead_fraction),
   'fourN_max_CIN_dead_fraction':qstats(vit.fourN_max_CIN_dead_fraction),
   'late_2N_predicted_mean_N':qstats(vit.late_2N_predicted_mean_N),
   'late_2N_observed_range':[float(vit.late_2N_observed_mean_min.min()),float(vit.late_2N_observed_mean_max.max())]},
 'invivo':{k:qstats(vivo[k]) for k in ['objective','burden_log_rmse','burden_tumor_balanced_log_rmse','terminal_mean_N_MAE','terminal_distribution_TV_mean','terminal_distribution_W1_mean','necrosis_fraction_MAE','necrosis_objective_reconstructed_no_half','necrosis_objective_reconstructed_half']},
 'joint':{
   'objective':qstats(joint.objective),'pair_best_objectives':pair_summary[['pair','objective_min']].to_dict('records'),
   'local_attempted_n':int(joint.local_attempted.sum()),'local_accepted_n':int(joint.local_accepted.sum()),
   'soft_penalty_fraction_of_cap':qstats(joint.soft_penalty_fraction_of_cap),
   'lam_max_ratio':qstats(soft[soft.parameter=='lam_max'].ratio_vivo_to_vitro),
   'p_misseg_ratio':qstats(soft[soft.parameter=='p_misseg'].ratio_vivo_to_vitro),
   'k_o_mis_ratio':qstats(soft[soft.parameter=='k_o_mis'].ratio_vivo_to_vitro),
   'survival':{c:qstats(sw[c]) for c in ['vivo_s44','vivo_s88','vitro_s44','vitro_s88','vivo_ratio_s88_s44','vitro_ratio_s88_s44','vivo_abs_gradient','vitro_abs_gradient']},
   'n_vivo_absolute_survival_higher_at_both':int(((sw.vivo_s44>sw.vitro_s44)&(sw.vivo_s88>sw.vitro_s88)).sum()),
   'n_vivo_ploidy_gradient_stronger':int((sw.vivo_abs_gradient>sw.vitro_abs_gradient).sum()),
   'misseg_ratio_nonzero_refs':qstats(mw[(mw.O2>0)&(mw.N==44)].ratio_vivo_vitro),
   'misseg_ratio_O2_0p1_N44':qstats(mw[(mw.O2==0.1)&(mw.N==44)].ratio_vivo_vitro),
   'misseg_ratio_O2_0p1_N88':qstats(mw[(mw.O2==0.1)&(mw.N==88)].ratio_vivo_vitro),
 },
 'fixed_o2':{
   'n_seed':int(len(by)),'n_o2':int(by.n_o2.iloc[0]),'n_rows_expected':int(len(by)*by.n_o2.iloc[0]),
   'class_counts':counts.to_dict('records'),'reliability_counts':rel.to_dict('records'),'final_interpretation_counts':finalc.to_dict('records'),
   'objective_source_counts':by.objective_source.value_counts().to_dict(),
   'spearman_raw_vs_MAP':float(pd.DataFrame(corr).query("variable_a=='objective' and variable_b=='separate_MAP_objective'").spearman_rho.iloc[0]),
   'cluster_class_cramers_v':cramer,
   'warm_start_reliable_n':int((pd.DataFrame(warm).monotonicity_reliability=='reliable').sum()),
 },
 'landscape':{
   'primary_clusters':cluster.to_dict('records'),'selected_primary_silhouette':sil[sil.selected==True].to_dict('records'),
   'selected_subcluster_silhouette':subsil[subsil.selected==True].to_dict('records'),
   'warm_start_invivo_seeds':[w['invivo_seed'] for w in warm],'shared_invitro_seed':10,
 },
 'necrosis_factor_two_mismatch':bool(np.allclose(vivo.necrosis_objective_reconstructed_no_half,vivo.necrosis_objective_reported,rtol=1e-10,atol=1e-10)),
}
with open(OUT/'analysis_summary.json','w') as fh: json.dump(summary,fh,indent=2)
print(json.dumps(summary,indent=2))
