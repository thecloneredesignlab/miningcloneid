"""Post-closure interpretation and independent table arithmetic; no forwards."""
from pathlib import Path
import math,sys
from common import read,table,sha,MECHANISMS,NUISANCES

def yes(x):return str(x).lower()=='true'
def num(x):return float(x)
def fmt(x):return format(float(x),'.7g')

def main(study):
 s=Path(study).resolve();c0=s/'verification/C0';out=s/'reports'
 assert (s/'controller/FINISHED').exists(),'Finite resident workflow must finish first'
 lock=read(out/'selection_lock.tsv')[0];assessment=read(c0/'assessment_gate.tsv')[0]
 assert assessment['status']=='PASS' and assessment['selection_lock_sha']==sha(out/'selection_lock.tsv')
 final=read(out/'final_gate.tsv')[0];assert final['active_owned_jobs']=='0'
 labels=read(c0/'exact_vector_labels.tsv');best={r['origin']:r['key'] for r in labels if r['label']=='best_found'}
 assert len(best)==6
 quality=read(out/'verified_candidates.tsv');params=read(c0/'all20_parameter_recovery.tsv')
 kernels={r['key']:r for r in read(c0/'sigma_kary_kernel_diagnostics.tsv')}
 roots=read(c0/'root_distribution_comparison.tsv');growth=read(c0/'growth_prediction_comparison.tsv')
 kary=read(c0/'full_karyotype_prediction_comparison.tsv');roles=read(s/'inputs/public/parameter_roles.tsv')
 bounds={r['parameter']:r for r in roles};metrics=[];distributions=[];arithmetic=[]
 for origin,key in best.items():
  q=next(r for r in quality if r['origin']==origin and r['key']==key);assert yes(q['verified'])
  pp=[r for r in params if r['key']==key]
  assert {r['parameter'] for r in pp}==set(MECHANISMS+NUISANCES) and len(pp)==20
  assert sum(r['role']=='unknown_mechanism' for r in pp)==14
  assert sum(r['role']=='unknown_nuisance' for r in pp)==5
  gg=[num(r['normalized_error']) for r in growth if r['key']==key];assert len(gg)==114
  rms=math.sqrt(sum(e*e for e in gg)/114);mg=max(abs(e) for e in gg)
  samples={r['passage_id'] for r in kary if r['key']==key};assert len(samples)==12
  h=[];tv=[];ltv=[]
  for passage in sorted(samples):
   rows=[r for r in kary if r['key']==key and r['passage_id']==passage];assert len(rows)==133
   a=[num(r['observed_probability']) for r in rows];b=[num(r['generating_observed_probability']) for r in rows]
   assert abs(sum(a)-1)<=1e-12 and abs(sum(b)-1)<=1e-12
   h.append(sum((math.sqrt(x/sum(a))-math.sqrt(y))**2 for x,y in zip(a,b))/2)
   tv.append(sum(abs(x/sum(a)-y) for x,y in zip(a,b))/2)
   ltv.append(sum(abs(num(r['latent_probability'])-num(r['generating_latent_probability'])) for r in rows)/2)
  dg=rms*rms/2;dk=sum(h)/12
  errors=[abs(rms-num(q['rms_growth'])),abs(mg-num(q['max_growth'])),abs(max(tv)-num(q['max_kary_tv'])),abs(dg-num(q['D_G'])),abs(dk-num(q['D_K'])),abs(dg+dk-num(q['D']))]
  assert max(errors)<=1e-10,'Independent arithmetic conflict '+origin
  arithmetic.append(dict(origin=origin,status='PASS',maximum_metric_discrepancy=max(errors),growth_entries=114,karyotype_samples=12,probabilities_per_sample=133))
  metrics.append(dict(origin=origin,key=key,D=q['D'],D_G=q['D_G'],D_K=q['D_K'],growth_RMS_over_reference_SD=q['rms_growth'],max_growth_error_over_reference_SD=q['max_growth'],max_sample_karyotype_TV=q['max_kary_tv'],prediction_pass=q['prediction_pass'],fresh_replay_verified=q['verified'],own_local_certificate=q['local_pass']))
  root_tvs={}
  for cohort in ['2N','4N']:
   a=[r for r in roots if r['key']==key and r['cohort']==cohort];b=[r for r in roots if r['key']=='generating_truth' and r['cohort']==cohort]
   assert len(a)==len(b)==133
   lookup={r['N']:num(r['probability']) for r in b}
   root_tvs[cohort]=sum(abs(num(r['probability'])-lookup[r['N']]) for r in a)/2
  distributions.append(dict(origin=origin,key=key,sigma_kary=kernels[key]['sigma_kary'],kernel_max_column_TV=kernels[key]['max_column_TV_from_generating'],kernel_binary64_identical=kernels[key]['binary64_identical_to_generating'],max_observed_sample_TV=max(tv),max_latent_sample_TV=max(ltv),root_2N_TV=root_tvs['2N'],root_4N_TV=root_tvs['4N']))
 table(out/'best_found_gate_metrics.tsv',metrics);table(out/'best_found_distribution_errors.tsv',distributions)
 table(out/'independent_summary_arithmetic.tsv',arithmetic)
 ranges=[]
 for name in MECHANISMS+NUISANCES:
  rows=[r for r in params if r['parameter']==name and best.get(r['origin'])==r['key']]
  assert len(rows)==6
  ranges.append(dict(parameter=name,role=rows[0]['role'],units=rows[0]['units'],generating=rows[0]['generating'],best_found_min=min(num(r['estimate']) for r in rows),best_found_max=max(num(r['estimate']) for r in rows),origins_at_exact_boundary=','.join(r['origin'] for r in rows if yes(r['lower_hit']) or yes(r['upper_hit'])) or 'NONE',scope='Six approximate finite-search endpoints; not a confidence or posterior interval'))
 table(out/'best_found_parameter_ranges.tsv',ranges)
 selected=lock['best_found_origin'];key=lock['best_found_key'];m=next(r for r in metrics if r['origin']==selected);p={r['parameter']:r for r in params if r['key']==key};dist=next(r for r in distributions if r['origin']==selected)
 budget=read(s/'accounting/final_budget.tsv')[0];attempts=read(s/'accounting/candidate_attempts.tsv')
 functions=read(c0/'function_error_summary.tsv')
 lines=['# Interpretation of the closed D19 restart','',
  '**'+lock['status']+'**. All four primary and two reserve search chains completed. Zero of six independent origins met the prediction gate; at least two were required. No N20 release follows.','',
  'This is one conditional deterministic expected-output recovery exercise at the frozen seed144 configuration. It reuses the original 114 growth means and 12 complete 133-state karyotype laws. It does not generate a finite-cell sample, a new noisy observation, a profile or a tumor fit.','',
  '## Prediction recovery and optimization quality','',
  '| Origin | Best verified D | Growth RMS / SD | Maximum growth error / SD | Maximum sample TV | Own local certificate |','|---|---:|---:|---:|---:|---|']
 for r in metrics:lines.append('| '+r['origin']+' | '+' | '.join(fmt(r[n]) for n in ['D','growth_RMS_over_reference_SD','max_growth_error_over_reference_SD','max_sample_karyotype_TV'])+' | '+r['own_local_certificate']+' |')
 lines += ['', 'Registered limits: D<=1e-4, growth RMS/SD<=0.01, maximum growth error/SD<=0.03 and maximum sample TV<=0.005. Each gate representative also requires two fresh-process replays. These are engineering tolerances, not statistical confidence limits.','',
  'The actual lowest point is '+selected+' / `'+key+'`. Its D_G='+fmt(m['D_G'])+' and D_K='+fmt(m['D_K'])+'. Growth RMS passes, while the total distance, maximum growth error and maximum sample TV fail. Its exact vector passed fresh replay but does not own a completed local certificate. A slightly worse checked anchor is kept separately and is not substituted for the lowest point. The finite local probes are not a global-optimum certificate.','',
  'Independent arithmetic from the saved full growth/karyotype tables reproduced all six endpoints within 1e-10. See `independent_summary_arithmetic.tsv`. Source/SIF, interface, native-scale, read-policy and replay acceptance are separate from optimizer/prediction acceptance.','',
  '## All 14 mechanisms and 6 nuisance roles','',
  'All 14 mechanisms, four shared root-distribution coordinates and sigma_kary remained unknown. Sigma_growth was the known native generating control injected without a transform roundtrip. It is not credited as an estimated parameter. All 19 free coordinates retained their declared bounds and transformations. No mechanism or nuisance was fixed because of a sensitivity result or a raw truth-error threshold.','',
  '| Parameter | Generating value | Best '+selected+' | Relative error (%) | Role | Boundary |','|---|---:|---:|---:|---|---|']
 for name in MECHANISMS+NUISANCES:
  r=p[name];b='lower' if yes(r['lower_hit']) else 'upper' if yes(r['upper_hit']) else 'interior'
  relative='NA' if r['relative_error']=='NA' else fmt(100*num(r['relative_error']))
  lines.append('| '+name+' | '+fmt(r['generating'])+' | '+fmt(r['estimate'])+' | '+relative+' | '+r['role']+' | '+b+' |')
 lines += ['', 'The selected point hits the declared lower bounds for gamma_growth and mu_hp; their generating values are also at those bounds. Some other endpoints hit gamma_mu at its generating lower bound. A boundary hit therefore needs context and is not, by itself, evidence for widening the domain or declaring identifiability. No selected mechanism-null inactivity was encountered. Full exact-vector labels, all20 coordinate errors and all endpoint boundary flags remain in the C0 tables.','',
  '## Sigma_kary and initial-distribution trade-offs','',
  'The generating sigma_kary is '+fmt(p['sigma_kary']['generating'])+' chromosomes, in the previously demonstrated local binary64-flat kernel region. In contrast, the six fitted widths span '+fmt(min(num(r['sigma_kary']) for r in distributions))+' to '+fmt(max(num(r['sigma_kary']) for r in distributions))+' chromosomes. None of their kernels is binary64-identical to the generating kernel. Their maximum-column kernel TV distances span '+fmt(min(num(r['kernel_max_column_TV']) for r in distributions))+' to '+fmt(max(num(r['kernel_max_column_TV']) for r in distributions))+'.','',
  'For '+selected+', sigma_kary='+fmt(dist['sigma_kary'])+', kernel maximum-column TV='+fmt(dist['kernel_max_column_TV'])+', maximum observed-sample TV='+fmt(dist['max_observed_sample_TV'])+' and maximum latent-sample TV='+fmt(dist['max_latent_sample_TV'])+'. Its 2N/4N root-law TV distances are '+fmt(dist['root_2N_TV'])+'/'+fmt(dist['root_4N_TV'])+'. The 4N initial SD is '+fmt(p['init_sd_4N']['estimate'])+' versus generating '+fmt(p['init_sd_4N']['generating'])+'.','',
  'Thus the observed width discrepancies cannot be dismissed as different raw widths within the same flat generating kernel. Observer width, unknown initial laws and propagated latent composition show partial trade-offs among approximate endpoints, while the strict observed-distribution gate still fails. These endpoints are not a validated family of equivalent predictions and do not prove structural non-identifiability. Sigma_kary was neither fixed nor judged by a 1% raw-width criterion.','',
  '## Effective response functions','',
  '| Response for '+selected+' | Declared grid | RMS absolute error | Maximum absolute error |','|---|---|---:|---:|']
 for r in functions:
  if r['key']==key and r['grid_scope'] in ['encountered_verified_union','N22_to154_error0_toN']:
   lines.append('| '+r['response']+' | '+r['grid_scope']+' | '+fmt(r['RMS_absolute_error'])+' | '+fmt(r['maximum_absolute_error'])+' |')
 lines += ['', 'Lambda and mu errors have units day^-1; effective missegregation, WGD and survival-modifier errors are dimensionless. These are unweighted diagnostic grids. The encountered mask records saved state fractions at least 1e-8; it is not an exposure distribution, a score threshold or a claim of equal experimental information. The explicit survival grid is N22:154 and integer errors0:N; it is not an event-flux calculation.','',
  'The response plots show comparatively similar survival curves at the displayed resolution, together with differences in growth/death and effective missegregation responses. For example, k_o_mis, mu_hp, gamma_mu and p_mis_base vary while effective induction is linked through mu/(k_o_mis+mu); buffer_smax, buffer_beta and buffer_n_exp change together. These are candidate trade-offs in approximate solutions, not identified combinations or proof of compensation at an accepted prediction optimum. Raw parameter agreement and response-function agreement must be reported separately.','',
  '## Resource use, stopping and scope','',
  'Allocated use: '+fmt(budget['settled_core_hours'])+' of 200 core-hours. Charged full-bundle attempts: '+str(sum(int(r['full_bundle_attempts']) for r in attempts))+'. One of one scientific reserve rounds was used; technical retries: '+budget['technical_retries_used']+' of one. All allocated jobs used xxlarge. There was no application concurrency cap or discretionary elapsed-time cutoff; site/budget scheduler envelopes and finite callback allowances remained in force.','',
  'The finite search/verification package and the single E/F round are complete. The 200-core-hour authority is a ceiling, not a spending target, and unspent retry/handoff/resource allowances do not authorize extra generations or families. DET_OPTIMIZER_UNRESOLVED follows the registered fewer-than-two-origin stop rule. It is not evidence that the full model or all 14 parameters are structurally non-identifiable. No N20, profiles, tumor analysis or new generating regime is released. The unresolved structural and experimental-contract questions keep their prior status.','',
  'The frozen read/source/phase audit is in C1; the full raw scientific archive remains on HPC. The local collection is an explicitly checksum-verified thin subset. `D19_recovery_diagnostics.pdf` contains 10 pages: normalized coordinate errors, raw width versus kernel difference, six growth comparisons, all 12 full karyotypes for the lowest point and response functions. Raw ranges across six optimization endpoints are descriptive and are not posterior samples or confidence intervals.']
 (out/'interpretation.md').write_text('\n'.join(lines)+'\n')
 inputs=['controller/FINISHED','controller/SELECTION_LOCKED','reports/selection_lock.tsv','reports/verified_candidates.tsv','reports/final_gate.tsv','accounting/final_budget.tsv','accounting/candidate_attempts.tsv','verification/C0/assessment_gate.tsv','verification/C0/exact_vector_labels.tsv','verification/C0/all20_parameter_recovery.tsv','verification/C0/sigma_kary_kernel_diagnostics.tsv','verification/C0/root_distribution_comparison.tsv','verification/C0/growth_prediction_comparison.tsv','verification/C0/full_karyotype_prediction_comparison.tsv','verification/C0/function_error_summary.tsv']
 outputs=['interpretation.md','best_found_gate_metrics.tsv','best_found_distribution_errors.tsv','best_found_parameter_ranges.tsv','independent_summary_arithmetic.tsv']
 rows=[dict(role='input',path=rel,bytes=(s/rel).stat().st_size,sha256=sha(s/rel)) for rel in inputs]
 rows += [dict(role='output',path='reports/'+n,bytes=(out/n).stat().st_size,sha256=sha(out/n)) for n in outputs]
 script=Path(__file__).resolve()
 for source in [script,script.with_name('common.py')]:
  rel=str(source.relative_to(s)) if s in source.parents else str(source)
  rows.append(dict(role='metadata_only_source',path=rel,bytes=source.stat().st_size,sha256=sha(source)))
 table(out/'interpretation_manifest.tsv',rows)
 print('INTERPRETATION_AND_INDEPENDENT_ARITHMETIC_PASS',len(best),'origins; zero model calls')

if __name__=='__main__':main(sys.argv[1])
