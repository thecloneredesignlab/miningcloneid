# Domain variants are audits of an explicit extrapolation, not new biological truth.
case<-rr_case(ctx,case_id);ref<-rr_cached(file.path(study,'cache',paste0(case$id,'_reference.rds')))
rows<-list();convergence<-list();quality<-list();outside<-list();baseline_contrasts<-list()
outside_buckets<-function(run,label){
 d<-rr_distribution(run);P<-rr_projection(d$grid,c(21,32,54,76,98,154));fraction<-d$fraction%*%t(P)
 data.frame(case_id=case_id,variant=label,time=rep(run$times,each=7),bucket=rep(c('below_production','class_1','class_2','class_3','class_4','class_5','above_production'),times=length(run$times)),fraction=as.numeric(t(fraction)),absolute_live=as.numeric(t(fraction*d$total)))
}
tail_metrics<-function(run){
 d<-rr_distribution(run);n<-length(run$op$grid);names<-c('direct',rr_counters);integral<-run$state[,n+7+seq_along(names),drop=FALSE];colnames(integral)<-names
 upper<-integral[,'upper_out'];generated<-integral[,'all_daughters'];per_generated<-ifelse(generated>0,upper/generated,0)
 edge<-rowSums(d$fraction[,d$grid>max(d$grid)-22,drop=FALSE])
 c(max_upper_per_live=max(upper/d$total),max_upper_per_generated=max(per_generated),max_upper_edge_fraction=max(edge),max_lower_per_live=max(integral[,'lower_out']/d$total),max_pruning_per_live=max(integral[,'pruning']/d$total))
}
tail_gate<-function(v)all(v[c('max_upper_per_live','max_upper_per_generated','max_upper_edge_fraction')]<=1e-6)
base_metrics<-tail_metrics(ref$tight)
rows[[1]]<-data.frame(case_id=case_id,variant='production_finite_domain',N_min=22,N_max=154,eps_tail=1e-8,t(base_metrics),numeric_reference_gate=TRUE,upper_flux_gate=tail_gate(base_metrics))
outside[[1]]<-outside_buckets(ref$tight,'production_finite_domain')
get_variant<-function(lo,hi,eps,label){
 grid<-seq.int(lo,hi);n0<-numeric(length(grid));n0[match(case$grid,grid)]<-case$n0;stopifnot(sum(n0)==sum(case$n0))
 op<-rr_operator(ctx,case,grid,eps);loose<-rr_propagate(ctx,op,n0,paste0(label,':loose'));tight<-rr_propagate(ctx,op,n0,paste0(label,':tight'),rtol=1e-12,atol=1e-14)
 alternate<-rr_propagate(ctx,op,n0,paste0(label,':independent'),rtol=1e-11,atol=1e-13,solver='radau')
 check<-rr_compare(tight,loose);independent<-rr_compare(tight,alternate);pass<-rr_gate_numeric(check)&&rr_gate_numeric(independent);metrics<-tail_metrics(tight)
 quality[[length(quality)+1L]]<<-data.frame(case_id=case_id,variant=label,comparison='loose_vs_tight',t(check),numeric_reference_gate=rr_gate_numeric(check))
 quality[[length(quality)+1L]]<<-data.frame(case_id=case_id,variant=label,comparison='tight_vs_radau',t(independent),numeric_reference_gate=rr_gate_numeric(independent))
 rows[[length(rows)+1L]]<<-data.frame(case_id=case_id,variant=label,N_min=lo,N_max=hi,eps_tail=eps,t(metrics),numeric_reference_gate=pass,upper_flux_gate=tail_gate(metrics))
 outside[[length(outside)+1L]]<<-outside_buckets(tight,label)
 baseline_contrasts[[length(baseline_contrasts)+1L]]<<-data.frame(case_id=case_id,variant=label,t(rr_compare(ref$tight,tight)))
 rr_save(list(contract=ctx$contract,case_id=case_id,variant=label,loose=loose,tight=tight,alternate=alternate,metrics=metrics,comparison=check,independent=independent),file.path(out,paste0(label,'.rds')))
 d0_emit(do.call(rbind,rows),file.path(out,'boundary_flux_ledger.tsv'));d0_emit(do.call(rbind,quality),file.path(out,'domain_solver_convergence.tsv'))
 d0_emit(do.call(rbind,outside),file.path(out,'absolute_and_normalized_outside_buckets.tsv'))
 d0_emit(do.call(rbind,baseline_contrasts),file.path(out,'production_domain_contrasts.tsv'))
 list(run=tight,metrics=metrics,valid=pass)
}
# Disable pruning on the original domain, keeping every other setting fixed.
unpruned<-get_variant(22,154,0,'original_domain_no_pruning');pruning_comparison<-rr_compare(ref$tight,unpruned$run)
d0_emit(data.frame(case_id=case_id,t(pruning_comparison)),file.path(out,'pruning_audit.tsv'))
lower<-get_variant(1,154,1e-8,'lower_only_1_154')
accepted<-list();numeric_valid<-unpruned$valid&&lower$valid
for(branch in c('upper_only','both_ends')){
 previous<-if(branch=='upper_only')list(run=ref$tight,metrics=base_metrics,valid=TRUE)else lower
 chosen<-NA_integer_
 for(upper in c(308L,616L,1232L)){
  label<-paste0(branch,'_',upper);current<-get_variant(if(branch=='upper_only')22L else 1L,upper,1e-8,label)
  v<-rr_compare(previous$run,current$run);good<-previous$valid&&current$valid&&rr_gate_numeric(v,1e-3,.1)&&tail_gate(previous$metrics)&&tail_gate(current$metrics)
  numeric_valid<-numeric_valid&&current$valid
  convergence[[length(convergence)+1L]]<-data.frame(case_id=case_id,branch=branch,previous_max=max(previous$run$op$grid),current_max=upper,t(v),adjacent_domains_gate=good)
  d0_emit(do.call(rbind,convergence),file.path(out,'domain_convergence.tsv'))
  cat('DOMAIN_CASE',case_id,label,'TV',v[['max_full_tv']],'log_gap',v[['max_log_live']],'gate',good,'\n');flush.console()
  if(good){chosen<-upper;break};previous<-current
 }
 accepted[[branch]]<-chosen
}
rr_result(list(domain_numeric_gate=numeric_valid,upper_only_converged=!is.na(accepted$upper_only),both_ends_converged=!is.na(accepted$both_ends),upper_only_selected_max=accepted$upper_only,both_ends_selected_max=accepted$both_ends,pruning_small_effect=rr_gate_numeric(pruning_comparison,1e-3,.1),production_domain_unchanged=TRUE,domain_variants_are_extrapolation_audits=TRUE))
