# Frozen-output postprocessing only. No production simulations are called here.
args<-commandArgs(TRUE);study<-normalizePath(args[1]);release<-Sys.getenv('RR_RELEASE');source(file.path(study,'inputs/d0_common.R'))
out<-file.path(study,'reports');cases<-d0_tab(file.path(study,'inputs/cases.tsv'));roles<-d0_tab(file.path(study,'inputs/parameter_roles.tsv'))
parameters<-functions<-identities<-time_rows<-list()
for(i in seq_len(nrow(cases))){
 row<-cases[i,];path<-file.path(study,'cache',paste0(row$case_id,'_reference.rds'));stopifnot(d0_hash(path)==readLines(paste0(path,'.sha256'),warn=FALSE)[1]);z<-readRDS(path);stopifnot(z$reference_gate)
 values<-vapply(roles$parameter,function(name)as.numeric(z$params[[name]]),numeric(1));stopifnot(length(values)==20L,all(is.finite(values)))
 parameters[[i]]<-data.frame(case_id=row$case_id,reference=row$parameter_reference,cohort=row$cohort,oxygen_pct=row$oxygen_pct,parameter=roles$parameter,role=roles$scientific_role,value=values,revision_role='Fixed archived value for numerical audit; no parameter estimation')
 op<-z$op;functions[[i]]<-data.frame(case_id=row$case_id,chromosome=op$grid,lambda=op$lambda,mu=op$mu,p_effective=op$probability,survival_per_error=op$survival,t(op$flows))
 if(row$parameter_reference=='induction_zero')stopifnot(z$params$p_misseg==0,max(abs(op$probability-z$params$p_mis_base))<1e-14)
 if(row$parameter_reference=='survival_zero')stopifnot(z$params$buffer_beta==0,max(abs(op$survival-z$params$buffer_smax))<1e-14)
 identities[[i]]<-data.frame(case_id=row$case_id,native_scaled_gap=z$identity$gap,old_seven_day_scaled_gap=z$identity$old_seven_day_gap,max_solver_full_tv=max(z$solver_checks$max_full_tv),max_solver_log_live=max(z$solver_checks$max_log_live),max_solver_mean_gap=max(z$solver_checks$max_mean_chromosome_gap))
 x<-d0_tab(file.path(study,'results',paste0('time_',row$case_id),'time_discretization_comparison.tsv'));x$case_id<-row$case_id;time_rows[[i]]<-x
}
d0_emit(do.call(rbind,parameters),file.path(out,'all_fixed_parameter_and_nuisance_values.tsv'));d0_emit(do.call(rbind,functions),file.path(out,'production_function_values.tsv'));d0_emit(do.call(rbind,identities),file.path(out,'all_reference_acceptance.tsv'));d0_emit(do.call(rbind,time_rows),file.path(out,'all_time_convergence.tsv'))
paired<-d0_tab(file.path(study,'results/within_summary/all_paired_comparisons.tsv'));groups<-split(paired,list(paired$candidate,paired$scheme,paired$contrast),drop=TRUE)
summary<-lapply(groups,function(x)data.frame(candidate=x$candidate[1],scheme=x$scheme[1],contrast=x$contrast[1],fixtures=nrow(x),state_pass=sum(x$state_gate),max_class_tv=max(x$max_class_tv),max_full_tv=max(x$max_full_tv),max_log_live=max(x$max_log_live),max_mean_chromosome_gap=max(x$max_mean_chromosome_gap),max_flux_tolerance_ratio=max(x$max_flux_tolerance_ratio)))
d0_emit(do.call(rbind,summary),file.path(out,'aggregation_contrast_summary.tsv'))
d0_guard(study,file.path(study,'results/report_postprocess'))
