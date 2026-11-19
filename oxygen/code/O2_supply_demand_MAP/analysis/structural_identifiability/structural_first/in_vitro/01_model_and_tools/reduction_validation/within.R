# Postprocess fixed candidate diagnostics. No trajectories or fits are generated.
metrics<-moments<-defects<-observations<-list()
for(id in d0_tab(file.path(study,'inputs/cases.tsv'))$case_id){
 path<-file.path(study,'results',paste0('aggregation_',id));stopifnot(file.exists(file.path(path,'COMPLETE')))
 metrics[[id]]<-d0_tab(file.path(path,'aggregation_only_comparison.tsv'));moments[[id]]<-d0_tab(file.path(path,'within_class_moments.tsv'));defects[[id]]<-d0_tab(file.path(path,'closure_defect.tsv'));observations[[id]]<-d0_tab(file.path(path,'observation_operator_comparison.tsv'))
}
all<-do.call(rbind,metrics);obs<-do.call(rbind,observations);defect<-do.call(rbind,defects);summary<-list()
for(candidate in c('C5_original','C5_uniform','C10_uniform','C20_uniform')){
 x<-all[all$candidate==candidate&all$contrast=='total',];o<-obs[obs$candidate==candidate,];d<-defect[defect$candidate==candidate,]
 summary[[candidate]]<-data.frame(candidate=candidate,scheme_cases=nrow(x),state_pass=sum(x$state_gate),state_total_gate=all(x$state_gate),observation_pass=sum(o$observation_gate),observation_total_gate=all(o$observation_gate),max_class_tv=max(x$max_class_tv),max_full_tv=max(x$max_full_tv),max_log_live=max(x$max_log_live),max_mean_chromosome_gap=max(x$max_mean_chromosome_gap),max_flux_tolerance_ratio=max(x$max_flux_tolerance_ratio),max_closure_defect_per_day=max(d$normalized_L1_closure_defect))
}
summary<-do.call(rbind,summary);selected<-summary$candidate[which(summary$state_total_gate)][1]
d0_emit(summary,file.path(out,'candidate_acceptance.tsv'));d0_emit(all,file.path(out,'all_paired_comparisons.tsv'));d0_emit(obs,file.path(out,'all_observation_operator_comparisons.tsv'));d0_emit(do.call(rbind,moments),file.path(out,'all_within_class_moments.tsv'));d0_emit(defect,file.path(out,'all_closure_defects.tsv'))
d0_write(list(release_history=!is.na(selected),selected_candidate=if(is.na(selected))'none'else selected,selection_rule='First registered smallest candidate passing every state/flux fixture; not selected by truth proximity',step05_release=FALSE),file.path(study,'reports/candidate_selection.tsv'))
rr_result(list(within_analysis_complete=TRUE,any_candidate_state_gate=any(summary$state_total_gate),selected_candidate=if(is.na(selected))'none'else selected,reference_calls=0,paired_contrasts_are_not_additive=TRUE))
