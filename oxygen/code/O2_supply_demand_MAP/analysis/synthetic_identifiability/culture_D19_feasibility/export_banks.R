# Metadata-only export of frozen bank objects. No simulator or optimizer loaded.
args<-commandArgs(TRUE);study<-normalizePath(args[1]);out<-file.path(study,'verification/C0')
source(file.path(study,'repair_release/production_bindings.R'))
suppressPackageStartupMessages(library(digest))
# Finish alias/historical/fit-object barrier tests without replaying any model.
status<-system2('gcc',c('-shared','-fPIC','-O2',shQuote(file.path(study,'closure_release/strict_read_guard.c')),'-ldl','-o',shQuote(file.path(out,'strict_read_guard.so'))));stopifnot(status==0)
private<-file.path(study,'inputs/private/truth.rds');foreign_path<-readLines(file.path(study,'controller/F_B_bank_path.txt'));foreign<-file.path(study,'verification',foreign_path,'bank.rds')
historical<-file.path(dirname(study),'20261116_culture_growth_karyotype_recovery_01/inputs/private/truth_parameters.rds')
aliases<-file.path(out,c('private_alias.rds','foreign_alias.rds'));stopifnot(all(file.symlink(c(private,foreign),aliases)),file.exists(historical))
dir.create(file.path(out,'canary'),showWarnings=FALSE);fixture<-file.path(out,'canary/fit_data.Rds');saveRDS(list(kind='non_scientific_forbidden_filename_canary_no_observations'),fixture)
paths<-c(aliases,historical,fixture)
real_fit<-file.path(dirname(dirname(dirname(study))),'ploidyOxygen/data/fit_objects/fit_data.Rds')
if(file.exists(real_fit))paths<-c(paths,real_fit)
writeLines(paths,file.path(out,'alias_canary_paths.txt'))
args<-c(paste0('LD_PRELOAD=',shQuote(paste(file.path(study,'runtime/guard.so'),file.path(out,'strict_read_guard.so'),sep=':'))),paste0('REC_STUDY_ROOT=',shQuote(study)),'REC_MODE=BUILD','REC_ORIGIN=F_A',paste0('REC_READ_DENIAL_LOG=',shQuote(file.path(out,'alias_read_denials.tsv'))),'/usr/local/bin/o2-hpc-exact-rscript','--vanilla',shQuote(file.path(study,'closure_release/alias_canaries.R')),shQuote(study))
status<-system2('env',args);unlink(aliases);stopifnot(status==0,!any(file.exists(aliases)))
sf_emit(data.frame(status='PASS',kind='late_metadata_only_alias_historical_and_fit_object_barrier_completeness',model_calls=0,existing_canaries=length(paths)),file.path(out,'alias_access_gate.tsv'))
contract<-readRDS(file.path(study,'verification/P0_R1/feasibility_contract.rds'));tab<-contract$design$bounds;free<-tab[tab$parameter!='sigma_growth',,drop=FALSE]
roles<-sf_tab(file.path(study,'inputs/parameter_roles.tsv'));edges<-contract$edges[contract$edges$kind=='explicit_parent',,drop=FALSE]
natural<-list();supply<-list();summary<-list();n<-0L;k<-0L
for(origin in paste0('F_',LETTERS[1:4])){
 path<-readLines(file.path(study,'controller',paste0(origin,'_bank_path.txt')));stopifnot(length(path)==1,startsWith(path,paste0(origin,'_')))
 bank<-readRDS(file.path(study,'verification',path,'bank.rds'));stopifnot(bank$origin==origin,length(bank$members)==200,bank$contract_sha==sf_hash(contract))
 for(i in seq_along(bank$members)){
  m<-bank$members[[i]];stopifnot(length(m$params)==20,setequal(names(m$params),tab$parameter),m$key==sf_hash(m$unit))
  p<-m$params;inactive<-character();if(p$p_misseg==0)inactive<-c(inactive,'k_o_mis');if(p$buffer_beta==0)inactive<-c(inactive,'buffer_n_exp');if(p$buffer_smax==0)inactive<-c(inactive,'buffer_beta','buffer_n_exp')
  values<-vapply(tab$parameter,function(a)as.numeric(p[[a]]),numeric(1));unit<-m$unit[match(tab$parameter,free$parameter)]
  stopifnot(all(is.finite(values)),all(values>=tab$lower),all(values<=tab$upper))
  n<-n+1L;natural[[n]]<-data.frame(origin=origin,member=i,key=m$key,parameter=tab$parameter,value=values,unit_coordinate=unit,role=roles$scientific_role[match(tab$parameter,roles$parameter)],status=ifelse(tab$parameter=='sigma_growth','neutral_placeholder_not_estimated',ifelse(tab$parameter=='sigma_kary','free_observer_not_in_physical_merit',ifelse(tab$parameter%in%inactive,'inactive_at_this_mechanism_null','physical_coordinate'))),lower_hit=values==tab$lower,upper_hit=values==tab$upper)
  selection<-m$pred$selection;totals<-numeric(nrow(edges));indices<-integer(nrow(edges));days<-numeric(nrow(edges))
  for(j in seq_len(nrow(edges))){e<-edges[j,];r<-selection[selection$cohort==e$cohort&selection$segment_id==e$child,,drop=FALSE];stopifnot(nrow(r)==1,is.finite(r$parent_index))
   states<-m$pred$states[[paste(e$cohort,e$parent)]];stopifnot(is.matrix(states),r$parent_index>=1,r$parent_index<=nrow(states));totals[j]<-sum(states[r$parent_index,]);indices[j]<-r$parent_index;days[j]<-r$parent_day
  }
  k<-k+1L;supply[[k]]<-data.frame(origin=origin,member=i,key=m$key,cohort=edges$cohort,parent=edges$parent,child=edges$child,parent_index=indices,parent_day=days,available_cells=totals,required_cells=edges$initial_cells_effective,supply_ratio=totals/edges$initial_cells_effective,required_count_source=edges$initial_source)
 }
 summary[[length(summary)+1L]]<-data.frame(origin=origin,bank_path=path,members=length(bank$members),construction_attempts=bank$construction_attempts,contract_sha=bank$contract_sha)
 rm(bank);gc(FALSE)
}
sf_emit(do.call(rbind,natural),file.path(out,'bank_parameters_natural_and_unit.tsv'));sf_emit(do.call(rbind,supply),file.path(out,'bank_supply_margins.tsv'));sf_emit(do.call(rbind,summary),file.path(out,'bank_manifest.tsv'))
sf_emit(data.frame(status='PASS',model_calls=0,optimizer_calls=0,bank_members=800,parameter_records=16000,parent_boundary_records=800*nrow(edges)),file.path(out,'export_receipt.tsv'));sf_inventory(out)
sf_emit(data.frame(status='PASS',task='C0',model_calls=0),file.path(out,'application_result.tsv'));sf_emit(data.frame(status='PASS',task='C0'),file.path(out,'gate.tsv'))
