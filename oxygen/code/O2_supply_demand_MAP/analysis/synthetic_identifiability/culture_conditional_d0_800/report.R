args<-commandArgs(TRUE);stopifnot(length(args)==3L);s<-normalizePath(args[1]);task<-args[2]
source(file.path(Sys.getenv('D0_RELEASE'),'common.R'))
out<-file.path(s,if(args[3]=='report_resume01')'reports/science_resume01'else'reports/science');stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
tables<-list();parameters<-list();functions<-list();jobs<-d0_tab(file.path(s,'accounting/jobs.tsv'));quality<-list();found<-list()
for(i in seq_len(nrow(jobs))) {
 j<-jobs[i,];mode<-j$mode
 if(!grepl('^(screen_|global_|globalcheck_|local_|check_|positive_)',mode))next
 phase<-if(grepl('^(globalcheck_|check_)',mode)||mode=='positive_check')'checks'else'fits'
 path<-file.path(s,phase,j$task);result_path<-file.path(path,'result.tsv');z<-if(file.exists(result_path))d0_kv(result_path)else list()
 logs<-list.files(path,pattern='^evaluations_[0-9]+[.]tsv$',full.names=TRUE)
 callback_rows<-sum(vapply(logs,function(p)nrow(d0_tab(p)),integer(1)))
 replay<-if(file.exists(file.path(path,'candidate_checks.tsv')))nrow(d0_tab(file.path(path,'candidate_checks.tsv')))else 0L
 native_path<-file.path(path,'native_deoptim.rds');native<-if(file.exists(native_path))readRDS(native_path)else NULL
 quality[[j$task]]<-data.frame(task=j$task,origin=sub('^.*_','',mode),slurm_state=j$state,slurm_exit=j$exit_code,release=j$release,actual_core_seconds=j$actual_core_seconds,recorded_callbacks=callback_rows,candidate_replay_calls=replay,total_recorded_forward_calls=callback_rows+replay,native_iterations=if(is.null(native))NA_integer_ else native$optim$iter,termination=if(is.null(z$termination))'See native/check receipt'else z$termination,objective=if(is.null(z$objective))NA_real_ else as.numeric(z$objective),sealed_payload=file.exists(file.path(path,'COMPLETE')))
 if(file.exists(file.path(path,'verified_candidate.rds'))) {
  checked<-readRDS(file.path(path,'verified_candidate.rds'))$loss$objective
  others<-list.files(path,pattern='^best_[0-9]+[.]rds$',full.names=TRUE)
  bestfound<-min(c(checked,vapply(others,function(p)readRDS(p)$loss$objective,numeric(1))))
  found[[j$task]]<-data.frame(task=j$task,checked_score=checked,best_found_score=bestfound,checked_minus_best_found=checked-bestfound,quality_evidence_belongs_to_checked_vector=TRUE)
 }
}
if(length(quality))d0_emit(do.call(rbind,quality),file.path(out,'optimizer_quality.tsv'))
if(length(found))d0_emit(do.call(rbind,found),file.path(out,'best_found_vs_checked.tsv'))
for(origin in c('positive','A','B','C','D','E','F')) {
 mode<-if(origin=='positive')'positive_check'else paste0('check_',origin);matching<-jobs[jobs$mode==mode,,drop=FALSE];if(!nrow(matching))next
 name<-tail(matching$task,1);path<-file.path(s,'checks',name)
 if(!file.exists(file.path(path,'COMPLETE')))next
 z<-d0_kv(file.path(path,'result.tsv'))
 tables[[origin]]<-data.frame(origin=origin,independent=origin!='positive',objective=as.numeric(z$objective),delta=as.numeric(z$objective_delta),growth_rms=as.numeric(z$rms_growth_error),growth_max=as.numeric(z$max_growth_error),kary_tv=as.numeric(z$max_kary_tv),flow_tv=as.numeric(z$max_flow_tv),prediction_pass=z$prediction_recovery,local_stability=z$local_directional_stability)
 p<-d0_tab(file.path(path,'parameter_recovery.tsv'));p$origin<-origin;parameters[[origin]]<-p
 f<-d0_tab(file.path(path,'function_recovery.tsv'));f$origin<-origin;functions[[origin]]<-f
}
if(length(tables)) {
 summary<-do.call(rbind,tables);allparams<-do.call(rbind,parameters);allfunc<-do.call(rbind,functions)
 d0_emit(summary,file.path(out,'prediction_recovery.tsv'));d0_emit(allparams,file.path(out,'parameter_recovery.tsv'));d0_emit(allfunc,file.path(out,'function_recovery.tsv'))
 blind<-summary[summary$independent,,drop=FALSE]
 if(nrow(blind)) {
  selected<-blind$origin[which.min(blind$objective)];p<-allparams[allparams$origin==selected&allparams$estimated,]
  d0_emit(p,file.path(out,'score_selected_14_mechanisms.tsv'))
  pdf(file.path(out,'parameter_recovery.pdf'),width=11.2,height=7.1);par(mar=c(10,5,3,1));barplot(abs(p$scaled_error),names.arg=p$parameter,las=2,col='#3d6da8',ylab='Absolute error / transformed parameter-domain width',main=paste('Known-nuisance deterministic recovery; score-selected origin',selected));abline(h=.01,lty=2,col='#aa4035');dev.off()
  near<-blind[blind$delta<=1e-6&blind$prediction_pass=='true',,drop=FALSE]
  d0_emit(near,file.path(out,'alternative_solutions.tsv'))
 }
} else d0_write(list(recovery='Not released',reason='Engineering or optimizer hold'),file.path(out,'recovery_not_released.tsv'))
d0_guard(s,out);d0_write(list(pass=TRUE,model_forward_calls=0L,noise_draws=0L),file.path(out,'result.tsv'));d0_seal(out)
