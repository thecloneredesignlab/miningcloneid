args<-commandArgs(TRUE);stopifnot(length(args)==4L)
study<-normalizePath(args[1]);task<-args[2];phase<-args[3];case_id<-args[4];release<-Sys.getenv('RR_RELEASE')
source(file.path(release,'numerics.R'));source(file.path(release,'production.R'))
ctx<-rr_load(study,task);out<-file.path(study,'results',task);stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
started<-proc.time()[3]
rr_save<-function(object,path){d0_rds(object,path);writeLines(d0_hash(path),paste0(path,'.sha256'))}
rr_cached<-function(path){stopifnot(file.exists(path),file.exists(paste0(path,'.sha256')),d0_hash(path)==readLines(paste0(path,'.sha256'),warn=FALSE)[1]);z<-readRDS(path);stopifnot(identical(z$contract,ctx$contract));z}
rr_result<-function(values){d0_write(c(list(execution_valid=TRUE,phase=phase,case_id=case_id,contract=ctx$contract,wall_seconds=proc.time()[3]-started,optimizations=0,synthetic_datasets=0,step05_release=FALSE),values),file.path(out,'result.tsv'))}
rr_reference<-function(case){
 path<-file.path(study,'cache',paste0(case$id,'_reference.rds'))
 if(file.exists(path))return(rr_cached(path))
 op<-rr_operator(ctx,case);native<-rr_native(ctx,case);rr_save(list(contract=ctx$contract,native=native),file.path(out,'native_base.rds'))
 replica<-rr_propagate(ctx,op,case$n0,paste0(case$id,':Euler_reference'),method='Euler')
 old_path<-file.path(ctx$paths$project,case$row$old_native_path);stopifnot(d0_hash(old_path)==case$row$old_native_sha256)
 identity<-rr_native_identity(native,replica,readRDS(old_path));d0_write(identity,file.path(out,'native_identity.tsv'))
 loose<-rr_propagate(ctx,op,case$n0,paste0(case$id,':continuous_loose'))
 tight<-rr_propagate(ctx,op,case$n0,paste0(case$id,':continuous_tight'),rtol=1e-12,atol=1e-14)
 alternate<-rr_propagate(ctx,op,case$n0,paste0(case$id,':continuous_independent'),rtol=1e-11,atol=1e-13,solver='radau')
 checks<-rbind(data.frame(comparison='loose_vs_tight',t(rr_compare(tight,loose))),data.frame(comparison='tight_vs_radau',t(rr_compare(tight,alternate))))
 d0_emit(checks,file.path(out,'solver_convergence.tsv'))
 valid<-isTRUE(identity$pass)&&all(apply(checks[,-1,drop=FALSE],1,function(v)rr_gate_numeric(v)))
 z<-list(contract=ctx$contract,case_id=case$id,params=case$params,n0=case$n0,op=op,native=native,replica=replica,loose=loose,tight=tight,alternate=alternate,identity=identity,solver_checks=checks,reference_gate=valid)
 rr_save(z,file.path(out,'reference_evidence.rds'))
 if(valid)rr_save(z,path)
 z
}
if(phase%in%c('engineering','time')){
 case<-rr_case(ctx,case_id);ref<-rr_reference(case)
 d0_emit(case$projection,file.path(out,'actual_parameters.tsv'))
 d0_write(list(Crowding=ctx$cfg$Crowding,crowding=ctx$cfg$crowding,capacity=ctx$cfg$K,production_DT=ctx$cfg$DT,N_MIN=ctx$cfg$N_MIN,N_MAX=ctx$cfg$N_MAX,N_UNIT=ctx$cfg$N_UNIT,oxygen_pct=case$oxygen,initial_cells=case$initial_cells,reference_scope='Unpassaged numerical fixture; not new experimental data'),file.path(out,'runtime_configuration.tsv'))
 if(phase=='engineering'){
  rr_result(list(reference_gate=ref$reference_gate,native_identity_gap=ref$identity$gap,old_seven_day_gap=ref$identity$old_seven_day_gap,continuous_gate=all(apply(ref$solver_checks[,-1,drop=FALSE],1,rr_gate_numeric)),native_calls=1,reference_calls=4))
 }else{
  stopifnot(ref$reference_gate);prior<-rr_native_run(ref$native,ref$op);prior_comparison<-rr_compare(ref$tight,prior);rows<-list(data.frame(level=0,DT=.05,t(prior_comparison),adjacent_full_tv=NA_real_,adjacent_log_live=NA_real_,adjacent_mean_gap=NA_real_,gate=FALSE));reached<-FALSE
  for(level in 1:10){
   dt<-.05/2^level;native<-rr_native(ctx,case,dt,label=paste0('refinement_',level));current<-rr_native_run(native,ref$op);v<-rr_compare(ref$tight,current);adjacent<-rr_compare(prior,current)
   good<-rr_gate_numeric(v,1e-3,.1)&&rr_gate_numeric(prior_comparison,1e-3,.1)&&rr_gate_numeric(adjacent,1e-3,.1)
   rr_save(list(contract=ctx$contract,native=native,metrics=v,adjacent=adjacent),file.path(out,paste0('native_level_',sprintf('%02d',level),'.rds')))
   rows[[length(rows)+1L]]<-data.frame(level=level,DT=dt,t(v),adjacent_full_tv=adjacent[['max_full_tv']],adjacent_log_live=adjacent[['max_log_live']],adjacent_mean_gap=adjacent[['max_mean_chromosome_gap']],gate=good)
   d0_emit(do.call(rbind,rows),file.path(out,'time_discretization_comparison.tsv'));cat('TIME_CASE',case_id,'level',level,'log_gap',v[['max_log_live']],'full_tv',v[['max_full_tv']],'converged',good,'\n');flush.console()
   if(good){reached<-TRUE;break};prior<-current;prior_comparison<-v
  }
  rr_result(list(reference_gate=ref$reference_gate,time_gate=reached,levels=length(rows)-1L,production_DT_unchanged=.05,production_dt_comparison_is_diagnostic=TRUE))
 }
}else if(phase=='aggregation'){
 case<-rr_case(ctx,case_id);ref<-rr_cached(file.path(study,'cache',paste0(case$id,'_reference.rds')));op<-ref$op;comparisons<-list();moments<-list();defects<-list();observations<-list()
 for(candidate in c('C5_original','C5_uniform','C10_uniform','C20_uniform')){
  part<-rr_partition(op$grid,candidate);coarse_op<-rr_aggregate(op,part);X0<-as.numeric(part$P%*%case$n0);lifted<-as.numeric(part$R%*%X0)
  lifted_D<-rr_propagate(ctx,op,lifted,paste(candidate,'lifted_Euler'),method='Euler');lifted_C<-rr_propagate(ctx,op,lifted,paste(candidate,'lifted_continuous'))
  coarse_D<-rr_propagate(ctx,coarse_op,X0,paste(candidate,'coarse_Euler'),method='Euler');coarse_C<-rr_propagate(ctx,coarse_op,X0,paste(candidate,'coarse_continuous'))
  for(scheme in c('D','C')){
   original<-if(scheme=='D')ref$replica else ref$tight;micro<-if(scheme=='D')lifted_D else lifted_C;coarse<-if(scheme=='D')coarse_D else coarse_C
   for(contrast in c('total','initial_lifting','closure')){
    v<-switch(contrast,total=rr_compare(original,coarse),initial_lifting=rr_compare(original,micro),closure=rr_compare(micro,coarse))
    comparisons[[length(comparisons)+1L]]<-data.frame(case_id=case_id,candidate=candidate,scheme=scheme,contrast=contrast,t(v),state_gate=rr_gate_state(v))
   }
   d<-rr_distribution(original);co<-rr_distribution(coarse);flow_sd<-max(2*(8/199),1/22);flow_grid<-seq(.5,8.5,length.out=200);kary_tv<-flow_tv<-numeric(length(rr_time))
   for(t in seq_along(rr_time)){
    k1<-ctx$backend$ivt_smooth_kary_distribution(op$grid,d$fraction[t,],case$params$sigma_kary);k2<-ctx$backend$ivt_smooth_kary_distribution(op$grid,co$fraction[t,],case$params$sigma_kary)
    f1<-ctx$backend$ivt_predicted_flow_density(flow_grid,op$grid,d$fraction[t,],flow_sd,22);f2<-ctx$backend$ivt_predicted_flow_density(flow_grid,op$grid,co$fraction[t,],flow_sd,22)
    kary_tv[t]<-sum(abs(k1-k2))/2;flow_tv[t]<-sum(abs(f1-f2))*(8/199)/2
   }
   observations[[length(observations)+1L]]<-data.frame(case_id=case_id,candidate=candidate,scheme=scheme,max_kary_tv=max(kary_tv),max_flow_mass_tv=max(flow_tv),observation_gate=max(kary_tv)<=.05&&max(flow_tv)<=.05)
   m<-rr_moments(original);m$candidate<-candidate;m$scheme<-scheme;m$model<-'full_original';moments[[length(moments)+1L]]<-m
   m<-rr_moments(coarse);m$candidate<-candidate;m$scheme<-scheme;m$model<-'coarse_lifted';moments[[length(moments)+1L]]<-m
   for(t in seq_along(rr_time)){
    y<-original$state[t,];n<-length(op$grid);xx<-y[seq_len(n)];yy<-y;yy[seq_len(n)]<-as.numeric(part$R%*%(part$P%*%xx))
    defect<-as.numeric(part$P%*%(rr_rhs(y,op)[seq_len(n)]-rr_rhs(yy,op)[seq_len(n)]));defects[[length(defects)+1L]]<-data.frame(case_id=case_id,candidate=candidate,scheme=scheme,time=rr_time[t],normalized_L1_closure_defect=sum(abs(defect))/sum(xx),crowding_total_preserved=abs(sum(xx)-sum(yy[seq_len(n)]))/sum(xx))
   }
  }
  rr_save(list(contract=ctx$contract,case_id=case_id,candidate=candidate,partition=part,lifted_D=lifted_D,lifted_C=lifted_C,coarse_D=coarse_D,coarse_C=coarse_C),file.path(out,paste0(candidate,'_trajectories.rds')))
  if(candidate=='C5_original'&&case$oxygen==20.5){native<-rr_native(ctx,case,n0=lifted,label='lifted_crosscheck');identity<-rr_native_identity(native,lifted_D);d0_write(identity,file.path(out,'lifted_native_identity.tsv'));stopifnot(identity$pass);rr_save(list(contract=ctx$contract,native=native),file.path(out,'lifted_native.rds'))}
 }
 d0_emit(do.call(rbind,comparisons),file.path(out,'aggregation_only_comparison.tsv'));d0_emit(do.call(rbind,moments),file.path(out,'within_class_moments.tsv'));d0_emit(do.call(rbind,defects),file.path(out,'closure_defect.tsv'));d0_emit(do.call(rbind,observations),file.path(out,'observation_operator_comparison.tsv'))
 rr_result(list(aggregation_executed=TRUE,candidates=4,reference_calls=16,no_truth_dependent_lifting=TRUE))
}else if(phase=='tail'){
 source(file.path(release,'tail.R'),local=TRUE)
}else if(phase=='within'){
 source(file.path(release,'within.R'),local=TRUE)
}else if(phase=='history'){
 stop('History execution requires a reviewed accepted candidate and its isolated adapter')
}else stop('Unregistered scientific phase')
d0_guard(study,out);d0_seal(out)
