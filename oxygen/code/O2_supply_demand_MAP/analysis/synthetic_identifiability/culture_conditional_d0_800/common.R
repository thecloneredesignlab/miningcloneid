# Conditional expected-observation adapter. Frozen biological functions are read-only.
d0_tab <- function(p) read.delim(p,stringsAsFactors=FALSE,check.names=FALSE)
d0_kv <- function(p) setNames(as.list(d0_tab(p)$value),d0_tab(p)$key)
d0_emit <- function(x,p) {
 dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);tmp<-tempfile('.pending_',dirname(p))
 old<-options(digits=17);on.exit(options(old));write.table(x,tmp,sep='\t',quote=TRUE,row.names=FALSE,na='NA');stopifnot(file.rename(tmp,p))
}
d0_write <- function(x,p) d0_emit(data.frame(key=names(x),value=vapply(x,function(v)if(is.logical(v))tolower(as.character(v)) else as.character(v),character(1))),p)
d0_rds <- function(x,p) {dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);q<-tempfile('.pending_',dirname(p));saveRDS(x,q,compress=TRUE);stopifnot(file.rename(q,p))}
d0_hash <- function(p) unname(digest::digest(file=p,algo='sha256'))
d0_guard <- function(s,out) {
 stopifnot(getRversion()=='4.4.2',Sys.getenv('D0_IMAGE_SHA256')=='dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372')
 stopifnot(!length(intersect(loadedNamespaces(),c('jsonlite','rjson','RJSONIO','rJava','V8'))))
 d0_emit(data.frame(namespace=sort(loadedNamespaces())),file.path(out,'loaded_namespaces.tsv'))
 packages<-intersect(c('Rcpp','Matrix','DEoptim','dplyr','digest','nloptr'),loadedNamespaces())
 d0_emit(data.frame(package=packages,version=vapply(packages,function(p)as.character(utils::packageVersion(p)),character(1)),loaded_path=vapply(packages,function(p)getNamespaceInfo(asNamespace(p),'path'),character(1))),file.path(out,'runtime_packages.tsv'))
 dlls<-lapply(getLoadedDLLs(),function(x)x[['path']]);dlls<-unlist(dlls,use.names=FALSE);dlls<-dlls[file.exists(dlls)]
 d0_emit(data.frame(path=dlls,sha256=vapply(dlls,d0_hash,character(1))),file.path(out,'loaded_native_libraries.tsv'))
}
d0_paths <- function(s) {
 p<-d0_kv(file.path(s,'study.tsv'));root<-p$project
 list(project=root,wp0=file.path(root,'oxygen/results/structural_identifiability/20261109_seed144_wp0'),
      archive=file.path(root,'oxygen/results/structural_identifiability/20261112_culture_v2_seven_point_review_nojson'))
}
d0_private_layout <- function(source_root,model) {
 dir.create(model,recursive=TRUE,showWarnings=FALSE)
 # Model bootstrap resolves ../util/common_semantics.R relative to model.
 # Preserve the entire frozen sibling layout through read-only source links;
 # only model/native-cache files live in the private writable process tree.
 parent<-dirname(model)
 siblings<-list.files(source_root,all.files=TRUE,full.names=TRUE,no..=TRUE)
 siblings<-siblings[basename(siblings)!='model'&basename(siblings)!='.DS_Store']
 for(f in siblings) {
  target<-file.path(parent,basename(f))
  if(!file.exists(target)&&!dir.exists(target))stopifnot(file.symlink(f,target))
 }
 stopifnot(file.exists(file.path(parent,'util/o2_supply_demand_map_common_semantics.R')))
 invisible(model)
}
d0_backend <- function(s,tag) {
 paths<-d0_paths(s);source(file.path(paths$archive,'implementation/helpers.R'),local=.GlobalEnv)
 source_root<-file.path(paths$wp0,'source/code/oxygen/code/O2_supply_demand_MAP')
 # Every R process uses a private byte-identical model tree and native cache.
 model<-file.path(Sys.getenv('D0_RUNTIME'),'processes',paste0(tag,'_',Sys.getpid()),'model')
 d0_private_layout(source_root,model)
 original<-file.path(source_root,'model');files<-list.files(original,recursive=TRUE,all.files=TRUE,full.names=TRUE,no..=TRUE)
 files<-files[!dir.exists(files)&!grepl('/[.]rcpp_cache|[.]DS_Store$',files)]
 proof<-lapply(files,function(f){relative<-substring(f,nchar(original)+2L);to<-file.path(model,relative);dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(f,to,overwrite=FALSE),d0_hash(f)==d0_hash(to));data.frame(path=relative,sha256=d0_hash(to))})
 d0_emit(do.call(rbind,proof),file.path(dirname(model),'model_source_manifest.tsv'))
 Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR=model)
 backend<-new.env(parent=.GlobalEnv)
 backend$commandArgs<-function(trailingOnly=FALSE){x<-base::commandArgs(trailingOnly);if(!trailingOnly)x<-x[!startsWith(x,'--file=')];x}
 source(file.path(source_root,'util/o2_supply_demand_map_fit_invitro_backend.R'),local=backend,chdir=TRUE)
 stopifnot(identical(backend$INVITRO_PASSAGE_IMPLEMENTATION,'v2'))
 list(backend=backend,paths=paths,source_root=source_root,model=model)
}
d0_context <- function(s,tag,engineering=FALSE) {
 c<-d0_backend(s,tag);b<-readRDS(file.path(s,'design/blind_bundle.rds'))
 c$bundle<-b;c$core<-c$backend$build_model_core(cfg=b$cfg)
 c$dataset<-readRDS(file.path(s,'data/expected_observations.rds'))
 c$data_sha<-d0_hash(file.path(s,'data/expected_observations.rds'))
 c$map<-b$map;stopifnot(sum(c$map$table$status=='active')==14L,length(c$map$fixed)==6L)
 c$spec<-o2prof_spec(c$map,c$map$base_params);c$engineering<-engineering;c$study<-s
 c
}
d0_parameters <- function(c,unit) {
 stopifnot(length(unit)==14L,all(is.finite(unit)),all(unit>=0),all(unit<=1))
 o2prof_from_free(c$map,setNames(c$spec$lower+unit*(c$spec$upper-c$spec$lower),c$spec$param_name))
}
d0_unit <- function(c,p) unname((o2prof_to_free(c$map,p)-c$spec$lower)/(c$spec$upper-c$spec$lower))
d0_forward <- function(c,p,task) {
 if(isTRUE(c$engineering)) {
  result<-system2('python3',c('-B',shQuote(file.path(Sys.getenv('D0_RELEASE'),'workflow.py')),'counter',shQuote(c$study),shQuote(task),shQuote(Sys.getenv('SLURM_JOB_ID')),Sys.getpid()),stdout=TRUE,stderr=FALSE)
  if(!is.null(attr(result,'status'))||length(result)!=1L||!grepl('^[0-9]+$',result))stop('Engineering forward reservation failed')
 }
 run<-lapply(c('2N','4N'),function(cohort)c$backend$ivt_run_lineage(c$bundle$adapters[[cohort]],cfg=c$bundle$cfg,run_params=p,model_core=c$core));names(run)<-c('2N','4N')
 run
}
d0_predict <- function(c,runs,p) {
 design<-c$bundle$observations;g<-numeric();k<-list();f<-list();selected<-list()
 for(cohort in c('2N','4N')) {
  run<-runs[[cohort]]
  for(seg in run$segment_results) {
   frac<-as.numeric(seg$selection$selected_frac);stopifnot(all(is.finite(frac)),all(frac>=0),abs(sum(frac)-1)<1e-12)
   rate<-(log(seg$selection$selected_live_cells)-log(seg$sim$Ntot_live_obs[1]))/seg$selection$selected_day
   kp<-c$backend$ivt_smooth_kary_distribution(run$grid_pre,frac,p$sigma_kary)
   for(id in seg$segment$data_ids) {
    o<-design[[id]]
    if(o$growth)g[id]<-rate
    if(o$kary)k[[id]]<-as.numeric(kp)
    if(o$flow){density<-c$backend$ivt_predicted_flow_density(o$flow_grid,run$grid_pre,frac,c$bundle$flow_sd,c$bundle$cfg$N_UNIT);f[[id]]<-as.numeric(density)}
    selected[[id]]<-data.frame(passage_id=id,cohort=cohort,segment_id=seg$segment$segment_id,selected_day=seg$selection$selected_day)
   }
  }
 }
 stopifnot(length(g)==114L,length(k)==12L,length(f)==20L,all(is.finite(g)))
 for(v in k)stopifnot(all(v>=0),abs(sum(v)-1)<1e-12)
 for(id in names(f))stopifnot(all(is.finite(f[[id]])),all(f[[id]]>=0),abs(sum(f[[id]])*design[[id]]$step-1)<1e-12)
 list(growth=g,kary=k,flow=f,selection=do.call(rbind,selected))
}
d0_loss <- function(c,pred,data=c$dataset) {
 if(!is.null(data$kind))stopifnot(data$kind=='deterministic_expected_not_finite_sample')
 stopifnot(identical(sort(names(pred$growth)),sort(names(data$growth))),identical(sort(names(pred$kary)),sort(names(data$kary))),identical(sort(names(pred$flow)),sort(names(data$flow))))
 sigma<-c$bundle$map$fixed$sigma_growth;gr<-pred$growth[names(data$growth)]
 a<-mean(stats::dnorm(data$growth,gr,sigma,log=TRUE))
 b<-mean(vapply(names(data$kary),function(id)sum(data$kary[[id]]*log(pmax(pred$kary[[id]],1e-12))),numeric(1)))
 z<-mean(vapply(names(data$flow),function(id){step<-c$bundle$observations[[id]]$step;mass<-data$flow[[id]]*step;mass<-mass/sum(mass);sum(mass*log(pmax(pred$flow[[id]],1e-12)))},numeric(1)))
 q<--sum(c(a,b,z));stopifnot(is.finite(q))
 # Independent expected-loss difference: KL terms plus Gaussian displacement.
 base_growth<--log(sigma)-log(2*pi)/2
 base_k<-mean(vapply(data$kary,function(v)sum(v*log(pmax(v,1e-12))),numeric(1)))
 base_f<-mean(vapply(names(data$flow),function(id){v<-data$flow[[id]];mass<-v*c$bundle$observations[[id]]$step;mass<-mass/sum(mass);sum(mass*log(pmax(v,1e-12)))},numeric(1)))
 diff_g<-mean((gr-data$growth)^2/(2*sigma^2))
 diff_k<-mean(vapply(names(data$kary),function(id){v<-data$kary[[id]];sum(v*log(pmax(v,1e-12)/pmax(pred$kary[[id]],1e-12)))},numeric(1)))
 diff_f<-mean(vapply(names(data$flow),function(id){v<-data$flow[[id]];mass<-v*c$bundle$observations[[id]]$step;mass<-mass/sum(mass);sum(mass*log(pmax(v,1e-12)/pmax(pred$flow[[id]],1e-12)))},numeric(1)))
 independent<--sum(c(base_growth,base_k,base_f))+sum(c(diff_g,diff_k,diff_f))
 stopifnot(abs(q-independent)<=1e-8)
 list(objective=q,growth=a,kary=b,flow=z,independent=independent,
      max_growth_error=max(abs(gr-data$growth))/sigma,
      rms_growth_error=sqrt(mean(((gr-data$growth)/sigma)^2)),
      max_kary_tv=max(vapply(names(data$kary),function(id)sum(abs(pred$kary[[id]]-data$kary[[id]]))/2,numeric(1))),
      max_flow_tv=max(vapply(names(data$flow),function(id)sum(abs(pred$flow[[id]]-data$flow[[id]]))*c$bundle$observations[[id]]$step/2,numeric(1))))
}
d0_evaluator <- function(c,out,task,deadline,max_calls=Inf,save_native=TRUE,infeasible_infinite=FALSE) {
 dir.create(out,recursive=TRUE,showWarnings=FALSE);state<-new.env(parent=emptyenv());state$n<-0L;state$best<-Inf
 log<-file.path(out,paste0('evaluations_',Sys.getpid(),'.tsv'))
 d0_emit(data.frame(call=integer(),pid=integer(),seconds=numeric(),feasible=logical(),objective=numeric(),error=character()),log)
 fn<-function(u) {
  if(state$n>=max_calls||as.numeric(Sys.time())>=deadline)stop(structure(list(message='Finite callback/time window ended',call=NULL),class=c('d0_budget_stop','error','condition')))
  state$n<-state$n+1L;t<-proc.time()[3];feasible<-TRUE;error<-'';penalty<-1e9;p<-d0_parameters(c,u)
  if(startsWith(task,'global_')&&!file.exists(file.path(out,'initial_counter_COMPLETE'))) {
   ordinal<-system2('python3',c('-B',shQuote(file.path(Sys.getenv('D0_RELEASE'),'workflow.py')),'initial',shQuote(c$study),task),stdout=TRUE)
   stopifnot(is.null(attr(ordinal,'status')),length(ordinal)==1L)
   if(as.integer(ordinal)>0L)d0_rds(u,file.path(out,'initial_inputs',paste0(sprintf('%04d',as.integer(ordinal)),'.rds')))
  }
  z<-tryCatch({runs<-d0_forward(c,p,task);pred<-d0_predict(c,runs,p);loss<-d0_loss(c,pred);list(params=p,unit=u,predictions=pred,loss=loss,runs=if(save_native)runs else NULL)},
      invitro_protocol_infeasible=function(e){feasible<<-FALSE;error<<-conditionMessage(e);penalty<<-c$backend$invitro_protocol_penalty_objective(e);NULL})
  q<-if(is.null(z))if(infeasible_infinite)Inf else penalty else z$loss$objective
  state$last_feasible<-feasible
  if(!is.null(z)&&q<state$best){state$best<-q;z$data_sha256<-c$data_sha;z$task<-task;z$pid<-Sys.getpid();z$call<-state$n
   d0_rds(z,file.path(out,'improvements',paste0(Sys.getpid(),'_',sprintf('%07d',state$n),'.rds')));d0_rds(z,file.path(out,paste0('best_',Sys.getpid(),'.rds')))}
  row<-data.frame(call=state$n,pid=Sys.getpid(),seconds=proc.time()[3]-t,feasible=feasible,objective=q,error=gsub('[\t\r\n]',' ',error))
  write.table(row,log,sep='\t',quote=TRUE,row.names=FALSE,col.names=FALSE,append=TRUE)
  q
 }
 list(fn=fn,state=state,context=c)
}
d0_best <- function(out) {
 paths<-list.files(out,pattern='^best_[0-9]+[.]rds$',full.names=TRUE);stopifnot(length(paths)>0L)
 z<-lapply(paths,readRDS);z[[which.min(vapply(z,function(v)v$loss$objective,numeric(1)))]]
}
d0_seal <- function(out) {
 paths<-list.files(out,recursive=TRUE,full.names=TRUE);paths<-paths[!dir.exists(paths)&!basename(paths)%in%c('payload_manifest.tsv','COMPLETE')]
 d0_emit(data.frame(path=substring(paths,nchar(out)+2L),sha256=vapply(paths,d0_hash,character(1))),file.path(out,'payload_manifest.tsv'))
 writeLines(d0_hash(file.path(out,'payload_manifest.tsv')),file.path(out,'COMPLETE'))
}
d0_cluster <- function(s,task,engineering,deadline,max_calls) {
  cl<-parallel::makePSOCKcluster(7L,outfile=file.path(s,'logs',paste0(task,'_workers.log')),rscript=Sys.getenv('RSCRIPT_BIN'),rscript_args='--vanilla')
 parallel::clusterCall(cl,function(s,task,engineering,deadline,max_calls,release){
  source(file.path(release,'common.R'),local=.GlobalEnv);cc<-d0_context(s,paste0(task,'_worker'),engineering)
  .GlobalEnv$d0_ev<-d0_evaluator(cc,file.path(s,if(engineering)'engineering' else 'fits',task),task,deadline,max_calls)
  d0_guard(s,file.path(Sys.getenv('D0_RUNTIME'),'workers',as.character(Sys.getpid())))
  list(pid=Sys.getpid(),data_sha256=cc$data_sha,free=sum(cc$map$table$status=='active'),known_nuisance=length(cc$map$fixed),nuisance_sha256=digest::digest(cc$map$fixed,algo='sha256'),model=cc$model)
 },s,task,engineering,deadline,max_calls,Sys.getenv('D0_RELEASE')) -> proof
 stopifnot(length(unique(vapply(proof,function(p)p$pid,integer(1))))==7L,
           length(unique(vapply(proof,function(p)p$data_sha256,character(1))))==1L,
           all(vapply(proof,function(p)p$free==14L&&p$known_nuisance==6L,logical(1))),length(unique(vapply(proof,function(p)p$nuisance_sha256,character(1))))==1L)
 d0_emit(do.call(rbind,lapply(proof,as.data.frame)),file.path(s,if(engineering)'engineering' else 'fits',task,'worker_identity.tsv'))
 cl
}
