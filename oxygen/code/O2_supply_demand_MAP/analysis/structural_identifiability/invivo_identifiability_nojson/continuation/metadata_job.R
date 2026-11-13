# Allocated, SIF-pinned archived-output validation. Model calls are forbidden.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==3L,args[2]=='metadata_audit')
study<-normalizePath(args[1]);task<-args[3];entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
Sys.setenv(IVC_METADATA_ONLY='true')
source(file.path(entry,'adapter.R'));source(file.path(entry,'helpers.R'))
source(file.path(entry,'metadata_stage.R'))
for(name in c('score','predictions')){
 rejection<-tryCatch({get(name)(NULL,numeric());'NOT_REJECTED'},error=function(error)conditionMessage(error))
 stopifnot(identical(rejection,'METADATA_MODEL_CALL_FORBIDDEN'))
}
stopifnot(as.character(getRversion())=='4.4.2',Sys.getenv('IVI_IMAGE_SHA256')=='dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372',as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))==1L)
dispatch<-file.path(study,'dispatch',task);d<-read_tsv(file.path(dispatch,'manifest.tsv'));stopifnot(all(vapply(file.path(dispatch,d$path),hash_file,character(1))==d$sha256))
spec<-read_kv(file.path(dispatch,'spec.tsv'));stopifnot(spec$mode=='metadata_audit',spec$release_manifest_sha256==hash_file(file.path(entry,'manifest.tsv')),spec$input_manifest_sha256==hash_file(file.path(study,'manifests/input_manifest.tsv')))
verify_inputs(study);collection<-file.path(study,'collections',spec$metadata_collection)
stopifnot(hash_file(file.path(collection,'manifest.tsv'))==spec$metadata_collection_sha256)
cm<-read_tsv(file.path(collection,'manifest.tsv'));stopifnot(all(vapply(file.path(collection,cm$path),hash_file,character(1))==cm$sha256))
out<-file.path(study,'tests',task);stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE);dir.create(file.path(out,'logs'))
view<-file.path(out,'analysis_view');dir.create(view);links<-character();started<-proc.time()[[3]]
cleanup<-function(){for(path in links){stopifnot(nzchar(Sys.readlink(path)));stopifnot(unlink(path,recursive=FALSE)==0L)};links<<-character();unlink(view,recursive=TRUE)}
tryCatch({
 for(name in c('source','manifests','reference','prior','releases','metadata_releases','cache_sources','fits','tests','dispatch','profile_epochs')){
  original<-file.path(study,name);if(!dir.exists(original))next
  target<-file.path(view,name);stopifnot(file.symlink(original,target));links<-c(links,target)
 }
 for(name in c('tables','reports','plots'))dir.create(file.path(view,name))
 original_tables<-list.files(file.path(study,'tables'),full.names=TRUE,all.files=TRUE,no..=TRUE)
 if(length(original_tables))stopifnot(all(file.copy(original_tables,file.path(view,'tables'),recursive=TRUE)))
 initial<-list.files(file.path(view,'tables'),recursive=TRUE,full.names=TRUE);initial<-initial[!file.info(initial)$isdir]
 initial_hash<-setNames(vapply(initial,hash_file,character(1)),substring(initial,nchar(view)+2L))
 stopifnot(file.copy(file.path(collection,'jobs.tsv'),file.path(view,'tables/jobs.tsv'),overwrite=TRUE),file.copy(file.path(collection,'raw_accounting.txt'),file.path(view,'tables/raw_accounting.txt'),overwrite=TRUE),file.copy(file.path(collection,'origin.tsv'),file.path(view,'tables/analysis_snapshot_origin.tsv'),overwrite=TRUE))
 jobs<-read.delim(file.path(view,'tables/jobs.tsv'),colClasses='character')
 wrapper<-Sys.getenv('RSCRIPT_BIN');stopifnot(nzchar(wrapper));sequence<-0L
 call<-function(script,...){
  stopifnot(as.numeric(Sys.time())<as.numeric(Sys.getenv('IVI_ALLOCATION_START'))+spec$wall_seconds-spec$margin_seconds)
  sequence<<-sequence+1L;log<-file.path(out,'logs',sprintf('%03d_%s.txt',sequence,script))
  status<-system2(wrapper,c('--vanilla',shQuote(file.path(entry,script)),shQuote(view),vapply(list(...),shQuote,character(1))),stdout=log,stderr=log)
  if(status!=0L)stop('Archived metadata child failed: ',script,'; see ',basename(log))
 }
 call('validate_audit_math.R');call('audit_predictions.R')
 eligible<-jobs$task[jobs$state=='COMPLETED'&jobs$exit_code=='0:0'&jobs$kind!='engineering'&jobs$mode%in%c('local','global','directions','profile_fresh')]
 # Batching is command-line transport, never an optimizer/attempt ceiling.
 if(length(eligible))for(chunk in split(eligible,ceiling(seq_along(eligible)/100L)))do.call(call,c(list('audit_records.R'),as.list(chunk)))
 if(any(jobs$mode=='global'&jobs$state=='COMPLETED'&jobs$exit_code=='0:0'))call('audit_populations.R')
 call('gates.R','engineering',read.delim(file.path(view,'tables/execution_release.tsv'),stringsAsFactors=FALSE)$value[1])
 baseline_review<-function(){
 call('gates.R','baseline')
 first_policy<-if(file.exists(file.path(view,'tables/baseline_policy_gates.tsv')))readBin(file.path(view,'tables/baseline_policy_gates.tsv'),'raw',n=file.info(file.path(view,'tables/baseline_policy_gates.tsv'))$size) else raw()
 first_chains<-if(file.exists(file.path(view,'tables/baseline_chain_summary.tsv')))readBin(file.path(view,'tables/baseline_chain_summary.tsv'),'raw',n=file.info(file.path(view,'tables/baseline_chain_summary.tsv'))$size) else raw()
 call('gates.R','baseline');usage<-read_kv(file.path(view,'tables/feasible_cache_usage.tsv'))
 if(length(first_policy))stopifnot(file.info(file.path(view,'tables/baseline_policy_gates.tsv'))$size==length(first_policy),identical(first_policy,readBin(file.path(view,'tables/baseline_policy_gates.tsv'),'raw',n=length(first_policy))))
 if(length(first_chains))stopifnot(file.info(file.path(view,'tables/baseline_chain_summary.tsv'))$size==length(first_chains),identical(first_chains,readBin(file.path(view,'tables/baseline_chain_summary.tsv'),'raw',n=length(first_chains))))
 stopifnot(usage$independent_audit_misses==0L)
 write_kv(list(pass=TRUE,identical_policy_and_recovery_tables=TRUE,repeat_cache_misses=0L,repeat_cache_hits=usage$cache_hits,model_evaluations=0L),file.path(view,'tables/cache_repeat_validation.tsv'))
 }
 epoch_exists<-file.exists(file.path(view,'tables/profile_epoch.tsv'))
 baseline_path<-file.path(view,'tables/baseline_gate.tsv');baseline_pass<-file.exists(baseline_path)&&isTRUE(read_kv(baseline_path)$pass)
 profile_path<-file.path(view,'tables/profile_gate.tsv');prior_invalid<-file.exists(profile_path)&&isTRUE(read_kv(profile_path)$free_baseline_invalidated)
 if(metadata_profile_phase(epoch_exists,baseline_pass,prior_invalid)){
  # The accepted reference belongs to its frozen baseline epoch. Preserve raw
  # small negative profile deltas; only the declared material counterexample
  # veto returns to baseline repair. Never re-invalidate a repaired baseline by
  # reviewing an already-invalidated old epoch.
  call('profile_gate.R','review')
  if(!isTRUE(read_kv(baseline_path)$pass))baseline_review()
 }else baseline_review()
 call('strategy_review.R')
 call('report.R');call('plot_results.R');call('export_points.R')
 snapshot<-file.path(out,'metadata_snapshot');dir.create(snapshot)
 for(name in c('tables','reports','plots'))stopifnot(file.copy(file.path(view,name),snapshot,recursive=TRUE))
 excluded<-c('jobs.tsv','raw_accounting.txt','execution_release.tsv','reviewed_failures.tsv','cache_source_manifest.tsv','metadata_release.tsv','metadata_execution_release.tsv','collection_files.txt','profile_epoch.tsv','search_batches.tsv','informed_roots.tsv','local_mirror_mode.tsv','local_replica_eviction.tsv')
 published<-list()
 for(name in c('tables','reports','plots')){
  paths<-list.files(file.path(view,name),recursive=TRUE,full.names=TRUE);paths<-paths[!file.info(paths)$isdir]
  for(path in paths){
   relative<-substring(path,nchar(view)+2L);if(name=='tables'&&basename(path)%in%excluded)next
   value<-hash_file(path);if(relative%in%names(initial_hash)&&identical(initial_hash[[relative]],value))next
   target<-file.path(study,relative);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE);temporary<-tempfile('.metadata_publish_',dirname(target));stopifnot(file.copy(path,temporary),file.rename(temporary,target))
   published[[length(published)+1L]]<-data.frame(path=relative,sha256=value)
  }
 }
 if(length(published))write_tsv(do.call(rbind,published),file.path(out,'published_files.tsv'))
 write_kv(list(source_collection=spec$metadata_collection,source_collection_sha256=spec$metadata_collection_sha256,mode='archived_output_only',actual_model_calls_forbidden=TRUE),file.path(out,'input_snapshot.tsv'))
 cleanup()
 write_kv(list(pass=TRUE,mode='metadata_audit',model_evaluations=0L,optimizer_callbacks=0L,wall_seconds=proc.time()[[3]]-started,published_files=length(published)),file.path(out,'result.tsv'));finalize_payload(out)
},error=function(error){write_kv(list(failed=TRUE,message=conditionMessage(error),model_evaluations=0L),file.path(out,'FAILED.tsv'));stop(error)})
