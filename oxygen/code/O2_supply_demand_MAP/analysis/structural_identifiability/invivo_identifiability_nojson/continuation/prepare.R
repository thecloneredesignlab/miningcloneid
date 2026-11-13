# Metadata-only bridge. Original numerical/input/contract bytes remain unchanged.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);old<-normalizePath(args[2]);entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));contract<-readRDS(file.path(study,'manifests/contract.rds'))
roots<-list();rows<-list()
for(policy in c('production','data_only'))for(kind in c('anchor','uniform_old','uniform_new')) {
 root<-if(kind=='anchor')paste0('archive_h25_',policy) else if(kind=='uniform_new')paste0('new_uniform_',policy) else paste0('formal_',policy,if(policy=='production')'_u1' else '_u2')
 path<-file.path(old,'tables',paste0('root_',root,'_selected.rds'));a<-readRDS(path)
 stopifnot(a$root_id==root,a$policy==policy,a$record$provenance$root_id==root,all(a$record$parameters>=contract$lower),all(a$record$parameters<=contract$upper))
 id<-paste(policy,kind,sep='_');z<-list(original_record=a$record,root_id=root,root_family=a$root_family,target_policy=policy,
  expected_objective=a$record$components$L,original_file_sha256=hash_file(path),original_path=paste0('tables/',basename(path)),
  original_manifest_sha256=hash_file(file.path(old,'manifests/input_manifest.tsv')),original_payload_sha256=a$payload_sha256)
 atomic_rds(z,file.path(study,'reference/imports',paste0(id,'.rds')))
 rows[[id]]<-data.frame(import_id=id,policy=policy,kind=kind,root_id=root,root_family=a$root_family,expected_objective=z$expected_objective,original_sha256=z$original_file_sha256)
}
write_tsv(do.call(rbind,rows),file.path(study,'manifests/roots.tsv'))
atomic_rds(list(prior_lifetime_core_seconds=301846L,prior_formal_core_seconds=271355L,
 aggregate_core_hour_limit=NA_real_,attempt_limit=NA_integer_,score_tolerance=1e-8,gain_epsilon=1e-10,
 baseline_recovery=.002,direction_improvement=.002,profile_recovery=.01,steps=c(.001,.0001)),file.path(study,'manifests/protocol.rds'))
write_kv(list(user_confirmed=TRUE,aggregate_core_hour_limit=NA_real_,culture_protected=TRUE,profile_requires_both_baselines=TRUE,local_window_seconds=7200L,local_callback_allowance=20000L),file.path(study,'manifests/authorization.tsv'))
paths<-sort(c(list.files(file.path(study,'source'),recursive=TRUE,full.names=TRUE,all.files=TRUE,no..=TRUE),list.files(file.path(study,'reference'),recursive=TRUE,full.names=TRUE,all.files=TRUE,no..=TRUE),list.files(file.path(study,'manifests'),recursive=TRUE,full.names=TRUE)))
paths<-paths[!file.info(paths)$isdir & basename(paths)!='input_manifest.tsv']
write_tsv(data.frame(path=substring(paths,nchar(study)+2L),sha256=vapply(paths,hash_file,character(1))),file.path(study,'manifests/input_manifest.tsv'))
cat('Six original root records retained; no aggregate/attempt ceiling; metadata preparation only.\n');assert_namespaces()
