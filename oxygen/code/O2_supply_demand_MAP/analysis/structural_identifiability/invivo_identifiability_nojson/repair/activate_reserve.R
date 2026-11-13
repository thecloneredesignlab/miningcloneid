# Finite reserve choices from closed, audited evidence only; no model calls.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2L)
study<-normalizePath(args[1]);action<-args[2]
entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));source(file.path(entry,'helpers.R'))
jobs<-read.delim(file.path(study,'tables/jobs.tsv'),colClasses='character')
stopifnot(all(jobs$state=='COMPLETED'),all(jobs$exit_code=='0:0'))
specs<-readRDS(file.path(study,'manifests/task_specs.rds'))
chain<-read_tsv(file.path(study,'tables/baseline_chain_summary.tsv'))
gates<-read_tsv(file.path(study,'tables/baseline_policy_gates.tsv'))
for(policy in gates$policy[!gates$pass])stopifnot(all(paste0('r_',policy,'_uniform_new_',c('global','local','fresh','directions'))%in%jobs$task))
contract<-readRDS(file.path(study,'manifests/contract.rds'))
activation_path<-file.path(study,'tables/reserve_activation.tsv')
activation<-if(file.exists(activation_path))read_tsv(activation_path) else data.frame(task=character(),selection_name=character(),reason=character())
decisions_path<-file.path(study,'tables/reserve_decision.tsv')
decisions<-if(file.exists(decisions_path))read_tsv(decisions_path) else data.frame(task=character(),policy=character(),root_id=character(),family=character(),selection_name=character(),selection_sha256=character(),selection_objective=numeric(),reason=character())
owned_point<-function(root_id) {
 path<-file.path(study,'tables',paste0('root_',root_id,'_selected.rds'))
 point<-readRDS(path);stopifnot(point$root_id==root_id,point$record$provenance$root_id==root_id)
 # The reserve global may improve its own root; it never supplies another root.
 for(task in jobs$task[jobs$role=='reserve_global']) {
  root<-file.path(study,'fits',task);verify_payload_R(root)
  sp<-readRDS(file.path(root,'dispatch_spec.rds'))
  if(sp$root_id!=root_id)next
  proof<-read_kv(file.path(study,'tables/record_audits',paste0(task,'_summary.tsv')))
  stopifnot(proof$pass,proof$payload_sha256==hash_file(file.path(root,'payload_manifest.tsv')))
  a<-readRDS(file.path(root,'accepted_candidate.rds'))
  stopifnot(a$policy==point$policy,a$record$provenance$root_id==root_id,
   abs(a$checked$L-a$record$components$L)<1e-8)
  if(a$checked$L<point$record$components$L)point<-list(record=a$record,policy=a$policy,root_id=root_id,
   root_family=sp$root_family,source_task=task,source_scope='fits',payload_sha256=hash_file(file.path(root,'payload_manifest.tsv')))
 }
 stopifnot(point$record$source_manifest_sha256==hash_file(file.path(study,'manifests/input_manifest.tsv')),
  point$record$contract_sha256==hash_file(file.path(study,'manifests/contract.rds')),
  all(point$record$parameters>=contract$lower),all(point$record$parameters<=contract$upper))
 point
}
record_choice<-function(task,point,reason) {
 stopifnot(!task%in%c(activation$task,jobs$task),task%in%names(specs),specs[[task]]$stage=='reserve',
  specs[[task]]$policy==point$policy,specs[[task]]$root_family==point$root_family,
  specs[[task]]$root_id%in%c('dispatch',point$root_id))
 name<-paste0('root_',point$root_id,'_',task,'_selected')
 path<-file.path(study,'tables',paste0(name,'.rds'));atomic_rds(point,path)
 write_kv(list(policy=point$policy,root_id=point$root_id,root_family=point$root_family,
  selection_sha256=hash_file(path)),file.path(study,'tables',paste0(name,'.tsv')))
 activation<<-rbind(activation,data.frame(task=task,selection_name=name,reason=reason))
 decisions<<-rbind(decisions,data.frame(task=task,policy=point$policy,root_id=point$root_id,family=point$root_family,
  selection_name=name,selection_sha256=hash_file(path),selection_objective=point$record$components$L,reason=reason))
}
if(action=='global') {
 stopifnot(!any(jobs$role=='reserve_global'),!any(grepl('_global$',activation$task)))
 failed<-gates[!gates$pass & gates$uniform_gap>.002,,drop=FALSE]
 if(nrow(failed)) {
  # Largest unresolved uniform gap / .002; exact ties favor data-only.
  failed<-failed[order(-failed$uniform_gap,failed$policy!='data_only'),,drop=FALSE]
  policy<-failed$policy[1];rows<-chain[chain$policy==policy & chain$family=='uniform',,drop=FALSE]
  row<-rows[which.min(rows$verified_best),,drop=FALSE];point<-owned_point(row$root_id)
  suffix<-if(point$root_id==paste0('new_uniform_',policy))'uniform_new' else 'uniform_old'
  task<-paste0('r_',policy,'_reserve_',suffix,'_global')
  record_choice(task,point,paste0('Largest failing independent-uniform gap / 0.002: ',sprintf('%.17g',failed$uniform_gap[1]/.002),'; own-root hybrid restart only'))
 }
} else if(action=='locals') {
 stopifnot(!any(jobs$role=='reserve_local'),!any(grepl('_local$',activation$task)))
 if(any(grepl('_global$',activation$task)))stopifnot(all(activation$task[grepl('_global$',activation$task)]%in%jobs$task))
 failed<-gates[!gates$pass,,drop=FALSE]
 stopifnot(nrow(failed)<=2L)
 for(policy in failed$policy) {
  rows<-chain[chain$policy==policy,,drop=FALSE]
  points<-lapply(unique(rows$root_id),owned_point)
  values<-vapply(points,function(p)p$record$components$L,numeric(1))
  families<-vapply(points,`[[`,character(1),'root_family')
  best<-min(values);family_gap<-vapply(c('archive','uniform'),function(f)min(values[families==f])-best,numeric(1))
  if(max(family_gap)>.002)family<-names(which.max(family_gap))
  else family<-families[which.min(values)] # Remaining directional/floor deficit.
  take<-which(families==family);point<-points[[take[which.min(values[take])]]]
  task<-paste0('r_',policy,'_reserve_',family,'_1_local')
  record_choice(task,point,paste0('One named existing root per still-failing policy; ',family,' family gap ',sprintf('%.17g',family_gap[family]),'; no winner transfer'))
  for(check in c('fresh','directions'))activation<-rbind(activation,data.frame(task=sub('_local$',paste0('_',check),task),selection_name='',reason='Separate reserve-local verification; original parent dependency'))
 }
} else stop('Expected global or locals')
stopifnot(sum(grepl('_global$',activation$task))<=1L,sum(grepl('_local$',activation$task))<=2L,
 sum(grepl('_fresh$',activation$task))<=2L,sum(grepl('_directions$',activation$task))<=2L,
 !anyDuplicated(activation$task))
write_tsv(activation,activation_path);write_tsv(decisions,decisions_path)
print(decisions,row.names=FALSE);assert_namespaces()
