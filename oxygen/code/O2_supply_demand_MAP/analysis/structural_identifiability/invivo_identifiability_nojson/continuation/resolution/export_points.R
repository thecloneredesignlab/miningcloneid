# Inspectable selected-parameter/component/native-stream exports; zero model calls.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));contract<-readRDS(file.path(study,'manifests/contract.rds'));tables<-file.path(study,'tables')
paths<-list.files(tables,pattern='^(baseline_(production|data_only)|baseline_stopping_rep_(production|data_only)|resolution_anchor_(production|data_only)|root_.*_selected|profile_e[0-9]+_.*_selected)[.]rds$',full.names=TRUE)
parameters<-components<-streams<-list()
for(path in paths){
 a<-readRDS(path);r<-a$record;p<-r$parameters;id<-sub('[.]rds$','',basename(path));kind<-if(startsWith(id,'baseline_stopping_rep_'))'checked_stopping_representative' else if(startsWith(id,'resolution_anchor_'))'completed_repair_anchor' else if(startsWith(id,'baseline_'))'free_candidate' else if(startsWith(id,'root_'))'root_candidate' else 'profile_branch_candidate'
 stopifnot(r$source_manifest_sha256==hash_file(file.path(study,'manifests/input_manifest.tsv')),identical(names(p),names(contract$endpoint)),all(p>=contract$lower),all(p<=contract$upper))
 natural<-p;islog<-startsWith(names(p),'log10_');natural[islog]<-10^p[islog];normalized<-(p-contract$lower)/(contract$upper-contract$lower)
 fixed<-names(r$fixed);parameters[[id]]<-data.frame(point=id,kind=kind,policy=r$policy,root_id=r$provenance$root_id,source_task=a$source_task,
  parameter=names(p),transformed=unname(p),natural=unname(natural),normalized_position=unname(normalized),fixed=names(p)%in%fixed,near_lower=unname(normalized<=1e-8),near_upper=unname(normalized>=1-1e-8),selected_record_sha256=hash_file(path))
 components[[id]]<-cbind(data.frame(point=id,kind=kind,policy=r$policy,root_id=r$provenance$root_id,source_task=a$source_task),as.data.frame(component_values(r$components)))
 for(i in seq_along(r$predictions$native)){
  n<-r$predictions$native[[i]]$result;scenario<-contract$scenarios[[i]]
  stopifnot(length(n$Vmm3_total_obs)==length(scenario$obs_days))
  streams[[length(streams)+1L]]<-data.frame(point=id,kind=kind,policy=r$policy,scenario_index=i,harvest=scenario$harvest,output_index=seq_along(n$Vmm3_total_obs),day=scenario$obs_days,
   Vmm3_total_obs=n$Vmm3_total_obs,Ntot_live_obs=n$Ntot_live_obs,Ntot_dead_hypoxia_obs=n$Ntot_dead_hypoxia_obs,Ntot_dead_buffer_obs=n$Ntot_dead_buffer_obs,Ntot_total_obs=n$Ntot_total_obs,O2_eff_obs=n$O2_eff_obs,O2_target_obs=n$O2_target_obs)
 }
}
if(length(parameters)){
 write_tsv(do.call(rbind,parameters),file.path(tables,'selected_point_parameters.tsv'));write_tsv(do.call(rbind,components),file.path(tables,'selected_point_components.tsv'));write_tsv(do.call(rbind,streams),file.path(tables,'selected_native_streams.tsv'))
}
assert_namespaces();cat(length(paths),'checked candidate selections exported; labels do not imply scientific acceptance; native stock/oxygen field names preserved; zero model calls.\n')
