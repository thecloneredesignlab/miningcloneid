# Snapshot/stage bookkeeping only. No objective or trajectory access.
metadata_profile_phase<-function(epoch_exists,baseline_pass,prior_invalid)epoch_exists&&baseline_pass&&!prior_invalid
freeze_baseline_evidence<-function(study,evidence){
 directory<-file.path(study,'tables/baseline_gate_epochs',basename(tempfile('g_')));dir.create(directory,recursive=TRUE)
 vapply(unique(evidence),function(path){
  if(!startsWith(path,'tables/'))return(path)
  target<-file.path(directory,path);dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  stopifnot(file.copy(file.path(study,path),target),file.exists(target))
  substring(target,nchar(study)+2L)
 },character(1),USE.NAMES=FALSE)
}
