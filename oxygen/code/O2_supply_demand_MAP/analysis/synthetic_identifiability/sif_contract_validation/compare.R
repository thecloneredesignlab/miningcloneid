source(file.path(dirname(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])),"common.R"))
load_defs();dir.create(file.path(study,"reports"),showWarnings=FALSE)
checks<-list()
for(case in tab(file.path(study,"cases.tsv"))$case) {
  folder<-file.path(study,"results/replays",paste0(case,"_rep"))
  a<-readRDS(file.path(paste0(folder,1),"replayed_components.rds"));b<-readRDS(file.path(paste0(folder,2),"replayed_components.rds"))
  assert_predictions(a,b)
  delta<-max(abs(score_components(a)-score_components(b)));stopifnot(delta<=1e-8)
  for(n in c("growth_df","ploidy_df","flow_df")) {
    x<-a[[n]];y<-b[[n]];stopifnot(identical(x$passage_id,y$passage_id))
    column<-if(n=="growth_df")"loglik" else "mean_loglik"
    stopifnot(max(abs(x[[column]]-y[[column]]))<=1e-8)
  }
  checks[[case]]<-data.frame(case=case,state_comparison_pass=TRUE,score_maximum_difference=delta,pass=TRUE)
}
emit(do.call(rbind,checks),file.path(study,"reports/repeat_state_comparison.tsv"))
record_environment(file.path(study,"reports"))
kv(list(pass=TRUE,compared_pairs=4L,model_propagations=0L,optimizer_calls=0L,synthetic_datasets=0L),file.path(study,"reports/independent_comparison.tsv"))
