source("oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step03/backend.R")
source(file.path(out,"search_helpers.R"))
cx <- branch_context()
stage <- read.delim(file.path(out,"search_ledger.tsv"),check.names=FALSE)
unit_names <- paste0("u_",cx$keep)
stage <- stage[stage$feasible & is.finite(stage$minimum_target_supply_ratio),]
source <- stage[order(abs(log(stage$minimum_target_supply_ratio))),]
unique_centers <- function(indices,n) {
  pick <- integer()
  for (i in indices) {
    u <- as.numeric(source[i,unit_names])
    if(length(pick)&&any(vapply(pick,function(j)max(abs(as.numeric(source[j,unit_names])-u))<1e-5,FALSE))) next
    pick<-c(pick,i)
    if(length(pick)>=n) break
  }
  pick
}
centers <- c(unique_centers(which(source$reset_mode>0),8L),unique_centers(which(source$reset_mode==0),8L))
ledger <- list(); cache <- list(); pairs <- list()
evaluate <- function(u,label) {
  ident <- paste0("reset_",length(cache)+1L)
  v<-cx$evaluate(cx$decode(u))
  cache[[ident]]<<-list(u=u,v=v)
  ledger[[ident]]<<-data.frame(id=ident,origin=label,feasible=v$ok,
    failure=if(v$ok)""else v$failure,minimum_log_supply=if(v$ok)min(log(v$supply),na.rm=TRUE)else NA_real_)
  ident
}
compare <- function(left,right,edge,kind) {
  a<-cache[[left]];b<-cache[[right]]
  if(!a$v$ok||!b$v$ok||a$v$fallback[edge]==b$v$fallback[edge]) return(NULL)
  d<-differences(a$v,b$v)
  data.frame(category="reset_mode",kind=kind,left=left,right=right,edge=names(a$v$days)[edge],
     unit_distance=max(abs(a$u[1:18]-b$u[1:18])),as.list(d),check.names=FALSE)
}

for(i in centers) {
  u<-as.numeric(source[i,unit_names]); left<-evaluate(u,paste0("center_",source$id[i]))
  y<-cache[[left]]$v; stopifnot(y$ok)
  edge<-which.min(y$supply)
  g<-log(y$supply[edge])
  for(j in c(1L,9L,10L,11L,12L)) {
    h<-if(u[j]>.9999)-1e-5 else 1e-5
    v<-u;v[j]<-v[j]+h
    probe_id<-evaluate(v,paste0("slope_",source$id[i],"_",cx$keep[j]))
    probe<-cache[[probe_id]]$v
    if(!probe$ok) next
    slope<-(log(probe$supply[edge])-g)/h
    if(!is.finite(slope)||abs(slope)<1e-6) next
    for(scale in c(1.05,1.8)) {
      step<- -scale*g/slope
      if(!is.finite(step)||abs(step)>.08||u[j]+step<0||u[j]+step>1||abs(step)<1e-8) next
      q<-u;q[j]<-q[j]+step
      right<-evaluate(q,paste0("cross_",source$id[i],"_",cx$keep[j],"_",scale))
      test<-compare(left,right,edge,"targeted_fallback")
      if(!is.null(test)) pairs[[length(pairs)+1L]]<-test
    }
  }
}
stopifnot(length(pairs)>0L)
candidates<-do.call(rbind,pairs)
score<-sqrt((candidates$growth_rms/cx$anchor$sigma_growth)^2+candidates$kary_tv_rms^2)
isolated<-which(candidates$endpoint==0 & candidates$mask==0 & candidates$reset_mode==1)
if(length(isolated)) selected<-isolated[order(score[isolated])][1L] else selected<-which.min(score)
seed<-candidates[selected,]
left<-seed$left;right<-seed$right
edge<-match(seed$edge,names(cache[[left]]$v$days))
start_mode<-cache[[left]]$v$fallback[edge]
for(i in seq_len(16L)) {
  mid<-evaluate((cache[[left]]$u+cache[[right]]$u)/2,paste0("isolated_bracket_",i))
  y<-cache[[mid]]$v
  if(!y$ok) break
  if(y$fallback[edge]==start_mode) left<-mid else right<-mid
}
tight<-compare(left,right,edge,"targeted_fallback_straddle")
stopifnot(!is.null(tight))
rows<-rbind(seed,tight)
ids<-unique(c(rows$left,rows$right))
pars<-setNames(lapply(ids,function(id)cx$decode(cache[[id]]$u)),ids)
emit_parameters(pars,file.path(out,"targeted_pair_parameters.tsv"),cx$keep)
outputs<-lapply(pars,function(p)cx$evaluate(p,exact_kary=TRUE))
stopifnot(all(vapply(outputs,function(r)r$ok,FALSE)))
emit_native_pair_payload(outputs,file.path(out,"targeted_native_outputs.tsv"),file.path(out,"targeted_native_branches.tsv"))
sf_emit(rows,file.path(out,"targeted_pairs.tsv"))
sf_emit(candidates,file.path(out,"targeted_fallback_candidates.tsv"))
sf_emit(do.call(rbind,ledger),file.path(out,"targeted_evaluations.tsv"))
schedule<-read.delim("oxygen/results/synthetic_identifiability/20261113_culture_contract_steps01_03/experimental_schedule.tsv")
export_changes(outputs,list(targeted_fallback=c(rows$left[1],rows$right[1]),
  targeted_fallback_straddle=c(rows$left[2],rows$right[2])),schedule,file.path(out,"targeted_branch_changes.tsv"))
print(rows)
cat("Callback count",length(cache),"isolated fallback candidates",length(isolated),"\n")
