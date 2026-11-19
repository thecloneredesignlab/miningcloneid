source("oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step03/backend.R")
source(file.path(out,"search_helpers.R"))
cx <- branch_context()
read_points <- function(path) {
  x <- read.delim(path)
  lapply(split(x,x$id),function(tab) {
    p <- cx$anchor
    p[tab$parameter] <- as.list(as.numeric(tab$binary64))
    p
  })
}
points <- read_points(file.path(out,"selected_pair_parameters.tsv"))
rows <- subset(read.delim(file.path(out,"selected_pairs.tsv")),kind=="switch_straddle")
if (file.exists(file.path(out,"targeted_pairs.tsv"))) {
  extra <- read.delim(file.path(out,"targeted_pairs.tsv"))
  extra <- extra[extra$kind=="targeted_fallback_straddle" &
    extra$endpoint==0 & extra$mask==0 & extra$reset_mode==1,]
  if (nrow(extra)) {
    rows <- rbind(rows[rows$category!="reset_mode",],extra[names(rows)])
    points <- c(points,read_points(file.path(out,"targeted_pair_parameters.tsv")))
  }
}
same <- function(a,b) a$ok && b$ok && identical(a$days,b$days) &&
  identical(a$parent_days,b$parent_days) && identical(a$fallback,b$fallback) && identical(a$mask,b$mask)
pars <- list(); outputs <- list(); summaries <- list()
for(i in seq_len(nrow(rows))) {
  row <- rows[i,]; category<-row$category
  left <- points[[as.character(row$left)]]; right <- points[[as.character(row$right)]]
  a <- cx$evaluate(left);b<-cx$evaluate(right)
  stopifnot(a$ok,b$ok,!same(a,b))
  first<-a;last<-b
  for(j in seq_len(18L)) {
    mid<-left
    mid[cx$keep]<-as.list((unlist(left[cx$keep])+unlist(right[cx$keep]))/2)
    z<-cx$evaluate(mid)
    if (same(z,first)) {left<-mid;a<-z} else if(same(z,last)) {right<-mid;b<-z} else break
  }
  d<-differences(a,b)
  key<-paste0("tight_",category,c("_left","_right"))
  pars[[key[1L]]]<-left; pars[[key[2L]]]<-right
  outputs[[key[1L]]]<-cx$evaluate(left,exact_kary=TRUE)
  outputs[[key[2L]]]<-cx$evaluate(right,exact_kary=TRUE)
  summaries[[i]]<-data.frame(category=category,kind="narrow_switch",left=key[1L],
        right=key[2L],unit_distance=max(abs(cx$encode(left)[1:18]-cx$encode(right)[1:18])),as.list(d))
}
emit_parameters(pars,file.path(out,"narrow_pair_parameters.tsv"),cx$keep)
emit_native_pair_payload(outputs,file.path(out,"narrow_native_outputs.tsv"),file.path(out,"narrow_native_branches.tsv"))
sf_emit(do.call(rbind,summaries),file.path(out,"narrow_pairs.tsv"))
print(do.call(rbind,summaries))
