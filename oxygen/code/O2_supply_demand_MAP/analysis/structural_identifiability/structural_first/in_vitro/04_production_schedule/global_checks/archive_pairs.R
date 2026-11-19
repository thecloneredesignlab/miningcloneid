source("oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step03/backend.R")
source("oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step03/search_helpers.R")
cx <- branch_context()
ctrl <- readRDS(file.path(out,"controls.rds"))
anchor <- ctrl$anchor
growth_ids <- names(anchor$growth)
kary_ids <- colnames(anchor$kary)
sg <- cx$anchor$sigma_growth
schedule <- read.delim("oxygen/results/synthetic_identifiability/20261113_culture_contract_steps01_03/experimental_schedule.tsv")
label <- paste(schedule$cohort,schedule$segment_id,sep="|")
parent <- match(paste(schedule$cohort,schedule$parent_segment_id,sep="|"),label)
bank <- "oxygen/results/synthetic_identifiability/20261117_culture_D19_recovery_restart_01/inputs/banks"

# All pairs are compared in a full-output Euclidean screen. Exact native
# growth-RMS/TV metrics are used to rerank the first 100 pairs in each class.
m <- expand.grid(N=22:154,O=sort(unique(schedule$oxygen_pct)))
events <- lapply(m$N,function(n)seq.int(0L,n))
eventm <- unlist(events); eventrow <- rep(seq_len(nrow(m)),lengths(events))
signature <- function(p) {
  h <- p$O2_crit^p$n_O/(p$O2_crit^p$n_O+m$O^p$n_O)
  mu <- p$mu_hp*h*(m$N/44)^p$gamma_mu
  prob <- p$p_mis_base+p$p_misseg*mu/(p$k_o_mis+mu)
  keep <- dbinom(eventm,m$N[eventrow],prob[eventrow])>=1e-8
  digest::digest(which(keep),algo="sha256")
}
rows <- list()
for (n in LETTERS[1:4]) {
  package <- readRDS(file.path(bank,paste0("D_",n,".rds")))
  for (i in seq_along(package$members)) {
    z <- package$members[[i]]
    s <- z$pred$selection
    ids <- match(label,paste(s$cohort,s$segment_id,sep="|"))
    stopifnot(!anyNA(ids))
    day <- s$day[ids]; parent_day <- s$parent_day[ids]
    fallback <- !is.na(parent_day)&parent_day!=day[parent]
    p <- cx$anchor
    p[names(z$params)] <- z$params
    rows[[paste0("bank_",n,"_",i)]] <- list(p=p,unit=cx$encode(p),day=day,fallback=fallback,mask=signature(p),
      growth=z$pred$growth[growth_ids],kary=do.call(cbind,z$pred$kary[kary_ids]))
  }
}
stopifnot(length(rows)==800L)
rows <- rows[vapply(rows,function(r) all(r$unit>=-1e-14&r$unit<=1+1e-14),FALSE)]
n_bank <- length(rows)
cat("Archived points within the nominal Step 02 box",n_bank,"of 800\n")
vectors <- do.call(rbind,lapply(rows,function(x)c(x$growth/(sg*sqrt(114)),as.numeric(x$kary)/sqrt(12))))
stopifnot(all(is.finite(vectors)))
gram <- tcrossprod(vectors)
norm <- diag(gram)
gap2 <- outer(norm,norm,"+")-2*gram
gap2[!upper.tri(gap2)] <- Inf
distance <- matrix(0,n_bank,n_bank)
units <- do.call(rbind,lapply(rows,function(r)r$unit))
for(j in 1:18) distance<-pmax(distance,abs(outer(units[,j],units[,j],"-")))
id_by <- function(field) vapply(rows,function(x)digest::digest(x[[field]],algo="sha256"),"")
categories <- c(endpoint="day",mask="mask",reset_mode="fallback")
ledger <- list(); pairs <- list(); ranks <- list()
for (category in names(categories)) {
  codes <- id_by(categories[[category]])
  qualifies <- is.finite(gap2)&distance>=.001&outer(codes,codes,"!=")
  wh <- which(qualifies)
  wh <- head(wh[order(gap2[wh])],100L)
  if(!length(wh)) {
    ledger[[category]] <- data.frame(category=category,screened=0L,reranked=0L)
    next
  }
  candidate <- arrayInd(wh,dim(gap2))
  score <- vapply(seq_len(nrow(candidate)),function(i) {
    a <- rows[[candidate[i,1]]]; b <- rows[[candidate[i,2]]]
    sqrt(mean(((a$growth-b$growth)/sg)^2)+mean((colSums(abs(a$kary-b$kary))/2)^2))
  },0)
  stopifnot(length(score)>0L)
  ordered <- order(score)
  i <- ordered[1L]
  sizes <- data.frame(category=category,screened=sum(qualifies),reranked=length(score))
  ledger[[category]] <- sizes
  ranks[[category]] <- data.frame(category=category,left=names(rows)[candidate[ordered,1]],
       right=names(rows)[candidate[ordered,2]],archived_score=score[ordered],
       archive_euclidean=sqrt(pmax(0,gap2[wh[ordered]])),unit_distance=distance[wh[ordered]])
  ids <- names(rows)[candidate[i,]]
  pairs[[paste0(category,"_archive_material_match")]] <- ids
}
sf_emit(do.call(rbind,ledger),file.path(out,"archive_pair_coverage.tsv"))
sf_emit(do.call(rbind,ranks),file.path(out,"archive_pair_candidates.tsv"))
ids <- sort(unique(unlist(pairs)))
pars <- setNames(lapply(ids,function(id)rows[[id]]$p),ids)
emit_parameters(pars,file.path(out,"additional_pair_parameters.tsv"),cx$keep)
outputs <- lapply(pars,function(p)cx$evaluate(p,exact_kary=TRUE))
for(id in ids) {
  a<-outputs[[id]];b<-rows[[id]]
  stopifnot(a$ok,max(abs(a$growth-b$growth))<1e-10,max(abs(a$kary-b$kary))<1e-10,
       identical(unname(a$days),b$day),identical(a$fallback,b$fallback))
}
selected <- do.call(rbind,lapply(names(pairs),function(name) {
  ids<-pairs[[name]]; l<-rows[[ids[1]]];r<-rows[[ids[2]]]
  d<-differences(outputs[[ids[1]]],outputs[[ids[2]]])
  data.frame(category=sub("_archive_material_match$","",name),kind="archive_material_match",left=ids[1],right=ids[2],
       unit_distance=max(abs(l$unit[1:18]-r$unit[1:18])),as.list(d),check.names=FALSE)
}))
sf_emit(selected,file.path(out,"additional_pairs.tsv"))
export_changes(outputs,pairs,schedule,file.path(out,"additional_branch_changes.tsv"))
emit_native_pair_payload(outputs,file.path(out,"additional_native_outputs.tsv"),file.path(out,"additional_native_branches.tsv"))
print(selected)
