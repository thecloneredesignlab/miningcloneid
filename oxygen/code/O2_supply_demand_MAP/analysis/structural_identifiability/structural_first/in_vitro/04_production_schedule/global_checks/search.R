source("oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step03/backend.R")
source("oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step03/search_helpers.R")
cx <- branch_context()
control <- readRDS(file.path(out,"controls.rds"))
anchor <- control$anchor
anchor$details <- NULL
sg <- cx$anchor$sigma_growth
schedule <- read.delim("oxygen/results/synthetic_identifiability/20261113_culture_contract_steps01_03/experimental_schedule.tsv")
children <- paste(schedule$cohort,schedule$parent_segment_id,sep="|")
design_labels <- paste(schedule$cohort,schedule$segment_id,sep="|")
stopifnot(identical(names(anchor$days),design_labels))
actual <- new.env(parent=emptyenv())
actual$seq <- list()
actual$records <- list()
actual$evaluations <- 0L
actual$label <- NA_character_
create <- function(u,label,b=NULL) {
  if(is.null(b)) {actual$evaluations<-actual$evaluations+1L; b<-cx$evaluate(cx$decode(u))}
  j <- length(actual$seq)+1L
  row <- data.frame(id=j,origin=label,feasible=b$ok,score=NA_real_,distance_to_anchor=max(abs(u[1:18]-cx$u[1:18])),
    failure=NA_character_,minimum_target_supply_ratio=NA_real_,growth_max=NA_real_,growth_rms=NA_real_,
    kary_tv_max=NA_real_,kary_tv_rms=NA_real_,endpoint=NA_integer_,mask=NA_integer_,reset_mode=NA_integer_,parent_day=NA_integer_)
  if(b$ok) {
    flags <- differences(anchor,b)
    row[names(flags)] <- as.list(flags)
    row$score <- sqrt((flags["growth_rms"]/sg)^2+flags["kary_tv_rms"]^2)
    row$minimum_target_supply_ratio <- if (all(is.na(b$supply))) NA_real_ else min(b$supply,na.rm=TRUE)
  } else row$failure <- b$failure
  actual$seq[[j]] <- list(id=j,u=u,value=b,row=row)
  actual$records[[j]] <- row
  j
}
ref <- function(i) actual$seq[[i]]
create(cx$u,"source_anchor",anchor)
set.seed(314159)
make_u <- function() {
  trials <- list()
  add <- function(u,val) {trials[[length(trials)+1L]] <<- list(u=u,tag=val)}
  for (j in 1:18) for (sign in c(-1,1)) for (h in c(0.003,0.015,0.075)) {
    u<-cx$u; u[j]<-u[j]+sign*h
    if (all(u>=0&u<=1)) add(u,paste("axis",cx$keep[j],sign,h,sep="_"))
  }
  for (band in c(.004,.02,.08,.22)) for (i in seq_len(24L)) {
    u <- pmin(1,pmax(0,cx$u+c(rnorm(18L)*band,0)))
    add(u,paste("mixed",band,i,sep="_"))
  }
  for (j in c(1,2,3,4,6,7,10,11,13,14)) {
    for (value in c(.1,.35,.65,.9)) {
      u<-cx$u;u[j]<-value
      add(u,paste("wide",cx$keep[j],value,sep="_"))
    }
  }
  trials
}
trials <- make_u()
for (t in trials) create(t$u,t$tag)
cat("New axis and mixed scan",actual$evaluations,"completed\n",flush=TRUE)

# Screen all archived feasible vectors with original predicted responses, then
# replay the most relevant ones in the current local runner for branch metrics.
bank <- "oxygen/results/synthetic_identifiability/20261117_culture_D19_recovery_restart_01/inputs/banks"
archived <- list()
for (origin in LETTERS[1:4]) {
  pkg <- readRDS(file.path(bank,paste0("D_",origin,".rds")))
  stopifnot(length(pkg$members)==200L)
  for (n in seq_along(pkg$members)) {
    m <- pkg$members[[n]]
    p <- m$params
    u <- pmin(1,pmax(0,cx$encode(p)))
    pred <- m$pred
    stopifnot(!is.null(pred),length(pred$growth)==114L,length(pred$kary)==12L)
    ss <- pred$selection
    rk <- match(design_labels,paste(ss$cohort,ss$segment_id,sep="|"))
    stopifnot(!anyNA(rk))
    day <- ss$day[rk]; parent <- ss$parent_day[rk]
    original_parent_day <- day[match(children,design_labels)]
    modes <- !is.na(parent)&parent!=original_parent_day
    growth <- pred$growth[names(anchor$growth)]
    k <- do.call(cbind,pred$kary[colnames(anchor$kary)])
    arch_score <- sqrt(mean(((growth-anchor$growth)/sg)^2)+mean((colSums(abs(k-anchor$kary))/2)^2))
    archived[[length(archived)+1L]] <- list(u=u,score=arch_score,fallback=sum(modes),endpoint=sum(day!=anchor$days),
                    id=paste(origin,n,sep="_"),parent_day=sum(parent!=anchor$parent_days,na.rm=TRUE))
  }
}
bank_table <- do.call(rbind,lapply(archived,function(a)as.data.frame(a[names(a)!="u"])))
sf_emit(bank_table,file.path(out,"archived_feasible_screen.tsv"))
order_bank <- order(bank_table$score)
top <- unique(c(head(order_bank,16L),head(order_bank[bank_table$fallback[order_bank]>0],16L)))
for (i in top) create(archived[[i]]$u,paste0("bank_",archived[[i]]$id))
cat("Archived feasible screen",nrow(bank_table),"archived vectors, fresh replays",length(top),"\n",flush=TRUE)

categories <- c("endpoint","mask","reset_mode")
eligible <- function(row,category,gap=1e-3) {
  isTRUE(row$feasible) && row$distance_to_anchor>=gap && row[[category]]>0
}
best <- function(category, k=1L) {
  rows <- do.call(rbind,actual$records)
  good <- vapply(seq_len(nrow(rows)),function(i)eligible(rows[i,,drop=FALSE],category),FALSE)
  ids <- rows$id[good]; head(ids[order(rows$score[good])],k)
}
before <- lapply(categories,best,k=5L); names(before)<-categories
print(before)
sf_emit(do.call(rbind,actual$records),file.path(out,"initial_scan.tsv"))

# A safeguarded independent optimization for the nearest anchor match in each
# branch category. The full 19-coordinate search uses the original feasible box.
for (category in categories) {
  selected <- head(before[[category]],2L)
  for (i in selected) {
    signal <- ref(i)
    f <- function(u) {
      id <- create(u,paste0("tune_",category,"_",i))
      r <- ref(id)$row
      if (eligible(r,category)) r$score else 10+if(is.finite(r$score)) r$score else 1
    }
    fit <- nloptr::nloptr(signal$u, f, lb=rep(0,19L),ub=rep(1,19L),
                         opts=list(algorithm="NLOPT_LN_SBPLX",maxeval=80L,xtol_rel=1e-7))
    cat("Optimization",category,"start",i,"status",fit$status,"calls",fit$iterations,"\n",flush=TRUE)
  }
}

# Approach both sides of an actual transition on rays from the anchor. The
# straddling pairs diagnose jumps or flatness, separately from material matches.
pairs <- list()
probe_summary <- list()
refined <- list()
for (category in categories) {
  final <- best(category,4L)
  if (!length(final)) next
  add_pair <- function(left,right,kind) {
    d <- differences(ref(left)$value,ref(right)$value)
    distance <- max(abs(ref(left)$u[1:18]-ref(right)$u[1:18]))
    probe_summary[[length(probe_summary)+1L]] <<- data.frame(category=category,kind=kind,left=left,right=right,
                           unit_distance=distance,as.list(d),check.names=FALSE)
    pairs[[paste(category,kind,sep="_")]] <<- c(left,right)
  }
  add_pair(1L,final[[1L]],"anchor_material_match")
  if(length(before[[category]])) add_pair(1L,before[[category]][[1L]],"anchor_material_before_tuning")
  jumps <- list()
  for (candidate in head(final,3L)) {
    left <- 1L; right <- candidate
    for (iter in seq_len(14L)) {
      new <- create((ref(left)$u+ref(right)$u)/2,paste0("bracket_",category,"_",candidate))
      val <- ref(new)
      if (!val$value$ok) break
      if (val$row[[category]]>0) right<-new else left<-new
    }
    diff <- differences(ref(left)$value,ref(right)$value)
    if(diff[[category]]==0) next
    jumps[[length(jumps)+1L]] <- data.frame(left=left,right=right,candidate=candidate,
     score=sqrt((diff["growth_rms"]/sg)^2+diff["kary_tv_rms"]^2),as.list(diff),check.names=FALSE)
  }
  if(length(jumps)) {
    all_jumps <- do.call(rbind,jumps)
    refined[[category]] <- data.frame(category=category,all_jumps)
    j <- order(all_jumps$score)[1L]
    add_pair(all_jumps$left[j],all_jumps$right[j],"switch_straddle")
  }
}
all <- do.call(rbind,actual$records)
units <- do.call(rbind,lapply(actual$seq,function(r) r$u))
colnames(units) <- paste0("u_",cx$keep)
sf_emit(cbind(all,units),file.path(out,"search_ledger.tsv"))
sf_emit(do.call(rbind,probe_summary),file.path(out,"selected_pairs.tsv"))
sf_emit(do.call(rbind,refined),file.path(out,"all_brackets.tsv"))

ids <- sort(unique(unlist(pairs)))
pars <- setNames(lapply(ids,function(id)if(id==1L)cx$anchor else cx$decode(ref(id)$u)),ids)
emit_parameters(pars,file.path(out,"selected_pair_parameters.tsv"),cx$keep)
outputs <- lapply(pars,function(p)cx$evaluate(p,keep_details=FALSE,exact_kary=TRUE))
stopifnot(all(vapply(outputs,function(a)a$ok,FALSE)))
observed <- do.call(rbind,lapply(names(outputs),function(id) {
  r<-outputs[[id]]
  data.frame(id=id,stream=c(rep("growth",114L),rep("kary",length(r$kary))),passage=c(names(r$growth),rep(colnames(r$kary),each=133L)),
      chr=c(rep("",114L),rep(rownames(r$kary),12L)),value=c(r$growth,as.numeric(r$kary)))
}))
sf_emit(observed,file.path(out,"paired_native_outputs.tsv"))
sf_rds(outputs,file.path(out,"paired_native.rds"))
sf_rds(pairs,file.path(out,"pairs.rds"))
export_changes(outputs,pairs,schedule,file.path(out,"paired_branch_changes.tsv"))
cat("Fresh production search callbacks",actual$evaluations,"plus",length(outputs),"selected exact-kary replays\n",flush=TRUE)
print(do.call(rbind,probe_summary))
