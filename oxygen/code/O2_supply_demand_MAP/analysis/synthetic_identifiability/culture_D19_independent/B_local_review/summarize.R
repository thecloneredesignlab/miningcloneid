# Independent offline closure and saved-distribution checks. No simulator is loaded.
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
root <- normalizePath(args[1L])
archive <- normalizePath(args[2L])
out <- normalizePath(args[3L])
code <- normalizePath(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))))
source(file.path(dirname(code),'closeout_replays.R'))
source(file.path(dirname(code),'closeout_polls.R'))
source(file.path(code,'review.R'))

jobs <- read_tsv(file.path(root,'accounting/jobs.tsv'))
stopifnot(!anyDuplicated(jobs$job_id),sum(jobs$attempt==2)<=1L,
          all(jobs$live_reserved_seconds==0),sum(jobs$settled_seconds)<=7200,
          sum(jobs$charged_calls)<=155)
phase <- function(stage) {
  r <- tail(jobs[jobs$stage==stage,,drop=FALSE],1)
  stopifnot(nrow(r)==1L,r$state=='COMPLETED',r$exit_code=='0:0')
  path <- file.path(root,'verification',r$task)
  app <- read_tsv(file.path(path,'application_result.tsv'))
  policy <- read_tsv(file.path(path,'policy_gate.tsv'))
  stopifnot(app$status=='PASS',as.character(app$job)==as.character(r$job_id),policy$status=='PASS')
  list(path=path,job=as.character(r$job_id))
}

case <- readRDS(file.path(root,'inputs/public/D19.rds'))
existing <- readRDS(file.path(root,'verification/D_B_anchor/replay_b.rds'))
a <- existing$candidate
stopifnot(a$key==anchor_key,digest::digest(as.numeric(a$unit),algo='sha256')==a$key)
gate <- read_tsv(file.path(archive,'reports/origin_verdicts.tsv'))
passing_origins <- sum(tolower(as.character(gate$prediction_pass))=='true')
stopifnot(passing_origins==3L)
admit <- read_tsv(file.path(phase('admit')$path,'anchor_admission.tsv'))
p <- phase('poll')$path
resolution <- read_tsv(file.path(p,'local_resolution_audit.tsv'))
coverage <- read_tsv(file.path(p,'poll_coverage.tsv'))
poll <- read_tsv(file.path(p,'poll_candidates.tsv'))
selected <- readRDS(file.path(p,'selected_point.rds'))
keys <- poll_keys(a,setdiff(case$design$bounds$parameter,'sigma_growth'))
candidate_key <- selected$token$key
check_key <- function(x) identical(digest::digest(as.numeric(x$token$unit),algo='sha256'),candidate_key)
stopifnot(nrow(coverage)==150L,identical(coverage$key[coverage$status=='RUN'],keys),
          identical(poll$key,keys),resolution$complete,
          resolution$evaluated==length(keys),admit$expected_unique_probes==length(keys),
          resolution$selected_key==candidate_key,check_key(selected),
          selected$score>=0,selected$score<=a$score,
          max(0,resolution$anchor_D-min(poll$D))<=1e-5,
          all(poll$D[poll$feasible]>=0),all(d19_metrics(case,existing$pred)$checks))

is_lower <- candidate_key!=anchor_key
stopifnot(identical(is_lower,resolution$lower_found))
receipt <- function(x,job_id,path) {
  m <- d19_metrics(case,x$prediction)
  stopifnot(check_key(x), all(m$checks), abs(m$metrics[['D']]-selected$score)<=1e-12)
  data.frame(key=candidate_key,slurm_job=job_id,path=path,as.list(m$metrics),
             all_four_prediction_gates=all(m$checks),worst_sample=m$worst_sample)
}
for_replay <- list()
if (is_lower) {
  final <- lapply(c('replayA','replayB'),phase)
  stopifnot(final[[1L]]$job!=final[[2L]]$job)
  for (i in seq_along(final)) {
    x <- readRDS(file.path(final[[i]]$path,'replayed_point.rds'))
    stopifnot(final[[i]]$job==x$allocation,identical(x$token$unit,selected$token$unit),
              identical(x$prediction$selection,selected$prediction$selection))
    compare <- function(left,right) max(vapply(names(left),function(k) sum(abs(left[[k]]-right[[k]]))/2,numeric(1)))
    stopifnot(compare(x$prediction$kary,selected$prediction$kary)<=1e-10,
              compare(x$prediction$latent_kary,selected$prediction$latent_kary)<=1e-10,
              max(abs(x$prediction$growth-selected$prediction$growth))<=1e-10,
              identical(names(x$prediction$states),names(selected$prediction$states)),
              all(vapply(names(x$prediction$states),function(k) {
                a <- x$prediction$states[[k]]
                b <- selected$prediction$states[[k]]
                identical(dim(a),dim(b)) && max(abs(a-b)/pmax(1,abs(a),abs(b)))<=1e-10
              },logical(1))))
    for_replay[[i]] <- receipt(x,final[[i]]$job,c('new_instrumented','new_independent')[i])
  }
} else {
  for (i in seq_len(2L)) {
    x <- readRDS(file.path(root,'verification/D_B_anchor',c('replay_a.rds','replay_b.rds')[i]))
    x$token <- list(unit=x$candidate$unit)
    x$prediction <- x$pred
    for_replay[[i]] <- receipt(x,x$slurm_job,'archived_verified_anchor')
  }
}
emit(do.call(rbind,for_replay),file.path(out,'selected_point_replays.tsv'))

# Assess the observation operator using the fixed saved latent laws.
archive_d <- file.path(archive,'verification/P')
parameters <- read_tsv(file.path(archive_d,'all20_parameter_recovery.tsv'))
best <- parameters[parameters$key %in% gate$best_found_key & parameters$parameter %in% width_names,]
stopifnot(nrow(best)==9L)
emit(best,file.path(out,'archived_width_and_kernel_comparison.tsv'))
prob <- read_tsv(file.path(archive_d,'full_karyotype_prediction_comparison.tsv'))
prob <- prob[prob$key==anchor_key,,drop=FALSE]
known <- best$generating[best$origin=='D_B'&best$parameter=='sigma_kary']
fit <- best$estimate[best$origin=='D_B'&best$parameter=='sigma_kary']
stopifnot(length(known)==1L,length(fit)==1L,nrow(prob)==133L*12L)
grid <- 22:154
kernel <- function(sigma) {
  weights <- outer(grid,grid,function(y,x)dnorm(y,x,sigma))
  sweep(weights,2,colSums(weights),'/')
}
kv <- list(generating=kernel(known),fitted_B=kernel(fit))
per_sample <- split(prob,prob$passage_id)
tbl <- list()
for (latent in names(kv)) for (observer in names(kv)) for (id in names(per_sample)) {
  x <- per_sample[[id]]
  stopifnot(identical(x$N,grid))
  distribution <- x[[if(latent=='generating') 'generating_latent_probability' else 'latent_probability']]
  q <- as.numeric(kv[[observer]] %*% distribution)
  target <- x$generating_observed_probability
  saved <- x[[if(latent=='generating') 'generating_observed_probability' else 'observed_probability']]
  tv_saved <- sum(abs(q-saved))/2
  if (latent==observer) stopifnot(tv_saved<=1e-12)
  stopifnot(abs(sum(q)-1)<=1e-12,all(q>=0))
  tbl[[length(tbl)+1L]] <- data.frame(sample=id,latent_law=latent,kernel=observer,
      TV_to_target=sum(abs(q-target))/2,TV_to_own_saved_observation=tv_saved,
      Hellinger_squared_to_target=sum((sqrt(q)-sqrt(target))^2)/2)
}
swap <- do.call(rbind,tbl)
emit(swap,file.path(out,'saved_distribution_kernel_swap.tsv'))
swap_max <- aggregate(TV_to_target~latent_law+kernel,swap,max)
swap_mean <- aggregate(TV_to_target~latent_law+kernel,swap,mean)
names(swap_max)[3]<-'max_TV'
names(swap_mean)[3]<-'mean_TV'
swap_s <- merge(swap_max,swap_mean)
emit(swap_s,file.path(out,'kernel_swap_summary.tsv'))
width <- read_tsv(file.path(p,'width_direction_diagnostics.tsv'))
width$maximum_kernel_column_TV_from_anchor <- vapply(width$sigma_kary,
  function(s)max(colSums(abs(kernel(s)-kv$fitted_B))/2),numeric(1))
emit(width,file.path(out,'width_direction_diagnostics.tsv'))

status <- if(is_lower) 'ELIGIBLE_FOR_SEPARATE_N20_REVIEW_AT_REGISTERED_D_RESOLUTION' else
  'ELIGIBLE_FOR_SEPARATE_N20_REVIEW_WITH_EXACT_MINIMUM_POLL'
emit(data.frame(status=status,prediction_origins=passing_origins,anchor_poll_complete=TRUE,
  selected_own_full_poll=!is_lower,absolute_D_improvement_bound=selected$score,registered_resolution=1e-5,
  N20_released=FALSE,profiles_released=FALSE,tumor_released=FALSE),file.path(out,'N20_review_gate.tsv'))
emit(data.frame(allocations=nrow(jobs),settled_core_seconds=sum(jobs$settled_seconds),
  settled_core_hours=sum(jobs$settled_seconds)/3600,hard_core_hours=2,
  live_reserved_core_seconds=sum(jobs$live_reserved_seconds),charged_forwards=sum(jobs$charged_calls),
  hard_forward_cap=155,technical_retries=sum(jobs$attempt==2)),file.path(out,'supplement_budget.tsv'))
cat('OFFLINE_B_LOCAL_REVIEW_PASS\n')
