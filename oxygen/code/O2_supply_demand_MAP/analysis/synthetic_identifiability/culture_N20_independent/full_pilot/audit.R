# Independent offline inspection of the completed pilot. No model forwards.
args <- commandArgs(TRUE)
r <- normalizePath(args[1L])
out <- normalizePath(args[2L])
tsv <- function(f)read.delim(file.path(r,f),check.names=FALSE,stringsAsFactors=FALSE)
emit <- function(x,f)write.table(x,file.path(out,f),sep='\t',row.names=FALSE,quote=TRUE)
source(file.path(r,'release/search_helpers.R'))
case <- readRDS(file.path(r,'inputs/public/N20.rds'))
quality <- tsv('reports/optimizer_quality.tsv')
verdict <- tsv('reports/origin_verdicts.tsv')
lock <- tsv('reports/selection_lock.tsv')
jobs <- tsv('accounting/jobs.tsv')
is_true <- function(x)tolower(as.character(x))=='true'
stopifnot(nrow(lock)==1L,case$kind=='N',is.null(case$known_sigma_growth),length(case$data$growth)==114L,
          length(case$data$kary)==12L,sum(lengths(case$data$kary))==220L,
          !anyDuplicated(quality[c('origin','key')]),all(is_true(quality$verified)))
loss <- function(pred,sigma) {
  g <- mean(dnorm(case$data$growth,pred$growth[names(case$data$growth)],sigma,log=TRUE))
  k <- mean(vapply(names(case$data$kary),function(id) {
    q <- pred$kary[[id]]
    stopifnot(length(q)==133L,all(is.finite(q)),all(q>=0),abs(sum(q)-1)<=1e-12)
    mean(log(pmax(q[match(case$data$kary[[id]],22:154)],1e-12)))
  },numeric(1)))
  c(Q=-g-k,growth=g,kary=k)
}
digest_key <- function(x)digest::digest(as.numeric(x),algo='sha256')
tv <- function(a,b)max(vapply(names(a),function(id)
  sum(abs(a[[id]]-b[[id]]))/2,numeric(1)))
rows <- poles <- list()
for (i in seq_len(nrow(quality))) {
  row <- quality[i,]
  task <- c(row$first_process_task,row$second_process_task)
  copies <- lapply(task,function(p) readRDS(file.path(r,'verification',p,paste0('point_',row$key,'.rds'))))
  scores <- list()
  for (j in seq_along(copies)) {
    current <- copies[[j]]; id <- jobs[jobs$task==task[j]&jobs$state=='COMPLETED'&jobs$exit_code=='0:0',,drop=FALSE]
    stopifnot(nrow(id)==1L,identical(as.character(id$job_id),current$slurm_job),
      current$task==task[j],current$candidate$origin==row$origin,
      current$candidate$key==row$key,digest_key(current$candidate$unit)==row$key,
      length(current$candidate$unit)==20L,length(current$candidate$params)==20L,
      identical(current$pred$selection,current$candidate$pred$selection))
    v <- loss(current$pred,current$candidate$params$sigma_growth)
    stopifnot(max(abs(v-c(row$Q,row$growth_log_mean,row$karyotype_log_mean)))<=1e-8)
    scores[[j]] <- v
  }
  x <- copies[[1L]]
  y <- copies[[2L]]
  stopifnot(x$slurm_job!=y$slurm_job,identical(x$candidate,y$candidate),
            identical(x$pred$selection,y$pred$selection),
            identical(names(x$pred$states),names(y$pred$states)),length(x$pred$states)>0L)
  max_state <- max(vapply(names(x$pred$states),function(id) {
    a<-x$pred$states[[id]];b<-y$pred$states[[id]]
    stopifnot(identical(dim(a),dim(b)))
    max(abs(a-b)/pmax(1,abs(a),abs(b)))
  },numeric(1)))
  max_growth <- max(abs(x$pred$growth-y$pred$growth))
  max_kary <- tv(x$pred$kary,y$pred$kary)
  stopifnot(max(max_state,max_growth,max_kary,tv(x$pred$latent_kary,y$pred$latent_kary))<=1e-10)
  rows[[i]] <- data.frame(origin=row$origin,key=row$key,job_a=x$slurm_job,job_b=y$slurm_job,
    Q_a=scores[[1L]][['Q']],Q_b=scores[[2L]][['Q']],max_growth_diff=max_growth,
    max_karyotype_TV=max_kary,max_scaled_state_diff=max_state,verified=TRUE)
}

for (origin in verdict$origin) {
  for (step in 1:2) {
    phase <- paste0('verification/',origin,'_verify',step)
    anchor <- readRDS(file.path(r,phase,'checked_representative.rds'))
    poll <- tsv(file.path(phase,'verification_polls.tsv'))
    variance <- tsv(file.path(phase,'variance_stationarity.tsv'))
    stopifnot(anchor$origin==origin,digest_key(anchor$unit)==anchor$key,variance$anchor_key==anchor$key,
              all(poll$anchor_key==anchor$key),!anyDuplicated(poll$key))
    dir <- r_directions(list(free=case$design$bounds))
    keys <- character()
    for (h in c(.001,.0001,.00002)) for (j in seq_len(nrow(dir))) for (sign in c(-1,1)) {
      u <- r_admissible(anchor$unit,sign*dir[j,],h)
      if (!is.null(u)) keys <- c(keys,digest_key(u))
    }
    keys <- unique(keys)
    complete <- identical(keys,poll$key)
    stopifnot(all(poll$key %in% keys))
    sigma <- subset(case$design$bounds,parameter=='sigma_growth')
    raw <- sqrt(mean((case$data$growth-anchor$pred$growth[names(case$data$growth)])^2))
    check_sigma <- min(sigma$upper,max(sigma$lower,raw))
    q <- loss(anchor$pred,check_sigma)[['Q']]
    max_gain <- max(0,anchor$score-poll$score)
    gain_from_sigma <- max(0,anchor$score-variance$replayed_Q)
    stopifnot(abs(variance$raw_residual_RMS-raw)<=1e-10,abs(variance$check_sigma-check_sigma)<=1e-10,
         max(abs(variance$analytic_Q-q),abs(variance$replayed_Q-q))<=1e-8,
         abs(variance$improvement-gain_from_sigma)<=1e-8)
    poles[[length(poles)+1L]] <- data.frame(origin=origin,phase=step,anchor_key=anchor$key,
       expected=length(keys),evaluated=nrow(poll),complete=complete,max_Q_gain=max_gain,
       sigma_gain=gain_from_sigma,local_pass=complete&&max(max_gain,gain_from_sigma)<=1e-5)
  }
}
local <- do.call(rbind,poles)
for (i in seq_len(nrow(quality))) {
  row <- quality[i,]
  x <- tail(local[local$origin==row$origin&local$anchor_key==row$key,,drop=FALSE],1)
  stopifnot(is_true(row$checked)==(nrow(x)==1),is_true(row$local_pass)==(nrow(x)==1&&x$local_pass))
}
lowest <- verdict[which.min(verdict$best_found_Q),]
best <- quality[which.min(quality$Q),]
selected_Q <- best$Q
agree <- sum(verdict$best_verified_Q-selected_Q<=1e-4)
qraw <- quality[quality$origin==lowest$origin&quality$key==lowest$best_found_key,,drop=FALSE]
exact <- nrow(qraw)==1L && is_true(qraw$local_pass)
stopifnot(agree==lock$independent_origins_within_1e_4,exact==is_true(lock$best_found_local_pass),
          lowest$origin==lock$best_found_origin,lowest$best_found_key==lock$best_found_key,
          best$origin==lock$selected_origin,best$key==lock$selected_key)
valid <- agree>=2L&&exact&&abs(lowest$best_found_Q-selected_Q)<=1e-8
stopifnot((lock$status=='N20_BOUNDED_OPTIMIZER_QUALIFIED')==valid)
emit(do.call(rbind,rows),'independent_replay_audit.tsv')
emit(local,'independent_exact_local_audit.tsv')
emit(data.frame(status=lock$status,agreeing_origins=agree,exact_lowest_local_pass=exact,
                lowest_found_Q=lowest$best_found_Q,selected_verified_Q=selected_Q),
     'independent_optimizer_audit.tsv')
cat('N20_REPLAY_LOCAL_AND_SELECTION_AUDIT_PASS\n')
