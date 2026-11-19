anchor_key <- '3691ca93d3c4d115dadf9afe0a699e76406149978680eb8056e646835f6bf716'
width_names <- c('init_sd_2N','init_sd_4N','sigma_kary')
call_limits <- c(admit=0L,poll=151L,replayA=1L,replayB=1L)
archive_jobs <- c('22051119', '22051124')

make_probes <- function(anchor, c) {
  d <- r_directions(c)
  u <- anchor$unit
  seen <- character()
  units <- list()
  rows <- list()
  for (mesh in c(.001,.0001,.00002)) for (k in seq_len(nrow(d))) for (sign in c(-1,1)) {
    y <- r_admissible(u, sign*d[k,], mesh)
    key <- if (is.null(y)) '' else sf_hash(y)
    status <- if (is.null(y)) 'ZERO_OR_BOUND' else if (key %in% seen) 'DUPLICATE' else 'RUN'
    if (status == 'RUN') {
      seen <- c(seen,key)
      units[[key]] <- y
    }
    rows[[length(rows)+1L]] <- data.frame(mesh=mesh,direction_id=k,
      direction=paste(colnames(d)[d[k,]!=0],collapse=' / '), sign=sign,
      requested_step=mesh,actual_step=if(is.null(y)) 0 else sqrt(sum((y-u)^2)),
      key=key,status=status)
  }
  independent <- poll_keys(anchor,c$free$parameter)
  stopifnot(length(rows)==150L, identical(independent, seen), !(anchor$key %in% seen))
  list(rows=do.call(rbind,rows), units=units)
}

load_records <- function(out, prior, units, study_sha) {
  results <- list()
  for (path in c(prior,out)) if (length(path) && dir.exists(file.path(path,'points'))) {
    files <- list.files(file.path(path,'points'),pattern='^[0-9a-f]{64}[.]rds$',full.names=TRUE)
    for (file in files) {
      key <- sub('[.]rds$','',basename(file))
      start <- readRDS(file.path(path,'starts',basename(file)))
      record <- readRDS(file)
      stopifnot(key %in% names(units),is.null(results[[key]]),identical(record$token,start),
                identical(start$study_sha,study_sha),identical(start$key,key),
                identical(start$unit,units[[key]]),identical(sf_hash(start$unit),key),
                identical(record$status,if(record$feasible) 'FEASIBLE' else 'PROTOCOL_INFEASIBLE'))
      results[[key]] <- record
    }
  }
  results
}

score_equals <- function(a,b) {
  fields <- c('D','growth_RMS','growth_max','kary_TV')
  all(abs(a$values[fields]-b$values[fields])<=1e-12) &&
    max(abs(a$tv-b$tv))<=1e-12 && identical(a$pass,b$pass)
}

run_review <- function() {
  args <- commandArgs(TRUE)
  study <- normalizePath(args[1L])
  task <- args[2L]
  stage <- sub('_retry$','',sub('^D_B_','',task))
  stopifnot(stage %in% names(call_limits))
  out <- file.path(study,'verification',task)
  for (f in c('production_bindings.R','interface.R','feasibility.R','search_helpers.R',
              'core.R','joint_numerics.R','constraints.R','closeout_polls.R'))
    source(file.path(study,'release',f))
  on.exit(sf_inventory(out), add=TRUE)
  rows <- list()
  check <- function(label,ok,detail='') {
    rows[[length(rows)+1L]] <<- data.frame(check=label,passed=isTRUE(ok),detail=detail)
    sf_emit(do.call(rbind,rows),file.path(out,'checks.tsv'))
    if (!isTRUE(ok)) stop('Check failed: ',label,call.=FALSE)
  }
  job <- Sys.getenv('SLURM_JOB_ID')
  check('native_allocation_receipt',grepl('^[0-9]+$',job))
  env <- sf_backend(study,out,native=stage!='admit')
  input <- d_read(study,'D_B')
  ctx <- d_context(env,input$case,input$contract,input$id)
  study_sha <- Sys.getenv('D19_REVIEW_SOURCE_SHA')
  check('supplement_manifest_identity',grepl('^[0-9a-f]{64}$',study_sha))
  archived <- lapply(c('replay_a.rds','replay_b.rds'),function(f)
    readRDS(file.path(study,'verification/D_B_anchor',f)))
  anchor <- archived[[2L]]$candidate
  check('exact_anchor_and_archive_processes',identical(anchor$key,anchor_key) &&
    identical(sf_hash(anchor$unit),anchor_key) && identical(anchor$origin,'D_B') &&
    identical(anchor,archived[[1L]]$candidate) &&
    identical(vapply(archived,`[[`,character(1),'slurm_job'),archive_jobs) &&
    identical(anchor$params,d_decode(ctx,anchor$unit)[cg_roles]) &&
    identical(anchor$identity,ctx$identity))
  archived_metrics <- exact_metrics(ctx,list(pred=archived[[2L]]$pred,loss=archived[[2L]]$loss))
  for (i in seq_along(archived)) {
    x <- archived[[i]]
    f_compare(archived[[2L]]$pred,x$pred,check,paste0('archived_',i))
    independent <- d_independent_loss(ctx,x$pred)
    check(paste0('archive_scores_',i),max(abs(independent-c(x$loss$growth,x$loss$kary,x$loss$value)))<=1e-12 &&
      abs(x$loss$value-anchor$score)<=1e-12 && x$loss$value>=0 && joint_pass(x$loss))
  }
  probes <- make_probes(anchor,ctx)
  if (stage=='admit') {
    check('registered_probe_set_reconstructed',length(probes$units)<=150L && length(probes$units)>0L)
    sf_emit(probes$rows,file.path(out,'poll_coverage_expected.tsv'))
    sf_emit(data.frame(key=anchor_key,D=anchor$score,archived_replay_a=archive_jobs[1],
      archived_replay_b=archive_jobs[2],expected_unique_probes=length(probes$units),status='PASS'),
      file.path(out,'anchor_admission.tsv'))
  }

  completed <- function(token,z,p,grid=NULL) {
    if (!z$feasible) {
      stopifnot(length(z$lineages)>0L)
      broken <- Filter(function(x)!x$feasible,z$lineages)
      stopifnot(length(broken)>0L,all(vapply(broken,function(x)length(x$error)>0L,logical(1))))
      return(list(token=token,feasible=FALSE,status='PROTOCOL_INFEASIBLE',score=Inf,
                  widths=p[width_names],error=lapply(broken,`[[`,'error')))
    }
    metrics <- exact_metrics(ctx,z)
    independent <- d_independent_loss(ctx,z$pred)
    stopifnot(z$loss$value>=0,max(abs(independent-c(z$loss$growth,z$loss$kary,z$loss$value)))<=1e-12,
              identical(metrics$pass,joint_pass(z$loss)))
    predicted <- z$pred
    selection <- if(is.null(z$lineages)) NULL else f_ledger(z,token$key,ctx$contract)
    reseed_modes <- if(is.null(z$lineages)) character() else unlist(lapply(z$lineages,function(a)
      vapply(a$edges,function(e) e$selection$reseed_mode,character(1))),use.names=FALSE)
    scored_days <- predicted$selection
    a <- archived[[2L]]$pred
    tv_between <- function(x,y) vapply(names(y),function(k) sum(abs(x[[k]]-y[[k]]))/2,numeric(1))
    record <- list(token=token,feasible=TRUE,status='FEASIBLE',score=z$loss$value,
                   components=c(D_G=z$loss$growth,D_K=z$loss$kary),metrics=metrics,
                   widths=p[width_names],branch=if(is.null(grid)) '' else joint_branch(z,joint_pruning(ctx,p,grid),p),
                   scored_days=scored_days,protocol=selection,reseed_modes=reseed_modes,allocation=job,process=Sys.getpid(),
                   observed_to_anchor=tv_between(predicted$kary,a$kary),
                   latent_to_anchor=tv_between(predicted$latent_kary,a$latent_kary),
                   bound_hits=paste(ctx$free$parameter[token$unit==0|token$unit==1],collapse=';'),
                   inactive=paste(r_activity(p),collapse=';'))
    if (stage!='poll' || z$loss$value<anchor$score || token$key==anchor_key) {
      record$prediction <- predicted
      record$loss <- z$loss
    }
    record
  }
  evaluate <- function(u,prior_cases=0L,reference=FALSE,grid=NULL) {
    key <- sf_hash(u)
    path <- file.path(out,'starts',paste0(key,'.rds'))
    stopifnot(!file.exists(path))
    count <- length(list.files(file.path(out,'starts'),pattern='^[0-9a-f]{64}[.]rds$'))
    if (prior_cases+count>=call_limits[[stage]]+(if(endsWith(task,'_retry')) 2L else 0L))
      stop('Forward-call allowance exhausted',call.=FALSE)
    token <- list(key=key,unit=u,study_sha=study_sha,job=job,stage=stage)
    sf_rds(token,path)
    p <- d_decode(ctx,u)
    z <- d_physics(ctx,p,reference=reference)
    record <- completed(token,z,p,grid)
    sf_rds(record,file.path(out,'points',basename(path)))
    record
  }

  original <- file.path(study,'verification',sub('_retry$','',task))
  prior <- if(endsWith(task,'_retry')) original else character()
  prior_calls <- if(length(prior)) length(list.files(file.path(prior,'starts'),pattern='^[0-9a-f]{64}[.]rds$')) else 0L

  if (stage=='poll') {
    units <- c(setNames(list(anchor$unit),anchor_key),probes$units)
    results <- load_records(out,prior,units,study_sha)
    grid <- joint_mask_grid(ctx)
    if (is.null(results[[anchor_key]]))
      results[[anchor_key]] <- evaluate(units[[anchor_key]],prior_calls,grid=grid)
    a <- results[[anchor_key]]
    check('anchor_feasible',a$feasible)
    f_compare(a$prediction,archived[[2L]]$pred,check,'current_anchor')
    check('anchor_component_and_gate_agreement',score_equals(a$metrics,archived_metrics) &&
      max(abs(a$components-c(archived[[2L]]$loss$growth,archived[[2L]]$loss$kary)))<=1e-12 &&
      abs(a$score-anchor$score)<=1e-12)
    for (key in names(probes$units)) if(is.null(results[[key]]))
      results[[key]] <- evaluate(units[[key]],prior_calls,grid=grid)
    check('complete_and_exact_probe_set',setequal(names(results),names(units)))
    candidates <- probes$rows[probes$rows$status=='RUN',,drop=FALSE]
    compare_protocol <- function(x,y,column) {
      stopifnot(identical(x$segment_id,y$segment_id),identical(x$cohort,y$cohort))
      lx <- x[[column]]; ly <- y[[column]]
      sum(xor(is.na(lx),is.na(ly)) | (!is.na(lx)&!is.na(ly)&lx!=ly))
    }
    summary <- lapply(seq_len(nrow(candidates)),function(i) {
      row <- candidates[i,]; x <- results[[row$key]]; ok <- x$feasible
      g <- if(ok) x$metrics$values else rep(NA_real_,4L)
      data.frame(row,feasible=ok,D=x$score,D_G=if(ok) x$components[1L] else NA_real_,
        D_K=if(ok) x$components[2L] else NA_real_,growth_RMS=g[2L],growth_max=g[3L],max_TV=g[4L],
        all_prediction_gates=if(ok) x$metrics$pass else FALSE,
        improvement=a$score-x$score,branch_same=ok&&identical(x$branch,a$branch),
        scored_day_changes=if(ok) compare_protocol(x$scored_days,a$scored_days,'day') else NA_integer_,
        parent_day_changes=if(ok) compare_protocol(x$scored_days,a$scored_days,'parent_day') else NA_integer_,
        reseeding_mode_changes=if(ok) sum(x$reseed_modes!=a$reseed_modes) else NA_integer_,
        bound_hits=if(ok) x$bound_hits else '',inactive=if(ok) x$inactive else '',
        init_sd_2N=x$widths$init_sd_2N,init_sd_4N=x$widths$init_sd_4N,sigma_kary=x$widths$sigma_kary)
    })
    tab <- do.call(rbind,summary)
    single_and_paired <- tab[tab$direction_id %in% c(match(width_names,ctx$free$parameter),24:25),,drop=FALSE]
    for (name in width_names) {
      column <- match(name,ctx$free$parameter)
      single_and_paired[[paste0('unit_change_',name)]] <- vapply(single_and_paired$key,
        function(key) results[[key]]$token$unit[column]-anchor$unit[column],numeric(1))
    }
    sf_emit(probes$rows,file.path(out,'poll_coverage.tsv'))
    full_tab <- cbind(tab,do.call(rbind,lapply(tab$key,function(key) {
      u <- as.list(results[[key]]$token$unit)
      names(u) <- paste0('u_',ctx$free$parameter)
      as.data.frame(u)
    })))
    sf_emit(full_tab,file.path(out,'poll_candidates.tsv'))
    sf_emit(single_and_paired,file.path(out,'width_direction_diagnostics.tsv'))
    tv <- lapply(results,function(x) if(x$feasible)
      data.frame(key=x$token$key,sample=names(x$metrics$tv),TV_to_target=unname(x$metrics$tv),
                 observed_TV_to_anchor=unname(x$observed_to_anchor),latent_TV_to_anchor=unname(x$latent_to_anchor)))
    sf_emit(do.call(rbind,tv),file.path(out,'per_sample_tv.tsv'))
    best <- if(min(tab$D)<a$score && min(tab$D)<anchor$score) results[[tab$key[which.min(tab$D)]]] else a
    gain <- max(0,a$score-min(tab$D))
    check('registered_resolution',gain<=1e-5&&gain<=a$score&&best$score>=0)
    check('selected_forward_prediction_gates',best$metrics$pass)
    sf_rds(best,file.path(out,'selected_point.rds'))
    sf_emit(data.frame(anchor_key=anchor_key,anchor_D=a$score,expected=nrow(candidates),
      evaluated=nrow(candidates),recognized_infeasible=sum(!tab$feasible),complete=TRUE,
      maximum_observed_improvement=gain,absolute_improvement_bound=a$score,
      local_resolution=1e-5,selected_key=best$token$key,selected_D=best$score,
      lower_found=best$token$key!=anchor_key,selected_exact_anchor_checked=best$token$key==anchor_key),
      file.path(out,'local_resolution_audit.tsv'))
  }

  if (stage %in% c('replayA','replayB')) {
    poll_folder <- file.path(study,'verification/D_B_poll')
    if(file.exists(file.path(study,'verification/D_B_poll_retry/selected_point.rds')))
      poll_folder <- file.path(study,'verification/D_B_poll_retry')
    selected <- readRDS(file.path(poll_folder,'selected_point.rds'))
    check('sealed_lower_point',selected$token$key!=anchor_key && selected$score<anchor$score &&
      identical(selected$token$study_sha,study_sha) && identical(sf_hash(selected$token$unit),selected$token$key))
    units <- setNames(list(selected$token$unit),selected$token$key)
    results <- load_records(out,character(),units,study_sha)
    item <- results[[selected$token$key]]
    if(is.null(item)) item <- evaluate(selected$token$unit,prior_calls,reference=stage=='replayB')
    check('fresh_lower_point_gate',item$feasible&&item$metrics$pass&&item$score<anchor$score)
    f_compare(selected$prediction,item$prediction,check,'selected_point')
    check('lower_point_component_identity',score_equals(selected$metrics,item$metrics) &&
      max(abs(selected$components-item$components))<=1e-12)
    sf_rds(item,file.path(out,'replayed_point.rds'))
    sf_emit(data.frame(key=item$token$key,D=item$score,job=item$allocation,R_process=item$process,
      path=if(stage=='replayA') 'instrumented' else 'independent',
      t(item$metrics$values[-1L]),strict_pass=item$metrics$pass),
      file.path(out,'selected_point_replays.tsv'))
  }
  sf_emit(data.frame(status='PASS',task=task,job=job),file.path(out,'application_result.tsv'))
}
if (sys.nframe()==0L) tryCatch(run_review(),error=function(e) {
  cat('D19_B_LOCAL_REVIEW_FAILED: ',conditionMessage(e),'\n',sep='')
  quit(status=1L,save='no')
})
