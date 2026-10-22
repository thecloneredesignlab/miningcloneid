# Figure8: full fitted-seed finite-time predictions, not steady-state estimates.
f8_label <- function(Q)ifelse(Q==0,"2N",ifelse(is.infinite(Q),"4N",paste0("2N + 4N (1:",Q,")")))
f8_tasks <- function(a) {
  rows <- list();ep <- fp_endpoint_spec(a)
  for(ctx in c("in vivo","in vitro"))for(mode in c("fitted","0.01","0.3")) {
    oxygen <- if(a$profile=="full")seq(0,if(ctx=="in vitro"&&mode=="fitted")20 else 5,by=a$cfg$oxygen_step) else
      a$cfg$representative$oxygen
    if(ctx=="in vitro")oxygen <- sort(unique(c(oxygen,20)))
    for(fam in c("C01","C02")) {
      seeds <- sort(unique(ep[pair_label==fam]$seed_number))
      blocks <- split(seeds,ceiling(seq_along(seeds)/a$cfg$seeds_per_task))
      for(o2 in oxygen)for(j in seq_along(blocks))rows[[length(rows)+1L]] <- data.table::data.table(
        model_context=ctx,pair_label=fam,O2_pct=o2,p_setting=mode,
        condition_mode=if(mode=="fitted")"fitted_p_misseg" else "fixed_p_misseg",
        p_misseg=if(mode=="fitted")NA_real_ else as.numeric(mode),seed_block=j,
        seed_numbers=paste(blocks[[j]],collapse=","))
    }
  }
  z <- data.table::rbindlist(rows);z[,task_id:=seq_len(.N)];z
}
f8_vivo <- function(pr,o2,p,days,qs) {
  l <- fp_ledger(pr,o2,p);shift <- max(Re(eigen(l$M,only.values=TRUE)$values))
  x <- vapply(qs,function(Q)fp_initial(pr,l,Q,FALSE),numeric(length(l$ngrid)))
  logmass <- rep(0,length(qs));rows <- list();previous <- 0;cache <- list()
  for(day in days) {
    delta <- day-previous;key <- as.character(delta)
    if(is.null(cache[[key]]))cache[[key]] <- as.matrix(Matrix::expm(delta*(l$M-diag(shift,nrow(l$M)))))
    x <- cache[[key]]%*%x
    if(any(!is.finite(x))||min(x)< -1e-10)stop("Invalid finite-time in vivo state")
    x <- pmax(x,0);mass <- colSums(x)
    if(any(mass<=0))stop("Finite-time in vivo state has no ploidy estimand")
    logmass <- logmass+log(mass)+shift*delta;x <- sweep(x,2,mass,"/")
    rows[[length(rows)+1L]] <- data.table::data.table(day=day,Q=qs,initial_condition=f8_label(qs),
      mean_ploidy=colSums(x*(l$ngrid/l$N_unit)),stochastic_mcse=0,executed_repeats=1L,
      n_stochastic_repeats=1L,deterministic_no_passage=TRUE,log_population_factor=logmass,
      negative_state_clipping_count=0L,passage_count=0L)
    previous <- day
  }
  data.table::rbindlist(rows)
}
f8_simulate_task <- function(a,b,t,dest) {
  seeds <- as.integer(strsplit(t$seed_numbers,",",fixed=TRUE)[[1]])
  ep <- b$endpoints[pair_label==t$pair_label&model_context==t$model_context&seed_number%in%seeds]
  if(nrow(ep)!=length(seeds))stop("Task seed block does not match prepared inputs")
  native <- t$p_setting=="fitted";culture <- t$model_context=="in vitro"
  qs <- if(culture&&native)c(0,1,Inf) else if(culture&&t$O2_pct==20)0 else c(0,1:5,Inf)
  rows <- repeats <- list();done <- integer();checkpoint <- sub(".rds$","_checkpoint.rds",dest)
  if(file.exists(checkpoint)) {
    old <- fp_read(checkpoint);if(!identical(old$fingerprint,b$fingerprint))stop("Incompatible seed checkpoint")
    rows <- old$rows;repeats <- old$repeats;done <- old$done
  }
  for(i in setdiff(seq_len(nrow(ep)),done)) {
    e <- ep[i];pr <- b$prepared[[e$prepared_index]]
    p <- if(native)pr$native_parameters[["p_misseg"]] else t$p_misseg
    meta <- data.table::data.table(task_id=t$task_id,pair_label=e$pair_label,seed_number=e$seed_number,
      model_context=t$model_context,O2_pct=t$O2_pct,p_setting=t$p_setting,
      condition_mode=t$condition_mode,p_misseg=t$p_misseg,p_misseg_used=p,
      p_misseg_fitted=pr$native_parameters[["p_misseg"]],n_native_parameters=length(pr$native_parameters),
      complete_input_signature=e$complete_input_signature)
    if(!culture) {
      result <- f8_vivo(pr,t$O2_pct,p,if(native)b$config$days else b$config$main_day,qs)
      # Extra mixtures need only the terminal day; three main starts also retain
      # all intermediate snapshots of this same finite-time trajectory.
      result <- result[day==b$config$main_day|Q%in%c(0,1,Inf)]
      rows[[length(rows)+1L]] <- cbind(meta[rep(1L,nrow(result))],result)
    } else {
      l <- fp_ledger(pr,t$O2_pct,p)
      for(Q in qs) {
        x0 <- fp_initial(pr,l,Q,TRUE)
        o <- fp_culture(b,pr,l,x0,e$pair_label,e$seed_number,t$O2_pct,p,Q,
          daily=native,events=a$profile=="representative")
        days <- if(native)b$config$days else b$config$main_day
        y <- o$mean_ploidy[,if(native)days+1L else 1L,drop=FALSE]
        result <- data.table::data.table(day=days,Q=Q,initial_condition=f8_label(Q),
          mean_ploidy=colMeans(y),stochastic_mcse=if(nrow(y)>1L)apply(y,2,sd)/sqrt(nrow(y)) else rep(0,ncol(y)),
          executed_repeats=o$executed_repeats,n_stochastic_repeats=b$replicates,
          deterministic_no_passage=o$deterministic_no_passage,log_population_factor=NA_real_,
          negative_state_clipping_count=o$negative_state_clipping_count,passage_count=sum(o$passage_count))
        rows[[length(rows)+1L]] <- cbind(meta[rep(1L,nrow(result))],result)
        repeats[[length(repeats)+1L]] <- list(seed_number=e$seed_number,Q=Q,days=days,ploidy=y,
          passage_count=o$passage_count,first_passage_day=o$first_passage_day,last_passage_day=o$last_passage_day,
          final_log_population=o$final_log_population,passage_events=o$passage_events,
          final_state=if(a$profile=="representative")o$final_state else NULL)
      }
    }
    done <- c(done,i)
    fp_save(list(rows=rows,repeats=repeats,done=done,fingerprint=b$fingerprint),checkpoint)
  }
  list(data=data.table::rbindlist(rows),repeats=repeats)
}
f8_aggregate <- function(a,b,tasks) {
  z <- data.table::rbindlist(lapply(tasks,`[[`,"data"))
  keys <- c("condition_mode","p_setting","p_misseg","model_context","O2_pct","day","Q","initial_condition")
  pooled <- fp_quantiles(z,keys);family <- fp_quantiles(z,c(keys,"pair_label"))
  terminal <- pooled[day==b$config$main_day];fterminal <- family[day==b$config$main_day]
  # A has no mixture series and uses the exact same anchor results as C/D.
  anchors <- terminal[Q==0&((model_context=="in vivo"&O2_pct==.5)|
    (model_context=="in vitro"&O2_pct%in%c(.5,20)))]
  anchors[,initial_ploidy:=2]
  write <- function(x,name)fp_write(x,file.path(a$run,"tables",name))
  write(anchors,"curve_summary.tsv")
  write(terminal[model_context=="in vivo"&p_setting=="fitted"&Q%in%c(0,1,Inf)&O2_pct<=5],"Figure8C_curve_summary.tsv")
  write(terminal[model_context=="in vitro"&p_setting=="fitted"&Q%in%c(0,1,Inf)&O2_pct<=5],"Figure8D_curve_summary.tsv")
  write(fterminal[model_context=="in vivo"&p_setting=="fitted"],"fitted_family_curve_summary.tsv")
  write(terminal[model_context=="in vivo"&p_setting=="0.01"],"supp_fig8_3_curve_summary.tsv")
  write(terminal[model_context=="in vivo"&p_setting=="0.3"],"supp_fig8_4_curve_summary.tsv")
  write(fterminal[model_context=="in vitro"&p_setting!="fitted"&O2_pct<=5],"supp_fig8_5_curve_summary.tsv")
  # Native-p snapshots are grouped directly from per-seed rows; family medians
  # or IQRs are never averaged to construct the pooled summaries.
  for(spec in list(c("in vitro","6","8"),c("in vivo","7","9"))) {
    keep <- pooled[model_context==spec[1]&p_setting=="fitted"&Q%in%c(0,1,Inf)]
    split <- family[model_context==spec[1]&p_setting=="fitted"&Q%in%c(0,1,Inf)]
    write(keep,paste0("supp_fig8_",spec[2],"_curve_summary.tsv"))
    write(split,paste0("supp_fig8_",spec[3],"_curve_summary.tsv"))
    name <- if(spec[1]=="in vitro")"invitro" else "invivo"
    fp_write(z[model_context==spec[1]&p_setting=="fitted"&Q%in%c(0,1,Inf)],
      file.path(a$run,"seed_data",paste0(name,"_seed_snapshots.tsv.gz")))
  }
  fp_write(z[day==b$config$main_day],file.path(a$run,"seed_data/seed_endpoint_predictions.tsv.gz"))
  fp_write(pooled,file.path(a$run,"tables/all_pooled_summaries.tsv"))
  fp_write(family,file.path(a$run,"tables/all_family_summaries.tsv"))
  # Derived display table is generated afresh, never copied from retained data.
  empirical <- data.table::fread(file.path(a$run,"tables/empirical_tissue_ploidy_summary.tsv"))
  model <- anchors[p_setting=="fitted"]
  model[,`:=`(site=paste0("Model: ",model_context,"; O2 = ",O2_pct,"%"),context=model_context,source="Model")]
  empirical[,`:=`(context=ifelse(cohort=="PCAWG","in vivo","in vitro"),source=cohort)]
  write(data.table::rbindlist(list(empirical[,.(site,context,source,median,q25,q75)],
    model[,.(site,context,source,median,q25,q75)])),"Figure8A_display_summary.tsv")
  invisible(z)
}
f8_validate <- function(a,b,check) {
  z <- data.table::fread(file.path(a$run,"seed_data/seed_endpoint_predictions.tsv.gz"))
  check("unique_seed_scenario_Q_endpoints",!anyDuplicated(z[,.(pair_label,seed_number,model_context,O2_pct,p_setting,Q)]))
  check("finite_live_weighted_mean_ploidy",all(is.finite(z$mean_ploidy)))
  check("no_negative_state_clipping",all(z$negative_state_clipping_count==0L))
  check("full_context_parameter_counts_in_predictions",all(z$n_native_parameters==ifelse(z$model_context=="in vitro",33L,26L)))
  check("native_p_not_replaced",all(z[p_setting=="fitted"]$p_misseg_used==z[p_setting=="fitted"]$p_misseg_fitted))
  check("fixed_native_p_is_explicit",all(abs(z[p_setting!="fitted"]$p_misseg_used-as.numeric(z[p_setting!="fitted"]$p_setting))<1e-12))
  keys <- c("model_context","p_setting","O2_pct","Q")
  counts <- z[,.(n=.N,unique_ids=data.table::uniqueN(paste(pair_label,seed_number))),by=keys]
  expected <- sum(b$endpoints$model_context=="in vivo")
  check("no_optimizer_seed_deduplication",all(counts$n==expected&counts$unique_ids==expected))
  for(ctx in c("in vivo","in vitro")) {
    name <- if(ctx=="in vivo")"invivo" else "invitro"
    raw <- data.table::fread(file.path(a$run,"seed_data",paste0(name,"_seed_snapshots.tsv.gz")))
    check(paste0(ctx,"_unique_snapshots"),!anyDuplicated(raw[,.(pair_label,seed_number,O2_pct,Q,day)]))
    nday <- raw[,.(n_days=.N),by=.(pair_label,seed_number,O2_pct,Q)]
    check(paste0(ctx,"_all_ten_snapshots"),all(nday$n_days==length(b$config$days)))
    cmp <- fp_quantiles(raw,c("model_context","p_setting","O2_pct","Q","day"))
    saved <- data.table::fread(file.path(a$run,"tables",paste0("supp_fig8_",if(ctx=="in vitro")6 else 7,"_curve_summary.tsv")))
    joined <- merge(cmp,saved,by=c("model_context","p_setting","O2_pct","Q","day"))
    error <- max(abs(joined$median.x-joined$median.y),abs(joined$q25.x-joined$q25.y),abs(joined$q75.x-joined$q75.y))
    check(paste0(ctx,"_independent_pooled_type7_quantiles"),error<1e-12,error)
    terminal <- raw[day==b$config$main_day&O2_pct<=5]
    summary <- fp_quantiles(terminal,c("model_context","O2_pct","Q"))
    main <- data.table::fread(file.path(a$run,"tables",paste0("Figure8",if(ctx=="in vivo")"C" else "D","_curve_summary.tsv")))
    joined <- merge(summary,main,by=c("model_context","O2_pct","Q"))
    error <- max(abs(joined$median.x-joined$median.y),abs(joined$q25.x-joined$q25.y),abs(joined$q75.x-joined$q75.y))
    check(paste0(ctx,"_main_equals_terminal_snapshot"),nrow(joined)==nrow(main)&&error<1e-12,error)
  }
  empirical <- data.table::fread(file.path(a$run,"tables/empirical_tissue_ploidy_summary.tsv"))
  author <- data.table::fread(file.path(a$run,"tables/author_figure5_points.tsv"))
  check("empirical_36_tissue_groups",nrow(empirical)==36L)
  check("author_figure5_19_original_rows",nrow(author)==19L)
}
