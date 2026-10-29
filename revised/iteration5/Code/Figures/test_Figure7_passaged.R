# Focused unit and independent numerical checks for the Figure7-only workflow.
f7p_unit_checks <- function(native,input,run) {
  checks <- list()
  add <- function(name,passed) {
    checks[[length(checks)+1L]] <<- data.table::data.table(check=name,passed=isTRUE(passed))
    if(!isTRUE(passed))stop("Figure7 passage unit check: ",name)
  }
  toy <- list(Ntot_live_obs=c(10,3,2),live_state_obs=rbind(c(10,0),c(2,1),c(1,1)))
  set.seed(20261024);before <- .Random.seed
  p <- native$ivt_extract_passage_end_state(toy,100,c(44,45),50,0:2)
  add("native_positive_day_nearest_target_not_threshold_trigger",p$selected_day==1)
  add("insufficient_stock_does_not_stop_or_cap_next_inoculum",
    p$selected_live_cells<100&&abs(sum(p$reseeded_state)-100)<1e-12)
  add("native_selection_consumes_no_rng",identical(before,.Random.seed))
  # Selection must see the full window, even if only one elapsed day remains.
  full <- native$ivt_extract_passage_end_state(list(Ntot_live_obs=c(1,2,6,10),
    live_state_obs=cbind(c(1,2,6,10),0)),100,c(44,45),10,0:3)
  add("final_window_selection_precedes_truncation",full$selected_day==3&&min(full$selected_day,1)==1)
  zero <- native$ivt_extract_passage_end_state(list(Ntot_live_obs=c(0,0),
    live_state_obs=matrix(0,2,2)),100,c(44,45),1,0:1)
  add("native_zero_mass_fallback_preserved",identical(zero$selected_frac,c(.5,.5)))
  protocol_rows <- list()
  for(cohort in c("2N","4N"))for(block in c(3L,5L,10L)) {
    p <- f7p_protocol(input$fit$best_components$invitro[[paste0("run_",cohort)]]$adapter,block)
    for(route in names(p$routes)) {
      z <- p$inventory[match(p$templates[[route]],p$inventory$segment_id)]
      add(paste(cohort,block,route,"raw_complete",sep="_"),nrow(z)==block&&all(z$complete))
      path_complete <- p$routes[[route]][p$inventory$complete[match(p$routes[[route]],p$inventory$segment_id)]]
      add(paste(cohort,block,route,"last_complete_original_order",sep="_"),
        identical(p$templates[[route]],tail(path_complete,block)))
      add(paste(cohort,block,route,"short_reoxygenation_not_long_term_leaf",sep="_"),
        !tail(p$routes[[route]],1)%in%p$short_terminal_ids)
      protocol_rows[[length(protocol_rows)+1L]] <- cbind(cohort=cohort,block=block,route=route,z)
    }
  }
  # Independent absolute-count Euler, using sparse native operators, not the
  # normalized C++ counter update. One full day for each blocked/unblocked case.
  independent <- list()
  for(beta in c(0,.5,1)) {
    op <- f7p_operator(input,.5,.3,7L,beta,66L);x <- fp_pure_initial(op$l,4)
    n0 <- 1e6;xx <- x*n0;dd <- xx*0;uu <- xx;H <- L <- 0
    sim <- f7p_segment_cpp(op$G,op$KH,op$GU,op$l$lambda,op$l$mu,op$c$hl,op$c$ll,
      op$loss$upper,op$loss$lower,op$l$ngrid/22,op$c$low,x,x*0,x,1L,.05,n0,TRUE,1e12)
    for(step in 1:20) {
      cf <- max(0,1-sum(xx)/1e12)
      H <- H+.05*cf*sum(op$c$hl*xx);L <- L+.05*cf*sum(op$c$ll*xx)
      dx <- cf*as.numeric(op$G%*%xx)-op$l$mu*xx
      dD <- cf*as.numeric(op$KH%*%xx)-(cf*op$l$lambda+op$l$mu)*dd
      dU <- cf*as.numeric((op$G-op$KH)%*%uu)-op$l$mu*uu
      xx <- xx+.05*dx;dd <- dd+.05*dD;uu <- uu+.05*dU
    }
    mass <- sum(xx)
    err <- c(X=max(abs(xx/mass-sim$X[,2])),D=max(abs(dd/mass-sim$D[,2])),
      U=max(abs(uu/mass-sim$U[,2])),log_count=abs(log(mass)-sim$metrics[2,"log_cells"]),
      HL=abs(H/n0-sim$metrics[2,"HL_per_start_cell"]),
      LL=abs(L/n0-sim$metrics[2,"LL_per_start_cell"]))
    independent[[length(independent)+1L]] <- data.table::data.table(beta=beta,
      quantity=names(err),absolute_error=as.numeric(err),passed=err<1e-10)
    add(paste0("independent_absolute_Euler_beta_",beta),all(err<1e-10))
  }
  f7p_write(data.table::rbindlist(checks),file.path(run,"results/unit_checks.tsv"))
  f7p_write(data.table::rbindlist(protocol_rows),file.path(run,"results/template_inventory.tsv"))
  f7p_write(data.table::rbindlist(independent),file.path(run,"results/independent_absolute_Euler.tsv"))
  invisible(TRUE)
}
f7p_native_checkpoint_checks <- function(op,native,routes) {
  stopifnot(op$meta$beta==0)
  cfg <- op$pr$config;rp <- op$pr$run_params
  core <- build_model_core(cfg=cfg);vol <- cell_volume_mm3_by_N(core$grid_pre,rp,cfg)
  rows <- list()
  # Original kernel at true root and after >500 adopted days; sample control/O1.
  for(route in c("control","O1"))for(cp in routes[[route]]$checkpoints) {
    ref <- native$ivt_run_segment_fixed_o2(cp$segment,cfg,rp,core,vol,
      cp$initial$x*cp$initial$n0,cp$initial$n0)
    rf <- ref$sim$live_state_obs/rowSums(ref$sim$live_state_obs)
    le <- f7p_distance(log(ref$sim$Ntot_live_obs),cp$metrics$log_cells)
    fe <- f7p_distance(as.numeric(t(rf)),as.numeric(cp$X))
    rs <- native$ivt_extract_passage_end_state(ref$sim,cp$initial$n0,op$l$ngrid,
      cp$row$target_cells,cp$segment$obs_days_local)
    rows[[length(rows)+1L]] <- cbind(op$meta,route=route,elapsed_start=cp$row$elapsed_start,
      segment_id=cp$segment$segment_id,max_log_count_error=le,max_state_fraction_error=fe,
      selection_identical=rs$selected_day==cp$selection$selected_day,
      passed=le<1e-8&&fe<1e-10&&rs$selected_day==cp$selection$selected_day)
  };data.table::rbindlist(rows)
}
