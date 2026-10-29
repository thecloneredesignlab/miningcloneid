#!/usr/bin/env Rscript
# Explicit numerical/unit acceptance; never changes formal data or images.
.script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))
.code <- dirname(.script)
for(.name in c("figure78_pipeline.R","figure78_model.R","figure7_simulation.R","figure7_current_pipeline.R"))
  source(file.path(.code,"util/analysis",.name),local=globalenv())
args <- fp_args();a <- fp_init(7L,file.path(.code,"run_Figure7.R"),args)
if(is.null(args$output_root))stop("An isolated --output-root is required")
output <- normalizePath(args$output_root,mustWork=FALSE)
stopifnot(startsWith(output,paste0(a$workspace,"/audit/")))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
checks <- list();check <- function(name,passed,value=NA_real_) {
  checks[[length(checks)+1L]] <<- data.table::data.table(test=name,value=value,passed=isTRUE(passed))
  if(!isTRUE(passed))stop("Figure7 formal-entry acceptance: ",name)
}
tasks <- f7c_tasks(a)
check("full_320_base_conditions",data.table::uniqueN(tasks$base_task_id)==320L)
check("660_manifest_workers",nrow(tasks)==660L)
check("170_culture_conditions_three_grids",nrow(tasks[model_context=="in vitro"])==510L&&
  all(tasks[model_context=="in vitro",.(n=.N),by=base_task_id]$n==3L))
check("150_unchanged_continuous_conditions",nrow(tasks[model_context=="in vivo"])==150L)
check("fixed_mechanistic_probabilities",setequal(tasks$p_misseg,c(.1,.3)))
check("pure_two_and_four_starts",setequal(tasks$initial_ploidy,c(2,4)))
check("only_approved_grid_triplet",setequal(tasks[model_context=="in vitro"]$maximum_ploidy,c(137,139,141)))
check("no_historical_plot_default",grepl('--analysis=instantaneous',paste(deparse(body(f7c_plot)),collapse=""),fixed=TRUE))
check("no_stochastic_culture_call",!grepl('fp_culture',paste(deparse(body(f7c_vitro_task)),collapse=""),fixed=TRUE))
b <- f7c_bundle(a)
runtime <- f7p_runtime(a$workspace,file.path(output,"runtime"))
input <- f7p_load_fit(a$workspace,"C01",a$fit_root)
source(file.path(.code,"test_Figure7_passaged.R"),local=globalenv())
f7p_unit_checks(runtime$native,input,output)
check("native_selection_protocol_and_independent_Euler_units",TRUE)
ids <- if(!is.null(args$task_ids))as.integer(strsplit(args$task_ids,",",fixed=TRUE)[[1]]) else integer()
for(id in ids) {
  t <- b$tasks[task_id==id];o <- fp_read(fp_task_path(a,id))
  check(paste0("task_",id,"_registered_fingerprint"),identical(o$fingerprint,b$fingerprint))
  check(paste0("task_",id,"_all_beta_QC"),all(vapply(o$outputs,function(v)all(v$qc$passed),TRUE)))
  if(t$model_context!="in vitro")next
  fi <- f7p_load_fit(a$workspace,t$pair_label,a$fit_root)
  protocol <- b$protocols[[paste(t$pair_label,t$initial_ploidy,sep="_")]]
  for(v in o$outputs) {
    op <- f7p_operator(fi,t$O2_pct,t$p_misseg,t$maximum_ploidy,v$metadata$q,t$threshold_N)
    for(path in v$route_paths) {
      r <- fp_read(path);key <- paste(id,v$metadata$q,r$meta$route,sep="_")
      H <- r$metrics$HL_per_start_cell;L <- r$metrics$LL_per_start_cell
      ratio <- ifelse(H+L>0,100*H/(H+L),NA_real_)
      check(paste0(key,"_daily_numerator_denominator"),f7p_distance(ratio,r$metrics$A)<1e-10)
      check(paste0(key,"_absolute_daughter_scale"),
        f7p_distance(r$metrics$expected_HL_daughters,H*exp(r$metrics$log_start_live_cells))<1e-6)
      # Independent native-kernel replay of the root and one late selected
      # passage. The pre-passage boundary composition is normalized explicitly;
      # at day>0 the saved count belongs to the previous segment, not the reseed.
      if(v$metadata$q==0) {
        events <- r$events
        at <- unique(c(1L,which(events$elapsed_start>=500)[1]))
        for(j in at[is.finite(at)]) {
          ev <- events[j];seg <- protocol$adapter$segments[[match(ev$template_source_id,protocol$inventory$segment_id)]]
          seg$oxygen_pct <- t$O2_pct
          x <- r$X[,ev$elapsed_start+1L];x <- x/sum(x)
          core <- build_model_core(cfg=op$pr$config)
          volume <- cell_volume_mm3_by_N(core$grid_pre,op$pr$run_params,op$pr$config)
          ref <- runtime$native$ivt_run_segment_fixed_o2(seg,op$pr$config,op$pr$run_params,core,volume,
            x*ev$initial_cells,ev$initial_cells)
          selected <- runtime$native$ivt_extract_passage_end_state(ref$sim,ev$initial_cells,op$l$ngrid,
            ev$target_cells,seg$obs_days_local)
          ix <- seq_len(ev$retained_days+1L)
          state <- t(ref$sim$live_state_obs[ix,,drop=FALSE]/ref$sim$Ntot_live_obs[ix])
          expected <- r$X[,ev$elapsed_start+ix,drop=FALSE]
          err <- max(abs(state-expected))
          check(paste0(key,"_native_checkpoint_",j),err<1e-9&&selected$selected_day==ev$selected_day,err)
        }
      }
    }
    A <- v$summary[metric=="flux_percent"];B <- v$summary[metric=="growth_percent"]
    check(paste0(id,"_",v$metadata$q,"_median_complement"),f7p_max(abs(A$median+B$median-100))<1e-10)
    check(paste0(id,"_",v$metadata$q,"_IQR_complement"),f7p_max(abs(A$q25+B$q75-100))<1e-10&&
      f7p_max(abs(A$q75+B$q25-100))<1e-10)
  }
}
fp_write(data.table::rbindlist(checks),file.path(output,"formal_entry_acceptance.tsv"))
message("PASS Figure7 formal-entry acceptance: ",length(checks)," checks plus native units")
