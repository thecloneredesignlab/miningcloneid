# Nine literal contract fixtures; no model simulation or culture data generation.
main <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  study <- args[[1]]; task <- args[[2]]
  stopifnot(task == "P03")
  out <- file.path(study,"verification",task)
  write_tsv <- function(x,name) write.table(x,file.path(out,name),sep="\t",quote=TRUE,row.names=FALSE,na="NA")
  write_tsv(data.frame(stage="before_allowlisted_expression_loading",task=task),"bootstrap_receipt.tsv")
  on.exit({
    ns <- sort(loadedNamespaces())
    write_tsv(data.frame(namespace=ns,version=vapply(ns,function(n)as.character(getNamespaceVersion(n)),character(1))),"loaded_namespaces.tsv")
  },add=TRUE)
  production <- file.path(study,"inputs","production")
  env <- new.env(parent=baseenv()); extracted <- list(); sources <- list()
  definitions <- list(
    c("o2_supply_demand_map_invitro_summary_utils.R","ivt_weighted_mean_kary_N"),
    c("o2_supply_demand_map_invitro_lineage_simulation_utils.R","ivt_select_segment_observation"),
    c("o2_supply_demand_map_invitro_lineage_simulation_utils.R","ivt_extract_passage_end_state"))
  for (definition in definitions) {
    expressions <- parse(file=file.path(production,definition[[1]]))
    matches <- Filter(function(e)is.call(e) && identical(e[[1]],as.name("<-")) && is.symbol(e[[2]]) && identical(as.character(e[[2]]),definition[[2]]),as.list(expressions))
    stopifnot(length(matches)==1L)
    eval(matches[[1]],envir=env)
    extracted[[length(extracted)+1L]] <- paste(deparse(matches[[1]],width.cutoff=500L),collapse="\n")
    sources[[length(sources)+1L]] <- data.frame(file=definition[[1]],function_name=definition[[2]],loading="unchanged allowlisted expression only")
  }
  writeLines(unlist(extracted),file.path(out,"extracted_function_source.R"))
  write_tsv(do.call(rbind,sources),"loaded_function_contract.tsv")
  select <- get("ivt_select_segment_observation",env)
  extract <- get("ivt_extract_passage_end_state",env)
  grid <- c(44,88)
  make_literal <- function(counts,states) list(Ntot_live_obs=counts,live_state_obs=matrix(states,ncol=2,byrow=TRUE))
  base <- make_literal(c(2,4,6),c(1,1,2,2,2,4))
  tests <- list()
  check <- function(id,condition,expectation) {
    tests[[length(tests)+1L]] <<- data.frame(fixture=id,passed=isTRUE(condition),independent_literal_expectation=expectation)
    write_tsv(do.call(rbind,tests),"passage_contract_tests.tsv")
    if(!isTRUE(condition)) stop("Passage contract fixture failed: ",id)
  }
  r <- select(base,1,grid,target_live_cells=2)
  check("F01_positive_day_subset",identical(r$selected_index,2L) && r$selected_day==1,"day0 is excluded when positive-day candidates exist")
  r <- select(base,1,grid,target_live_cells=5)
  check("F02_nearest_target_first_index_tie",identical(r$selected_index,2L),"counts4 and6 tie around5; index2 wins")
  r <- select(base,1,grid,target_live_cells=NA_real_)
  check("F03_missing_target_last_candidate",identical(r$selected_index,3L),"last positive-day index3 wins")
  fixture <- make_literal(c(1,3,7,5),c(1,0,1,2,3,4,1,4))
  r <- extract(fixture,4,grid,target_live_cells=3)
  check("F04_inoculum_smallest_excess_fallback",identical(r$selected_index,4L) && isTRUE(r$passage_executed) && r$available_cells==5,"target count3 is insufficient; eligible5 beats7 for inoculum4")
  fixture <- make_literal(c(1,2,3),c(1,0,1,1,1,2))
  r <- extract(fixture,4,grid,target_live_cells=2)
  check("F05_protocol_infeasible",identical(r$passage_executed,FALSE) && is.null(r$reseeded_state) && is.na(r$selected_index) && grepl("required_inoculum_not_reached",r$passage_failure_reason,fixed=TRUE),"no candidate supplies4; no artificial rescaling passage")
  r <- extract(base,3,grid,target_live_cells=6)
  check("F06_fraction_preserving_reseed",isTRUE(all.equal(r$reseeded_state,c(1,2))) && r$boundary_scale==0.5 && sum(r$reseeded_state)==3,"selected state(2,4) becomes(1,2), preserving fractions")
  fixture <- base;fixture$live_state_obs <- fixture$live_state_obs[1:2,,drop=FALSE]
  error <- tryCatch({select(fixture,1,grid);NULL},error=function(e)e)
  check("F07_incompatible_shapes_rejected",inherits(error,"error") && grepl("compatible daily trajectories",conditionMessage(error),fixed=TRUE),"3 totals versus2 state rows must fail")
  error <- tryCatch({select(base,0,grid);NULL},error=function(e)e)
  check("F08_invalid_inoculum_rejected",inherits(error,"error") && grepl("reseed_live_cells must be positive",conditionMessage(error),fixed=TRUE),"zero requested inoculum must fail")
  fixture <- make_literal(c(1,2,5,5),c(1,0,1,1,5,0,0,5))
  r <- extract(fixture,4,grid,target_live_cells=2)
  check("F09_fallback_tie_state_bookkeeping",identical(r$selected_index,3L) && r$selected_day==2 && r$selected_live_cells==5 && r$predicted_mean_kary_N==44 && isTRUE(all.equal(r$reseeded_state,c(4,0))),"equal eligible totals5 use index3 and its state, not index4")
  stopifnot(length(tests)==9L)
  write_tsv(data.frame(status="PASS",fixtures=length(tests),unchanged_functions=length(definitions),simulator_calls=0,kernel_calls=0,generated_culture_datasets=0),"passage_gate.tsv")
  write_tsv(data.frame(status="PASS",task=task,phase="pure_production_passage_contract"),"application_result.tsv")
  cat("NINE_PURE_PASSAGE_CONTRACT_FIXTURES_PASS\n")
}
main()
