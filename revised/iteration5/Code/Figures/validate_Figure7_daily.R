#!/usr/bin/env Rscript
# Independent two-state numerical examples for the passive daily counters.
script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))
code <- dirname(script)
source(file.path(code,"util/analysis/figure78_pipeline.R"))
source(file.path(code,"util/analysis/figure7_simulation.R"))
suppressPackageStartupMessages({library(Rcpp);library(Matrix);library(data.table)})
args <- fp_args();root <- args$output_root
if(is.null(root))stop("Specify --output-root for isolated numerical acceptance")
dir.create(root,recursive=TRUE,showWarnings=FALSE)
Rcpp::sourceCpp(file.path(code,"util/analysis/figure7_daily.cpp"),cacheDir=file.path(root,"cache"))
checks <- list();check <- function(n,ok,v=NA_real_) {
  checks[[length(checks)+1L]] <<- data.frame(test=n,passed=isTRUE(ok),value=v)
  if(!isTRUE(ok))stop("FAIL: ",n)
}
# J already counts both ordinary daughters. High-mother H->L=0.4/day;
# low-mother L->L=2/day. Both total populations grow at net rate 1/day.
l <- list(ngrid=c(44,88),M=matrix(c(1,0,.4,.6),2),G=matrix(c(1,0,.4,.6),2),
  lambda=c(1,1),mu=c(0,0))
c <- list(KH=matrix(c(0,0,.4,0),2),hl=c(0,.4),ll=c(2,0))
initial <- c(0,1);empty <- matrix(numeric(),2,0)
events <- data.frame(day=integer(),repeat_id=integer())
op <- f7_direct_operator(l,c,FALSE)
X <- as.numeric(op$op$P%*%initial);U <- as.numeric(op$op$PU%*%initial)
ref <- list(final_state=matrix(X/sum(X),2),final_untagged_state=matrix(U/sum(X),2),passage_events=events,
  pre_passage_state=empty,post_passage_state=empty,pre_passage_untagged_state=empty,post_passage_untagged_state=empty)
v <- f7_daily_cpp(op$op,ref,initial,c(TRUE,FALSE),c(2,4),c$hl,c$ll,op$count_rows,
  FALSE,1,FALSE,0L,1e12,10,1L,TRUE,c(0L,1L))
H <- .4*expm1(.6)/.6;L <- 2*(expm1(1)-expm1(.6)/.6)
check("exact_vivo_HL_integral",abs(v$HL_per_start_live_cell[1,2]-H)<1e-12)
check("exact_vivo_LL_integral",abs(v$LL_per_start_live_cell[1,2]-L)<1e-12)
check("exact_vivo_complement",abs(v$high_parent_percent[1,2]+v$low_parent_percent[1,2]-100)<1e-12)
instant <- 100*.4*v$snapshot_X[2,2]/(.4*v$snapshot_X[2,2]+2*v$snapshot_X[1,2])
check("daily_integral_not_endpoint_rate_fraction",abs(instant-v$high_parent_percent[1,2])>1)
check("no_interval_at_day_zero",all(is.na(v$high_parent_percent[,1]))&&all(is.na(v$low_parent_percent[,1])))
check("exact_replay_state",v$final_state_error<1e-12)
# Independent hand-Euler example, dt=0.5, with a saved integer-day passage.
# X_abs: (0,1)->(.2,1.3)->(.56,1.69). Counts HL=.46, LL=.20.
# Passage chooses X=(.5,.5); it contributes no births of either source.
preX <- c(.56,1.69)/2.25;preU <- c(0,1.69)/2.25
ref <- list(final_state=matrix(c(.5,.5),2),final_untagged_state=matrix(c(0,.5),2),
  passage_events=data.frame(day=1L,repeat_id=1L),pre_passage_state=matrix(preX,2),
  post_passage_state=matrix(c(.5,.5),2),pre_passage_untagged_state=matrix(preU,2),
  post_passage_untagged_state=matrix(c(0,.5),2))
v <- f7_daily_cpp(op$op,ref,initial,c(TRUE,FALSE),c(2,4),c$hl,c$ll,op$count_rows,
  TRUE,.5,FALSE,0L,1e12,10,1L,TRUE,c(0L,1L))
check("culture_hand_HL_count",abs(v$HL_per_start_live_cell[1,2]-.46)<1e-12)
check("culture_hand_LL_count",abs(v$LL_per_start_live_cell[1,2]-.20)<1e-12)
check("passage_does_not_add_births",abs(v$high_parent_percent[1,2]-100*.46/.66)<1e-12)
check("passage_anchors_preserved",v$passage_error<1e-12&&v$final_state_error<1e-12)
# At a single low chromosome state, origin changes neither rate nor ploidy.
check("same_ploidy_same_division_rate",v$metrics[[4]][1,2]==v$metrics[[9]][1,2])
check("same_ploidy_same_low_birth_rate",v$metrics[[3]][1,2]==v$metrics[[8]][1,2])
check("direct_and_history_are_distinct",v$metrics[[7]][1,2]<v$metrics[[17]][1,2])
check("initial_empty_group_mean_NA",is.na(v$metrics[[10]][1,1]))
# Full blockade with a pure high start has no low production at all.
zero <- l;zero$M <- zero$G <- diag(c(1,1))
cz <- c;cz$KH[,] <- 0;cz$hl[] <- 0
oz <- f7_direct_operator(zero,cz,FALSE)
rz <- ref;rz$passage_events <- events
rz$pre_passage_state <- rz$post_passage_state <- empty
rz$pre_passage_untagged_state <- rz$post_passage_untagged_state <- empty
rz$final_state <- rz$final_untagged_state <- matrix(initial,2)
vz <- f7_daily_cpp(oz$op,rz,initial,c(TRUE,FALSE),c(2,4),cz$hl,cz$ll,oz$count_rows,
  FALSE,1,FALSE,0L,1e12,10,1L,FALSE,integer())
check("zero_production_both_shares_NA",all(is.na(vz$high_parent_percent))&&all(is.na(vz$low_parent_percent)))
check("no_artificial_counter_under_full_blockade",vz$HL_per_start_live_cell[1,2]==0&&vz$LL_per_start_live_cell[1,2]==0)
if(!is.null(args$data_root)) {
  a <- fp_init(7L,file.path(code,"run_Figure7.R"),args);b <- fp_bundle(a)
  numerical <- data.table::fread(file.path(a$run,"provenance/numerical_validation.tsv"))
  check("campaign_numerical_gates_passed",nrow(numerical)>0&&all(numerical$passed))
  if(identical(b$engine,"Figure7_native_fit_last5_current")) {
    source(file.path(code,"util/analysis/figure7_current_pipeline.R"))
    f7c_bundle(a)
    grid <- data.table::fread(file.path(a$run,"provenance/grid_certificates.tsv"))
    qc <- data.table::fread(file.path(a$run,"provenance/task_numerical_validation.tsv"))
    check("native_last5_all_grid_certificates",nrow(grid)>0&&all(grid$passed))
    check("native_last5_all_route_counter_QC",all(qc$passed)&&
      all(qc$daily_complement_error<1e-10))
    fp_write(data.table::rbindlist(checks),file.path(root,"daily_source_acceptance.tsv"))
    cat("PASS: independent daily examples and current native-last5 campaign receipts\n")
    quit(save="no",status=0L)
  }
  n <- 0L;complement <- recompute <- 0
  for(id in b$tasks$task_id) {
    task <- fp_read(fp_task_path(a,id))
    for(o in task$outputs) {
      d <- o$daily_counts;H <- d$HL_per_start_live_cell;L <- d$LL_per_start_live_cell
      ok <- is.finite(H)&is.finite(L)&H+L>0
      check(paste0("paired_NA_mask_task_",id,"_beta_",o$metadata$q),
        identical(is.na(d$high_parent_percent),is.na(d$low_parent_percent)))
      if(any(ok)) {
        complement <- max(complement,max(abs(d$high_parent_percent[ok]+d$low_parent_percent[ok]-100)))
        recompute <- max(recompute,max(abs(100*H[ok]/(H[ok]+L[ok])-d$high_parent_percent[ok])),
          max(abs(100*L[ok]/(H[ok]+L[ok])-d$low_parent_percent[ok])))
      }
      n <- n+sum(ok)
    }
  }
  check("raw_repeat_counters_independent_recalculation",recompute<1e-10,recompute)
  check("raw_repeat_complement",complement<1e-10,complement)
  check("nonempty_daily_windows_exercised",n>0,n)
  qc <- data.table::fread(file.path(a$run,"provenance/task_numerical_validation.tsv"))
  check("actual_stochastic_passages_exercised",sum(qc[model_context=="in vitro"]$n_passages)>0,
    sum(qc[model_context=="in vitro"]$n_passages))
  # The preserved history result is a post-computation regression reference,
  # never an input to trajectory generation. Repeat count differs in culture;
  # exact historical reference parity is required for deterministic in vivo.
  keys <- c("pair_label","model_context","O2_pct","p_misseg","threshold_N","initial_ploidy","q","metric","day")
  old <- data.table::fread(file.path(a$workspace,"data/Figures/Figure7/tables/history_time_summary.tsv"))
  now <- data.table::fread(file.path(a$run,"tables/history_time_summary.tsv"))
  old <- old[model_context=="in vivo"];now <- now[model_context=="in vivo"]
  joined <- merge(now,old,by=keys);error <- 0
  for(m in c("median","q25","q75")) {
    x <- joined[[paste0(m,".x")]];y <- joined[[paste0(m,".y")]]
    if(!identical(is.na(x),is.na(y)))error <- Inf
    good <- is.finite(x)&is.finite(y);if(any(good))error <- max(error,max(abs(x[good]-y[good])))
  }
  check("unchanged_history_vivo_reference_parity",nrow(joined)==nrow(now)&&error<1e-6,error)
}
fp_write(data.table::rbindlist(checks),file.path(root,"daily_toy_validation.tsv"))
fp_write(data.frame(window="[0,1)",HL_per_initial_live_cell=.46,LL_per_initial_live_cell=.20,
  denominator=.66,A_percent=100*.46/.66,B_percent=100*.20/.66),file.path(root,"meeting_numerical_example.tsv"))
cat("PASS:",length(checks),"independent daily-counter tests\n")
