#!/usr/bin/env Rscript
# Reporting-policy boundary tests; no simulation or graphics file is generated.
suppressPackageStartupMessages(library(data.table))
script <- normalizePath(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1]))
source(file.path(dirname(script),"util/analysis/figure7_population_validity.R"))
z <- data.table(q=0,initial_ploidy=2,day=100,group="direct_first_generation",mass_fraction=1e-8,
  log_total_cells=log(1e9),ancestry_identifiable=TRUE,
  mean_low_birth_rate=0,mean_division_rate=0,mean_ploidy=2.5)
tests <- list()
check <- function(name,ok) {stopifnot(isTRUE(ok));tests[[length(tests)+1L]]<<-data.table(test=name,passed=TRUE)}
check("relative_boundary_inclusive_and_zero_rates_valid",f7pv_flags(z)$valid)
x <- copy(z);x[,mass_fraction:=.999e-8]
check("below_relative_threshold_excluded_despite_many_cells",!f7pv_flags(x)$valid)
x <- copy(z);x[,`:=`(mass_fraction=1,log_total_cells=0)]
check("one_expected_cell_inclusive",f7pv_flags(x)$valid)
check("one_expected_cell_not_flagged_lt1",!f7pv_flags(x)$expected_cells_lt1)
x[,log_total_cells:=log(.999)]
check("below_one_cell_does_not_filter_mean",f7pv_flags(x)$valid)
check("below_one_cell_flagged_separately",f7pv_flags(x)$expected_cells_lt1)
x[,log_total_cells:= -1000]
check("absolute_count_underflow_does_not_filter_mean",f7pv_flags(x)$valid)
check("underflow_retains_log_count_and_lt1_flag",f7pv_flags(x)$log_expected_cells== -1000&&f7pv_flags(x)$expected_cells_lt1)
x[,log_total_cells:=NA_real_]
check("unavailable_count_does_not_filter_finite_composition_mean",f7pv_flags(x)$valid)
check("unavailable_count_flagged_without_inventing_zero",is.na(f7pv_flags(x)$expected_cells_lt1)&&f7pv_flags(x)$expected_cells_status=="unavailable")
x <- copy(z);x[,`:=`(mass_fraction=1,log_total_cells=1000)]
check("large_cell_counts_checked_without_exponential_overflow",f7pv_flags(x)$valid)
x[,`:=`(mass_fraction=1e-320,log_total_cells=1000)]
check("tiny_relative_mass_excluded_even_when_absolute_count_is_large",!f7pv_flags(x)$valid)
x <- copy(z);x[,mass_fraction:=0]
check("zero_mass_excluded",!f7pv_flags(x)$valid)
x <- copy(z);x[,q:=1]
check("complete_block_direct_is_structurally_empty",f7pv_flags(x)$reason=="structurally_empty")
x[,group:="historical_lineage"]
check("complete_block_history_is_structurally_empty",f7pv_flags(x)$reason=="structurally_empty")
x[,group:="never_HL_lineage"]
check("complete_block_never_lineage_not_assumed_empty",f7pv_flags(x)$valid)
x[,initial_ploidy:=4]
check("pure4N_never_low_lineage_is_theoretically_empty",!f7pv_flags(x)$valid)
x[,`:=`(group="all_low",q=1)]
check("pure4N_complete_block_has_no_low_group",!f7pv_flags(x)$valid)
x <- copy(z);x[,day:=0]
check("initial_direct_label_is_theoretically_empty",!f7pv_flags(x)$valid)
x[,group:="all_low"]
check("initial_pure2N_low_group_not_assumed_empty",f7pv_flags(x)$valid)
x <- copy(z);x[,ancestry_identifiable:=FALSE]
check("unknown_ancestry_excluded",!f7pv_flags(x)$valid)
x <- copy(z);x[,mean_ploidy:=NaN]
check("nonfinite_mean_excluded",!f7pv_flags(x)$valid)
x <- copy(z);x[,mass_fraction:=-1e-15]
check("negative_mass_not_replaced_with_pseudocount",!f7pv_flags(x)$valid)
x <- rbind(copy(z),copy(z),copy(z));x[2,mass_fraction:=0]
f <- f7pv_flags(x)
check("three_routes_not_complete_when_one_fails",sum(f$valid)!=3L)
print(rbindlist(tests));cat("PASS:",length(tests),"reporting-policy checks; no propagation\n")
