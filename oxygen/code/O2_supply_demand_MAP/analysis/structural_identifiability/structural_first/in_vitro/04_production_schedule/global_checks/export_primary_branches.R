source("oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step03/backend.R")
source(file.path(out,"search_helpers.R"))
emit_native_pair_payload(readRDS(file.path(out,"paired_native.rds")),
                         file.path(out,"paired_native_outputs.tsv"),file.path(out,"paired_native_branches.tsv"))
