# Export the exact binary64 anchor encodings and the snapshot configuration.
bridge <- "oxygen/results/structural_identifiability/20261116_culture_step05_bridge_design_01"
anchor <- readRDS(file.path(bridge, "inputs/evidence/d0_truth_parameters.rds"))
roles <- read.delim(file.path(bridge, "inputs/parameter_roles.tsv"), check.names = FALSE)
table <- read.csv(file.path(bridge, "inputs/data/parameter_table_input.csv"))
i <- match(roles$parameter, table$param_symbol)
stopifnot(!anyNA(i), all(roles$parameter %in% names(anchor)))
result <- data.frame(parameter = roles$parameter,
                     binary64 = vapply(anchor[roles$parameter], function(x) sprintf("%a", as.numeric(x)), ""),
                     lower = table$lower_bound[i], upper = table$upper_bound[i])
write.table(result, "oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step02/anchor.tsv",
            quote = FALSE, row.names = FALSE, sep = "\t")
