#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", arg[1L])))
source(file.path(here, "util/analysis/supp_figure4_4.R"))
fails <- function(expr) inherits(try(force(expr), silent = TRUE), "try-error")
stopifnot(identical(supp4_4_quantiles(1:4), c(1.75, 2.5, 3.25)))
dt <- CJ(seed_number = 1:4, harvest = c("a", "b"), day = 0:1)
dt[, cohort := "2N"]
dt[, ploidy := seed_number + ifelse(harvest == "a", 0, 10)]
dt[, effective_o2_pct := seed_number / 10]
dt[, tumor_size_mm3 := seed_number * 100 + ifelse(harvest == "a", 0, 50)]
ans <- supp4_4_aggregate(dt, 1:4)
stopifnot(all(ans$cohort$n_tumors == 2),
          all(ans$summary[metric == "ploidy", median] == 7.5),
          all(ans$summary[metric == "ploidy", q25] == 6.75),
          all(ans$summary[metric == "ploidy", q75] == 8.25),
          fails(supp4_4_aggregate(rbind(dt, dt[1L]), 1:4)),
          fails(supp4_4_aggregate(dt[seed_number != 4L], 1:4)),
          fails(supp4_4_aggregate(dt[-1L], 1:4)),
          fails(supp4_4_quantiles(c(1, NA_real_))))
changed <- dt[!(harvest == "b" & day == 1L)]
seg <- supp4_4_aggregate(changed, 1:4)
stopifnot(all(seg$cohort[day == 0L, segment] == 1L),
          all(seg$cohort[day == 1L, segment] == 2L),
          all(seg$summary[day == 1L, n_seeds] == 4L))
cat("supp_figure4_4_unit_tests_ok\n")
