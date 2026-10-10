#!/usr/bin/env Rscript
# Export canonical matrices for numerical diagnostics; no analysis outputs are changed.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
root <- normalizePath(args[[1L]])
out <- args[[2L]]
tasks <- as.integer(strsplit(args[[3L]], ",", fixed=TRUE)[[1L]])
dir.create(out, recursive=TRUE, showWarnings=FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1L]])
here <- dirname(normalizePath(script))
source(file.path(here, "../../simulation/o2/fixed_o2/run_fixed_o2_simulation.R"), local=TRUE)
suppressPackageStartupMessages(library(jsonlite))
plan <- fromJSON(file.path(root, "slurm/submission_plan.json"))
active <- c("lam_max", "p_mis_base", "p_misseg", "k_o_mis", "buffer_smax",
            "buffer_beta", "buffer_n_exp", "p_wgd", "alpha_o2", "gamma_growth",
            "mu_hp", "gamma_mu", "O2_crit", "n_O")
reports <- list()
options(digits=17)
for (task in tasks) {
  seed <- (task-1L) %/% 70L + 1L
  repeat_id <- ((task-1L) %% 70L) %/% 14L + 1L
  pindex <- (task-1L) %% 14L + 1L
  part <- file.path(root, "runs", paste0("seed", seed), "runs", paste0("N513_R", repeat_id),
                    "trajectories", sprintf("P%02d_%s", pindex, active[[pindex]]))
  raw <- read.delim(gzfile(file.path(part, "outputs.tsv.gz.part01.gz")))
  bad <- raw[raw$status != "ok" | is.na(raw$eigenvector_nonnegative) |
               !raw$eigenvector_nonnegative | !is.finite(raw$dominant_mean_ploidy) |
               !is.finite(raw$dominant_growth_rate), ]
  selected <- rbind(head(bad, 2L), raw[raw$sample_id == 1L & raw$O2_pct %in% c(0,2.5,5), ])
  samples <- read.delim(gzfile(file.path(part, "samples.tsv.gz")), check.names=FALSE)
  fit <- file.path(plan$fit_root, paste0("seed", seed))
  params <- read.delim(file.path(fit, "best_params.tsv"))
  base <- as.list(setNames(as.numeric(params$value), params$parameter))
  config <- readRDS(file.path(fit, "fit_config.rds"))
  cfg <- prepare_sim_cfg(config, list(), 0, prepare_run_params(base, "invivo", config, 0))
  for (j in seq_len(nrow(selected))) {
    row <- selected[j, ]
    values <- base
    for (name in active) values[[name]] <- samples[[name]][[row$sample_id]]
    rp <- prepare_run_params(values, "invivo", cfg, 0)
    fm <- fixo2_fixed_matrix(globalenv(), cfg, rp, row$O2_pct)
    eig <- eigen(fm$M)
    idx <- which.max(Re(eig$values))
    v <- Re(eig$vectors[,idx]); if (sum(v) < 0) v <- -v
    tag <- sprintf("task%05d_sample%03d_o2_%g", task, row$sample_id, row$O2_pct)
    saveRDS(fm$M, file.path(out, paste0(tag, ".matrix.rds")))
    con <- gzfile(file.path(out, paste0(tag, ".matrix.tsv.gz")), "wt")
    write.table(fm$M, con, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE); close(con)
    con <- gzfile(file.path(out, paste0(tag, ".vector.tsv.gz")), "wt")
    write.table(data.frame(N=fm$ngrid, raw_R_vector=v), con, sep="\t", quote=FALSE, row.names=FALSE); close(con)
    off <- fm$M; diag(off) <- 0
    reports[[length(reports)+1L]] <- list(tag=tag, task_id=task, sample_id=row$sample_id,
      O2_pct=row$O2_pct, original=as.list(row), min_off_diagonal=min(off),
      raw_vector_min=min(v), negative_l1_mass=sum(abs(v[v<0]))/sum(abs(v)),
      lambda=Re(eig$values[[idx]]), lambda_imaginary=Im(eig$values[[idx]]),
      raw_residual=max(abs(fm$M %*% v-Re(eig$values[[idx]])*v)),
      eigenvalues_real=Re(eig$values), eigenvalues_imaginary=Im(eig$values))
  }
}
write_json(reports, file.path(out, "matrix_manifest.json"), pretty=TRUE, auto_unbox=TRUE, digits=17)
cat("Exported", length(reports), "canonical matrices\n")
