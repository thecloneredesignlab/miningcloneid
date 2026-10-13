#!/usr/bin/env Rscript
# Canonical fixed-O2 matrix, checked dominant eigenpair, unchanged-matrix recovery.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
part <- normalizePath(args[[1L]])
fit_root <- normalizePath(args[[2L]])
reference_file <- normalizePath(args[[3L]])
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[[1L]])
here <- dirname(normalizePath(script))
source(file.path(here, '../../../simulation/o2/fixed_o2/run_fixed_o2_simulation.R'), local=TRUE)
suppressPackageStartupMessages(library(jsonlite))
hash_file <- function(p) strsplit(system2('sha256sum', shQuote(p), stdout=TRUE), ' ')[[1L]][[1L]]
meta <- fromJSON(file.path(part, 'metadata.json'))
stopifnot(meta$method %in% c('Sobol_second_order', 'Sobol_response_surface'))
sample_path <- file.path(part, 'samples.tsv.gz')
stopifnot(hash_file(sample_path) == meta$samples_sha256)
samples <- read.delim(gzfile(sample_path), check.names=FALSE)
stopifnot(nrow(samples) == meta$n_samples, identical(as.integer(samples$sample_id), seq_len(nrow(samples))))
stopifnot(identical(names(samples)[-1L], as.character(meta$parameters)))
fit <- file.path(fit_root, meta$fit_seed)
stopifnot(hash_file(file.path(fit,'best_params.tsv')) == meta$best_params_sha256,
          hash_file(file.path(fit,'fit_config.rds')) == meta$fit_config_sha256)
raw <- readRDS(file.path(fit,'fit_config.rds'))
bp <- read.delim(file.path(fit,'best_params.tsv'))
base <- as.list(setNames(as.numeric(bp$value), bp$parameter))
cfg <- prepare_sim_cfg(raw, list(), 0, prepare_run_params(base,'invivo',raw,0))
stopifnot(cfg$N_MIN == 22L, cfg$N_MAX == 154L, cfg$N_UNIT == 22L,
          isTRUE(cfg$O2_growth), identical(cfg$ploidy_O2_death,'ploidy_related'))

solve_one <- function(p, o2, id, save_recovery=TRUE) {
  rp <- prepare_run_params(p,'invivo',cfg,0)
  fm <- fixo2_fixed_matrix(globalenv(),cfg,rp,o2)
  eg <- eigen(fm$M, only.values=FALSE)
  k <- which.max(Re(eg$values)); lam <- Re(eg$values[[k]])
  v <- Re(eg$vectors[,k]); if (sum(v) < 0) v <- -v
  nonneg <- all(is.finite(v)) && all(v >= -1e-8)
  v <- fixo2_normalize_state(v)
  residual <- max(abs(fm$M %*% v - lam*v))
  solver <- 'double'; original <- NULL
  if (!nonneg || any(!is.finite(v)) || !is.finite(residual) ||
      residual > 1e-10 * max(1,max(abs(fm$M))) || abs(Im(eg$values[[k]])) > 1e-8) {
    original <- list(ploidy=sum(fm$ngrid*v)/cfg$N_UNIT, growth=lam,
                     nonnegative=nonneg, residual=residual)
    if (!exists('perron_high_precision', mode='function')) {
      Rcpp::sourceCpp(file.path(here,'../perron_high_precision_safe.cpp'),
                      cacheDir=Sys.getenv('SOBOL_PERRON_CACHE'), rebuild=FALSE)
    }
    p50 <- perron_high_precision(fm$M,50L)
    p100 <- perron_high_precision(fm$M,100L)
    stopifnot(all(p50$vector >= 0), all(p100$vector >= 0),
              abs(p50$lambda-p100$lambda) < 1e-12,
              sum(abs(p50$vector-p100$vector)) < 1e-10,
              abs(sum(fm$ngrid*(p50$vector-p100$vector))/cfg$N_UNIT) < 1e-12,
              p50$upper_minus_lambda < 1e-20, p100$upper_minus_lambda < 1e-20)
    v <- p50$vector; lam <- p50$lambda; nonneg <- TRUE; solver <- 'decimal50_100'
    residual <- max(abs(fm$M %*% v-lam*v))
    stopifnot(residual < 1e-12)
    if (save_recovery) {
      recovery_dir <- file.path(part,'recovery'); dir.create(recovery_dir, showWarnings=FALSE)
      stem <- sprintf('sample%05d_o2_%s',id,format(o2,scientific=FALSE,trim=TRUE))
      matrix_path <- file.path(recovery_dir,paste0(stem,'_matrix.rds'))
      saveRDS(fm$M,matrix_path)
      write_json(list(original=original, precision50=p50, precision100=p100,
        matrix_sha256=hash_file(matrix_path), matrix_file=basename(matrix_path),
        eigen_residual=residual),file.path(recovery_dir,paste0(stem,'.json')),
        pretty=TRUE,auto_unbox=TRUE,digits=17)
    }
  }
  stopifnot(nonneg, all(is.finite(v)), abs(sum(v)-1) < 1e-12, is.finite(lam))
  data.frame(sample_id=id,O2_pct=o2,dominant_mean_ploidy=sum(fm$ngrid*v)/cfg$N_UNIT,
    dominant_growth_rate=lam,spectral_gap=diff(sort(Re(eg$values),decreasing=TRUE)[2:1]),
    eigenvector_nonnegative=TRUE,eigen_residual=residual,solver=solver,status='ok')
}

validation <- list()
if (isTRUE(meta$validate_reference)) {
  reference <- read.delim(reference_file)
  for (o2 in c(0,2.5,5)) {
    expected <- reference[reference$seed_id == meta$fit_seed & abs(reference$O2_pct-o2)<1e-12,]
    actual <- solve_one(base,o2,0,FALSE)
    canonical <- fixo2_dominant_attractor_one(meta$fit_seed,
      prepare_run_params(base,'invivo',cfg,0),globalenv(),cfg,o2)
    stopifnot(nrow(expected)==1L, canonical$status=='ok',
      abs(actual$dominant_mean_ploidy-expected$dominant_mean_ploidy)<1e-8,
      abs(actual$dominant_growth_rate-expected$dominant_growth_rate)<1e-8,
      abs(actual$dominant_mean_ploidy-canonical$dominant_mean_ploidy)<1e-8,
      abs(actual$dominant_growth_rate-canonical$dominant_growth_rate)<1e-8)
    validation[[length(validation)+1L]] <- actual
  }
}
common <- list()
if (isTRUE(meta$validate_reference) && !is.null(meta$common_reference_params)) {
  p <- base
  for (name in names(meta$common_reference_params)) p[[name]] <- meta$common_reference_params[[name]]
  for (o2 in meta$oxygen_pct) common[[length(common)+1L]] <- solve_one(p,o2,0,FALSE)
}
start <- Sys.time()
out <- file.path(part,'outputs.tsv.gz')
temporary <- paste0(out,'.incomplete')
con <- gzfile(temporary,'wt')
for (i in seq_len(nrow(samples))) {
  p <- base
  for (name in names(samples)[-1L]) p[[name]] <- samples[[name]][[i]]
  rows <- do.call(rbind,lapply(meta$oxygen_pct,function(o2) solve_one(p,o2,i)))
  write.table(rows,con,sep='\t',quote=FALSE,row.names=FALSE,col.names=i==1L,na='NA')
  if (i %% 128L == 0L) message('Evaluated ',i,'/',nrow(samples),' parameter vectors')
}
close(con)
stopifnot(file.rename(temporary,out))
write_json(list(status='complete',n_rows=nrow(samples)*length(meta$oxygen_pct),
  metadata_sha256=hash_file(file.path(part,'metadata.json')),samples_sha256=hash_file(sample_path),
  outputs_sha256=hash_file(out),reference_checks=if(length(validation)) do.call(rbind,validation) else NULL,
  common_reference=if(length(common)) do.call(rbind,common) else NULL,
  elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs')),
  host=Sys.info()[['nodename']],job_id=Sys.getenv('SLURM_JOB_ID')),
  paste0(out,'.receipt.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
