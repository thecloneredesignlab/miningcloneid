#!/usr/bin/env Rscript
# Read saved telemetry only. No model source, propagator, fit, or RNG is loaded.
suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite)})
setDTthreads(2L)
args <- commandArgs(TRUE)
opts <- as.list(setNames(sub("^[^=]*=", "", args), sub("^--([^=]*)=.*$", "\\1", args)))
stopifnot(all(c("current-cache", "vivo-cache", "output-root") %in% names(opts)))
cur <- normalizePath(opts[["current-cache"]], mustWork=TRUE)
vivo <- normalizePath(opts[["vivo-cache"]], mustWork=TRUE)
out <- opts[["output-root"]]
if (dir.exists(out)) stop("Use a new output directory; existing results are never overwritten")
dir.create(out, recursive=TRUE)
sha <- function(p) vapply(p, digest, "", file=TRUE, algo="sha256", USE.NAMES=FALSE)
sources <- list(); checks <- list(); rows <- list()
record <- function(path, role) {
  sources[[length(sources)+1L]] <<- data.table(path=path, role=role, sha256=sha(path))
}
cp <- file.path(cur, "prepared/analysis_bundle.rds")
vp <- file.path(vivo, "prepared/analysis_bundle.rds")
cb <- readRDS(cp); vb <- readRDS(vp); record(cp,"current_bundle"); record(vp,"saved_vivo_bundle")
stopifnot(cb$engine=="Figure7_native_fit_last5_current",cb$profile=="full")
tasks <- as.data.table(cb$tasks)
tasks <- tasks[threshold_N==66L & O2_pct %in% c(.5,5) &
  (model_context=="in vivo" | maximum_ploidy==141L)]
key <- c("pair_label","model_context","O2_pct","p_misseg","threshold_N","initial_ploidy")
vt <- as.data.table(vb$tasks)
groups <- c("all_low","direct_first_generation","other_low","historical_lineage","never_HL_lineage")
metric_names <- c("mean_low_birth_rate","mean_division_rate","mean_ploidy")
metadata <- function(t,beta) as.list(cbind(t,q=beta))
for (i in seq_len(nrow(tasks))) {
  t <- tasks[i]
  tp <- file.path(cur,"tasks",sprintf("task_%05d.rds",t$task_id))
  task <- readRDS(tp);record(tp,"current_task")
  receipt <- fread(file.path(cur,"receipts",sprintf("task_%05d.tsv",t$task_id)))
  stopifnot(task$fingerprint==cb$fingerprint,receipt$sha256==sha(tp),receipt$status=="task_validated")
  old <- NULL
  if (t$model_context=="in vivo") {
    ix <- which(Reduce(`&`,lapply(key,function(k)vt[[k]]==t[[k]])))
    stopifnot(length(ix)==1L)
    op <- file.path(vivo,"tasks",sprintf("task_%05d.rds",vt$task_id[ix]))
    old <- readRDS(op);record(op,"saved_matching_vivo_telemetry")
  }
  for (o in task$outputs) {
    beta <- o$metadata$q
    if (t$model_context=="in vitro") {
      stopifnot(length(o$route_paths)==3L)
      for (j in seq_along(o$route_paths)) {
        path <- file.path(cur,"routes",basename(o$route_paths[j]))
        stopifnot(sha(path)==o$route_sha256[j])
        r <- readRDS(path);record(path,"current_culture_route")
        stopifnot(r$fingerprint==cb$fingerprint,r$meta$beta==beta)
        m <- as.data.table(r$metrics)
        # Low chromosome states are the leading rows of the ordered grid;
        # D/U in the saved route contain only these rows. Verify against the
        # independently saved total low fraction before using their masses.
        X <- r$X[seq_len(nrow(r$D)),,drop=FALSE]
        stopifnot(max(abs(colSums(X)-m$low_fraction))<1e-12)
        G <- list(X,r$D,pmax(X-r$D,0),pmax(X-r$U,0),r$U)
        prefixes <- c("all","direct","other","historical","never")
        for (g in seq_along(groups)) {
          mass <- colSums(G[[g]])
          z <- data.table(day=seq_len(nrow(m))-1L,group=groups[g],route=r$meta$route,
            mass_fraction=mass,log_total_cells=m$log_cells,
            ancestry_identifiable=if("ancestry_identifiable"%in%names(m))m$ancestry_identifiable else TRUE,
            mean_low_birth_rate=m[[paste0(prefixes[g],"_low_birth_rate")]],
            mean_division_rate=m[[paste0(prefixes[g],"_division_rate")]],
            mean_ploidy=m[[paste0(prefixes[g],"_ploidy")]])
          rows[[length(rows)+1L]] <- cbind(t[,..key],q=beta,z)
        }
      }
    } else {
      j <- which(vapply(old$outputs,function(x)x$metadata$q,0.0)==beta)
      stopifnot(length(j)==1L)
      v <- old$outputs[[j]]
      fields <- c("group","metric","day","median","q25","q75")
      a <- as.data.table(v$population_summary)[,..fields]
      b <- as.data.table(o$population_summary)[,..fields]
      setorderv(a,c("group","metric","day"));setorderv(b,c("group","metric","day"))
      stopifnot(isTRUE(all.equal(a,b,tolerance=0,check.attributes=FALSE)))
      counts_a <- as.data.table(v$source_counts_summary)[,.(counter,day,median,q25,q75)]
      counts_b <- as.data.table(o$source_counts_summary)[,.(counter,day,median,q25,q75)]
      stopifnot(isTRUE(all.equal(counts_a,counts_b,tolerance=0,check.attributes=FALSE)))
      checks[[length(checks)+1L]] <- data.table(task_id=t$task_id,beta=beta,
        saved_vivo_task_id=vt$task_id[ix],population_and_count_parity=TRUE)
      low <- as.numeric(v$history$low_fraction[1,])
      start <- as.numeric(v$daily_counts$log_start_live_cells[1,])
      logn <- c(start[-1L],v$trajectories$final_log_population[1])
      stopifnot(length(low)==1001L,length(logn)==1001L,all(is.finite(logn)))
      for (g in groups) {
        pop <- as.data.table(o$population_summary)[group==g]
        cell <- pop[metric=="cell_share_percent"][order(day)]
        values <- dcast(pop[metric%in%metric_names],day~metric,value.var="median")
        z <- data.table(day=cell$day,group=g,route="continuous",mass_fraction=low*cell$median/100,
          log_total_cells=logn,ancestry_identifiable=TRUE)
        for (nm in metric_names) z[,(nm):=values[[nm]]]
        rows[[length(rows)+1L]] <- cbind(t[,..key],q=beta,z)
      }
    }
  }
  message("Read saved telemetry ",i,"/",nrow(tasks),": ",t$task_id)
}
z <- rbindlist(rows,fill=TRUE)
stopifnot(!anyDuplicated(z[,c(key,"q","route","day","group"),with=FALSE]))
saveRDS(z,file.path(out,"route_metrics.rds"),compress="gzip")
fwrite(unique(rbindlist(sources)),file.path(out,"cache_sources.tsv"),sep="\t")
fwrite(rbindlist(checks),file.path(out,"saved_vivo_parity.tsv"),sep="\t")
write_json(list(status="saved_telemetry_only",fingerprint=cb$fingerprint,
  route_metrics_sha256=sha(file.path(out,"route_metrics.rds")),
  rows=nrow(z),conditions=nrow(unique(z[,c(key,"q"),with=FALSE])),
  simulation_calls=0L,job_submissions=0L),file.path(out,"export_receipt.json"),pretty=TRUE,auto_unbox=TRUE)
