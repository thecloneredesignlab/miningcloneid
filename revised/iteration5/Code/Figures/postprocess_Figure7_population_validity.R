#!/usr/bin/env Rscript
# Derive reporting masks exclusively from saved route telemetry. Never source
# numerical engines, update populations, refit, sample, or submit a job.
suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite)})
setDTthreads(2L)
script <- normalizePath(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))][1]))
workspace <- normalizePath(file.path(dirname(script),"../.."))
args <- commandArgs(TRUE)
opts <- as.list(setNames(sub("^[^=]*=","",args),sub("^--([^=]*)=.*$","\\1",args)))
`%||%` <- function(x,y)if(is.null(x))y else x
if(length(setdiff(names(opts),c("workspace","input-root","telemetry-root","output-root"))))stop("Unknown option")
if("workspace"%in%names(opts))workspace <- normalizePath(opts[["workspace"]])
root <- normalizePath(opts[["input-root"]] %||% file.path(workspace,"data/Figures/Figure7"),mustWork=TRUE)
telemetry <- normalizePath(opts[["telemetry-root"]],mustWork=TRUE)
out <- opts[["output-root"]] %||% file.path(root,"postprocess/population_validity")
if(dir.exists(out))stop("Use an isolated new output root; promote only after acceptance")
dir.create(out,recursive=TRUE)
helper <- file.path(workspace,"Code/Figures/util/analysis/figure7_population_validity.R")
policy_path <- file.path(workspace,"Code/config/figure7_population_validity.json")
source(helper)
policy <- fromJSON(policy_path)
stopifnot(policy$schema==2L,policy$relative_mass_min==1e-8,
  isFALSE(policy$absolute_count_filter),isTRUE(policy$expected_cells_lt1_flag_only),
  isFALSE(policy$count_availability_filters_means),
  isTRUE(policy$require_all_protocol_routes))
sha <- function(p)vapply(p,digest,"",file=TRUE,algo="sha256",USE.NAMES=FALSE)
write <- function(z,name)fwrite(z,file.path(out,name),sep="\t",na="NA")
checks <- list()
check <- function(name,ok,detail="") {
  checks[[length(checks)+1L]] <<- data.table(test=name,passed=isTRUE(ok),detail=as.character(detail))
  if(!isTRUE(ok)) {write(rbindlist(checks),"checks.tsv");stop("Failed: ",name,": ",detail)}
}
cache_files <- c("route_metrics.rds","cache_sources.tsv","saved_vivo_parity.tsv","export_receipt.json")
stopifnot(all(file.exists(file.path(telemetry,cache_files))))
export <- fromJSON(file.path(telemetry,"export_receipt.json"))
check("telemetry_receipt",export$status=="saved_telemetry_only"&&export$simulation_calls==0&&
  export$job_submissions==0&&export$route_metrics_sha256==sha(file.path(telemetry,"route_metrics.rds")))
bundle_path <- file.path(root,"prepared/analysis_bundle.rds")
completion_path <- file.path(root,"provenance/completion.tsv")
done <- fread(completion_path)
check("current_scientific_fingerprint",done$fingerprint==export$fingerprint&&done$status=="validated_complete")
for(n in cache_files)stopifnot(file.copy(file.path(telemetry,n),file.path(out,n),overwrite=FALSE))
z <- as.data.table(readRDS(file.path(out,"route_metrics.rds")))
id <- c("pair_label","model_context","O2_pct","p_misseg","threshold_N","initial_ploidy","q","group","day")
check("unique_route_points",!anyDuplicated(z[,c(id,"route"),with=FALSE]))
check("known_source_groups",setequal(unique(z$group),c("all_low","direct_first_generation","other_low","historical_lineage","never_HL_lineage")))
check("no_negative_or_excess_group_mass",all(is.na(z$mass_fraction)|
  (z$mass_fraction>=0&z$mass_fraction<=1+1e-10)))
z[,required_routes:=ifelse(model_context=="in vitro",3L,1L)]
counts <- z[,.(n_routes=.N,required_routes=required_routes[1]),by=id]
check("complete_route_coverage",all(counts$n_routes==counts$required_routes))
raw_path <- file.path(root,"tables/low_population_time_summary.tsv")
raw <- fread(raw_path)
means <- policy$mean_metrics
check("unique_summary_points",!anyDuplicated(raw[,c(id,"metric"),with=FALSE]))
check("complete_telemetry_keys",setequal(do.call(paste,counts[,..id]),
  do.call(paste,unique(raw[metric%in%means,..id]))))
# Before masking, independently aggregate every saved route mean and compare
# it with the authenticated current table. This also verifies reuse of the
# matching saved in-vivo telemetry, rather than trusting a filename or fit ID.
for(metric_name in means) {
  z[,value:=get(metric_name)]
  # The saved design has one or three routes. Vectorized closed forms for
  # type-7 quartiles avoid millions of tiny quantile calls during a cache-only
  # check; all results are still compared against the original independent
  # summaries below, and no retained mean is replaced by this calculation.
  wide <- dcast(z,as.formula(paste(paste(id,collapse="+"),"~route")),value.var="value")
  v <- as.matrix(wide[,setdiff(names(wide),id),with=FALSE]);finite <- is.finite(v)
  n <- rowSums(finite);check(paste0("at_most_three_routes_",metric_name),all(n<=3L))
  vn <- v;vn[!finite] <- Inf;lo <- do.call(pmin,as.data.frame(vn))
  vn[!finite] <- -Inf;hi <- do.call(pmax,as.data.frame(vn))
  vn[!finite] <- 0;total <- rowSums(vn)
  med <- ifelse(n==3,total-lo-hi,ifelse(n==2,total/2,ifelse(n==1,total,NA_real_)))
  q25 <- ifelse(n==3,lo+(med-lo)*.5,ifelse(n==2,lo+(hi-lo)*.25,med))
  q75 <- ifelse(n==3,med+(hi-med)*.5,ifelse(n==2,lo+(hi-lo)*.75,med))
  q <- wide[,..id];q[,`:=`(median=med,q25=q25,q75=q75)]
  r <- raw[metric==metric_name,c(id,"median","q25","q75"),with=FALSE]
  setorderv(q,id);setorderv(r,id)
  # Original RDS roots store 2/4 as doubles; TSV readers store them as integers.
  # Require exact key values, without treating their storage type as a change.
  check(paste0("raw_route_summary_keys_",metric_name),
    isTRUE(all.equal(q[,..id],r[,..id],tolerance=0,check.attributes=FALSE)))
  deltas <- unlist(lapply(c("median","q25","q75"),function(n) {
    check(paste0("raw_NA_parity_",metric_name,"_",n),identical(is.na(q[[n]]),is.na(r[[n]])))
    abs(q[[n]]-r[[n]])
  }))
  error <- if(any(is.finite(deltas)))max(deltas[is.finite(deltas)]) else 0
  check(paste0("raw_route_summary_parity_",metric_name),error<1e-12,format(error,digits=17))
}
z[,value:=NULL]
flags <- f7pv_flags(z,policy$relative_mass_min)
z <- cbind(z,flags)
primary <- z[,.(n_eligible=sum(valid),n_routes=.N,required_routes=required_routes[1],
  n_expected_cells_lt1=sum(expected_cells_lt1,na.rm=TRUE),
  n_count_unavailable=sum(is.na(expected_cells_lt1)),structurally_empty=all(structurally_empty)),by=id]
primary[,summary_valid:=n_eligible==required_routes]
check("no_partial_route_summary",all(primary[summary_valid==TRUE]$n_eligible==primary[summary_valid==TRUE]$required_routes))
derived <- copy(raw)
raw_means <- raw[metric%in%means]
mask <- merge(raw_means[,c(id,"metric"),with=FALSE],primary,by=id,all.x=TRUE,sort=FALSE)
check("all_mean_points_have_validity",nrow(mask)==nrow(raw_means)&&!anyNA(mask$summary_valid))
setkeyv(derived,c(id,"metric"));setkeyv(mask,c(id,"metric"))
ix <- derived[mask,which=TRUE]
bad <- ix[!mask$summary_valid]
good <- ix[mask$summary_valid]
derived[bad,`:=`(median=NA_real_,q25=NA_real_,q75=NA_real_)]
derived[ix,n_valid:=mask$n_eligible]
# Restore original row order; keying is only for exact joins.
setkeyv(raw,c(id,"metric"))
check("masked_means_are_NA",all(is.na(derived[bad,median]))&&all(is.na(derived[bad,q25]))&&all(is.na(derived[bad,q75])))
check("valid_means_unchanged",identical(derived[good,.(median,q25,q75)],raw[good,.(median,q25,q75)]))
check("nonmean_values_unchanged",identical(derived[!metric%in%means],raw[!metric%in%means]))
check("structural_zero_means_masked",all(!primary[q==1&group%in%c("direct_first_generation","historical_lineage")]$summary_valid))
check("all_theoretical_empty_groups_masked",all(!primary[structurally_empty==TRUE]$summary_valid))
check("subcell_expectations_retained",nrow(z[valid&expected_cells_lt1==TRUE])>0L)
reference <- primary[pair_label=="C01"&model_context=="in vitro"&O2_pct==.5&
  p_misseg==.1&initial_ploidy==2&q==0&day==1000&
  group%in%c("all_low","direct_first_generation","other_low")]
check("known_day1000_subcell_means_retained",nrow(reference)==3L&&
  all(reference$summary_valid)&&all(reference$n_expected_cells_lt1==3L))
check("zero_rates_not_excluded",all(z[valid&(mean_division_rate==0|mean_low_birth_rate==0)]$valid))
summary_cols <- c("pair_label","model_context","O2_pct","p_misseg","initial_ploidy","q","group")
reason_summary <- z[,.(route_points=.N),by=c(summary_cols,"reason")]
write(reason_summary,"reason_summary.tsv")
coverage <- primary[,.(summary_points=.N,valid_points=sum(summary_valid),
  first_valid_day=if(any(summary_valid))min(day[summary_valid]) else NA_integer_,
  last_valid_day=if(any(summary_valid))max(day[summary_valid]) else NA_integer_),by=summary_cols]
write(coverage,"primary_coverage.tsv")
count_flags <- z[,.(route_points=.N,mean_valid=sum(valid),
  expected_cells_lt1=sum(expected_cells_lt1,na.rm=TRUE),
  valid_means_with_expected_cells_lt1=sum(valid&expected_cells_lt1,na.rm=TRUE),
  expected_count_unavailable=sum(is.na(expected_cells_lt1))),by=summary_cols]
write(count_flags,"expected_cell_count_flags_summary.tsv")
# Only relative thresholds are sensitivity filters. Absolute expected counts
# remain reporting fields and cannot alter any primary or sensitivity mask.
sens <- list();support <- list()
for(rel in policy$sensitivity_relative_mass) {
  f <- f7pv_flags(z,rel)
  z[,candidate_valid:=f$valid]
  s <- z[,.(candidate_valid=sum(candidate_valid)==required_routes[1]),by=id]
  s <- merge(s,primary[,c(id,"summary_valid"),with=FALSE],by=id,sort=FALSE)
  sens[[length(sens)+1L]] <- s[,.(summary_points=.N,valid_points=sum(candidate_valid),
    changed_from_primary=sum(candidate_valid!=summary_valid),
    valid_in_both=sum(candidate_valid&summary_valid)),by=.(model_context,group)][,
      `:=`(relative_mass_min=rel,absolute_count_filter=FALSE)]
  support[[length(support)+1L]] <- s[day%in%c(20L,100L,500L,1000L),
    c(id,"candidate_valid","summary_valid"),with=FALSE][,
      `:=`(relative_mass_min=rel,absolute_count_filter=FALSE)]
}
z[,candidate_valid:=NULL]
write(rbindlist(sens),"sensitivity_summary.tsv")
write(rbindlist(support),"sensitivity_selected_days.tsv.gz")
write(z,"route_validity.tsv.gz")
write(primary,"summary_validity.tsv.gz")
# Keep each retained scalar exactly as supplied, including zero rates. Only
# conditional means and their n_valid field can differ from the raw table.
raw_original <- fread(raw_path)
setorderv(derived,c("task_id","q","group","metric","day"))
setorderv(raw_original,c("task_id","q","group","metric","day"))
check("same_rows_and_schema",identical(names(derived),names(raw_original))&&nrow(derived)==nrow(raw_original))
check("unchanged_metadata",identical(derived[,setdiff(names(derived),c("median","q25","q75","n_valid")),with=FALSE],
  raw_original[,setdiff(names(raw_original),c("median","q25","q75","n_valid")),with=FALSE]))
write(derived,"low_population_time_summary.tsv")
check("three_relative_sensitivity_cases",nrow(unique(rbindlist(sens)[,.(relative_mass_min)]))==3L)
check("no_absolute_sensitivity_filter",all(!rbindlist(sens)$absolute_count_filter))
write(rbindlist(checks),"checks.tsv")
workspace_files <- c(script,helper,policy_path,file.path(workspace,"Code/Figures/export_Figure7_population_cache.R"))
overlay_files <- list.files(out,full.names=TRUE)
manifest <- rbindlist(list(
  data.table(scope="workspace",relative_path=substring(workspace_files,nchar(workspace)+2L),sha256=sha(workspace_files)),
  data.table(scope="figure7",relative_path=c("tables/low_population_time_summary.tsv","prepared/analysis_bundle.rds","provenance/completion.tsv"),
    sha256=sha(c(raw_path,bundle_path,completion_path))),
  data.table(scope="overlay",relative_path=basename(overlay_files),sha256=sha(overlay_files))))
write(manifest,"manifest.tsv")
write_json(list(schema=2L,status="postprocessed_validated",fingerprint=done$fingerprint,
  relative_mass_min=policy$relative_mass_min,absolute_count_filter=FALSE,
  expected_cells_lt1_flag_only=TRUE,count_availability_filters_means=FALSE,
  all_protocol_routes_required=TRUE,simulation_calls=0L,job_submissions=0L,
  manifest_sha256=sha(file.path(out,"manifest.tsv")),
  summary_rows=nrow(derived),mean_rows=nrow(raw_means),valid_group_days=sum(primary$summary_valid),
  invalid_group_days=sum(!primary$summary_valid),checks=length(checks),
  sensitivity_cases=3L,valid_subcell_route_points=nrow(z[valid&expected_cells_lt1==TRUE]),
  zero_rate_valid_route_points=nrow(z[valid&(mean_division_rate==0|mean_low_birth_rate==0)])),
  file.path(out,"completion.json"),pretty=TRUE,auto_unbox=TRUE)
message("PASS: relative-only cached postprocessing and three sensitivity cases; expected count <1 is flag only; output: ",out)
