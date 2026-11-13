# Isolated tumor adapter. Numerical functions remain in the frozen source.
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  x[] <- lapply(x, function(v) if (is.double(v)) { z <- sprintf("%.17g",v); z[is.na(v)] <- NA_character_; z } else v)
  tmp <- tempfile(".pending_",dirname(path)); on.exit(unlink(tmp))
  write.table(x,tmp,sep="\t",quote=TRUE,row.names=FALSE,na="NA",qmethod="double")
  stopifnot(file.rename(tmp,path))
}
write_kv <- function(x,path) {
  types <- vapply(x,function(v) if(is.logical(v)) "boolean" else if(is.integer(v)) "integer" else if(is.numeric(v)) "double" else "string",character(1))
  vals <- vapply(x,function(v) if(length(v)!=1L) stop("Non-scalar setting") else if(is.double(v) && is.nan(v)) "NaN" else if(is.na(v)) NA_character_ else if(is.logical(v)) tolower(as.character(v)) else if(is.numeric(v)) sprintf("%.17g",v) else as.character(v),character(1))
  missing <- vapply(x,function(v) is.na(v) && !is.nan(v),logical(1))
  write_tsv(data.frame(key=names(x),type=types,value=vals,missing=missing),path)
}
read_tsv <- function(path) read.delim(path,stringsAsFactors=FALSE,check.names=FALSE,na.strings="NA")
parse_exported_double <- function(s) {
  # Some base-R decimal readers differ by one ULP across platforms. A %.17g
  # export uniquely identifies its IEEE double; verify the formatted value and
  # inspect adjacent bit patterns rather than relaxing numeric tolerances.
  v <- as.numeric(s)
  if(is.na(s)) return(NA_real_)
  if(s=="NaN") return(NaN)
  if(is.na(v)) stop("Invalid numeric export")
  if(is.infinite(v) && s%in%c("Inf","-Inf")) return(v)
  if(is.infinite(v)) v <- sign(v)*.Machine$double.xmax
  if(identical(sprintf("%.17g",v),s)) return(v)
  initial <- as.integer(writeBin(v,raw(),size=8L,endian="big"))
  candidates <- numeric()
  for(direction in c(-1L,1L)) {
    bits <- initial
    for(step in 1:4) {
      carry <- direction
      for(i in 8:1) {
        n <- bits[i]+carry
        bits[i] <- n%%256L
        carry <- if(n<0) -1L else if(n>255) 1L else 0L
        if(carry==0L) break
      }
      x <- readBin(as.raw(bits),numeric(),n=1L,size=8L,endian="big")
      if(is.finite(x) && identical(sprintf("%.17g",x),s)) candidates <- c(candidates,x)
    }
  }
  if(length(candidates)==1L) return(candidates[[1]])
  stop("Noncanonical or ambiguous full-precision decimal export: ",s)
}
read_kv <- function(path) {
  x <- read.delim(path,colClasses="character",na.strings=NULL,check.names=FALSE)
  stopifnot(identical(names(x),c("key","type","value")) || identical(names(x),c("key","type","value","missing")),!anyDuplicated(x$key))
  if("missing"%in%names(x)) {
    stopifnot(all(tolower(x$missing)%in%c("true","false")))
    x$value[tolower(x$missing)=="true"] <- NA_character_
  } else x$value[x$value=="NA"] <- NA_character_
  setNames(lapply(seq_len(nrow(x)),function(i) switch(x$type[i],boolean={stopifnot(is.na(x$value[i]) || x$value[i]%in%c("true","false"));as.logical(x$value[i])},integer=as.integer(x$value[i]),double=parse_exported_double(x$value[i]),string=x$value[i],stop("Unknown setting type"))),x$key)
}
verify_inputs <- function(study) {
  for(path in c(file.path(study,"manifests/input_manifest.tsv"))) {
    m <- read_tsv(path)
    stopifnot(all(vapply(file.path(study,m$path),hash_file,character(1))==m$sha256))
  }
}
hash_file <- function(path) digest::digest(file=path,algo="sha256")
atomic_rds <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(".pending_",dirname(path)); on.exit(unlink(tmp))
  saveRDS(x,tmp,version=3); stopifnot(file.rename(tmp,path))
}
assert_namespaces <- function() {
  bad <- intersect(loadedNamespaces(),c("rJava","V8","jsonlite","rjson","RJSONIO"))
  if(length(bad)) stop("Forbidden executed namespace: ",paste(bad,collapse=","))
}
load_backend <- function(study, numerical=FALSE) {
  root <- file.path(study,"source/numerical")
  e <- new.env(parent=.GlobalEnv)
  e$commandArgs <- function(trailingOnly=FALSE) {
    x <- base::commandArgs(trailingOnly=trailingOnly)
    if(!trailingOnly) x <- x[!startsWith(x,"--file=")]
    x
  }
  Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR=file.path(root,"model"))
  source(file.path(root,"util/o2_supply_demand_map_fit_invivo_backend.R"),local=e,chdir=TRUE)
  if(numerical) source(file.path(root,"model/model_O2_supply_demand_MAP.R"),local=.GlobalEnv)
  assert_namespaces(); e
}
load_contract <- function(study,numerical=TRUE,policy="production") {
  verify_inputs(study)
  e <- load_backend(study,numerical)
  a <- readRDS(file.path(study,"manifests/contract.rds"))
  a$study <- study
  a$backend <- e
  a$evaluation_mode <- Sys.getenv("IVI_EVALUATION_MODE","legacy")
  if(numerical && a$evaluation_mode=="fresh_fork") {
    stopifnot(.Platform$OS.type=="unix",as.character(getRversion())=="4.4.2",
      Sys.getenv("IVI_IMAGE_SHA256")=="dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372")
  }
  a$cfg$model_path <- file.path(study,"source/numerical/model/model_O2_supply_demand_MAP.R")
  a$cfg$parameter_table <- file.path(study,"source/seed25/parameter_table_input.csv")
  a$cfg$necrosis_mapping_csv <- file.path(study,"source/inputs/histology_to_dt_Gem_VT_20260209_v5_mapping.csv")
  if(numerical) {
    d <- o2simps_cpp_dll_info()
    a$cfg$cpp_dll_name <- d$name; a$cfg$cpp_dll_path <- d$path; a$cfg$cpp_wrapper_path <- d$wrapper_path
  }
  stopifnot(policy %in% c("production","data_only"))
  if(policy=="data_only") a$cfg$use_soft_prior <- FALSE
  a$policy <- policy; a
}
component_values <- function(x) x[c("L","L_data","L_prior","L_prior_raw","L_b","L_p","L_n","L_n_raw","n_burden","n_burden_obs_total","n_ploidy","n_ploidy_obs_total","n_necrosis","n_necrosis_obs_total","n_prior_terms")]
score <- function(ctx,p) {
  if(identical(ctx$evaluation_mode,"fresh_fork")) return(fork_model_value(function() {
    raw <- ctx;raw$evaluation_mode <- "legacy"
    score(raw,p)
  },kind="score",trajectories=length(ctx$scenarios)))
  stopifnot(identical(names(p),names(ctx$endpoint)),all(is.finite(p)),all(p>=ctx$lower-1e-12),all(p<=ctx$upper+1e-12))
  z <- ctx$backend$evaluate_objective_components(p,ctx$scenarios,ctx$cfg)
  stopifnot(is.finite(z$L),z$n_burden_obs_total==ctx$expected_counts$burden,
    z$n_ploidy_obs_total==ctx$expected_counts$ploidy,z$n_necrosis_obs_total==ctx$expected_counts$necrosis)
  z
}
predictions <- function(ctx,p,full_states=FALSE) {
  if(identical(ctx$evaluation_mode,"fresh_fork")) return(fork_model_value(function() {
    raw <- ctx;raw$evaluation_mode <- "legacy"
    invisible(score(raw,p))
    predictions(raw,p,full_states)
  },kind="predictions",trajectories=2L*length(ctx$scenarios)))
  native <- get("cpp_o2simps_simulate_one",envir=.GlobalEnv)
  capture <- list(); counter <- 0L
  ctx$backend$cpp_o2simps_simulate_one <- function(sim_args) {
    sim_args$return_full_trajectory <- full_states
    z <- native(sim_args); counter <<- counter+1L
    capture[[counter]] <<- list(arguments=sim_args,result=z)
    z
  }
  on.exit(rm("cpp_o2simps_simulate_one",envir=ctx$backend))
  rp <- ctx$backend$decode_params(p,FALSE,FALSE,ctx$cfg)
  out <- ctx$backend$collect_predictions(rp,ctx$scenarios,ctx$cfg)
  stopifnot(length(capture)==length(ctx$scenarios))
  for(i in seq_along(capture)) {
    s <- ctx$scenarios[[i]]; z <- capture[[i]]$result
    stopifnot(all(is.finite(z$Vmm3_total_obs)),all(z$Vmm3_total_obs>0),
      all(is.finite(z$frac_N_live)),abs(sum(z$frac_N_live)-1)<1e-10,
      all(is.finite(z$O2_eff_obs)),is.finite(z$Vmm3_total_terminal),z$Vmm3_total_terminal>0)
    k <- match(s$harvest,out$necrosis$harvest)
    if(!is.na(k)) {
      out$necrosis$pred_dead_total_volume[k] <- z$Vmm3_dead_total_terminal
      out$necrosis$pred_dead_hypoxia_volume[k] <- z$Vmm3_dead_hypoxia_terminal
      out$necrosis$pred_dead_buffer_volume[k] <- z$Vmm3_dead_buffer_terminal
      out$necrosis$pred_total_volume[k] <- z$Vmm3_total_terminal
      f <- z$Vmm3_dead_total_terminal/z$Vmm3_total_terminal
      out$necrosis$pred_necrosis_fraction[k] <- f
      eps <- ctx$cfg$necrosis_fraction_eps
      out$necrosis$residual_logit[k] <- qlogis(pmin(1-eps,pmax(eps,f)))-qlogis(pmin(1-eps,pmax(eps,s$obs_necrosis_fraction)))
    }
  }
  list(tables=out,native=capture)
}
meter_values <- function() {
 if(!exists(".ivi_runtime_meter",envir=.GlobalEnv,inherits=FALSE)) return(list(score_children=0L,prediction_children=0L,completed_trajectory_calls=0L,score_seconds=0,prediction_seconds=0))
 as.list(get(".ivi_runtime_meter",envir=.GlobalEnv))
}
fork_model_value <- function(fn,kind="score",trajectories=8L) {
  if(.Platform$OS.type!="unix") stop("Process isolation requires the pinned Unix runtime")
 if(!exists(".ivi_runtime_meter",envir=.GlobalEnv,inherits=FALSE)) {
  e<-list2env(meter_values(),parent=emptyenv());assign(".ivi_runtime_meter",e,envir=.GlobalEnv)
 }
 t<-proc.time()[[3]]
  job <- parallel::mcparallel(list(value=fn(),pid=Sys.getpid()),mc.set.seed=FALSE,silent=TRUE)
  result <- parallel::mccollect(job,wait=TRUE)[[1]]
  if(inherits(result,"try-error") || !is.list(result) || is.null(result$pid)) stop("Isolated model evaluation failed: ",paste(result,collapse=" "))
  assign(".ivi_last_model_pid",result$pid,envir=.GlobalEnv)
 e<-get(".ivi_runtime_meter",envir=.GlobalEnv)
 if(kind=="score") {e$score_children<-e$score_children+1L;e$score_seconds<-e$score_seconds+proc.time()[[3]]-t}
 else {e$prediction_children<-e$prediction_children+1L;e$prediction_seconds<-e$prediction_seconds+proc.time()[[3]]-t}
 e$completed_trajectory_calls<-e$completed_trajectory_calls+as.integer(trajectories)
  result$value
}
independent_components <- function(ctx,p,pred) {
  rp <- ctx$backend$decode_params(p,FALSE,FALSE,ctx$cfg)
  b <- pp <- nn <- numeric(length(ctx$scenarios)); pp2 <- pp4 <- numeric()
  nn[] <- NA_real_
  for(i in seq_along(ctx$scenarios)) {
    s <- ctx$scenarios[[i]]; z <- pred$native[[i]]$result
    keep <- ctx$cfg$scenario_cpp$keep_burden[[i]] & is.finite(s$obs_burden) & s$obs_burden>0
    b[i] <- mean(log(rp$sigma_burden)+.5*log(2*pi)+.5*((log(s$obs_burden[keep])-log(z$Vmm3_total_obs[keep]))/rp$sigma_burden)^2)
    means <- ctx$cfg$model_core$grid_pre
    fractions <- as.numeric(z$frac_N_live); fractions <- fractions/sum(fractions)
    v <- vapply(s$endpoint_obs_z,function(obs) {
      v <- log(pmax(fractions,1e-300))+dnorm(obs,means,ctx$cfg$sigma_ploidy,log=TRUE)
      m <- max(v); -(m+log(sum(exp(v-m))))
    },numeric(1))
    pp[i] <- mean(v)
    if(s$cohort=="2N") pp2 <- c(pp2,pp[i]) else pp4 <- c(pp4,pp[i])
    if(ctx$cfg$scenario_cpp$keep_necrosis[i]) {
      row <- match(s$harvest,pred$tables$necrosis$harvest)
      nn[i] <- (pred$tables$necrosis$residual_logit[row]/ctx$cfg$sigma_necrosis_logit)^2
    }
  }
  c(L_b=mean(b),L_p=.5*mean(pp2)+.5*mean(pp4),L_n=ctx$cfg$lambda_necrosis*mean(nn,na.rm=TRUE))
}
finalize_payload <- function(out) {
  stopifnot(requireNamespace("digest",quietly=TRUE))
  assert_namespaces()
  write_kv(meter_values(),file.path(out,"runtime_meter_master.tsv"))
  write_tsv(data.frame(namespace=sort(loadedNamespaces()),version=vapply(sort(loadedNamespaces()),function(n) as.character(packageVersion(n)),character(1))),file.path(out,"loaded_namespaces.tsv"))
  paths <- sort(list.files(out,recursive=TRUE,full.names=TRUE,all.files=TRUE,no..=TRUE)); paths <- paths[!file.info(paths)$isdir & !basename(paths)%in%c("payload_manifest.tsv","COMPLETE")]
  write_tsv(data.frame(path=substring(paths,nchar(out)+2L),sha256=vapply(paths,hash_file,character(1))),file.path(out,"payload_manifest.tsv"))
  writeLines(hash_file(file.path(out,"payload_manifest.tsv")),file.path(out,"COMPLETE"))
}
