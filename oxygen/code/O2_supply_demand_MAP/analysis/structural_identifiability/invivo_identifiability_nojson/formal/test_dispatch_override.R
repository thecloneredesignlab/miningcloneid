# Pure metadata tests; no model loading, scoring, propagation or allocation.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
entry <- normalizePath(args[1]);source(file.path(entry,"dispatch_override.R"))
base <- list(role="reserve_baseline_global",mode="global",fixed=NULL,policy="data_only",
 init_mode="uniform",family="uniform",start_seed=0L,seed=31002L,NP=256L,itermax=1000L,reltol=1e-6,steptol=100L)
stopifnot(identical(base,apply_reserve_override(base,list())))
hybrid <- apply_reserve_override(base,list(override_init_mode="hybrid",override_start_seed="292"))
unchanged <- setdiff(names(base),c("init_mode","family","start_seed"))
stopifnot(identical(base[unchanged],hybrid[unchanged]),hybrid$init_mode=="hybrid",hybrid$start_seed==292L)
uniform <- apply_reserve_override(base,list(override_init_mode="uniform",override_start_seed="0"))
stopifnot(identical(base,uniform))
bad <- list(list(override_init_mode="hybrid",override_start_seed="0"),
 list(override_init_mode="uniform",override_start_seed="25"),
 list(override_init_mode="hybrid",override_start_seed="293"),
 list(override_init_mode="hybrid",override_start_seed="292",override_NP="512"))
for(context in bad) stopifnot(inherits(try(apply_reserve_override(base,context),silent=TRUE),"try-error"))
primary <- base;primary$role <- "baseline_global"
stopifnot(inherits(try(apply_reserve_override(primary,list(override_init_mode="hybrid",override_start_seed="292")),silent=TRUE),"try-error"))
cat("Reserve-dispatch metadata tests passed; only family and archived center can change; zero model calls.\n")
