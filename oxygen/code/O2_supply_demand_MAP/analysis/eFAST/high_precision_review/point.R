#!/usr/bin/env Rscript
# Preserve the canonical matrix and every normalization stage. No repairs.
args <- commandArgs(TRUE);stopifnot(length(args)==2L)
root <- normalizePath(args[[1L]]);id <- as.integer(args[[2L]])
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
here <- dirname(normalizePath(script));efast <- dirname(here)
options(digits=17)
suppressPackageStartupMessages(library(jsonlite));suppressPackageStartupMessages(library(Rcpp))
hash <- function(p) strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' ')[[1L]][[1L]]
selection <- fromJSON(file.path(root,'selection.json'),simplifyVector=FALSE)
stopifnot(hash(file.path(root,'point_manifest.tsv'))==selection$point_manifest_sha256)
for(name in names(selection$code_sha256))stopifnot(hash(file.path(here,name))==selection$code_sha256[[name]])
manifest <- read.delim(file.path(root,'point_manifest.tsv'),check.names=FALSE,colClasses='character')
row <- as.list(manifest[id,,drop=FALSE]);stopifnot(as.integer(row$point_id)==id)
config <- if(row$fit_seed=='seed138')'N8193_M4' else 'N4097_M8'
plan <- selection$source_plans[[config]]
for(name in names(plan$code_sha256))stopifnot(hash(file.path(efast,name))==plan$code_sha256[[name]])
source(file.path(efast,'../../simulation/o2/fixed_o2/run_fixed_o2_simulation.R'),local=TRUE)
sourceCpp(file.path(here,'perron_review.cpp'),cacheDir=Sys.getenv('EFAST_REVIEW_CACHE'),rebuild=FALSE)
repo <- normalizePath(file.path(efast,'../../../../..'))
sample_path <- file.path(repo,row$source_part,'samples.tsv.gz')
stopifnot(hash(sample_path)==row$samples_sha256,hash(file.path(repo,row$source_part,'metadata.json'))==row$metadata_sha256)
samples <- read.delim(gzfile(sample_path),check.names=FALSE)
active <- c('lam_max','p_mis_base','p_misseg','k_o_mis','buffer_smax','buffer_beta','buffer_n_exp','p_wgd','alpha_o2','gamma_growth','mu_hp','gamma_mu','O2_crit','n_O')
sid <- as.integer(row$sample_id);stopifnot(samples$sample_id[[sid]]==sid)
for(name in active)stopifnot(identical(as.numeric(row[[name]]),as.numeric(samples[[name]][[sid]])))
fit <- file.path(plan$fit_root,row$fit_seed);proof <- selection$fit_proofs[[row$fit_seed]]
stopifnot(hash(file.path(fit,'best_params.tsv'))==proof$best_params_sha256,hash(file.path(fit,'fit_config.rds'))==proof$fit_config_sha256)
params <- read.delim(file.path(fit,'best_params.tsv'));base <- as.list(setNames(as.numeric(params$value),params$parameter))
raw_cfg <- readRDS(file.path(fit,'fit_config.rds'))
cfg <- prepare_sim_cfg(raw_cfg,list(),0,prepare_run_params(base,'invivo',raw_cfg,0))
stopifnot(cfg$N_MIN==22L,cfg$N_MAX==154L,cfg$N_UNIT==22L,isTRUE(cfg$O2_growth),cfg$ploidy_O2_death=='ploidy_related')
values <- base;for(name in active)values[[name]] <- as.numeric(row[[name]])
rp <- prepare_run_params(values,'invivo',cfg,0)
fm <- fixo2_fixed_matrix(globalenv(),cfg,rp,as.numeric(row$O2_pct));M <- as.matrix(fm$M);ngrid <- fm$ngrid
out <- file.path(root,'points',sprintf('point_%04d',id));dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,'result.json')))stop('Refuse to overwrite completed point')
saveRDS(M,file.path(out,'matrix.rds'));saveRDS(list(run_params=rp,cfg=cfg,values=values,ngrid=ngrid),file.path(out,'context.rds'))
con <- gzfile(file.path(out,'matrix.tsv.gz'),'wt');write.table(format(M,digits=17,scientific=TRUE),con,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE);close(con)
eig <- eigen(M);saveRDS(eig,file.path(out,'eigensystem_double.rds'))
index <- which.max(Re(eig$values));lambda <- eig$values[[index]];raw <- eig$vectors[,index]
v <- Re(raw);sign_flip <- sum(v)<0;if(sign_flip)v <- -v
clipped <- pmax(v,0);denom <- sum(clipped);normalized <- fixo2_normalize_state(v)
stopifnot(all(is.finite(normalized)),abs(sum(normalized)-1)<1e-12)
norm_inf <- max(rowSums(abs(M)));relative <- function(x,l)max(Mod(M%*%x-l*x))/((norm_inf+Mod(l))*max(Mod(x)))
signed_ploidy <- if(sum(v)!=0)sum(ngrid*v)/sum(v)/cfg$N_UNIT else NA_real_
double <- list(lambda_real=Re(lambda),lambda_imaginary=Im(lambda),selected_eigen_index=index,sign_flip=sign_flip,
    raw_sum_real=sum(Re(raw)),raw_sum_imaginary=sum(Im(raw)),signed_sum=sum(v),clipped_sum=denom,
    raw_min=min(v),negative_components=sum(v<0),negative_l1_mass=sum(abs(v[v<0]))/sum(abs(v)),
    signed_ploidy=signed_ploidy,clipped_ploidy=sum(ngrid*clipped)/denom/cfg$N_UNIT,
    normalized_ploidy=sum(ngrid*normalized)/cfg$N_UNIT,raw_relative_residual=relative(raw,lambda),
    signed_relative_residual=relative(v,Re(lambda)),normalized_relative_residual=relative(normalized,Re(lambda)),
    reproduced_ploidy_difference=sum(ngrid*normalized)/cfg$N_UNIT-as.numeric(row$original_ploidy),
    reproduced_growth_difference=Re(lambda)-as.numeric(row$original_growth),
    spectral_gap=Re(lambda)-sort(Re(eig$values),decreasing=TRUE)[[2L]])
vectors <- data.frame(N=ngrid,raw_real=Re(raw),raw_imaginary=Im(raw),signed_real=v,clipped=clipped,normalized=normalized)
con <- gzfile(file.path(out,'vectors_double.tsv.gz'),'wt');write.table(vectors,con,sep='\t',quote=FALSE,row.names=FALSE);close(con)
write_json(data.frame(real=Re(eig$values),imaginary=Im(eig$values)),file.path(out,'eigenvalues_double.json'),digits=17,pretty=TRUE)
reach <- M>0;diag(reach) <- TRUE
for(k in seq_len(nrow(M)))reach <- reach|outer(reach[,k],reach[k,],'&')
classes <- reach&t(reach);components <- length(unique(apply(classes,1L,paste,collapse='')))
elapsed <- list()
for(digits in c(50L,100L)){
    start <- proc.time()[['elapsed']]
    safely <- function(A)tryCatch(review_perron(A,ngrid,cfg$N_UNIT,digits),error=function(e)list(error=conditionMessage(e),converged=FALSE))
    result <- list(right=safely(M),left=safely(t(M)))
    write_json(result,file.path(out,sprintf('precision_%d.json',digits)),digits=17,pretty=TRUE,auto_unbox=TRUE)
    elapsed[[as.character(digits)]] <- proc.time()[['elapsed']]-start
}
files <- list.files(out,full.names=TRUE);names(files) <- basename(files)
report <- list(status='evaluated',identity=row,double=double,
    matrix=list(n=nrow(M),norm_inf=norm_inf,min_off_diagonal=min(M[row(M)!=col(M)]),
                strongly_connected=all(reach),strong_components=components,precision_scope='Multiprecision solve of the unchanged canonical double matrix'),
    elapsed_seconds=elapsed,files_sha256=as.list(vapply(files,hash,'')),
    runtime=list(R=R.version.string,host=Sys.info()[['nodename']],slurm_job_id=Sys.getenv('SLURM_JOB_ID'),
                 code_sha256=selection$code_sha256,selection_sha256=hash(file.path(root,'selection.json'))))
write_json(report,file.path(out,'result.incomplete.json'),auto_unbox=TRUE,digits=17,pretty=TRUE,na='null')
file.rename(file.path(out,'result.incomplete.json'),file.path(out,'result.json'))
cat('Completed point',id,'ploidy double',double$normalized_ploidy,'elapsed',unlist(elapsed),'\n')
