#!/usr/bin/env Rscript
# Change only the chromosome ceiling; retain exact input vectors and all stages.
args <- commandArgs(TRUE);stopifnot(length(args)==2L)
root <- normalizePath(args[[1L]]);id <- as.integer(args[[2L]])
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
here <- dirname(normalizePath(script));efast <- dirname(here);prior_code <- file.path(efast,'high_precision_review')
options(digits=17)
suppressPackageStartupMessages(library(jsonlite));suppressPackageStartupMessages(library(Rcpp))
hash <- function(p) strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' ')[[1L]][[1L]]
selection <- fromJSON(file.path(root,'selection.json'),simplifyVector=FALSE)
stopifnot(hash(file.path(root,'task_manifest.tsv'))==selection$task_manifest_sha256)
for(name in names(selection$code_sha256))stopifnot(hash(file.path(here,name))==selection$code_sha256[[name]])
manifest <- read.delim(file.path(root,'task_manifest.tsv'),check.names=FALSE,colClasses='character')
row <- as.list(manifest[id,,drop=FALSE]);stopifnot(as.integer(row$task_id)==id)
config <- if(row$fit_seed=='seed138')'N8193_M4' else 'N4097_M8'
plan <- selection$source_plans[[config]]
for(name in names(plan$code_sha256))stopifnot(hash(file.path(efast,name))==plan$code_sha256[[name]])
source(file.path(efast,'../../simulation/o2/fixed_o2/run_fixed_o2_simulation.R'),local=TRUE)
sourceCpp(file.path(prior_code,'perron_review.cpp'),cacheDir=Sys.getenv('EFAST_REVIEW_CACHE'),rebuild=FALSE)
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
repo <- normalizePath(file.path(efast,'../../../../..'))
for(name in names(selection$inherited_code_sha256))stopifnot(hash(file.path(repo,name))==selection$inherited_code_sha256[[name]])
prior <- file.path(repo,selection$prior_root,'points',sprintf('point_%04d',as.integer(row$point_id)))
proof <- selection$baseline_proofs[[row$point_id]]
for(name in names(proof))stopifnot(hash(file.path(prior,name))==proof[[name]])
original <- readRDS(file.path(prior,'context.rds'))
values <- base;for(name in active)values[[name]] <- as.numeric(row[[name]])
stopifnot(identical(values,original$values),identical(cfg,original$cfg))
rp <- prepare_run_params(values,'invivo',cfg,0)
stopifnot(identical(rp,original$run_params),identical(rp$boundary,'drop'))
old_cfg <- cfg;cfg$N_MAX <- as.integer(row$N_MAX)
restore_cfg <- cfg;restore_cfg$N_MAX <- old_cfg$N_MAX;stopifnot(identical(restore_cfg,old_cfg))
fm <- fixo2_fixed_matrix(globalenv(),cfg,rp,as.numeric(row$O2_pct));M <- as.matrix(fm$M);ngrid <- fm$ngrid
old_M <- readRDS(file.path(prior,'matrix.rds'));old_n <- nrow(old_M)
common_diff <- max(abs(M[seq_len(old_n),seq_len(old_n)]-old_M))
# The diagnostic calls below do not feed back into the canonical operator.
old_fm <- if(cfg$N_MAX==154L)fm else fixo2_fixed_matrix(globalenv(),old_cfg,rp,as.numeric(row$O2_pct))
stopifnot(identical(old_fm$M,old_M))
shared_rates_exact <- identical(fm$mu_all[seq_len(old_n)],old_fm$mu_all)
stopifnot(shared_rates_exact,identical(fm$ngrid[seq_len(old_n)],original$ngrid))
wide_cfg <- cfg;wide_cfg$N_MAX <- 2L*cfg$N_MAX
wide_G <- o2pr_build_G(globalenv(),wide_cfg,rp,as.numeric(row$O2_pct))
tri <- attr(fm$G,'triplet');wide_tri <- attr(wide_G,'triplet');k <- length(ngrid)
# All possible daughters of an included parent lie at <=2*N_MAX. The wide
# kernel therefore separates upper losses from lower losses without changing
# the eigensolve grid. It is used for transition diagnostics only.
upper <- as.numeric(Matrix::colSums(wide_G[(k+1L):nrow(wide_G),seq_len(k),drop=FALSE]))
lower <- as.numeric(wide_tri$boundary_dropped_rate[seq_len(k)])
stopifnot(max(abs(upper+lower-tri$boundary_dropped_rate))<1e-12)
lam <- .lambda_eff_of_O2(as.numeric(row$O2_pct),rp,N=ngrid,O2_growth=rp$O2_growth)
wgd_upper <- lam*rp$p_wgd*as.numeric(2L*ngrid>cfg$N_MAX)
mis_upper <- upper-wgd_upper;stopifnot(min(mis_upper)> -1e-12);mis_upper <- pmax(mis_upper,0)
flux <- data.frame(N=ngrid,upper_dropped_rate=upper,upper_misseg_dropped_rate=mis_upper,
    upper_wgd_dropped_rate=wgd_upper,lower_dropped_rate=lower,
    total_boundary_dropped_rate=as.numeric(tri$boundary_dropped_rate),
    misseg_nonviable_rate=as.numeric(tri$misseg_nonviable_rate),mu=fm$mu_all,lambda_division=lam)
out <- file.path(root,'points',sprintf('task_%04d',id));dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,'result.json')))stop('Refuse to overwrite completed point')
write.table(flux,file.path(out,'transition_rates.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
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
    if(cfg$N_MAX==154L){
        stopifnot(identical(M,old_M))
        stopifnot(file.copy(file.path(prior,sprintf('precision_%d.json',digits)),file.path(out,sprintf('precision_%d.json',digits))))
    }else{
        result <- list(right=safely(M),left=safely(t(M)))
        write_json(result,file.path(out,sprintf('precision_%d.json',digits)),digits=17,pretty=TRUE,auto_unbox=TRUE)
    }
    elapsed[[as.character(digits)]] <- proc.time()[['elapsed']]-start
}
files <- list.files(out,full.names=TRUE);names(files) <- basename(files)
report <- list(status='evaluated',identity=row,double=double,
    domain=list(N_MAX=cfg$N_MAX,baseline_matrix_exact=identical(old_fm$M,old_M),
        common_block_max_difference=common_diff,rates_exact_on_shared_states=shared_rates_exact,
        precision_reused=cfg$N_MAX==154L,flux_diagnostic_target_ceiling=wide_cfg$N_MAX,
        scope='Wide transition kernel only decomposes losses; eigensolve uses N_MAX from manifest'),
    matrix=list(n=nrow(M),norm_inf=norm_inf,min_off_diagonal=min(M[row(M)!=col(M)]),
                strongly_connected=all(reach),strong_components=components,precision_scope='Multiprecision solve of the unchanged canonical double matrix'),
    elapsed_seconds=elapsed,files_sha256=as.list(vapply(files,hash,'')),
    runtime=list(R=R.version.string,host=Sys.info()[['nodename']],slurm_job_id=Sys.getenv('SLURM_JOB_ID'),
                 code_sha256=selection$code_sha256,selection_sha256=hash(file.path(root,'selection.json'))))
write_json(report,file.path(out,'result.incomplete.json'),auto_unbox=TRUE,digits=17,pretty=TRUE,na='null')
file.rename(file.path(out,'result.incomplete.json'),file.path(out,'result.json'))
cat('Completed boundary task',id,'ploidy double',double$normalized_ploidy,'elapsed',unlist(elapsed),'\n')
