# Meaningful type, constraint, observation and native scoring oracles.
unit <- file.path(study,"tests",paste0(task,"_unit")); dir.create(unit,recursive=TRUE,showWarnings=FALSE)
x <- list(zero=0,negative=-1.2345678901234567,missing=NA_real_,infinite=Inf,false=FALSE,index=12L,label="A quoted field\twith a tab\" and newline\n")
write_kv(x,file.path(unit,"typed_roundtrip.tsv")); stopifnot(identical(read_kv(file.path(unit,"typed_roundtrip.tsv")),x))
atomic_rds(ctx$endpoint,file.path(unit,"parameter_roundtrip.rds")); stopifnot(identical(ctx$endpoint,readRDS(file.path(unit,"parameter_roundtrip.rds"))))
write_tsv(data.frame(parameter=names(ctx$endpoint),value=ctx$endpoint),file.path(unit,"parameter_roundtrip.tsv"))
stopifnot(identical(unname(ctx$endpoint),read_tsv(file.path(unit,"parameter_roundtrip.tsv"))$value))
fixed <- setNames(5.198970004336019,"log10_rho_2N")
free <- setdiff(names(ctx$endpoint),names(fixed)); p <- ctx$endpoint; p[names(fixed)] <- fixed
stopifnot(length(free)==19L,p[[names(fixed)]]==unname(fixed))
data <- load_contract(study,TRUE,"data_only")
prod <- load_contract(study,TRUE,"production")
changed <- names(prod$cfg)[!vapply(names(prod$cfg),function(k) identical(prod$cfg[[k]],data$cfg[[k]]),logical(1))]
stopifnot(identical(changed,"use_soft_prior"))
base_pred <- predictions(prod,prod$endpoint,TRUE)
archived <- readRDS(file.path(study,"tests/repair1_rep1/replay.rds"))
stopifnot(isTRUE(all.equal(base_pred$tables,archived$predictions$tables,tolerance=1e-8)))
base1 <- score(prod,prod$endpoint)
for(i in 1:3) stopifnot(identical(component_values(base1),component_values(score(prod,prod$endpoint))))
a <- score(prod,p); b <- score(data,p)
stopifnot(abs(a$L_data-b$L)<1e-8,b$L_prior==0,identical(a[c("L_b","L_p","L_n")],b[c("L_b","L_p","L_n")]))
pred0 <- predictions(prod,p,FALSE); pred1 <- predictions(prod,p,TRUE)
stopifnot(isTRUE(all.equal(pred0$tables,pred1$tables,tolerance=1e-12)))
ind <- independent_components(prod,p,pred1)
stopifnot(max(abs(ind-unlist(a[names(ind)])))<1e-8)
for(n in c("log10_o2_S0","log10_kappa_O","log10_eta_o2")) {
 q <- p;q[n] <- min(prod$upper[n],p[n]+.002)
 if(q[n]==p[n]) q[n] <- p[n]-.002
 invisible(score(prod,q))
 stopifnot(identical(component_values(a),component_values(score(prod,p))))
}
for(i in seq_along(pred1$native)) {
 z <- pred1$native[[i]]$result
 stopifnot(ncol(z$live_state_obs)==133L,all(is.finite(z$live_state_obs)),all(is.finite(z$dead_hypoxia_state_obs)),all(is.finite(z$dead_buffer_state_obs)))
}
# A simple nuisance-refit oracle distinguishes refitting from fixing nuisances.
for(t in c(0,1,2,3,4)) {
 opt <- optim(1,function(u) (t+u-3)^2+(u-1)^2,method="L-BFGS-B",lower=-5,upper=5)
 stopifnot(abs(opt$par-(2-t/2))<1e-6,abs(opt$value-(t-2)^2/2)<1e-8)
}
write_kv(list(pass=TRUE,typed_roundtrip=TRUE,exact_constraint=TRUE,prior_and_cap_invariance=TRUE,nuisance_refit_oracle=TRUE,native_full_state_output_equivalent=TRUE,independent_score=TRUE,coordinates_full=20L,coordinates_free_fixed_density=19L),file.path(unit,"result.tsv"))
finalize_payload(unit)
