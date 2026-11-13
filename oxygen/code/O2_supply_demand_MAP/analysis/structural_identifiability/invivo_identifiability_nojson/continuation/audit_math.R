# Vectorized independent arithmetic for archived outputs only, never the scorer.
mixture_nll_vectorized<-function(obs,means,fractions,sigma){
 n<-length(obs);k<-length(means)
 if(!n)return(numeric())
 values<-matrix(dnorm(rep(obs,times=k),rep(means,each=n),sigma,log=TRUE),nrow=n,ncol=k)
 values<-sweep(values,2,log(pmax(fractions,1e-300)),'+')
 maximum<-do.call(pmax,unname(as.data.frame(values)))
 -(maximum+log(rowSums(exp(sweep(values,1,maximum,'-')))))
}
independent_components_vectorized<-function(ctx,p,pred){
 rp<-ctx$backend$decode_params(p,FALSE,FALSE,ctx$cfg)
 b<-pp<-nn<-numeric(length(ctx$scenarios));pp2<-pp4<-numeric();nn[]<-NA_real_
 for(i in seq_along(ctx$scenarios)){
  s<-ctx$scenarios[[i]];z<-pred$native[[i]]$result
  keep<-ctx$cfg$scenario_cpp$keep_burden[[i]]&is.finite(s$obs_burden)&s$obs_burden>0
  b[i]<-mean(log(rp$sigma_burden)+.5*log(2*pi)+.5*((log(s$obs_burden[keep])-log(z$Vmm3_total_obs[keep]))/rp$sigma_burden)^2)
  fractions<-as.numeric(z$frac_N_live);fractions<-fractions/sum(fractions)
  pp[i]<-mean(mixture_nll_vectorized(s$endpoint_obs_z,ctx$cfg$model_core$grid_pre,fractions,ctx$cfg$sigma_ploidy))
  if(s$cohort=='2N')pp2<-c(pp2,pp[i]) else pp4<-c(pp4,pp[i])
  if(ctx$cfg$scenario_cpp$keep_necrosis[i]){row<-match(s$harvest,pred$tables$necrosis$harvest);nn[i]<-(pred$tables$necrosis$residual_logit[row]/ctx$cfg$sigma_necrosis_logit)^2}
 }
 c(L_b=mean(b),L_p=.5*mean(pp2)+.5*mean(pp4),L_n=ctx$cfg$lambda_necrosis*mean(nn,na.rm=TRUE))
}
install_validated_audit_math<-function(study,entry){
 path<-file.path(study,'tables/audit_math_validation.tsv')
 if(!file.exists(path))return(FALSE)
 v<-read_kv(path)
 if(!isTRUE(v$pass)||v$input_manifest_sha256!=hash_file(file.path(study,'manifests/input_manifest.tsv'))||v$helper_sha256!=hash_file(file.path(entry,'audit_math.R'))||v$original_adapter_sha256!=hash_file(file.path(entry,'adapter.R')))return(FALSE)
 assign('independent_components',independent_components_vectorized,envir=.GlobalEnv);TRUE
}
