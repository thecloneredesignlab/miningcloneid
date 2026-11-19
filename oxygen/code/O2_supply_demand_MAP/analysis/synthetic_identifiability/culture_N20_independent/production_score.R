# Reuse the production score cross-check from the frozen growth/karyotype recovery adapter.
cg_score <- function(env,design,runs,p,pred,data) {
  entries<-cg_entries(design,data)
  summary<-dplyr::bind_rows(lapply(runs,function(r)env$ivt_collect_lineage_summary(r,entries)))
  gd<-env$ivt_growth_loglik_df(summary,p$sigma_growth)
  kd<-dplyr::bind_rows(lapply(runs,function(r)env$ivt_ploidy_loglik_df(r,entries,p$sigma_kary)))
  a<-mean(gd$loglik);v<-kd$mean_loglik[is.finite(kd$mean_loglik)];k<-mean(v)
  residual<-data$growth[names(pred$growth)]-pred$growth
  ai<-mean(-log(p$sigma_growth)-log(2*pi)/2-residual^2/(2*p$sigma_growth^2))
  ki<-mean(vapply(names(pred$kary),function(id)mean(log(pmax(pred$kary[[id]][match(data$kary[[id]],22:154)],1e-12))),numeric(1)))
  stopifnot(nrow(gd)==114,length(v)==12,is.finite(a+k),abs(a-ai)<=1e-8,abs(k-ki)<=1e-8)
  list(objective=-(a+k),independent=-(ai+ki),growth=a,kary=k,difference=abs(a+k-ai-ki),growth_rows=nrow(gd),kary_rows=length(v))
}
