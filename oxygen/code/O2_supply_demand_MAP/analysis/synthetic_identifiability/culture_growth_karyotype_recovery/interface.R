# Pure conditional design, explicit candidate injection and masked observation law.
cg_roles <- c('lam_max','p_misseg','k_o_mis','buffer_smax','buffer_beta','buffer_n_exp','p_wgd','alpha_o2','gamma_growth','mu_hp','gamma_mu','O2_crit','n_O','p_mis_base','sigma_growth','sigma_kary','init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N')
cg_bounds <- function(b)b$map$table[,c('parameter','lower','upper','transform','scale','coordinate'),drop=FALSE]
cg_neutral <- function(env,tab) {
  setNames(lapply(seq_len(nrow(tab)),function(i){r<-tab[i,,drop=FALSE];env$o2prof_transform((env$o2prof_transform(r$lower,r)+env$o2prof_transform(r$upper,r))/2,r,inverse=TRUE)}),tab$parameter)
}
cg_candidate <- function(design,candidate) {
  stopifnot(is.list(candidate),length(candidate)==20,!anyDuplicated(names(candidate)),setequal(names(candidate),cg_roles))
  for(i in seq_len(nrow(design$bounds))) {
    r<-design$bounds[i,];x<-candidate[[r$parameter]]
    stopifnot(is.numeric(x),length(x)==1,is.finite(x),x>=r$lower,x<=r$upper)
  }
  c(design$fixed,candidate[cg_roles])
}
cg_from_free <- function(env,design,z) {
  t<-design$bounds;stopifnot(length(z)==20,!anyDuplicated(names(z)),setequal(names(z),t$coordinate))
  # The original helper owns endpoint-roundoff semantics; its base is independent neutral.
  t$status<-'active';t$fixed_value<-NA_real_
  map<-list(table=t,base_params=c(design$fixed,cg_neutral(env,t)))
  values<-env$o2prof_from_free(map,z[t$coordinate])
  cg_candidate(design,values[cg_roles])
}
cg_conditioning <- function(objects) {
  do.call(rbind,lapply(names(objects$fit_data),function(id){o<-objects$fit_data[[id]];data.frame(passage_id=id,passage_duration=as.numeric(o$passage_duration)[1],initial_cells=as.numeric(o$initial_cells)[1],final_cells=as.numeric(o$final_cells)[1])}))
}
cg_design <- function(env,b) {
  t<-cg_bounds(b);neutral<-cg_neutral(env,t);cfg<-b$cfg
  for(n in cg_roles) {
    if(n%in%names(cfg))cfg[[n]]<-neutral[[n]]
    key<-paste0(n,'_init');if(key%in%names(cfg))cfg[[key]]<-neutral[[n]]
  }
  cfg$o2_crit_init<-neutral$O2_crit
  mask<-do.call(rbind,lapply(names(b$objects$fit_data),function(id){o<-env$ivt_observed_passage_summary(b$objects$fit_data[[id]]);data.frame(passage_id=id,growth=is.finite(o$observed_growth),kary=length(o$observed_kary)>0,n_kary=length(o$observed_kary))}))
  adapters<-lapply(b$adapters,function(a)list(cohort=a$cohort,segments=lapply(a$segments,function(s)s[setdiff(names(s),'observed')])))
  list(bounds=t,fixed=b$p[setdiff(names(b$p),cg_roles)],cfg=cfg,adapters=adapters,conditioning=cg_conditioning(b$objects),mask=mask,model='production v2;133 live states;ploidy_related',scope='artificial fixed-input endpoint growth and finite integer karyotype',weights=c(growth=1,karyotype=1,flow=0,death=0,prior=0))
}
cg_historical <- function(env,b) {
  obs<-lapply(b$objects$fit_data,env$ivt_observed_passage_summary)
  g<-vapply(obs,function(x)x$observed_growth,numeric(1));g<-g[is.finite(g)]
  k<-lapply(obs,function(x)x$observed_kary);k<-k[lengths(k)>0]
  list(growth=g,kary=k)
}
cg_entries <- function(design,data) {
  setNames(lapply(seq_len(nrow(design$conditioning)),function(i){r<-design$conditioning[i,];id<-r$passage_id;list(g=if(id%in%names(data$growth))unname(data$growth[id])else NA_real_,kary=if(id%in%names(data$kary))data$kary[[id]]else numeric(),flow=NULL,passage_duration=r$passage_duration,initial_cells=r$initial_cells,final_cells=r$final_cells)}),design$conditioning$passage_id)
}
cg_predict <- function(env,design,runs,p) {
  g<-numeric();k<-list();latent<-list();selection<-list();states<-list()
  for(cohort in names(runs))for(s in runs[[cohort]]$segment_results) {
    pi<-as.numeric(s$selection$selected_frac);stopifnot(all(is.finite(pi)),min(pi)>=0,abs(sum(pi)-1)<=1e-12)
    rate<-env$ivt_log_growth_rate(s$sim$Ntot_live_obs[1],s$selection$selected_live_cells,s$selection$selected_day)
    for(id in s$segment$data_ids) {
      m<-design$mask[match(id,design$mask$passage_id),];stopifnot(nrow(m)==1,!is.na(m$growth))
      if(m$growth)g[id]<-rate
      if(m$kary){q<-env$ivt_smooth_kary_distribution(runs[[cohort]]$grid_pre,pi,p$sigma_kary);stopifnot(all(is.finite(q)),min(q)>=0,abs(sum(q)-1)<=1e-12);k[[id]]<-q;latent[[id]]<-pi}
    }
    parent<-s$parent_boundary_selection
    selection[[length(selection)+1L]]<-data.frame(cohort=cohort,segment_id=s$segment$segment_id,index=s$selection$selected_index,day=s$selection$selected_day,parent_index=if(is.null(parent))NA_integer_ else parent$selected_index,parent_day=if(is.null(parent))NA_real_ else parent$selected_day)
    states[[paste(cohort,s$segment$segment_id)]]<-s$sim$live_state_obs
  }
  stopifnot(length(g)==114,length(k)==12)
  list(growth=g,kary=k,latent_kary=latent,selection=do.call(rbind,selection),states=states)
}
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
cg_vector <- function(pred,sd) {
  v<-pred$growth/sd;names(v)<-paste0('growth:',names(pred$growth))
  for(id in names(pred$kary))v<-c(v,setNames(pred$kary[[id]][-133],paste0('kary:',id,':',1:132)))
  stopifnot(length(v)==1698);v
}
