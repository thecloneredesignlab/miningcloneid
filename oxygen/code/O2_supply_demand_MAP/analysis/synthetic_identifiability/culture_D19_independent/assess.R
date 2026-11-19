# Post-selection descriptive assessment only. No simulation, optimization or sampling.
main<-function(){
 args<-commandArgs(TRUE);study<-normalizePath(args[1]);out<-normalizePath(args[2]);release<-file.path(study,'release')
 for(n in c('production_bindings.R','interface.R','feasibility.R','search_helpers.R','core.R'))source(file.path(release,n))
 on.exit(sf_inventory(out),add=TRUE)
 stopifnot(file.exists(file.path(study,'controller/SELECTION_LOCKED')),file.exists(file.path(study,'reports/selection_lock.tsv')))
 lock<-sf_tab(file.path(study,'reports/selection_lock.tsv'));stopifnot(nrow(lock)==1,lock$N20==0,lock$profiles==0,lock$tumor==0)
 lock_sha<-digest::digest(file=file.path(study,'reports/selection_lock.tsv'),algo='sha256')
 env<-sf_backend(study,out,FALSE);case<-readRDS(file.path(study,'inputs/public/D19.rds'));contract<-readRDS(file.path(study,'inputs/public/contract.rds'))
 # Truth is opened only after the exact objective/gate/quality lock above.
 role_table<-sf_tab(file.path(study,'inputs/public/parameter_roles.tsv'));truth<-readRDS(file.path(study,'inputs/private/truth.rds'));ref<-readRDS(file.path(study,'inputs/private/reference.rds'))
 quality<-sf_tab(file.path(study,'reports/verified_candidates.tsv'));selected<-list();labels<-list();summary<-list()
 add<-function(p,label,origin){key<-p$key;if(is.null(selected[[key]]))selected[[key]]<<-p;labels[[length(labels)+1L]]<<-data.frame(origin=origin,label=label,key=key,D=p$score)}
 for(letter in LETTERS[1:6]){
  origin<-paste0('D_',letter);path<-file.path(study,'verification',paste0(origin,'_replay2b'),'state.rds');if(!file.exists(path))next
  s<-readRDS(path);add(s$best,'best_found',origin)
  rows<-quality[quality$origin==origin,,drop=FALSE]
  if(nrow(rows)){
   b<-rows[which.min(rows$D),,drop=FALSE];p<-s$quality_records[[b$key]]$candidate;add(p,'best_verified',origin)
   g<-rows[rows$prediction_pass,,drop=FALSE];if(nrow(g)){r<-g[which.min(g$D),,drop=FALSE];add(s$quality_records[[r$key]]$candidate,'gate_representative',origin)}
  }
  # Keep each exact checked anchor separately from any lower subsequent probe.
  for(pass in 1:2){p<-readRDS(file.path(study,'verification',paste0(origin,'_verify',pass),'checked_representative.rds'));add(p,paste0('checked_representative_pass',pass),origin)}
  summary[[origin]]<-data.frame(origin=origin,total_calls=sum(s$counts),D=s$best$score,best_found_key=s$best$key,global=s$counts['global'],handoff=s$counts['handoff'],local=s$counts['local'],pattern=s$counts['pattern'],verify=s$counts['verify'])
 }
 stopifnot(length(selected)>0);sf_emit(do.call(rbind,labels),file.path(out,'exact_vector_labels.tsv'));sf_emit(do.call(rbind,summary),file.path(out,'origin_search_summary.tsv'))
 parameters<-list();roots<-list();kernels<-list();qrows<-list();grows<-list();functions<-list();survivals<-list()
 grid<-22:154
 kernel<-function(sd){x<-outer(grid,grid,function(o,t)dnorm(o,t,sd));sweep(x,2,colSums(x),'/')}
 K0<-kernel(truth$sigma_kary);stopifnot(max(abs(K0%*%ref$pred$latent_kary[[1]]-env$ivt_smooth_kary_distribution(grid,ref$pred$latent_kary[[1]],truth$sigma_kary)))<=1e-12)
 # Read saved, already verified full states; no forward is performed here.
 for(key in names(selected)){
  q<-quality[quality$key==key&quality$origin==selected[[key]]$origin,,drop=FALSE]
  if(nrow(q))selected[[key]]$pred<-readRDS(file.path(study,q$state_file[1]))$pred
 }
 encountered<-function(pred){
  keys<-character();if(is.null(pred$states))return(keys)
  for(cohort in names(contract$design$adapters))for(segment in contract$design$adapters[[cohort]]$segments){
   state<-pred$states[[paste(cohort,segment$segment_id)]];stopifnot(ncol(state)==133)
   fractions<-state/pmax(rowSums(state),.Machine$double.xmin);hit<-apply(fractions,2,max)>=1e-8
   keys<-c(keys,paste(grid[hit],as.numeric(segment$oxygen_pct),sep='@'))
  };unique(keys)
 }
 et<-encountered(ref$pred);eu<-unique(c(et,unlist(lapply(selected,function(p)encountered(p$pred)),use.names=FALSE)))
 O<-sort(unique(contract$edges$oxygen_pct));common<-sort(unique(c(seq(0,21,.25),O)));fg<-expand.grid(N=grid,oxygen_pct=common);fg$assigned_historical_oxygen<-fg$oxygen_pct%in%O
 fg$encountered_generating<-paste(fg$N,fg$oxygen_pct,sep='@')%in%et;fg$encountered_verified_union<-paste(fg$N,fg$oxygen_pct,sep='@')%in%eu
 sf_emit(data.frame(definition='Any saved trajectory day with state fraction at least1e-8; diagnostic support mask, not a score threshold',generating_pairs=length(et),generating_and_verified_pairs=length(eu)),file.path(out,'encountered_grid_definition.tsv'))
 definitions<-c(selected,list(generating_truth=list(params=truth,origin='generating_truth',key='generating_truth',score=0,pred=ref$pred)))
 for(key in names(definitions)){
  candidate<-definitions[[key]];p<-c(contract$design$fixed,candidate$params[cg_roles]);if(key=='generating_truth')p<-truth
  if(key!='generating_truth'){
   origin<-candidate$origin;id<-readRDS(file.path(study,'verification/T00/identities',paste0(origin,'.rds')));c<-d_context(env,case,contract,id)
   stopifnot(identical(d_decode(c,candidate$unit)[cg_roles],candidate$params),identical(p$sigma_growth,case$known_sigma_growth))
   for(i in seq_len(nrow(contract$design$bounds))){
    r<-contract$design$bounds[i,];name<-r$parameter;v<-p[[name]];t<-truth[[name]];uu<-(env$o2prof_transform(v,r)-env$o2prof_transform(r$lower,r))/(env$o2prof_transform(r$upper,r)-env$o2prof_transform(r$lower,r));tu<-(env$o2prof_transform(t,r)-env$o2prof_transform(r$lower,r))/(env$o2prof_transform(r$upper,r)-env$o2prof_transform(r$lower,r))
    parameters[[length(parameters)+1L]]<-data.frame(origin=origin,key=key,parameter=name,units=role_table$units[match(name,role_table$parameter)],role=if(name=='sigma_growth')'known_control'else if(name%in%c('sigma_kary','init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N'))'unknown_nuisance'else'unknown_mechanism',estimate=v,generating=t,absolute_error=abs(v-t),relative_error=if(t==0)NA_real_ else (v-t)/t,log10_error=if(v>0&&t>0)log10(v/t)else NA_real_,normalized_coordinate=uu,normalized_error=uu-tu,lower_hit=v==r$lower,upper_hit=v==r$upper,near_boundary=uu<=1e-4||uu>=1-1e-4,inactive=name%in%r_activity(p),generating_inactive=name%in%r_activity(truth))
   }
  }
  K<-kernel(p$sigma_kary);kernels[[key]]<-data.frame(key=key,origin=candidate$origin,sigma_kary=p$sigma_kary,max_column_TV_from_identity=max(colSums(abs(K-diag(133)))/2),max_column_TV_from_generating=max(colSums(abs(K-K0))/2),max_element_error_from_generating=max(abs(K-K0)),binary64_identical_to_generating=identical(K,K0))
  for(cohort in c('2N','4N')){
   pi<-env$ivt_gaussian_initial_state(grid,p[[paste0('init_mean_',cohort)]],p[[paste0('init_sd_',cohort)]]);roots[[length(roots)+1L]]<-data.frame(key=key,cohort=cohort,N=grid,probability=pi,mean_parameter=p[[paste0('init_mean_',cohort)]],sd_parameter=p[[paste0('init_sd_',cohort)]])
  }
  f<-fg;f$key<-key;f$lambda<-env$.lambda_eff_of_O2(f$oxygen_pct,p,f$N,O2_growth=contract$design$cfg$O2_growth);f$mu<-env$.mu_eff_of_O2(f$oxygen_pct,p,f$N);f$p_mis_effective<-env$.pmisseg_of_O2(f$oxygen_pct,p,f$N);f$p_wgd<-p$p_wgd;functions[[key]]<-f
  surv<-expand.grid(N=grid,error_size=0:154);surv<-surv[surv$error_size<=surv$N,,drop=FALSE];sN<-p$buffer_smax*exp(-p$buffer_beta*(44/surv$N)^p$buffer_n_exp);surv$survival<-ifelse(surv$error_size==0,1,sN^surv$error_size);surv$key<-key;survivals[[key]]<-surv
  pred<-candidate$pred
  # Prefer the fresh-process full record whenever the selected point was verified.
  if(key!='generating_truth'){
   q<-quality[quality$key==key&quality$origin==candidate$origin,,drop=FALSE];if(nrow(q))pred<-readRDS(file.path(study,q$state_file[1]))$pred
  }
  if(!is.null(pred)){
   grows[[key]]<-data.frame(key=key,passage_id=names(case$data$growth),prediction=unname(pred$growth[names(case$data$growth)]),generating=unname(case$data$growth),normalized_error=unname((pred$growth[names(case$data$growth)]-case$data$growth)/case$known_sigma_growth))
   for(id in names(case$data$prob))qrows[[length(qrows)+1L]]<-data.frame(key=key,passage_id=id,N=grid,observed_probability=pred$kary[[id]],generating_observed_probability=case$data$prob[[id]],latent_probability=if(is.null(pred$latent_kary[[id]]))NA_real_ else pred$latent_kary[[id]],generating_latent_probability=ref$pred$latent_kary[[id]])
  }
 }
 sf_emit(do.call(rbind,parameters),file.path(out,'all20_parameter_recovery.tsv'));sf_emit(do.call(rbind,kernels),file.path(out,'sigma_kary_kernel_diagnostics.tsv'));sf_emit(do.call(rbind,roots),file.path(out,'root_distribution_comparison.tsv'));sf_emit(do.call(rbind,functions),file.path(out,'effective_function_comparison.tsv'));sf_emit(do.call(rbind,survivals),file.path(out,'survival_function_comparison.tsv'));sf_emit(do.call(rbind,grows),file.path(out,'growth_prediction_comparison.tsv'));sf_emit(do.call(rbind,qrows),file.path(out,'full_karyotype_prediction_comparison.tsv'))
 # Summaries separate assigned oxygen conditions from the common extrapolation grid.
 fe<-list();ft<-functions[['generating_truth']]
 for(key in names(selected))for(name in c('lambda','mu','p_mis_effective','p_wgd'))for(scope in c('assigned_oxygen','encountered_verified_union','common_grid')){
  take<-if(scope=='assigned_oxygen')ft$assigned_historical_oxygen else if(scope=='encountered_verified_union')ft$encountered_verified_union else rep(TRUE,nrow(ft));delta<-functions[[key]][[name]][take]-ft[[name]][take]
  fe[[length(fe)+1L]]<-data.frame(key=key,origin=selected[[key]]$origin,response=name,grid_scope=scope,points=sum(take),RMS_absolute_error=sqrt(mean(delta^2)),maximum_absolute_error=max(abs(delta)))
 }
 for(key in names(selected)){delta<-survivals[[key]]$survival-survivals[['generating_truth']]$survival;fe[[length(fe)+1L]]<-data.frame(key=key,origin=selected[[key]]$origin,response='survival_modifier',grid_scope='N22_to154_error0_toN',points=length(delta),RMS_absolute_error=sqrt(mean(delta^2)),maximum_absolute_error=max(abs(delta)))}
 sf_emit(do.call(rbind,fe),file.path(out,'function_error_summary.tsv'));sf_emit(data.frame(scope='Unweighted diagnostic response-function grids; not an oxygen exposure distribution or offspring flux',chromosome_arguments='22:154',error_arguments='all integers0:N',survival_source='buffering_survival_modifier in the pinned production C++ source'),file.path(out,'response_grid_definition.tsv'))
 pdf(file.path(out,'D19_recovery_diagnostics.pdf'),width=11,height=8)
 ps<-do.call(rbind,parameters);unknown<-unique(ps$parameter[ps$role!='known_control']);lab<-do.call(rbind,labels);bk<-unique(lab$key[lab$label=='best_found']);colors<-c('#0072B2','#D55E00','#009E73','#CC79A7','#56B4E9','#8F6B00')[seq_along(bk)]
 mat<-sapply(bk,function(k){a<-ps[ps$key==k,];a$normalized_error[match(unknown,a$parameter)]});if(is.null(dim(mat)))mat<-matrix(mat,ncol=1)
 par(mar=c(10,5,3,1));matplot(seq_along(unknown),mat,type='p',pch=16,col=colors,xaxt='n',xlab='',ylab='Estimate minus truth in normalized coordinate',main='Best found per origin: descriptive coordinate recovery');axis(1,seq_along(unknown),unknown,las=2,cex.axis=.65);abline(h=0,col='grey50');legend('topright',vapply(selected[bk],`[[`,character(1),'origin'),pch=16,col=colors,cex=.8)
 par(mar=c(5,5,3,1));kd<-do.call(rbind,kernels);z<-kd[match(bk,kd$key),];plot(z$sigma_kary,z$max_column_TV_from_generating,pch=16,col=colors,xlab='Fitted sigma_kary (chromosomes)',ylab='Maximum column TV from generating observation kernel',main='Raw width and effective kernel are distinct recovery targets');text(z$sigma_kary,z$max_column_TV_from_generating,vapply(selected[bk],`[[`,character(1),'origin'),pos=3,cex=.8)
 for(key in bk){
  a<-grows[[key]];plot(a$generating,a$prediction,xlab='Generating expected growth (day^-1)',ylab='Fitted expected growth (day^-1)',main=paste('D19',selected[[key]]$origin,substr(key,1,10)));abline(0,1,col='grey50')
 }
 # Full133-state observed laws for the objective-selected global best.
 global_key<-as.character(lock$best_found_key);qq<-do.call(rbind,qrows);par(mfrow=c(3,4),mar=c(3,3,2,1),oma=c(0,0,3,0))
 for(id in names(case$data$prob)){
  a<-qq[qq$key==global_key&qq$passage_id==id,];plot(a$N,a$generating_observed_probability,type='l',col='black',ylim=range(c(a$generating_observed_probability,a$observed_probability)),xlab='Chromosomes',ylab='Probability',main=id,cex.main=.75);lines(a$N,a$observed_probability,col='red',lty=2);if(id==names(case$data$prob)[1])legend('topright',c('Generating','Fitted'),col=c('black','red'),lty=c(1,2),cex=.6,bty='n')
 }
 mtext(paste('Full karyotypes: objective-selected best',selected[[global_key]]$origin,substr(global_key,1,10)),outer=TRUE,line=1,cex=.9)
 par(mfrow=c(2,2),mar=c(4,4,3,1),oma=c(0,0,0,0))
 for(name in c('lambda','mu','p_mis_effective')){
  g<-ft[ft$N==88,];yr<-range(c(g[[name]],unlist(lapply(functions[bk],function(a)a[[name]][a$N==88]),use.names=FALSE)));plot(g$oxygen_pct,g[[name]],type='l',ylim=yr,xlab='Oxygen (%)',ylab=name,main=paste(name,'at N=88'))
  for(j in seq_along(bk)){a<-functions[[bk[j]]];a<-a[a$N==88,];lines(a$oxygen_pct,a[[name]],col=colors[j])}
 }
 st<-survivals[['generating_truth']];st<-st[st$error_size==1,];plot(st$N,st$survival,type='l',ylim=c(0,1),xlab='Survival-function chromosome argument',ylab='Survival modifier',main='One-error survival response')
 for(j in seq_along(bk)){a<-survivals[[bk[j]]];a<-a[a$error_size==1,];lines(a$N,a$survival,col=colors[j])};legend('bottomright',c('Generating',vapply(selected[bk],`[[`,character(1),'origin')),col=c('black',colors),lty=1,cex=.6)
 dev.off()
 stopifnot(identical(lock_sha,digest::digest(file=file.path(study,'reports/selection_lock.tsv'),algo='sha256')))
 sf_emit(data.frame(status='PASS',selection_lock_sha=lock_sha,selected_vectors=length(selected),full_forward_calls=0,truth_use='descriptive only after selection lock',N20=0,profiles=0,tumor=0),file.path(out,'assessment_gate.tsv'));sf_emit(data.frame(status='PASS',full_forward_calls=0),file.path(out,'application_result.tsv'));sf_emit(data.frame(status='PASS'),file.path(out,'gate.tsv'))
}
tryCatch(main(),error=function(e){cat('ASSESSMENT_FAILED: ',conditionMessage(e),'\n',sep='');quit(status=1,save='no')})
