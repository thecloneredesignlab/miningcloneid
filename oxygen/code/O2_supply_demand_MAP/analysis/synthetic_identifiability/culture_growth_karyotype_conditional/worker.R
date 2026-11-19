# Finite source/interface, replay and observation-generation phases. No fitting.
main<-function(){
 args<-commandArgs(TRUE);study<-normalizePath(args[1]);task_id<-args[2];task<-sub('_R1$','',task_id);out<-file.path(study,'verification',task_id)
 source(file.path(study,'repair_release/production_bindings.R'));source(file.path(study,'repair_release/interface.R'))
 sf_emit(data.frame(stage='before_optional_imports',task=task),file.path(out,'bootstrap_receipt.tsv'));on.exit(sf_inventory(out),add=TRUE)
 env<-sf_backend(study,out,native=task%in%c('R01','R02'))
 checks<-list();check<-function(n,ok,d=''){checks[[length(checks)+1L]]<<-data.frame(check=n,passed=isTRUE(ok),detail=d);sf_emit(do.call(rbind,checks),file.path(out,'checks.tsv'));if(!isTRUE(ok))stop('Check failed: ',n)}
 failed<-function(f)isTRUE(tryCatch({f();FALSE},error=function(e)TRUE))
 finish<-function(status='PASS'){
   sf_emit(data.frame(status=status,task=task,fits=0,profiles=0,tumor=0),file.path(out,'gate.tsv'))
   sf_emit(data.frame(status='PASS',task=task,optimizer_calls=0),file.path(out,'application_result.tsv'));cat('CGK_PHASE_PASS ',task,'\n',sep='')
 }
 if(task=='B01'){
   bundle<-readRDS(file.path(study,'verification/D01/evaluation_bundle.rds'));design<-bundle$design
   check('public_schema',identical(sort(names(bundle)),sort(c('design','observations','scope'))))
   neutral<-cg_neutral(env,design$bounds);candidate<-cg_candidate(design,neutral)
   check('all20_explicit_candidate',all(cg_roles%in%names(candidate)))
   check('reject_missing_unknown',failed(function()cg_candidate(design,neutral[-1])))
   private<-file.path(study,'inputs/private/anchor_bundle.rds')
   check('OS_truth_read_denied',failed(function()readRDS(private)))
   denied<-readLines(file.path(out,'read_denials.tsv'));check('only_registered_truth_probe',length(denied)>=1&&all(denied==private))
   check('dataset_has_no_latent_truth',setequal(names(bundle$observations),c('growth','kary')))
   check('no_active_fixed_values',!any(cg_roles%in%names(design$fixed)))
   for(n in cg_roles){x<-neutral;x[[n]]<-design$bounds$lower[match(n,design$bounds$parameter)];check(paste0('explicit_injection_',n),cg_candidate(design,x)[[n]]==x[[n]])}
   writeLines('The scoped evaluator read the public bundle and injected20 explicit coordinates. Its direct readRDS attempt on the private generating object was denied by the separately loaded native read guard. Only the registered denied probe was accepted. The guard also denies earlier private phase outputs and the upstream Step5 tree. No optimizer or model propagation ran in this phase.',file.path(out,'truth_access_tests.md'))
   finish();return(invisible(NULL))
 }
 if(task=='I01'){
   b<-readRDS(file.path(study,'inputs/private/anchor_bundle.rds'));design<-cg_design(env,b);t<-b$map$table
   check('20_active_original_roles',setequal(t$parameter,cg_roles)&&nrow(t)==20&&all(t$status=='active')&&length(b$map$fixed)==0)
   check('schema_counts',nrow(design$mask)==131&&sum(design$mask$growth)==114&&sum(design$mask$kary)==12&&sum(design$mask$n_kary)==220)
   check('all122_segments',sum(vapply(design$adapters,function(a)length(a$segments),integer(1)))==122)
   expected<-sort(c(12,15,16,18,19,rep(20,7)));check('sample_specific_sizes',identical(sort(as.numeric(design$mask$n_kary[design$mask$kary])),expected))
   z<-setNames(vapply(seq_len(nrow(t)),function(i)env$o2prof_transform(b$p[[t$parameter[i]]],t[i,,drop=FALSE]),numeric(1)),t$coordinate)
   p<-cg_from_free(env,design,z);check('20_roundtrip',max(abs(vapply(cg_roles,function(n)p[[n]]-b$p[[n]],numeric(1))))<=1e-12)
   for(which_bound in c('lower','upper')) {
     endpoint<-setNames(env$o2prof_transform(t[[which_bound]],t),t$coordinate)
     back<-cg_from_free(env,design,endpoint)
     check(paste0('all20_endpoint_roundtrip_',which_bound),all(vapply(seq_len(nrow(t)),function(i){v<-back[[t$parameter[i]]];v>=t$lower[i]&&v<=t$upper[i]&&abs(v-t[[which_bound]][i])<=64*.Machine$double.eps*max(1,abs(t$lower[i]),abs(t$upper[i]))},logical(1))))
   }
   neutral<-cg_neutral(env,design$bounds)
   for(n in cg_roles){x<-neutral;i<-match(n,t$parameter);x[[n]]<-t$lower[i];got<-cg_candidate(design,x);check(paste0('injection_',n),got[[n]]==t$lower[i])}
   check('missing_coordinate_rejected',failed(function()cg_candidate(design,neutral[-1])))
   check('duplicate_coordinate_rejected',failed(function(){x<-neutral;names(x)[1]<-names(x)[2];cg_candidate(design,x)}))
   for(n in c('sigma_growth','sigma_kary'))check(paste0('nonpositive_',n),failed(function(){x<-neutral;x[[n]]<-0;cg_candidate(design,x)}))
   for(n in c('p_misseg','buffer_beta')){x<-neutral;x[[n]]<-0;check(paste0('exact_zero_',n),cg_candidate(design,x)[[n]]==0)}
   for(cfg in c('N0_known6','Ninit','Nerror','Nall')){
     known<-switch(cfg,N0_known6=setdiff(cg_roles,cg_roles[!grepl('^init_|^sigma_',cg_roles)]),Ninit=c('sigma_growth','sigma_kary'),Nerror=grep('^init_',cg_roles,value=TRUE),Nall=character())
     map<-env$o2prof_map(file.path(study,'inputs/private/parameter_table_input.csv'),b$p,'zero_extended',fixed=b$p[known]);want<-c(N0_known6=14,Ninit=18,Nerror=16,Nall=20)[cfg]
     check(paste0('configuration_',cfg),sum(map$table$status=='active')==want)
   }
   m<-env$o2prof_map(file.path(study,'inputs/private/parameter_table_input.csv'),b$p,'zero_extended',fixed=list(p_misseg=0));check('null_induction_inactive_k',m$table$status[m$table$parameter=='k_o_mis']=='inactive')
   m<-env$o2prof_map(file.path(study,'inputs/private/parameter_table_input.csv'),b$p,'zero_extended',fixed=list(buffer_beta=0));check('null_survival_inactive_exponent',m$table$status[m$table$parameter=='buffer_n_exp']=='inactive')
   check('mixed_run_settings_filter20',length(b$p)==36&&length(b$p[cg_roles])==20&&all(vapply(b$p[cg_roles],is.numeric,logical(1))))
   check('no_observed_response_in_design',!any(c('objects','p','map','truth','best_params','fit_data')%in%names(design)))
   check('no_observed_slots_in_adapters',all(vapply(design$adapters,function(a)all(vapply(a$segments,function(s)!'observed'%in%names(s),logical(1))),logical(1))))
   hist<-cg_historical(env,b);entries<-cg_entries(design,hist)
   check('masked_entries_no_flow',all(vapply(entries,function(x)is.null(x$flow),logical(1))))
   check('conditioning_separate_from_response',identical(design$conditioning,cg_design(env,b)$conditioning))
   sf_rds(design,file.path(out,'public_design.rds'));sf_emit(design$mask,file.path(out,'observation_mask.tsv'));sf_emit(design$conditioning,file.path(out,'conditioning_manifest.tsv'));sf_emit(design$bounds,file.path(out,'public_bounds.tsv'))
   sf_rds(b$p,file.path(out,'truth_parameters.rds'));sf_emit(data.frame(parameter=cg_roles,value=vapply(b$p[cg_roles],as.character,character(1))),file.path(out,'private_truth_registry.tsv'))
   code<-system2('gcc',c('-shared','-fPIC','-O2',shQuote(file.path(study,'repair_release/read_guard.c')),'-ldl','-o',shQuote(file.path(out,'read_guard.so'))))
   check('read_guard_compiled',code==0&&file.exists(file.path(out,'read_guard.so')))
   sf_emit(data.frame(configuration=c('N0_known6','Ninit','Nerror','Nall'),free=c(14,18,16,20),fits=0),file.path(out,'nuisance_configurations.tsv'))
   finish();return(invisible(NULL))
 }
 design<-readRDS(file.path(study,'verification/I01_R1/public_design.rds'))
 if(task=='P01'){
   base<-readRDS(file.path(study,'inputs/private/reference/baseline.rds'));all<-readRDS(file.path(study,'inputs/private/reference/bridge_results.rds'));op<-readRDS(file.path(study,'inputs/private/reference/observer_perturbations.rds'))
   keep<-!startsWith(names(base$pred$vector),'flow:');v<-base$pred$vector[keep];check('1698_projection',length(v)==1698)
   emit<-function(v,path)sf_emit(data.frame(output=names(v),value=unname(v)),path)
   emit(v,file.path(out,'baseline_vector.tsv'));rows<-list();bf<-sf_hash(lapply(base$pred$kary,function(q)q<1e-12))
   for(id in names(all)){r<-all[[id]];check(paste0('returned_',id),r$status=='returned');a<-r$pred$vector[keep];check(paste0('projection_',id),identical(names(a),names(v)));emit(a,file.path(out,paste0(id,'_vector.tsv')));rows[[length(rows)+1L]]<-data.frame(eval_id=id,kary_floor_switch=sf_hash(lapply(r$pred$kary,function(q)q<1e-12))!=bf)}
   for(i in seq_along(op)){a<-op[[i]]$pred$vector[keep];emit(a,file.path(out,paste0('sigma_',i,'_vector.tsv')));check(paste0('sigma_flat_',i),max(abs(a-v))==0)}
   sf_emit(do.call(rbind,rows),file.path(out,'kary_floor_switches.tsv'));finish();return(invisible(NULL))
 }
 if(task%in%c('R01','R02')){
   p<-readRDS(file.path(study,'verification/I01_R1/truth_parameters.rds'));b<-readRDS(file.path(study,'inputs/private/anchor_bundle.rds'));hist<-cg_historical(env,b)
   r<-cg_forward(study,out,env,design,p[cg_roles],task,hist);old<-readRDS(file.path(study,'inputs/private/reference/baseline.rds'))
   check('original_selection',identical(r$pred$selection,old$pred$selection))
   gap<-max(vapply(names(old$pred$states),function(n)max(abs(r$pred$states[[n]]-old$pred$states[[n]])/pmax(1,abs(r$pred$states[[n]]),abs(old$pred$states[[n]]))),numeric(1)))
   ktv<-max(vapply(names(old$pred$kary),function(n)sum(abs(r$pred$kary[[n]]-old$pred$kary[[n]]))/2,numeric(1)))
   ggap<-max(abs(r$pred$growth[names(old$pred$growth)]-old$pred$growth))
   check('source_state_agreement',gap<=1e-10,as.character(gap));check('complete_karyotype_agreement',ktv<=1e-10);check('growth_agreement',ggap<=1e-10)
   check('masked_reference_score',abs(r$score$objective+old$score$growth+old$score$kary)<=1e-8)
   if(task=='R02'){
     first<-readRDS(file.path(study,'verification/R01/result.rds'));check('fresh_process_selection',identical(r$pred$selection,first$pred$selection));check('fresh_process_predictions',isTRUE(all.equal(r$pred,first$pred,tolerance=1e-10)));check('fresh_process_score',abs(r$score$objective-first$score$objective)<=1e-8)
   }
   sf_rds(r,file.path(out,'result.rds'));sf_emit(r$pred$selection,file.path(out,'selection.tsv'));sf_emit(data.frame(state_scaled_gap=gap,kary_tv=ktv,growth_gap=ggap,score_difference=r$score$difference),file.path(out,'replay_checks.tsv'));sf_emit(data.frame(forward_seconds=r$seconds,includes_compilation=FALSE),file.path(out,'cost_measurement.tsv'));finish();return(invisible(NULL))
 }
 r<-readRDS(file.path(study,'verification/R01/result.rds'));p<-readRDS(file.path(study,'verification/I01_R1/truth_parameters.rds'))
 if(task%in%c('Q01','Q02')){
   if(task=='Q01'){
     states<-cg_states();sf_rds(states,file.path(out,'RNG_states.rds'));sf_emit(data.frame(ordinal=1:1025,kind=c(rep('QA',1024),'primary'),state_sha=vapply(states,sf_hash,character(1))),file.path(out,'RNG_manifest.tsv'))
     # Eight registered standalone fixture passes, never new propagation.
     fixtures<-data.frame(id=c('baseline','sigma_growth_half','sigma_growth_twice','sigma_kary_width1',paste0('sigma_local_',1:4)),status='registered');sf_emit(fixtures,file.path(out,'fixture_attempts.tsv'))
     hist<-cg_historical(env,readRDS(file.path(study,'inputs/private/anchor_bundle.rds')))
     for(i in seq_len(nrow(fixtures))){fixtures$status[i]<-'started';sf_emit(fixtures,file.path(out,'fixture_attempts.tsv'));pp<-p;pred<-r$pred
       if(i%in%2:3)pp$sigma_growth<-p$sigma_growth*c(0.5,2)[i-1]
       if(i>=4){pp$sigma_kary<-if(i==4)1 else sf_tab(file.path(study,'inputs/private/reference/observer_points.tsv'))$new_value[i-4];pred$kary<-lapply(pred$latent_kary,function(pi)env$ivt_smooth_kary_distribution(22:154,pi,pp$sigma_kary))}
       sc<-cg_score(env,design,r$runs,pp,pred,hist);check(paste0('fixture_score_',i),sc$difference<=1e-8)
       if(i>=5)check(paste0('local_sigma_flat_',i),identical(pred$kary,r$pred$kary))
       if(i==4)check('nonflat_width_control',max(abs(unlist(pred$kary)-unlist(r$pred$kary)))>1e-8)
       fixtures$status[i]<-'returned';sf_emit(fixtures,file.path(out,'fixture_attempts.tsv'))
     }
     residual<-p$sigma_growth*c(-1,0,1);for(sd in p$sigma_growth*c(.5,1,2)){a<-env$ivt_growth_loglik_df(data.frame(observed_growth=residual,predicted_growth_rate=0),sd)$loglik;z<--log(sd)-log(2*pi)/2-residual^2/(2*sd^2);check(paste0('literal_variance_',sd),max(abs(a-z))<1e-12)}
     acc<-list(n=0L,g_sum=setNames(numeric(114),names(r$pred$growth)),g_sumsq=setNames(numeric(114),names(r$pred$growth)),k=lapply(r$pred$kary,function(q)integer(133)));ord<-1:64
   }else{states<-readRDS(file.path(study,'verification/Q01/RNG_states.rds'));acc<-readRDS(file.path(study,'verification/Q01/QA_accumulator.rds'));ord<-65:1024}
   start<-proc.time()[3]
   for(i in ord){id<-sprintf('QA%04d',i);ledger<-cg_observation_start(study,out,id);data<-cg_draw(r$pred,design,p$sigma_growth,states[[i]]);cg_data_valid(data,design)
     g<-data$growth[names(acc$g_sum)];acc$g_sum<-acc$g_sum+g;acc$g_sumsq<-acc$g_sumsq+g*g
     for(n in names(acc$k))acc$k[[n]]<-acc$k[[n]]+tabulate(match(data$kary[[n]],22:154),nbins=133)
     acc$n<-acc$n+1L;sf_rds(acc,file.path(out,'QA_accumulator.rds'));ledger$status[nrow(ledger)]<-'returned';ledger$payload_hash[nrow(ledger)]<-sf_hash(data);sf_emit(ledger,file.path(out,'observation_attempts.tsv'))
   }
   sf_emit(data.frame(realizations=length(ord),seconds=proc.time()[3]-start,model_calls=0),file.path(out,'cost_measurement.tsv'))
   if(task=='Q02'){
     Rn<-acc$n;check('1024_QA_distinct',Rn==1024);alpha<-.01/1824;rows<-list()
     means<-acc$g_sum/Rn;vars<-(acc$g_sumsq-acc$g_sum^2/Rn)/(Rn-1)
     for(id in names(r$pred$growth)){
       z<-(means[id]-r$pred$growth[id])/(p$sigma_growth/sqrt(Rn));v<-(Rn-1)*vars[id]/p$sigma_growth^2
       rows[[length(rows)+1L]]<-data.frame(check=paste0('growth_mean:',id),value=unname(z),lower=qnorm(alpha/2),upper=qnorm(1-alpha/2))
       rows[[length(rows)+1L]]<-data.frame(check=paste0('growth_variance:',id),value=unname(v),lower=qchisq(alpha/2,Rn-1),upper=qchisq(1-alpha/2,Rn-1))
     }
     for(id in names(acc$k)){n<-design$mask$n_kary[match(id,design$mask$passage_id)]*Rn;q<-r$pred$kary[[id]]
       # Source probabilities normalize to floating precision; use the same normalized sampling law.
       q<-q/sum(q)
       for(j in 1:133)rows[[length(rows)+1L]]<-data.frame(check=paste0('kary:',id,':',j+21),value=acc$k[[id]][j],lower=qbinom(alpha/2,n,q[j]),upper=qbinom(1-alpha/2,n,q[j]))
       check(paste0('multinomial_total_',id),sum(acc$k[[id]])==n)
     }
     qa<-do.call(rbind,rows);qa$passed<-qa$value>=qa$lower&qa$value<=qa$upper;sf_emit(qa,file.path(out,'generator_QA.tsv'));check('all1824_predeclared_bands',nrow(qa)==1824&&all(qa$passed),paste(sum(!qa$passed),'failed'))
   }
   finish();return(invisible(NULL))
 }
 if(task%in%c('D01','D02')){
   states<-readRDS(file.path(study,'verification/Q01/RNG_states.rds'));ledger<-cg_observation_start(study,out,task)
   data<-cg_draw(r$pred,design,p$sigma_growth,states[[1025]]);check('finite_exact_sample_schema',cg_data_valid(data,design));hash<-sf_hash(data)
   ledger$status[nrow(ledger)]<-'returned';ledger$payload_hash[nrow(ledger)]<-hash;sf_emit(ledger,file.path(out,'observation_attempts.tsv'))
   if(task=='D02')check('same_stream_independent_reproduction',identical(data,readRDS(file.path(study,'verification/D01/synthetic_dataset.rds'))))
   sc<-cg_score(env,design,r$runs,p,r$pred,data);check('synthetic_same_vector_score',sc$difference<=1e-8)
   sf_rds(data,file.path(out,'synthetic_dataset.rds'));sf_emit(data.frame(payload_sha=hash,negative_growth=sum(data$growth<0),growth=length(data$growth),kary_samples=length(data$kary),kary_cells=sum(lengths(data$kary)),fits=0),file.path(out,'dataset_receipt.tsv'))
   if(task=='D01'){
     sf_rds(list(design=design,observations=data,scope='conditional artificial response study; no identified-parameter claim'),file.path(out,'evaluation_bundle.rds'))
     sf_emit(data.frame(passage_id=names(data$growth),growth=data$growth),file.path(out,'synthetic_growth.tsv'))
     sf_emit(do.call(rbind,lapply(names(data$kary),function(id)data.frame(passage_id=id,cell_index=seq_along(data$kary[[id]]),observed_chromosomes=data$kary[[id]]))),file.path(out,'synthetic_karyotype_cells.tsv'))
     sf_emit(design$mask,file.path(out,'synthetic_observation_manifest.tsv'))
   }
   finish();return(invisible(NULL))
 }
 stop('Unknown bounded task')
}
main()
