# Finite production diagnostics only. Never generates observations or fits.
main <- function() {
  args<-commandArgs(TRUE);study<-normalizePath(args[1]);task<-args[2]
  out<-file.path(study,'verification',task)
  source(file.path(study,'release/api.R'))
  sf_emit(data.frame(stage='before_optional_imports',task=task),file.path(out,'bootstrap_receipt.tsv'))
  on.exit(sf_inventory(out),add=TRUE)
  env<-sf_backend(study,out,native=task%in%c('M01','R01','R02','R03'))
  checks<-list()
  check<-function(name,condition,detail='') {
    checks[[length(checks)+1L]]<<-data.frame(check=name,passed=isTRUE(condition),detail=detail)
    sf_emit(do.call(rbind,checks),file.path(out,'checks.tsv'));if(!isTRUE(condition))stop('Check failed: ',name)
  }
  empty_ledger<-function()data.frame(eval_id=character(),status=character(),parameter_sha=character(),seconds=double())
  if(task=='M01') {
    b<-sf_context(study,env);tab<-b$map$table
    check('all20_unknown_mapping',nrow(tab)==20 && sum(tab$status=='active')==20 && length(b$map$fixed)==0)
    roles<-sf_tab(file.path(study,'inputs/parameter_roles.tsv'));check('role_identity',setequal(tab$parameter,roles$parameter))
    check('full_precision_bounds_match_metadata',max(abs(tab$lower-as.numeric(roles$lower[match(tab$parameter,roles$parameter)])))<=1e-6 && max(abs(tab$upper-as.numeric(roles$upper[match(tab$parameter,roles$parameter)])))<=1e-6)
    u<-sf_unit(env,b,b$p);back<-sf_params(env,b,u)
    check('20_coordinate_roundtrip',max(abs(vapply(tab$parameter,function(n)back[[n]]-b$p[[n]],numeric(1))))<=1e-12)
    for(n in tab$parameter) {
      v<-u;v[n]<-if(u[n]+0.001<=1)u[n]+0.001 else u[n]-0.001
      changed<-sf_params(env,b,v)
      check(paste0('injection_',n),changed[[n]]!=b$p[[n]] && max(abs(vapply(setdiff(tab$parameter,n),function(k)changed[[k]]-b$p[[k]],numeric(1))))<=1e-12)
    }
    nullmap<-env$o2prof_map(file.path(study,'inputs/data/parameter_table_input.csv'),b$p,'zero_extended',fixed=list(p_misseg=0))
    check('exact_zero_induction_and_inactive_k',nullmap$fixed$p_misseg==0 && nullmap$table$status[nullmap$table$parameter=='k_o_mis']=='inactive')
    nullmap<-env$o2prof_map(file.path(study,'inputs/data/parameter_table_input.csv'),b$p,'zero_extended',fixed=list(buffer_beta=0))
    check('exact_zero_survival_and_inactive_nexp',nullmap$fixed$buffer_beta==0 && nullmap$table$status[nullmap$table$parameter=='buffer_n_exp']=='inactive')
    check('native_simulator_interface',identical(names(formals(env$cpp_o2simps_simulate_one)),'sim_args'))
    descriptors<-lapply(b$objects$fit_data,env$ivt_observed_passage_summary)
    check('actual_observation_counts',sum(vapply(descriptors,function(o)is.finite(o$observed_growth),logical(1)))==114 && sum(vapply(descriptors,function(o)length(o$observed_kary)>0,logical(1)))==12 && sum(vapply(descriptors,function(o)length(o$observed_kary),integer(1)))==220 && sum(vapply(descriptors,function(o)!is.null(o$observed_flow)&&nrow(o$observed_flow)>0,logical(1)))==20)
    sf_emit(tab,file.path(out,'parameter_map20.tsv'));sf_emit(b$projection,file.path(out,'requested_actual_anchor.tsv'))
    sf_emit(data.frame(parameter=names(b$p),value=vapply(b$p,as.character,character(1))),file.path(out,'anchor_parameters.tsv'))
    points<-list();observer<-list();ordinal<-0L
    for(n in tab$parameter[!tab$parameter%in%c('sigma_growth','sigma_kary')])for(h in c(-0.01,0.01,-0.005,0.005)) {
      ordinal<-ordinal+1L;i<-match(n,tab$parameter);unit<-u[n]+h;admissible<-unit>=0&&unit<=1
      natural<-if(admissible)env$o2prof_transform(env$o2prof_transform(tab$lower[i],tab[i,,drop=FALSE])+unit*(env$o2prof_transform(tab$upper[i],tab[i,,drop=FALSE])-env$o2prof_transform(tab$lower[i],tab[i,,drop=FALSE])),tab[i,,drop=FALSE],inverse=TRUE)else NA_real_
      points[[ordinal]]<-data.frame(eval_id=sprintf('E%03d',ordinal),parameter=n,offset=h,base_unit=unname(u[n]),new_unit=unname(unit),new_value=natural,status=if(admissible)'registered' else 'SKIP_OUTSIDE_DOMAIN')
    }
    i<-match('sigma_kary',tab$parameter)
    for(h in c(-0.01,0.01,-0.005,0.005)) {
      unit<-u['sigma_kary']+h;z<-env$o2prof_transform(tab$lower[i],tab[i,,drop=FALSE])+unit*(env$o2prof_transform(tab$upper[i],tab[i,,drop=FALSE])-env$o2prof_transform(tab$lower[i],tab[i,,drop=FALSE]))
      observer[[length(observer)+1L]]<-data.frame(offset=h,new_value=env$o2prof_transform(z,tab[i,,drop=FALSE],inverse=TRUE),status=if(unit>=0&&unit<=1)'registered' else 'SKIP_OUTSIDE_DOMAIN')
    }
    sf_emit(do.call(rbind,points),file.path(out,'bridge_points.tsv'));sf_emit(do.call(rbind,observer),file.path(out,'observer_points.tsv'))
    sf_rds(b,file.path(out,'anchor_bundle.rds'))
    sf_emit(data.frame(status='PASS',coordinates=20,full_forward_calls=0,synthetic_observations=0,fit_calls=0),file.path(out,'gate.tsv'))
  } else if(task%in%c('R01','R02','R03')) {
    b<-readRDS(file.path(study,'verification/M01/anchor_bundle.rds'));ledger<-empty_ledger()
    if(task!='R03') {
      result<-sf_run(env,b,b$p,out,task,ledger);check('anchor_protocol_feasible',result$result$status=='returned')
      sf_rds(result$result,file.path(out,'result.rds'));sf_emit(result$result$pred$selection,file.path(out,'selection.tsv'))
      sf_emit(as.data.frame(result$result$branch,stringsAsFactors=FALSE),file.path(out,'branch.tsv'))
      sf_emit(data.frame(output=names(result$result$pred$vector),value=unname(result$result$pred$vector)),file.path(out,'mean_vector.tsv'))
      sf_emit(data.frame(quantity=names(result$result$score),value=unlist(result$result$score)),file.path(out,'score.tsv'))
      sf_emit(data.frame(forward_seconds=result$result$seconds,phase=task,includes_compilation=FALSE),file.path(out,'cost_measurement.tsv'))
      if(task=='R02') {
        baseline<-readRDS(file.path(study,'verification/R01/result.rds'));r<-result$result
        check('same_vector_score',abs(r$score$objective-baseline$score$objective)<=1e-8)
        check('selection_and_boundary_identity',identical(r$pred$selection,baseline$pred$selection))
        state_gap<-max(vapply(names(r$pred$states),function(n)max(abs(r$pred$states[[n]]-baseline$pred$states[[n]])/pmax(1,abs(baseline$pred$states[[n]]))),numeric(1)))
        check('state_repeat',state_gap<=1e-10,as.character(state_gap))
        k_tv<-max(vapply(names(r$pred$kary),function(n)sum(abs(r$pred$kary[[n]]-baseline$pred$kary[[n]]))/2,numeric(1)))
        f_tv<-max(vapply(names(r$pred$flow),function(n){o<-env$ivt_observed_passage_summary(b$objects$fit_data[[n]])$observed_flow;sum(abs(r$pred$flow[[n]]-baseline$pred$flow[[n]]))*median(diff(sort(unique(o$ploidy))))/2},numeric(1)))
        check('distribution_repeat',k_tv<=1e-10&&f_tv<=1e-10)
        discrepancy<-max(abs(r$pred$vector-baseline$pred$vector));sf_emit(data.frame(state_scaled_gap=state_gap,kary_tv=k_tv,flow_tv=f_tv,mean_vector_discrepancy=discrepancy),file.path(out,'anchor_replay_checks.tsv'))
      }
    } else {
      points<-sf_tab(file.path(study,'verification/M01/bridge_points.tsv'));results<-list();branch_rows<-list()
      for(i in seq_len(nrow(points)))if(points$status[i]=='registered') {
        p<-b$p;p[[points$parameter[i]]]<-points$new_value[i]
        # P01/P02 are static/observer tasks; only R01/R02 previously propagated.
        stopifnot(nrow(ledger)+2L<74L)
        x<-sf_run(env,b,p,out,points$eval_id[i],ledger);ledger<-x$ledger
        results[[points$eval_id[i]]]<-x$result
        if(x$result$status=='returned') {
          slim<-x$result;slim$runs<-NULL;slim$pred$states<-NULL
          sf_rds(slim,file.path(out,paste0(points$eval_id[i],'.rds')))
          sf_emit(data.frame(output=names(slim$pred$vector),value=unname(slim$pred$vector)),file.path(out,paste0(points$eval_id[i],'_vector.tsv')))
          branch_rows[[length(branch_rows)+1L]]<-data.frame(eval_id=points$eval_id[i],parameter=points$parameter[i],offset=points$offset[i],objective=x$result$score$objective,seconds=x$result$seconds,as.data.frame(x$result$branch,stringsAsFactors=FALSE))
        }
      }
      sf_emit(if(length(branch_rows))do.call(rbind,branch_rows) else data.frame(eval_id=character()),file.path(out,'selection_branch_audit.tsv'))
      sf_rds(lapply(results,function(r){r$runs<-NULL;r$pred$states<-NULL;r}),file.path(out,'bridge_results.rds'))
    }
    sf_emit(data.frame(status='PASS',forward_bundles=if(task=='R03')nrow(ledger) else 1,optimizer_calls=0),file.path(out,'gate.tsv'))
  } else if(task=='O02') {
    b<-readRDS(file.path(study,'verification/M01/anchor_bundle.rds'));baseline<-readRDS(file.path(study,'verification/R01/result.rds'))
    widths<-unique(c(b$map$table$lower[b$map$table$parameter=='sigma_kary'],b$p$sigma_kary,2*b$p$sigma_kary,0.05,0.1,0.5,1,2,3))
    bounds<-b$map$table[b$map$table$parameter=='sigma_kary',];inputs<-list()
    for(cohort in c('2N','4N'))inputs[[cohort]]<-env$ivt_gaussian_initial_state(b$core$grid_pre,b$p[[paste0('init_mean_',cohort)]],b$p[[paste0('init_sd_',cohort)]])
    inputs$uniform<-rep(1/133,133);rows<-list();calls<-0L;grid<-b$core$grid_pre
    for(width in widths)if(width>=bounds$lower&&width<=bounds$upper) {
      K<-exp(-outer(grid,grid,'-')^2/(2*width^2));K<-sweep(K,2,colSums(K),'/')
      for(name in names(inputs)) {
        calls<-calls+1L;stopifnot(calls<=27L)
        sf_emit(data.frame(ordinal=seq_len(calls),status=c(rep('returned',calls-1),'started')),file.path(out,'observer_call_ledger.tsv'))
        input<-inputs[[name]]/sum(inputs[[name]]);q<-env$ivt_smooth_kary_distribution(grid,input,width);ref<-as.numeric(K%*%input)
        check(paste0('observer_reference_',calls),max(abs(q-ref))<=1e-12&&abs(sum(q)-1)<=1e-12)
        rows[[calls]]<-data.frame(width=width,input=name,identity_max_gap=max(abs(q-input)),reference_max_gap=max(abs(q-ref)),all_reference_offdiagonal_zero=all(K[row(K)!=col(K)]==0),normalization_gap=abs(sum(q)-1),floor_cells=sum(q<1e-12))
        sf_emit(data.frame(ordinal=seq_len(calls),status='returned'),file.path(out,'observer_call_ledger.tsv'))
      }
    }
    sf_emit(do.call(rbind,rows),file.path(out,'floating_observer_audit.tsv'))
    residual<-b$p$sigma_growth*c(0,1,2);variance<-list()
    for(factor in c(0.5,1,2)) {
      sigma<-b$p$sigma_growth*factor;native<-env$ivt_growth_loglik_df(data.frame(observed_growth=residual,predicted_growth_rate=0),sigma)$loglik
      reference<--log(sigma)-log(2*pi)/2-residual^2/(2*sigma^2)
      check(paste0('variance_score_',factor),max(abs(native-reference))<=1e-12)
      variance[[length(variance)+1L]]<-data.frame(SD=sigma,variance=sigma^2,max_score_gap=max(abs(native-reference)),mean_rank_role='outside mean map',control='literal residual fixture; no generated observations')
    }
    sf_emit(do.call(rbind,variance),file.path(out,'variance_role_audit.tsv'))
    observer<-sf_tab(file.path(study,'verification/M01/observer_points.tsv'));op<-list()
    for(i in seq_len(nrow(observer)))if(observer$status[i]=='registered') {
      p<-b$p;p$sigma_kary<-observer$new_value[i];pred<-sf_predictions(env,b,baseline$runs,p)
      op[[as.character(observer$offset[i])]]<-list(offset=observer$offset[i],params=p,pred=pred)
      op[[as.character(observer$offset[i])]]$pred$states<-NULL
      sf_emit(data.frame(output=names(pred$vector),value=unname(pred$vector)),file.path(out,paste0('sigma_',i,'_vector.tsv')))
    }
    sf_rds(op,file.path(out,'observer_perturbations.rds'))
    sf_emit(data.frame(status='PASS',observer_calls=calls,observer_perturbation_passes=length(op),forward_bundles=0,generated_observations=0),file.path(out,'gate.tsv'))
  } else stop('Unsupported finite Step5 task')
  sf_emit(data.frame(status='PASS',task=task,synthetic_observations=0,optimizer_calls=0),file.path(out,'application_result.tsv'))
  cat('STEP5_PHASE_PASS ',task,'\n',sep='')
}
main()
