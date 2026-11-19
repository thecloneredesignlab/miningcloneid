# Fixed-parameter reduction validation inside the pinned production SIF.
args<-commandArgs(TRUE);stopifnot(length(args)==2L);s<-normalizePath(args[1]);task<-args[2]
source(file.path(s,'inputs/d0_common.R'));out<-file.path(s,'results',task);dir.create(out,recursive=TRUE,showWarnings=FALSE)
c<-d0_backend(s,task);source(file.path(c$paths$archive,'implementation/diagnostics.R'))
fit0<-readRDS(file.path(c$paths$wp0,'source/seed144/fit_result.rds'));cfg<-fit0$cfg
cfg$parameter_table<-file.path(c$paths$wp0,'source/seed144/parameter_table_input.csv')
stopifnot(isTRUE(cfg$Crowding),cfg$crowding=='logistic',cfg$K==1e12,!isTRUE(cfg$o2_burden_feedback),isTRUE(cfg$O2_growth),cfg$N_MIN==22,cfg$N_MAX==154)
core<-c$backend$build_model_core(cfg=cfg);grid<-cfg$N_MIN:cfg$N_MAX;nodes<-c(22,44,66,88,110);fine<-1:440;times<-0:7
# Explicit right-closed classes; all positive upper-tail destinations are routed.
P<-function(N,shift=0){edges<-c(32,54,76,98)+shift;idx<-vapply(N,function(n)1L+sum(n>edges),integer(1));outer(1:5,idx,'==')*1}
lift<-function(N,projection,weighted=FALSE){R<-matrix(0,length(N),5);for(j in 1:5){if(!weighted)R[match(nodes[j],N),j]<-1 else{R[,j]<-dnorm(N,nodes[j],2)*projection[j,];R[,j]<-R[,j]/sum(R[,j])}};stopifnot(max(abs(projection%*%R-diag(5)))<1e-12);R}
aug<-function(A,P,R){n<-nrow(A)-3L;L<-matrix(0,7,n+3);U<-matrix(0,n+3,7);L[1:5,1:n]<-P;U[1:n,1:5]<-R;L[6:7,(n+1):(n+2)]<-diag(2);U[(n+1):(n+2),6:7]<-diag(2);L%*%A%*%U}
derivative<-function(A,M,x,n){crowd<-max(0,1-sum(x[seq_len(n)])/cfg$K);as.numeric(A%*%x)-(1-crowd)*as.numeric(M%*%x)}
integrate<-function(A,M,x,n){z<-deSolve::ode(y=x,times=times,func=function(t,y,p)list(derivative(A,M,y,n)),parms=NULL,method='lsoda',rtol=1e-10,atol=1e-12);ans<-as.matrix(z[,-1,drop=FALSE]);stopifnot(min(ans)>-1e-9,all(is.finite(ans)));ans[ans<0]<-0;ans}
summary5<-function(state,P,n){live<-state[,1:n,drop=FALSE];list(total=rowSums(live),fraction=t(P%*%t(live))/rowSums(live),dead=rowSums(state[,(n+1):(n+2),drop=FALSE]))}
compare<-function(a,b){c(max_class_tv=max(rowSums(abs(a$fraction-b$fraction))/2),max_log_live_error=max(abs(log(a$total/b$total))),max_dead_absolute_gap=max(abs(a$dead-b$dead)))}
operator<-function(N,oxygen,p){
 n<-length(N);A<-matrix(0,n+3,n+3);M<-A;bio<-numeric(n);overflow<-numeric(n);maximum_ledger<-0
 h<-p$O2_crit/(p$O2_crit+oxygen);b<-p$lam_max/(1+p$alpha_o2*h*N/44);d<-p$mu_hp*h*N/44;prob<-p$p_mis_base+p$p_misseg*d/(p$k_o_mis+d)
 s2<-p$buffer_smax*exp(-p$buffer_beta);s4<-p$buffer_smax*exp(-p$buffer_beta*(.5)^p$buffer_n_exp)
 weight<-pmin(1,pmax(0,(N-44)/44));survival<-(1-weight)*s2+weight*s4
 stopifnot(all(prob>=0&prob<=1),all(survival>=0&survival<=1))
 for(j in seq_along(N)){
  m<-0:N[j];pn<-dbinom(m,N[j],prob[j]);v<-pn*survival[j]^m;offspring<-numeric(n);valid<-0
  for(sign in c(-1,1)){
   destination<-N[j]+sign*m;keep<-destination>0;valid<-valid+sum(v[keep]);capped<-pmin(max(N),destination[keep]);idx<-match(capped,N)
   added<-rowsum(matrix(v[keep],ncol=1),idx,reorder=FALSE);where<-as.integer(rownames(added));offspring[where]<-offspring[where]+added[,1];overflow[j]<-overflow[j]+sum(v[keep&destination>max(N)])
  }
  bio[j]<-2-valid;maximum_ledger<-max(maximum_ledger,abs(sum(offspring)+bio[j]-2))
  A[1:n,j]<-b[j]*(1-p$p_wgd)*offspring;wgd<-min(max(N),2*N[j]);A[match(wgd,N),j]<-A[match(wgd,N),j]+b[j]*p$p_wgd
  A[j,j]<-A[j,j]-b[j]-d[j];A[n+1,j]<-d[j];A[n+2,j]<-b[j]*(1-p$p_wgd)*bio[j]
  A[n+3,j]<-b[j]*((1-p$p_wgd)*overflow[j]+p$p_wgd*as.numeric(2*N[j]>max(N)))
  M[,j]<-A[,j];M[j,j]<-M[j,j]+d[j];M[n+1,j]<-0
 }
 A[n+1,n+1]<-A[n+2,n+2]<--1e-4
 stopifnot(maximum_ledger<1e-12)
 list(A=A,M=M,ledger=maximum_ledger,s2=s2,s4=s4)
}
cases<-d0_tab(file.path(s,'inputs/production_cases.tsv'));rows<-list();curves<-list();ledgers<-list();native_calls<-0L;reused<-0L;kernel_roundoff<-list()
reuse<-d0_tab(file.path(Sys.getenv('SI4_RELEASE'),'native_reuse_manifest.tsv'))
for(i in seq_len(nrow(cases))){
 row<-cases[i,];path<-file.path(c$paths$project,row$path);stopifnot(d0_hash(path)==row$sha256);fit<-readRDS(path)
 map<-o2prof_map(cfg$parameter_table,fit$best_params,'zero_extended');projected<-o2prof_project_start(map,fit$best_params);stopifnot(!any(projected$projection$reason=='projected_to_bound'));p<-projected$params
 d0_emit(projected$projection,file.path(out,paste0(row$case,'_projection.tsv')))
 volumes<-c$backend$cell_volume_mm3_by_N(grid,run_params=p,cfg=cfg);stopifnot(max(abs(volumes-1/p$rho_2N))<1e-15)
 for(oxygen in c(.2,1,20.5)){
  rates<-pilotdiag_rates(c,p,cfg,oxygen);prod<-as.matrix(rates$generator)-diag(rates$rates$mu);n<-length(grid);Apr<-matrix(0,n+2,n+2);Apr[1:n,1:n]<-prod;Apr[n+1,1:n]<-rates$rates$mu
  decomposed_buffer<-rates$rates$ordinary_biological_nonviability+rates$rates$numerical_tail_pruning+rates$rates$ordinary_boundary_loss+rates$rates$ordinary_boundary_gain+rates$rates$wgd_boundary
  tri<-get('cpp_o2simps_build_G_for_o2_triplet',envir=c$backend,inherits=TRUE)(O2=oxygen,O2_crit=p$O2_crit,N0min=cfg$N_MIN,N0max=cfg$N_MAX,N1min=cfg$N_MIN,N1max=cfg$N_MAX,lam_max=p$lam_max,p_mis_base=p$p_mis_base,p_misseg=p$p_misseg,k_o_mis=p$k_o_mis,p_wgd=p$p_wgd,boundary='drop',eps_tail=1e-8,buffer_smax=p$buffer_smax,buffer_beta=p$buffer_beta,buffer_n_exp=p$buffer_n_exp,N_unit=cfg$N_UNIT,beta_size=p$beta_size,O2_growth=cfg$O2_growth,alpha_o2=p$alpha_o2,gamma_growth=p$gamma_growth,mu_hp=p$mu_hp,gamma_mu=p$gamma_mu,n_O=p$n_O,ploidy_O2_death=cfg$ploidy_O2_death)
  Apr[n+2,1:n]<-tri$dead_buffer_rate
  kernel_roundoff[[length(kernel_roundoff)+1L]]<-data.frame(case=row$case,oxygen_pct=oxygen,maximum_native_minus_decomposed_dead_rate=max(abs(tri$dead_buffer_rate-decomposed_buffer)),uses_native_rate_for_propagation=TRUE)
  d0_emit(do.call(rbind,kernel_roundoff),file.path(out,'native_vs_decomposed_dead_rate.tsv'))
  Apr[n+1,n+1]<-Apr[n+2,n+2]<--p$k_clear
  Mpr<-matrix(0,n+2,n+2);Mpr[1:n,1:n]<-as.matrix(rates$generator);Mpr[n+2,1:n]<-Apr[n+2,1:n]
  surrogate<-operator(fine,oxygen,p);PP<-P(fine);RR<-lift(fine,PP);Ac<-aug(surrogate$A,PP,RR);Mc<-aug(surrogate$M,PP,RR)
  # Independent five-parent operator constructed on a destination-complete grid.
  expected<-matrix(0,7,7)
  for(j in 1:5){N<-nodes[j];m<-0:N;h<-p$O2_crit/(p$O2_crit+oxygen);b<-p$lam_max/(1+p$alpha_o2*h*N/44);d<-p$mu_hp*h*N/44;q<-p$p_mis_base+p$p_misseg*d/(p$k_o_mis+d);u<-pmin(1,pmax(0,(N-44)/44));ss<-(1-u)*surrogate$s2+u*surrogate$s4;v<-dbinom(m,N,q)*ss^m;B<-numeric(5)
   for(sign in c(-1,1)){dst<-N+sign*m;ok<-dst>0;B<-B+as.numeric(P(dst[ok])%*%v[ok])}
   expected[1:5,j]<-b*(1-p$p_wgd)*B;dest<-which(P(2*N)==1);expected[dest,j]<-expected[dest,j]+b*p$p_wgd;expected[j,j]<-expected[j,j]-b-d;expected[6,j]<-d;expected[7,j]<-b*(1-p$p_wgd)*(2-sum(B))
  };expected[6,6]<-expected[7,7]<--1e-4;projection_gap<-max(abs(expected-Ac));stopifnot(projection_gap<1e-12)
  for(cohort in c('2N','4N')){
   f0<-c$backend$ivt_gaussian_initial_state(grid,p[[paste0('init_mean_',cohort)]],p[[paste0('init_sd_',cohort)]]);initial_cells<-fit$best_components[[paste0('run_',cohort)]]$segment_results[[1]]$sim$Ntot_live_obs[1];stopifnot(is.finite(initial_cells),initial_cells>0,initial_cells<cfg$K)
   segment<-list(cohort=cohort,oxygen_pct=oxygen,initial_cells=initial_cells,obs_days_local=times,duration_days=7)
   label<-paste(row$case,cohort,oxygen,sep='_');cached<-reuse[reuse$label==label,,drop=FALSE]
   if(nrow(cached)){
    stopifnot(nrow(cached)==1L,d0_hash(file.path(s,cached$path))==cached$sha256);native<-readRDS(file.path(s,cached$path));reused<-reused+1L
   }else{
    slot<-system2('python3',c('-B',shQuote(file.path(Sys.getenv('SI4_RELEASE'),'workflow.py')),shQuote(s),'native','--task',shQuote(label)),stdout=TRUE);stopifnot(is.null(attr(slot,'status')),length(slot)==1L,as.integer(slot)<=48L)
    native_calls<-native_calls+1L;native<-c$backend$ivt_run_segment_fixed_o2(segment,cfg,p,core,volumes,init_state_override=f0,init_cells_override=initial_cells)$sim
   }
   saveRDS(native,file.path(out,paste0(label,'_native.rds')))
   xx<-c(f0*initial_cells,0,0);euler<-matrix(0,8,n+2);euler[1,]<-xx;step<-0L
   for(k in 1:round(7/cfg$DT)){xx<-pmax(0,xx+cfg$DT*derivative(Apr,Mpr,xx,n));if(abs(k*cfg$DT-round(k*cfg$DT))<1e-10){step<-step+1L;euler[step+1L,]<-xx}}
   native_state<-cbind(native$live_state_obs,native$Ntot_dead_hypoxia_obs,native$Ntot_dead_buffer_obs)
   native_gap<-max(abs(euler-native_state)/pmax(1,abs(native_state)));d0_write(list(case=label,maximum_scaled_gap=native_gap,threshold=1e-12,reused=nrow(cached)==1L),file.path(out,paste0(label,'_native_identity.tsv')));stopifnot(native_gap<1e-12)
   continuous<-integrate(Apr,Mpr,c(f0*initial_cells,0,0),n);fine0<-numeric(length(fine));fine0[match(grid,fine)]<-f0*initial_cells;fine_run<-integrate(surrogate$A,surrogate$M,c(fine0,0,0,0),length(fine));coarse<-integrate(Ac,Mc,c(as.numeric(PP%*%fine0),0,0),5)
   a<-summary5(euler,P(grid),length(grid));b<-summary5(continuous,P(grid),length(grid));f<-summary5(fine_run,PP,length(fine));co<-summary5(coarse,diag(5),5)
   for(pair in c('production_euler_vs_continuous','production_continuous_vs_rational_fine','rational_fine_vs_five_class','production_euler_vs_five_class')){
    comparison<-switch(pair,production_euler_vs_continuous=compare(a,b),production_continuous_vs_rational_fine=compare(b,f),rational_fine_vs_five_class=compare(f,co),production_euler_vs_five_class=compare(a,co))
    rows[[length(rows)+1L]]<-data.frame(case=row$case,cohort=cohort,oxygen_pct=oxygen,comparison=pair,t(comparison),engineering_resolution_pass=comparison[1]<=.05&&comparison[2]<=.1)
   }
   for(alternative in c('weighted_lifting','shifted_edges')){pa<-P(fine,if(alternative=='shifted_edges')2 else 0);ra<-lift(fine,pa,alternative=='weighted_lifting');aa<-aug(surrogate$A,pa,ra);mm<-aug(surrogate$M,pa,ra);alt<-integrate(aa,mm,c(as.numeric(pa%*%fine0),0,0),5);comparison<-compare(summary5(fine_run,pa,length(fine)),summary5(alt,diag(5),5));rows[[length(rows)+1L]]<-data.frame(case=row$case,cohort=cohort,oxygen_pct=oxygen,comparison=alternative,t(comparison),engineering_resolution_pass=comparison[1]<=.05&&comparison[2]<=.1)}
   ledgers[[length(ledgers)+1L]]<-data.frame(case=row$case,cohort=cohort,oxygen_pct=oxygen,initial_cells=initial_cells,carrying_capacity=cfg$K,maximum_native_crowding_fraction=max(rowSums(native$live_state_obs))/cfg$K,native_euler_scaled_gap=native_gap,instantaneous_PAR_gap=projection_gap,binomial_ledger_gap=surrogate$ledger,upper_tail_routed_equivalents=fine_run[8,length(fine)+3],upper_tail_route_per_final_live=fine_run[8,length(fine)+3]/sum(fine_run[8,seq_along(fine)]))
   for(name in c('a','b','f','co')){z<-get(name);curves[[length(curves)+1L]]<-data.frame(case=row$case,cohort=cohort,oxygen_pct=oxygen,model=name,day=times,live=z$total,retained_dead=z$dead,z$fraction)}
   d0_emit(do.call(rbind,rows),file.path(out,'reduction_validation.tsv'));d0_emit(do.call(rbind,ledgers),file.path(out,'accounting_and_native_checks.tsv'))
   cat('REDUCTION_CASE',label,'native_gap',native_gap,'\n');flush.console()
  }
 }
}
d0_emit(do.call(rbind,curves),file.path(out,'validation_trajectories.tsv'));d0_emit(do.call(rbind,kernel_roundoff),file.path(out,'native_vs_decomposed_dead_rate.tsv'));tab<-do.call(rbind,rows);ledger<-do.call(rbind,ledgers)
d0_write(list(pass=TRUE,source_and_kernel_integrity=TRUE,production_to_reduction_approximation_pass=all(tab$engineering_resolution_pass[tab$comparison=='production_euler_vs_five_class']),coarse_projection_approximation_pass=all(tab$engineering_resolution_pass[tab$comparison=='rational_fine_vs_five_class']),fine_tail_routing_pass=all(ledger$upper_tail_route_per_final_live<=1e-6),native_segment_calls=native_calls,reused_native_segments=reused,cumulative_native_segment_calls=nrow(d0_tab(file.path(s,'state/native_calls.tsv'))),model_optimizations=0L,synthetic_datasets=0L,step05_release=FALSE,scope='Conditional rational reduced model only; approximation failures remain transfer holds'),file.path(out,'result.tsv'))
d0_guard(s,out);d0_seal(out)
