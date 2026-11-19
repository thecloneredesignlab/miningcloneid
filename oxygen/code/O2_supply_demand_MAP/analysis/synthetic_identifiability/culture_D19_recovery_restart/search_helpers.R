# Reused finite direction and polling rules; no randomized screening.
r_activity<-function(p){
 inactive<-character()
 if(p$p_misseg==0)inactive<-c(inactive,'k_o_mis')
 if(p$buffer_beta==0)inactive<-c(inactive,'buffer_n_exp')
 if(p$buffer_smax==0)inactive<-c(inactive,'buffer_beta','buffer_n_exp')
 unique(inactive)
}
r_directions<-function(c){
 n<-nrow(c$free);D<-diag(n);rownames(D)<-c$free$parameter
 pairs<-list(c('lam_max','alpha_o2',1),c('mu_hp','k_o_mis',1),c('p_misseg','p_mis_base',-1),c('buffer_smax','buffer_beta',1),c('init_sd_2N','sigma_kary',-1),c('init_sd_4N','sigma_kary',-1))
 for(pair in pairs){v<-numeric(n);v[match(pair[1],c$free$parameter)]<-1/sqrt(2);v[match(pair[2],c$free$parameter)]<-as.numeric(pair[3])/sqrt(2);D<-rbind(D,v)}
 colnames(D)<-c$free$parameter;D
}
r_admissible<-function(x,d,h){
 dist<-c((1-x[d>0])/d[d>0],-x[d<0]/d[d<0]);step<-min(c(h,dist))
 if(!is.finite(step)||step<=0)return(NULL)
 y<-x+step*d;tol<-64*.Machine$double.eps
 if(any(y < -tol|y>1+tol))stop('Material boundary error')
 y<-pmin(1,pmax(0,y));if(identical(as.numeric(y),as.numeric(x)))NULL else as.numeric(y)
}
r_poll<-function(fn,x,directions,meshes,tolerance,max_calls){
 calls<-0L;q<-fn(x);calls<-calls+1L;rows<-list();exhausted<-FALSE
 for(h in meshes){repeat{
  oldq<-q;winner<-x;winq<-q;seen<-character();complete<-TRUE;maximum<-0
  for(j in seq_len(nrow(directions)))for(sign in c(-1,1)){
   y<-r_admissible(x,sign*directions[j,],h);if(is.null(y))next;key<-sf_hash(y);if(key%in%seen)next;seen<-c(seen,key)
   if(calls>=max_calls){complete<-FALSE;exhausted<-TRUE;break}
   z<-tryCatch(fn(y),recovery_budget_stop=function(e){exhausted<<-TRUE;NA_real_});if(is.na(z)){complete<-FALSE;break};calls<-calls+1L
   maximum<-max(maximum,q-z);if(z<winq){winner<-y;winq<-z}
  }
  x<-winner;q<-winq;rows[[length(rows)+1L]]<-data.frame(mesh=h,complete=complete,before=oldq,after=q,maximum_improvement=maximum,calls=calls)
  if(exhausted||oldq-q<=tolerance)break
 }
 if(exhausted)break}
 list(unit=x,score=q,rows=do.call(rbind,rows),budget_limited=exhausted,calls=calls)
}
r_result<-function(out,x){sf_rds(x,file.path(out,'state.rds'));sf_emit(data.frame(key=names(x)[vapply(x,function(v)is.atomic(v)&&length(v)==1,logical(1))],value=vapply(x[vapply(x,function(v)is.atomic(v)&&length(v)==1,logical(1))],as.character,character(1))),file.path(out,'result.tsv'))}
