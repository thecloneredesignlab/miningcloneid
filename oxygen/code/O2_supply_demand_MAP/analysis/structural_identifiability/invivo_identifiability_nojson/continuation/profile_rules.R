# Arithmetic over verified records only, with no simulator access.
summarize_profile_anchor<-function(x){
 if(is.null(x)||!nrow(x))return(list(branches=0L,objective=NA_real_,raw_delta=NA_real_,branch_gap=NA_real_,directional_gain=NA_real_,pass=FALSE))
 # Multiple windows from one history are one branch, regardless of task count.
 x<-do.call(rbind,lapply(split(x,x$branch),function(z)z[which.min(z$objective),,drop=FALSE]));x<-x[order(x$objective),,drop=FALSE]
 gap<-if(nrow(x)>=2L)x$objective[2L]-x$objective[1L] else Inf
 gain<-min(x$center_objective)-min(x$objective)
 list(branches=nrow(x),objective=x$objective[1L],raw_delta=x$raw_delta[1L],branch_gap=gap,directional_gain=gain,pass=gap<=.01&&gain<=.002&&min(x$raw_delta)>=-.002)
}
