# Prospective resolution01 profile arithmetic. No simulator/model calls.
resolution_profile_case<-function(x,actual_floor,free_reference){
 empty<-list(branches=0L,objective=actual_floor,raw_delta=actual_floor-free_reference,
   branch_gap=Inf,directional_gain=Inf,representative_available=FALSE,pass=FALSE,fine_shape_unresolved=TRUE)
 if(is.null(x)||!nrow(x))return(empty)
 required<-c('branch','objective','representative_objective','probe_gain','complete','has_local_history')
 stopifnot(all(required%in%names(x)),is.finite(actual_floor),is.finite(free_reference))
 # A branch supplies its own checked representative, never another branch's
 # center or the smallest probe improvement across unrelated histories.
 valid<-x$complete&x$has_local_history&is.finite(x$probe_gain)&x$probe_gain<=.002&
    x$representative_objective-x$objective<=1e-4
 z<-x[valid,,drop=FALSE]
 if(!nrow(z))return(empty)
 z<-do.call(rbind,lapply(split(z,z$branch),function(a)a[which.min(a$objective),,drop=FALSE]))
 z<-z[order(z$objective),,drop=FALSE]
 # Preserve an unchecked lower canonical point as a deficit, not a discarded
 # outlier that two worse agreeing histories can conceal.
 near<-z$objective<=actual_floor+.01
 pair<-z[near,,drop=FALSE]
 gap<-if(nrow(pair)>=2L)pair$objective[2L]-pair$objective[1L] else Inf
 rise<-actual_floor-free_reference
 representative<-any(z$representative_objective<=actual_floor+1e-4)
 list(branches=nrow(z),objective=actual_floor,raw_delta=rise,branch_gap=gap,
   directional_gain=if(nrow(pair))max(pair$probe_gain[seq_len(min(2L,nrow(pair)))]) else Inf,
   representative_available=representative,
   pass=nrow(pair)>=2L&&representative&&gap<=.01&&rise>=-.002,
   fine_shape_unresolved=rise<.01||(!is.infinite(gap)&&rise<=gap))
}
