# Static, source-backed diagnostics from closed records only; zero model calls.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));figures<-file.path(study,'plots');dir.create(figures,showWarnings=FALSE)
draw<-function(name,fn,width=12,height=5.5){png(file.path(figures,paste0(name,'.png')),width=width,height=height,units='in',res=180,type='cairo');tryCatch(fn(),finally=dev.off())}
if(file.exists(file.path(study,'tables/baseline_chain_summary.tsv'))){
 x<-read_tsv(file.path(study,'tables/baseline_chain_summary.tsv'));gates<-read_tsv(file.path(study,'tables/baseline_policy_gates.tsv'));labels<-list()
 root_counts<-vapply(split(x$root_id,x$policy),function(v)length(unique(v)),integer(1))
 draw('baseline_recovery',function(){
  par(mfrow=c(1,2),mar=c(5.5,4.5,3,1),oma=c(2,0,0,0),las=1)
  for(policy in c('production','data_only')){
   z<-x[x$policy==policy,,drop=FALSE];g<-gates[gates$policy==policy,,drop=FALSE]
   if(!nrow(z)){plot.new();title(policy);next}
   roots<-unique(z$root_id);z<-do.call(rbind,lapply(roots,function(root){a<-z[z$root_id==root,,drop=FALSE];a[which.min(a$verified_best),,drop=FALSE]}))
   key<-ifelse(z$family=='archive','Archive',paste('Uniform',cumsum(z$family=='uniform')));if(any(z$family=='informed'))key[z$family=='informed']<-paste('Informed',seq_len(sum(z$family=='informed')))
   labels[[policy]]<<-data.frame(policy=policy,label=key,root_id=z$root_id,task=z$task,objective=z$verified_best,gap=z$verified_best-g$best_objective)
   yrange<-range(c(z$verified_best,g$known_feasible_bound,g$best_objective+.002));padding<-max(.004,diff(yrange)*.12)
   plot(seq_len(nrow(z)),z$verified_best,type='n',xaxt='n',xlab='',ylab='Verified objective',ylim=yrange+c(-padding,padding),main=paste(if(policy=='production')'Production penalties' else 'Data only',if(g$pass)'(accepted)' else '(unresolved)'))
   axis(1,seq_len(nrow(z)),key,las=2,cex.axis=.8);abline(h=g$best_objective,col='#326aa3',lwd=1.4);abline(h=g$best_objective+.002,col='#326aa3',lty=3);abline(h=g$known_feasible_bound,col='#b85b25',lty=2)
   points(seq_len(nrow(z)),z$verified_best,pch=ifelse(z$family=='uniform',17,ifelse(z$family=='archive',16,18)),cex=1.2,col=ifelse(z$verified_best-g$best_objective<=.002,'#26714d','#9d392e'))
   legend('topleft',c('Lowest verified feasible','Recovery +0.002','All-known feasible floor'),col=c('#326aa3','#326aa3','#b85b25'),lty=c(1,3,2),bty='n',cex=.7)
  }
  mtext('One point per original root; continuations and forks count once. Finite recovery diagnostics, no calibrated confidence intervals.',outer=TRUE,side=1,line=.3,cex=.78)
 },width=max(12,2+.65*max(root_counts)))
 write_tsv(do.call(rbind,labels),file.path(study,'tables/baseline_plot_roots.tsv'))
}
if(file.exists(file.path(study,'tables/fit_summary.tsv'))){
 x<-read_tsv(file.path(study,'tables/fit_summary.tsv'));x<-x[x$calls>0,,drop=FALSE]
 if(nrow(x))draw('window_timing',function(){layout(matrix(c(1,2),ncol=1),heights=c(1,.16));par(mar=c(4.8,5,2,1),las=1);plot(x$calls,x$parent_seconds/60,log='xy',pch=ifelse(x$mode=='global',17,16),col=ifelse(x$policy=='production','#326aa3','#b85b25'),xlab='Optimizer callbacks in the closed window',ylab='Slurm parent minutes',main='Measured window timing (all successful fit windows)');par(mar=c(0,0,0,0));plot.new();legend('center',c('Production','Data only','Serial local','8-core global'),col=c('#326aa3','#b85b25','black','black'),pch=c(16,16,16,17),bty='n',cex=.85,ncol=4)})
}
if(file.exists(file.path(study,'tables/profile_branch_history.tsv'))){
 x<-read_tsv(file.path(study,'tables/profile_branch_history.tsv'));ep<-read_kv(file.path(study,'tables/profile_epoch.tsv'));x<-x[x$epoch==ep$epoch,,drop=FALSE];gates<-read_tsv(file.path(study,'tables/profile_anchor_gates.tsv'))
 draw('density_profile',function(){par(mfrow=c(1,2),mar=c(4.8,4.8,3,1),oma=c(1.7,0,0,0),las=1)
  for(policy in c('production','data_only')){
   z<-x[x$policy==policy,,drop=FALSE];if(!nrow(z)){plot.new();title(policy);next}
   yrange<-range(c(z$raw_delta,0));padding<-max(.01,diff(yrange)*.12);plot(c(50000,500000),yrange,type='n',log='x',xaxt='n',xlim=c(50000,500000),ylim=yrange+c(-padding,padding),xlab='Density (cell equivalents / mm^3)',ylab='Objective minus frozen free baseline',main=if(policy=='production')'Production penalties' else 'Data only')
   axis(1,at=c(50000,100000,200000,500000),labels=c('50k','100k','200k','500k'));abline(h=0,col='#777777',lty=3)
   branches<-c('global_local','ascending','descending');colors<-c('#326aa3','#26714d','#b85b25')
   for(i in seq_along(branches)){a<-z[z$branch==branches[i],,drop=FALSE];if(!nrow(a))next;a<-do.call(rbind,lapply(split(a,a$index),function(b)b[which.min(b$objective),,drop=FALSE]));a<-a[order(a$index),,drop=FALSE];lines(a$density,a$raw_delta,col=colors[i],lty=i);points(a$density,a$raw_delta,col=colors[i],pch=ifelse(gates$pass[match(paste(a$policy,a$index),paste(gates$policy,gates$index))],16,1))}
   legend('topleft',c('Global + local','Ascending','Descending'),col=colors,lty=1:3,bty='n',cex=.75)
  };mtext('Raw deltas are not clipped; open symbols are unresolved anchors. No confidence interval calibration.',outer=TRUE,side=1,line=.2,cex=.78)
 })
}
assert_namespaces();cat('Available closed-record PNG diagnostics exported; zero model calls.\n')
