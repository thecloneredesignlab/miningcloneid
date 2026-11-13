# Static scientific QA figure from previously validated saved predictions only.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study <- normalizePath(args[1]);z <- readRDS(file.path(study,"tests/repair1_rep1/replay.rds"))
ctx <- readRDS(file.path(study,"manifests/contract.rds"));b <- z$predictions$tables$burden
output <- file.path(study,"figures/production_replay_burden.png")
png(output,width=2200,height=1350,res=180,bg="white")
par(mfrow=c(2,4),mar=c(3.4,3.9,3,1),oma=c(3.2,0,4.6,0),mgp=c(2.2,.65,0),tcl=-.25,family="sans",col.axis="#334155",col.lab="#17212B")
ylim <- c(0,max(c(b$obs_burden,b$pred_burden_volume_mm3),na.rm=TRUE)*1.06)
for(s in ctx$scenarios) {
 d <- b[b$harvest==s$harvest,];keep <- d$day>0 & is.finite(d$obs_burden) & d$obs_burden>0
 short <- sub("SUM159-","",sub("_harvest","",s$harvest))
 plot(d$day,d$obs_burden,type="n",xlim=c(0,100),ylim=ylim,xlab="Day",ylab="Volume (mm3)",main=short,bty="l",cex.main=1.05)
 abline(h=pretty(ylim),col="#E5E7EB",lwd=.7)
 points(d$day[!keep],d$obs_burden[!keep],pch=1,col="#64748B",cex=1.05,lwd=1.2)
 points(d$day[keep],d$obs_burden[keep],pch=16,col="#2463A6",cex=.95)
 points(d$day,d$pred_burden_volume_mm3,pch=2,col="#C65D17",cex=1.05,lwd=1.2)
}
mtext("Production endpoint replay: volume predictions at recorded times",side=3,outer=TRUE,line=2.3,cex=1.3,font=2,col="#17212B")
mtext("Eight harvest scenarios; 92 scored positive burden values. Archived and replayed predictions agree within declared tolerances.",side=3,outer=TRUE,line=.8,cex=.9,col="#475569")
par(mfrow=c(1,1),oma=c(0,0,0,0),fig=c(0,1,0,.045),new=TRUE,mar=c(0,0,0,0));plot.new()
legend("center",legend=c("Scored observation","Native replay prediction","Excluded by production scoring"),pch=c(16,2,1),col=c("#2463A6","#C65D17","#64748B"),horiz=TRUE,bty="n",cex=1.05)
invisible(dev.off());cat(output,"\n")
