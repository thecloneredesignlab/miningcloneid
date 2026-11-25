args <- commandArgs(trailingOnly=TRUE)
out <- args[1]
rows <- read.delim(file.path(out,"profile_grid.tsv"),stringsAsFactors=FALSE)
summary <- read.delim(file.path(out,"mechanism_summary.tsv"),stringsAsFactors=FALSE)
stopifnot(nrow(summary)==14)
png(file.path(out,"joint_profile_objectives.png"),width=1840,height=2820,res=145,bg="white",type="cairo")
par(mfrow=c(5,3),mar=c(3.65,3.55,2.3,.65),oma=c(4.8,1.1,3.9,.5),las=1,bty="n",
    family="sans",fg="#20394f",col.axis="#40586c",col.lab="#344f67")
verified <- "#28648c"
provisional <- "#a85b2e"
zero <- "#195b4f"
for(i in seq_len(nrow(summary))) {
  p <- summary$parameter[i]
  a <- subset(rows,parameter==p)
  ok <- a$status=="verified"
  open <- a$status=="provisional"
  envelope <- !is.na(a$feasibility_evidence) & a$feasibility_evidence!=""
  no <- a$status=="feasibility_unresolved" & !envelope
  fixed_zero <- a$exact_zero=="True" & ok & p %in% c("p_misseg","buffer_beta")
  # Show score rises on a log1p scale to retain detail near the estimate.
  y <- log1p(pmax(0,a$rise))
  top <- max(c(1,y),na.rm=TRUE)
  plot(NA,xlim=c(-.025,1.025),ylim=c(-.08,top*1.06),xlab="",ylab="",main=p,
       cex.main=1.1,xaxt="n",yaxt="n")
  ticks <- seq(0,floor(top)+1L)
  visible <- ticks <= top*1.06
  ticks <- ticks[visible]
  axis(1,at=c(0,.25,.5,.75,1),labels=c("0",".25",".5",".75","1"),cex.axis=.83,col="#92a5b6")
  axis(2,at=ticks,labels=formatC(expm1(ticks),digits=2,format="fg"),cex.axis=.76,col="#92a5b6")
  abline(h=ticks,col="#e6edf4")
  abline(v=summary$estimate_normalized[i],col="#7890a3",lty=3)
  points(a$normalized[ok & !fixed_zero],y[ok & !fixed_zero],pch=19,cex=1.02,col=verified)
  points(a$normalized[open],y[open],pch=1,cex=1.1,lwd=1.35,col=provisional)
  points(a$normalized[fixed_zero],y[fixed_zero],pch=18,cex=1.3,col=zero)
  points(a$normalized[no],rep(0,sum(no)),pch=4,cex=.96,col="#6d7d87")
  points(a$normalized[envelope],rep(0,sum(envelope)),pch=15,cex=.80,col="#866b68")
  if(i %% 3L == 1L) mtext("Objective rise",2,line=2.7,las=0,cex=.80)
}
plot.new()
legend("center",legend=c("Verified objective rise","Provisional candidate","Verified exact zero",
       "Feasibility unresolved","Ruled out by donor supply","Unrestricted J estimate"),pch=c(19,1,18,4,15,NA),
       lty=c(NA,NA,NA,NA,NA,3),lwd=1.2,cex=.90,
       col=c(verified,provisional,zero,"#6d7d87","#866b68","#7890a3"),bty="n")
mtext("Mechanism profiles from corrected counts, death and chromosomes",outer=TRUE,side=3,line=2.1,cex=1.38,font=2,col="#172f46")
mtext("Strict donor inheritance | Observed data | Each panel has its own vertical range",outer=TRUE,side=3,line=.7,cex=1.0,col="#456176")
mtext("Normalized fitting coordinate",outer=TRUE,side=1,line=2.6,cex=1.03,col="#223e55")
mtext("Vertical position uses log1p of the rise. Crosses and squares mark fixed values without a fitted score.",outer=TRUE,
      side=1,line=1.15,cex=.86,col="#486277")
invisible(dev.off())
