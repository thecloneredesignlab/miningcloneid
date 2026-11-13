# Static source-backed diagnostics; no model evaluations or web rendering.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study <- normalizePath(args[1]);entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"))
figures <- file.path(study,"figures");dir.create(figures,showWarnings=FALSE)
draw_both <- function(name,fn,width=10,height=4.8) {
 png(file.path(figures,paste0(name,".png")),width=width,height=height,units="in",res=180,type="cairo")
 fn();dev.off()
 pdf(file.path(figures,paste0(name,".pdf")),width=width,height=height,family="Helvetica",useDingbats=FALSE)
 fn();dev.off()
}
if(file.exists(file.path(study,"tables/baseline_chain_summary.tsv"))) {
 x <- read_tsv(file.path(study,"tables/baseline_chain_summary.tsv"))
 gates <- read_tsv(file.path(study,"tables/baseline_policy_gates.tsv"))
 draw_both("baseline_recovery",function() {
  layout(matrix(c(1,2,3,3),nrow=2,byrow=TRUE),heights=c(1,.13))
  par(mar=c(5.4,4.4,3,1),oma=c(1.4,0,0,0),las=1)
  for(policy in c("production","data_only")) {
   z <- x[x$policy==policy,,drop=FALSE];g <- gates[gates$policy==policy,,drop=FALSE]
   labels <- sub(paste0("^b_",policy,"_"),"",sub("_local$","",z$task))
   yrange <- range(c(z$objective,g$known_feasible_bound,g$best_objective+.002));padding <- max(.004,diff(yrange)*.12)
   plot(seq_len(nrow(z)),z$objective,type="n",xaxt="n",xlab="",ylab="Objective",ylim=yrange+c(-padding,padding),
    main=paste(if(policy=="production") "Production penalties" else "Data only",if(g$pass) "(accepted)" else "(unresolved)"))
   axis(1,seq_len(nrow(z)),labels,las=2,cex.axis=.85)
   abline(h=g$best_objective,col="#3b6ea8",lwd=1.4)
   abline(h=g$best_objective+.002,col="#3b6ea8",lty=3)
   abline(h=g$known_feasible_bound,col="#bb5a21",lty=2)
   points(seq_len(nrow(z)),z$objective,pch=ifelse(!z$independent,5,ifelse(z$family=="uniform",17,16)),cex=1.2,
    col=ifelse(z$gap<=.002,"#26714d","#8d3030"))
  }
  par(mar=c(0,0,0,0));plot.new()
  legend("center",c("Best verified chain","Archive hybrid","Recovery tolerance +0.002","Full-box uniform","Known feasible bound","Propagated restart"),
   col=c("#3b6ea8","#555555","#3b6ea8","#555555","#bb5a21","#555555"),
   lty=c(1,NA,3,NA,2,NA),pch=c(NA,16,NA,17,NA,5),bty="n",cex=.72,ncol=3)
  mtext("Finite optimization diagnostics; the tolerance is not a confidence threshold.",outer=TRUE,side=1,line=.1,cex=.8)
 })
}
if(file.exists(file.path(study,"tables/profile_point_review.tsv"))) {
 x <- read_tsv(file.path(study,"tables/profile_point_review.tsv"));parameter_rows <- list();prediction_rows <- list()
 contract <- readRDS(file.path(study,"manifests/contract.rds"))
 for(j in seq_len(nrow(x))) {
  z <- x[j,];a <- readRDS(file.path(study,"tables",paste0("profile_",z$policy,"_",z$density_index,"_verified.rds")))
  p <- a$record$parameters;nat <- p;islog <- startsWith(names(p),"log10_");nat[islog] <- 10^p[islog]
  parameter_rows[[j]] <- data.frame(policy=z$policy,density_index=z$density_index,parameter=names(p),transformed=unname(p),natural=unname(nat),
   normalized_position=unname((p-contract$lower)/(contract$upper-contract$lower)),point_checks_passed=z$pass)
  for(scenario_index in seq_along(a$record$predictions$native)) {
   native <- a$record$predictions$native[[scenario_index]]
   n <- native$result;arg <- native$arguments
   scenario <- contract$scenarios[[scenario_index]]
   prediction_rows[[length(prediction_rows)+1L]] <- data.frame(policy=z$policy,density_index=z$density_index,
    harvest=scenario$harvest,output_index=seq_along(n$Vmm3_total_obs),day=scenario$obs_days,
    volume=n$Vmm3_total_obs,live=n$Ntot_live_obs,dead_hypoxia=n$Ntot_dead_hypoxia_obs,
    dead_buffer=n$Ntot_dead_buffer_obs,oxygen=n$O2_eff_obs,oxygen_target=n$O2_target_obs,point_checks_passed=z$pass)
  }
 }
 write_tsv(do.call(rbind,parameter_rows),file.path(study,"tables/profile_parameters.tsv"))
 write_tsv(do.call(rbind,prediction_rows),file.path(study,"tables/profile_prediction_streams.tsv"))
 draw_both("density_profile_diagnostic",function() {
  par(mfrow=c(1,2),mar=c(4.7,4.4,3,1),oma=c(1.6,0,0,0),las=1)
  for(policy in c("production","data_only")) {
   z <- x[x$policy==policy,,drop=FALSE];z <- z[order(z$density),,drop=FALSE]
   if(!nrow(z)) {plot.new();title(policy);next}
   ylim <- range(c(z$raw_delta,0));padding <- max(.01,diff(ylim)*.12)
   plot(z$density,z$raw_delta,log="x",type="l",lwd=1.4,col="#3b6ea8",xaxt="n",xlim=c(50000,500000),
    ylim=ylim+c(-padding,padding),xlab="Density (cell equivalents / mm^3)",ylab="Objective minus matching free baseline",
    main=if(policy=="production") "Production penalties" else "Data only")
   axis(1,at=c(50000,100000,200000,500000),labels=c("50k","100k","200k","500k"))
   abline(h=0,col="#777777",lty=3)
   points(z$density,z$raw_delta,pch=ifelse(z$pass,16,1),cex=1.05,col=ifelse(z$pass,"#26714d","#8d3030"))
   baseline <- readRDS(file.path(study,"tables",paste0("baseline_",policy,".rds")))
   points(10^baseline$record$parameters[['log10_rho_2N']],0,pch=8,cex=1.2)
   legend("topleft",c("Passed anchor checks","Unresolved anchor","Free candidate reference"),pch=c(16,1,8),
    col=c("#26714d","#8d3030","black"),bty="n",cex=.75)
  }
  mtext("Conditional on fixed initial stocks and demand reference; no confidence interval calibration.",outer=TRUE,side=1,line=.2,cex=.78)
 })
}
assert_namespaces();cat("Available static diagnostics exported; zero model evaluations.\n")
