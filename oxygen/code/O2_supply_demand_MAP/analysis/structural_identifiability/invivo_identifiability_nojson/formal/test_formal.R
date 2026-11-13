# Analytic direction/constraint oracles and exact scalar transport checks.
unit <- file.path(study,"tests",paste0(task,"_unit"));dir.create(unit,recursive=TRUE,showWarnings=FALSE)
p <- c(x=.5,y=.5);lo <- c(x=0,y=0);hi <- c(x=1,y=1)
probe <- make_direction_probes(p,lo,hi,NULL,c(.25,.1))
stopifnot(length(probe)==8L)
quad <- function(z) (z[['x']]-.5)^2+(z[['y']]-.25)^2
stopifnot(min(vapply(probe,function(z) quad(z$parameters),numeric(1)))==0)
fixed <- make_direction_probes(p,lo,hi,c(x=.5),c(.25,.1))
stopifnot(length(fixed)==4L,all(vapply(fixed,function(z) z$parameters[['x']]==.5,logical(1))))
boundary <- make_direction_probes(c(x=0,y=1),lo,hi,c(y=1),c(.25,.1))
stopifnot(length(boundary)==2L,all(vapply(boundary,function(z) z$parameters[['x']]>0 && z$parameters[['y']]==1,logical(1))))
for(t in 0:4) {
 z <- optim(1,function(u) (t+u-3)^2+(u-1)^2,method="L-BFGS-B",lower=-5,upper=5)
 stopifnot(abs(z$par-(2-t/2))<1e-6,abs(z$value-(t-2)^2/2)<1e-8)
}
contract <- readRDS(file.path(study,"manifests/contract.rds"))
set.seed(52717);values <- c(contract$endpoint,setNames(rnorm(1000),paste0("random",1:1000)),
 edge=as.numeric("0x1.ef6ab1ac5bf8fp+3"),small=.Machine$double.xmin,large=.Machine$double.xmax,
 positive_zero=0,negative_zero=-0,infinite=Inf,missing=NA_real_,not_a_number=NaN)
write_kv(as.list(values),file.path(unit,"scalars.tsv"));stopifnot(identical(values,unlist(read_kv(file.path(unit,"scalars.tsv")))))
misc <- list(literal_na="NA",empty="",missing=NA_character_,true=TRUE,false=FALSE,integer=7L,quoted="tab\tquote\" and newline\n")
write_kv(misc,file.path(unit,"strings.tsv"));stopifnot(identical(misc,read_kv(file.path(unit,"strings.tsv"))))
atomic_rds(contract$endpoint,file.path(unit,"parameters.rds"));stopifnot(identical(contract$endpoint,readRDS(file.path(unit,"parameters.rds"))))
write_kv(list(pass=TRUE,analytic_box_probe=TRUE,exact_typed_transport=TRUE,nuisance_refit_oracle=TRUE,
 roundtrip_values=length(values),model_evaluations=0L),file.path(unit,"result.tsv"));finalize_payload(unit)
