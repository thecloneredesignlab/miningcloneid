# Analytic aggregation fixtures only; no biological propagation or random data.
d0_unit_loss <- function(out) {
 g<-setNames(rep(0,114),paste0('g',1:114));k<-setNames(rep(list(rep(1/133,133)),12),paste0('k',1:12))
 f<-setNames(rep(list(rep(1/8,200)),20),paste0('f',1:20));data<-list(growth=g,kary=k,flow=f)
 c<-list(bundle=list(map=list(fixed=list(sigma_growth=1)),observations=setNames(rep(list(list(step=.04)),20),names(f))),dataset=data)
 z<-d0_loss(c,data);expected<-.5*log(2*pi)+log(133)+log(8)
 stopifnot(abs(z$objective-expected)<1e-10,z$max_growth_error==0,z$max_kary_tv==0,z$max_flow_tv==0)
 changed<-data;changed$growth[1]<-1;q<-d0_loss(c,data,changed)$objective
 stopifnot(abs(q-expected-1/(2*114))<1e-10)
 d0_emit(data.frame(fixture=c('uniform_expected_entropy','single_growth_displacement'),actual=c(z$objective,q-z$objective),expected=c(expected,1/(2*114)),pass=TRUE,model_calls=0L),file.path(out,'analytic_loss_fixtures.tsv'))
}
