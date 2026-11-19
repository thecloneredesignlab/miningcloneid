# Model-free regression: TSV values are strings, including numeric controls.
source(commandArgs(TRUE)[1],local=.GlobalEnv)
p<-tempfile();writeLines(c('key\tvalue','local_max_callbacks\t11183','truth_selection_by_distance\tfalse'),p)
controls<-d0_kv(p);stopifnot(is.character(controls$local_max_callbacks))
opts<-d0_local_options(controls,6900)
stopifnot(opts$maxeval==7828L,opts$maxtime==4830,opts$algorithm=='NLOPT_LN_SBPLX')
stopifnot(inherits(try(d0_local_options(list(local_max_callbacks='not_a_number'),6900),silent=TRUE),'try-error'))
unlink(p);cat('Local TSV numeric-control regression PASS; zero model calls\n')
