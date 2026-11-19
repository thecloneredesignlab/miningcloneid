# Filesystem-only regression. No model source is loaded or propagated.
source('common.R')
root<-tempfile('d0_layout_');dir.create(root)
original<-file.path(root,'frozen');dir.create(file.path(original,'model'),recursive=TRUE);dir.create(file.path(original,'util'))
writeLines('fixture only; no model code',file.path(original,'util/o2_supply_demand_map_common_semantics.R'))
private<-file.path(root,'runtime/worker/model');d0_private_layout(original,private)
stopifnot(nzchar(Sys.readlink(file.path(dirname(private),'util'))),
 identical(readLines(file.path(dirname(private),'util/o2_supply_demand_map_common_semantics.R')),'fixture only; no model code'),
 normalizePath(file.path(dirname(private),'util'))==normalizePath(file.path(original,'util')))
unlink(root,recursive=TRUE)
cat('Private sibling-layout filesystem fixture PASS; zero model loads/calls\n')
