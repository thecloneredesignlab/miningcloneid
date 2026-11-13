# Parse all formal R code and exercise metadata-only checks without model loading.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study <- normalizePath(args[1]);entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
for(path in list.files(entry,pattern="[.]R$",full.names=TRUE)) parse(path)
source(file.path(entry,"adapter.R"));source(file.path(entry,"helpers.R"))
task <- "formal_metadata_local"
source(file.path(entry,"test_formal.R"))
cat("All formal R files parsed; analytic and exact-transport metadata checks passed; zero tumor evaluations.\n")
