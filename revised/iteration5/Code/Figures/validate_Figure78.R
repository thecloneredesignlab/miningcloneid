#!/usr/bin/env Rscript
# Validate current retained/generated full-data receipts without numerical work.
# This metadata CLI may write an audit receipt; static checking uses Python instead.
script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),mustWork=TRUE)
source(file.path(dirname(script),"util/analysis/figure78_render_only.R"))
figure <- sub("^--figure=","",grep("^--figure=",commandArgs(TRUE),value=TRUE))
if(length(figure)>1L||any(!figure%in%c("7","8")))stop("Specify --figure=7 or --figure=8 once")
scope <- if(!length(figure))"all" else if(figure=="7")"all7" else "figure8"
f78_render_cli(script,default_scope=scope,validate_only=TRUE)
