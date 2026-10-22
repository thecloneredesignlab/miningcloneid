# Thin compatibility routing only. No old workers, configs, or task-ID mapping.
f78_compat_cli <- function(script,figure=NULL,role="generation") {
  args <- commandArgs(TRUE)
  code <- dirname(script)
  if("--help=yes"%in%args) {
    cat("Legacy filename; only the explicit current Figure7/8 interface is supported.\n",
      "Generation: run_Figure7.R / run_Figure8.R --stage=...\n",
      "Metadata validation: validate_Figure78.R --figure=7|8\n",
      "Static check: python3 Code/Figures/figure78_entry_tools.py static\n")
    return(invisible(NULL))
  }
  if(any(grepl("^--(action|task-index|task-ids|version|config|run-id)=",args)))
    stop("Legacy actions/task IDs/configs cannot be mapped safely. Use the current run_Figure7.R or run_Figure8.R interface")
  f <- grep("^--figure=",args,value=TRUE)
  if(length(f)) {
    if(length(f)!=1L||!sub("^--figure=","",f)%in%c("7","8"))stop("--figure must be 7 or 8")
    given <- as.integer(sub("^--figure=","",f))
    if(!is.null(figure)&&figure!=given)stop("Figure conflicts with compatibility filename")
    figure <- given;args <- args[!grepl("^--figure=",args)]
  }
  if(role=="generation") {
    if(is.null(figure)||!any(grepl("^--stage=",args)))
      stop("Specify --figure=7|8 where needed and an explicit --stage; no implicit full run is allowed through legacy filenames")
    target <- file.path(code,paste0("run_Figure",figure,".R"))
  } else {
    if(!"--mode=retained"%in%args)
      stop("Use validate_Figure78.R for retained metadata, or figure78_entry_tools.py static. This alias requires --mode=retained and never runs simulations")
    args <- args[args!="--mode=retained"]
    if(!is.null(figure))args <- c(args,paste0("--figure=",figure))
    target <- file.path(code,"validate_Figure78.R")
  }
  status <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(target),vapply(args,shQuote,"")))
  if(status!=0L)stop("Unified Figure7/8 entry returned status ",status)
  invisible(NULL)
}
