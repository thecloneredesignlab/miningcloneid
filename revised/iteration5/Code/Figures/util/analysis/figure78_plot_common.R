# Canonical plotting implementation extracted from the accepted figures.
# Only names and input-path resolution changed; no scientific recomputation.
f78_require_files <- function (paths, label = "required input") 
{
    missing <- paths[!file.exists(paths)]
    if (length(missing)) {
        stop("Missing ", label, " file(s):\n", paste(missing, collapse = "\n"))
    }
    invisible(normalizePath(paths, mustWork = TRUE))
}
f78_read_tsv <- function (path) 
{
    f78_require_files(path)
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
f78_args <- function () 
{
    x <- commandArgs(TRUE)
    out <- list(action = "all", workspace = "", task_index = NA_integer_)
    for (s in x) {
        z <- strsplit(sub("^--", "", s), "=", fixed = TRUE)[[1]]
        if (length(z) == 2) 
            out[[gsub("-", "_", z[1])]] <- z[2]
    }
    out
}
f78_normalize_annotation_theme <- function (plot) 
{
    if (inherits(plot, "patchwork")) {
        th <- plot$patches$annotation$theme
        if (!inherits(th, "ggplot2::theme")) 
            plot$patches$annotation$theme <- do.call(ggplot2::theme, as.list(th))
        plot$patches$plots <- lapply(plot$patches$plots, f78_normalize_annotation_theme)
    }
    plot
}
f78_plot_theme <- function () 
ggplot2::theme_classic(base_size = 9, base_family = "Helvetica") + ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(), 
    plot.title = ggplot2::element_text(face = "bold", size = 11), strip.background = ggplot2::element_rect(fill = "#F2F2F2", 
        colour = "#BEBEBE"), panel.border = ggplot2::element_rect(fill = NA, colour = "#555555", linewidth = 0.22))
f78_context_colors <- c("in vivo" = "#0072B2", "in vitro" = "#CC79A7")
`%||%` <- function(x,y) if(is.null(x)) y else x
f78_publish <- function(a,plot,width,height,stem) {
  if(is.null(a$render_ctx)) stop("Use draw_Figure7.R or draw_Figure8.R")
  f78ro_publish(a,plot,width,height,stem)
}
