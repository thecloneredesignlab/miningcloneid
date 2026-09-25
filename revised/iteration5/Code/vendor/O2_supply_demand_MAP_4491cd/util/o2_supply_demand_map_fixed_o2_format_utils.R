# Shared fixed-O2 labels and stable path slugs.

num_path_tag <- function(x) {
  val <- suppressWarnings(as.numeric(x))
  if (length(val) != 1L || !is.finite(val)) return("NA")
  txt <- format(val, scientific = FALSE, trim = TRUE)
  txt <- sub("^-", "m", txt)
  txt <- gsub("\\.", "p", txt)
  txt <- gsub("[^A-Za-z0-9]+", "", txt)
  if (!nzchar(txt)) "NA" else txt
}

format_o2_label <- function(x) {
  paste0("O2 = ", format(as.numeric(x), scientific = FALSE, trim = TRUE))
}

initial_condition_from_ploidy <- function(x) {
  paste0("init_", format(as.numeric(x), scientific = FALSE, trim = TRUE), "N")
}

method_slug <- function(x) {
  x <- tolower(trimws(as.character(x[[1L]])))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "analytical" else x
}

method_label <- function(x) {
  x <- method_slug(x)
  if (identical(x, "eigen")) return("Eigen analytical")
  if (identical(x, "expm")) return("Expm analytical")
  paste0(x, " analytical")
}
