#!/usr/bin/env Rscript

main <- function() {
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=","",arg)))
workspace <- normalizePath(file.path(script_dir,"..",".."))
out <- Sys.getenv("FIGURE2_ANIMATION_OUTPUT_DIR",file.path(workspace,"Figures","figure2_animation"))
for (package in c("png","jsonlite","xml2")) {
  if (!requireNamespace(package,quietly=TRUE)) stop("Missing R package: ",package)
}
for (command in c("ffmpeg","ffprobe","shasum")) {
  if (!nzchar(Sys.which(command))) stop("Missing executable: ",command)
}
tmp <- tempfile("animation-validation-")
dir.create(tmp)
on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
stages <- read.delim(file.path(out,"storyboard.tsv"),check.names=FALSE)
checks <- data.frame(check=character(),passed=logical(),detail=character())
record <- function(name,ok,detail) {
  checks <<- rbind(checks,data.frame(check=name,passed=isTRUE(ok),detail=detail))
}
sha <- function(path) {
  sub(" .*","",system2("shasum",c("-a","256",shQuote(path)),stdout=TRUE))
}
source_paths <- file.path(script_dir,c("draw_Figure2.R","draw_Figure2_animation.R","validate_Figure2_animation.R"))
Sys.setenv(FIGURE2_DRAW_WORKER="1",FIGURE_WORKSPACE_ROOT=workspace,FIGURE2_OUTPUT_DIR=tmp)
src <- new.env(parent=globalenv())
sys.source(source_paths[1],envir=src)
png_original <- file.path(workspace,"Figures","assembled_fig2.png")
png_reference <- file.path(tmp,"assembled_fig2.png")
original <- png::readPNG(png_original)
reference <- png::readPNG(png_reference)
same_size <- identical(dim(original),dim(reference))
error <- if (same_size) mean(abs(original-reference)) else Inf
record("Full static figure regeneration",same_size && error<.01,
       sprintf("Mean absolute pixel error %.5f; tolerance allows platform font rendering differences.",error))

stage_files <- file.path(out,"stages",stages$png)
for (i in seq_along(stage_files)) {
  a <- png::readPNG(stage_files[i])
  nonwhite <- mean(apply(a[,,1:3,drop=FALSE],c(1,2),min)<.93)
  record(paste("Stage",i,"dimensions and content"),
    identical(dim(a)[1:2],c(1080L,1920L)) && nonwhite>.025,
    sprintf("1920 x 1080; %.1f%% nonwhite pixels",100*nonwhite))
}

probe <- system2("ffprobe",c("-v","error","-show_entries",
  "format=duration,size:stream=codec_name,width,height,r_frame_rate,nb_frames","-of","json",
  shQuote(file.path(out,"figure2_progressive_build.mp4"))),stdout=TRUE)
writeLines(probe,file.path(out,"video_metadata.json"))
meta <- jsonlite::fromJSON(paste(probe,collapse="\n"))
record("MP4 format and duration",meta$streams$codec_name=="h264" &&
  meta$streams$width==1920 && meta$streams$height==1080 &&
  meta$streams$r_frame_rate=="25/1" && abs(as.numeric(meta$format$duration)-tail(stages$end_seconds,1))<.045,
  paste("H.264; 1920 x 1080; 25 fps;",meta$format$duration,"seconds"))

for (i in seq_along(stage_files)) {
  file <- file.path(tmp,paste0("decoded_",i,".png"))
  code <- system2("ffmpeg",c("-hide_banner","-loglevel","error","-y",
    "-ss",as.character(stages$start_seconds[i]+2),"-i",shQuote(file.path(out,"figure2_progressive_build.mp4")),
    "-frames:v","1","-threads","2","-update","1",shQuote(file)))
  if (code!=0 || !file.exists(file)) stop("Could not decode MP4 frame for stage ",i)
  a <- png::readPNG(stage_files[i])[,,1:3]
  b <- png::readPNG(file)[,,1:3]
  error <- mean(abs(a-b))
  record(paste("Decoded hold",i,"matches stage"),code==0 && error<.008,
    sprintf("Mean absolute RGB error %.5f on [0,1] scale",error))
}
unzip(file.path(out,"figure2_progressive_build.pptx"),exdir=tmp)
ns <- c(p="http://schemas.openxmlformats.org/presentationml/2006/main")
presentation <- xml2::read_xml(file.path(tmp,"ppt/presentation.xml"))
size <- xml2::xml_find_first(presentation,"//p:sldSz",ns)
record("PowerPoint 16:9 canvas",xml2::xml_attr(size,"cx")=="12192000" &&
  xml2::xml_attr(size,"cy")=="6858000","13.3333 x 7.5 inches")
slide_files <- list.files(file.path(tmp,"ppt/slides"),pattern="^slide[0-9]+.xml$",full.names=TRUE)
note_files <- list.files(file.path(tmp,"ppt/notesSlides"),pattern="^notesSlide[0-9]+.xml$",full.names=TRUE)
record("Seven slides and seven speaker notes",length(slide_files)==7 && length(note_files)==7,
       "Each stage has a dedicated slide and notes page.")
for (file in slide_files) {
  slide <- xml2::read_xml(file)
  tr <- xml2::xml_find_first(slide,"//p:transition",ns)
  record(paste(basename(file),"click-only fade"),xml2::xml_attr(tr,"advClick")=="1" &&
    is.na(xml2::xml_attr(tr,"advTm")) && length(xml2::xml_find_all(tr,"p:fade",ns))==1,
    "Presenter controls advance; no automatic timing.")
}
media <- list.files(file.path(tmp,"ppt/media"),full.names=TRUE)
record("PowerPoint stage image fidelity",all(unname(tools::md5sum(stage_files)) %in% unname(tools::md5sum(media))),
       "All seven generated PNGs are embedded byte for byte.")
inputs <- c(source_paths,png_original)
write.csv(data.frame(path=substring(inputs,nchar(workspace)+2L),sha256=vapply(inputs,sha,character(1))),
  file.path(out,"input_hashes.csv"),row.names=FALSE)
outputs <- c(stage_files,file.path(out,c("figure2_progressive_build.mp4","figure2_progressive_build.pptx",
  "figure2_progressive_build.pdf","storyboard.tsv","storyboard.json","speaker_notes.md")))
write.csv(data.frame(path=substring(outputs,nchar(out)+2L),bytes=file.info(outputs)$size,
  sha256=vapply(outputs,sha,character(1))),file.path(out,"output_hashes.csv"),row.names=FALSE)
write.csv(checks,file.path(out,"validation_checks.csv"),row.names=FALSE)
writeLines(c("# Automated Validation","",sprintf("%d / %d checks passed.",sum(checks$passed),nrow(checks)),"",
  paste0("- ",ifelse(checks$passed,"PASS","FAIL"),": ",checks$check,". ",checks$detail),"",
  "Visual QA and interpretation boundaries are recorded in Code/Figures/README_Figure2.md."),file.path(out,"validation.md"))
if (!all(checks$passed)) stop("Validation failed; see validation_checks.csv")
message("All ",nrow(checks)," checks passed.")
}

if (sys.nframe() == 0L) main()
