#!/usr/bin/env Rscript

main <- function() {
args <- commandArgs(TRUE)
if (any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript --vanilla draw_Figure2_animation.R [--stills]\n",
      "FIGURE2_ANIMATION_OUTPUT_DIR optionally overrides the output directory.\n")
  return(invisible(NULL))
}
if (length(setdiff(args, "--stills"))) stop("Unknown argument: ", paste(setdiff(args, "--stills"), collapse=", "))
packages <- if ("--stills" %in% args) "jsonlite" else c("jsonlite", "officer", "xml2", "zip")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing,collapse=", "))
if (!"--stills" %in% args && !nzchar(Sys.which("ffmpeg"))) stop("ffmpeg is required to encode the MP4")
if (!capabilities("cairo")) stop("R requires Cairo graphics support")
suppressPackageStartupMessages(library(grid))
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
workspace <- normalizePath(file.path(script_dir, "..", ".."))
out <- Sys.getenv("FIGURE2_ANIMATION_OUTPUT_DIR", file.path(workspace, "Figures", "figure2_animation"))
stills <- file.path(out, "stages")
scratch <- tempfile("figure2-animation-")
frames <- file.path(scratch, "frames")
on.exit(unlink(scratch, recursive=TRUE), add=TRUE)
for (d in c(out, stills, frames, file.path(scratch, "reference"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
Sys.setenv(FIGURE2_DRAW_WORKER = "1", FIGURE_WORKSPACE_ROOT = workspace,
           FIGURE2_OUTPUT_DIR = file.path(scratch, "reference"))
src <- new.env(parent = globalenv())
sys.source(file.path(script_dir, "draw_Figure2.R"), envir = src)

stages <- data.frame(
  stage = 0:6,
  name = c("three_nodes", "death_hazard", "variation_and_adaptation",
           "stress_induced_mutagenesis", "proliferation", "cin_linked_loss", "full_model"),
  hold = c(4, 6, 7, 12, 7, 6, 8),
  title = c(
    "Resource limitation pushes ploidy in two opposing directions",
    "Experienced resource stress links the environment to CIN",
    "Generating variation and retaining adapted descendants are different processes",
    "Stress-induced mutagenesis: adaptation can turn CIN back down",
    "Adapted composition also changes population growth",
    "Missegregation can also produce nonviable daughters",
    "Death-hazard-linked CIN and ploidy selection create a population-composition feedback"
  ),
  takeaway = c(
    "Resource costs can oppose high ploidy; buffering can favor its retention under CIN.",
    "CIN responds to the modeled death hazard, not to oxygen alone.",
    "CIN supplies karyotype variation; selection and retention reshape the population.",
    "At unchanged oxygen, adaptation can reduce experienced stress and population-average CIN.",
    "Resource limitation suppresses proliferation; better-performing states can expand.",
    "Stress-associated death and division-linked daughter loss are distinct sources of cell loss.",
    "A constant-probability WGD branch supplies an additional route into high-ploidy states."
  ),
  notes = c(
    paste("Start with the central tension from slide 6. Extra chromosomes can buffer segregation",
          "errors, but maintaining them can be costly under resource limitation. The simple CIN-to-ploidy",
          "edge summarizes variation and retention, not a claim that each error increases ploidy.",
          "Ask the audience whether resource limitation should raise or lower ploidy.",
          "The original slide's WGD shorthand is deliberately deferred to a separate branch."),
    paste("First unpack the resource-to-CIN arrow. Resource limitation changes the modeled",
          "death hazard, which represents experienced stress. The hazard increases inducible",
          "per-chromosome missegregation. A cell need not die in order to missegregate.",
          "The coarse direct resource-cost edge is set aside here and resolved through growth",
          "and death later in the build."),
    paste("Now unpack ploidy. Missegregation generates chromosome-number variation, and selection",
          "and survival filtering determine which descendants persist. Adapted composition is a",
          "population-level outcome, not a claim that every new variant is adaptive.",
          "Ploidy-dependent buffering is represented within selection and retention;",
          "the approved final diagram has no separate post-missegregation-survival node."),
    paste("Pause here. This is the core stress-induced-mutagenesis architecture.",
          "At the same environmental oxygen, a shift toward better-performing chromosome states",
          "can lower the population-average resource-stress hazard. Because inducible CIN depends",
          "on that hazard, average CIN can then decrease without oxygen changing.",
          "This is a model possibility, not a claim that all fitted populations must adapt."),
    paste("Complete the direct fitness effects. Resources inhibit proliferation, while changes",
          "in composition can favor faster-growing states. This restores the direct resource-cost",
          "side of the original shorthand in its explicit growth and death components.",
          "The newly drawn green growth connection does not mean every altered state grows faster."),
    paste("Missegregation can also produce daughters that fail to survive. Add the yellow",
          "upward arrow and broaden the aggregate label from Death hazard to Cell death.",
          "Crucially, this daughter-loss flux is separate from the continuous resource-stress hazard.",
          "The model does not feed total CIN-associated loss back into the inducible-CIN equation."),
    paste("Finally add WGD as a distinct constant-probability branch per division.",
          "Stress does not directly increase its per-division probability in this model.",
          "WGD changes the states supplied to selection, whereas ploidy-dependent buffering affects",
          "their retention. This is the unchanged current Figure 2A diagram, not a new model.",
          "Re-emphasize the core takeaway: changing population composition can change CIN at fixed oxygen.")
  ), stringsAsFactors = FALSE
)
stages$png <- sprintf("stage_%02d_%s.png", stages$stage + 1L, stages$name)
transition <- 1.6
fps <- 25L
stages$start_seconds <- c(0, cumsum(head(stages$hold, -1) + transition))
stages$end_seconds <- stages$start_seconds + stages$hold
write.table(stages, file.path(out, "storyboard.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
jsonlite::write_json(stages, file.path(out, "storyboard.json"), pretty = TRUE, auto_unbox = TRUE)

clamp <- function(x) pmin(1, pmax(0, x))
ease <- function(x) { x <- clamp(x); x * x * (3 - 2 * x) }
with_alpha <- function(a, code) {
  if (a <= 0) return(invisible(NULL))
  pushViewport(viewport(gp = gpar(alpha = clamp(a))))
  force(code)
  popViewport()
}

# Reveal paths by physical arc length, keeping arrowheads attached to their tips.
path <- function(x, y, col, progress = 1, inhibitory = FALSE, lwd = 1.45) {
  progress <- clamp(progress)
  if (progress <= 0) return(invisible(NULL))
  lengths <- c(0, cumsum(sqrt(diff(x * 7.1)^2 + diff(y * 6.35)^2)))
  target <- progress * tail(lengths, 1)
  j <- min(which(lengths >= target)[1], length(x))
  if (j <= 1) return(invisible(NULL))
  f <- (target - lengths[j - 1]) / (lengths[j] - lengths[j - 1])
  xx <- c(x[seq_len(j - 1)], x[j - 1] + f * (x[j] - x[j - 1]))
  yy <- c(y[seq_len(j - 1)], y[j - 1] + f * (y[j] - y[j - 1]))
  grid.lines(xx, yy, arrow = if (!inhibitory) arrow(length = unit(2.1, "mm"), type = "closed") else NULL,
             gp = gpar(col = col, fill = col, lwd = lwd, lineend = "round", linejoin = "round"))
  if (inhibitory && progress >= 0.995) {
    dx <- diff(tail(xx, 2)) * 7.1
    dy <- diff(tail(yy, 2)) * 6.35
    n <- sqrt(dx^2 + dy^2)
    grid.lines(unit(rep(tail(xx, 1), 2), "npc") + unit(c(-1, 1) * (-dy / n) * 1.5, "mm"),
               unit(rep(tail(yy, 1), 2), "npc") + unit(c(-1, 1) * (dx / n) * 1.5, "mm"),
               gp = gpar(col = col, lwd = lwd, lineend = "butt"))
  }
}
bezier_points <- function(x, y) {
  t <- seq(0, 1, length.out = 70)
  b <- cbind((1-t)^3, 3*(1-t)^2*t, 3*(1-t)*t^2, t^3)
  cbind(drop(b %*% x), drop(b %*% y))
}
feedback_curve <- rbind(
  bezier_points(c(.290,.330,.490,.625), c(.580,.710,.745,.755)),
  bezier_points(c(.625,.649,.640,.667), c(.755,.755,.790,.790))[-1,]
)

draw_partial <- function(q) {
  s <- floor(q + .5) + 1L
  a <- ease(q)
  split <- ease(q - 1)
  feedback <- ease(q - 2)
  growth <- ease(q - 3)
  loss <- ease(q - 4)
  wgd <- ease(q - 5)
  src$draw_panel(.5, .735, .96, .50, "", stages$title[s], fill = "#FCFDFE", title_size = 8.8)
  grid.text("Environmental input", .50, .875, gp = src$gp_text(6.8, src$muted))

  cin_x <- .19 + .62 * a
  cin_y <- .545 + .135 * sqrt(a)
  # The original triangle gives way to the explicitly resolved hazard path.
  with_alpha(1 - ease(q / .20), {
    path(c(.42,.25), c(.783,.592), src$magenta)
    grid.text("Stress can increase CIN", .19, .72, gp = src$gp_text(7.4, src$magenta_dark))
    path(c(.585,.740,.780), c(.783,.613,.582), src$blue, inhibitory = TRUE)
    grid.text("Resource costs", .75, .715, gp = src$gp_text(7.4, src$blue_dark))
    path(c(.325,.675), c(.545,.545), src$magenta)
    grid.text("Variation and retention", .50, .57, gp = src$gp_text(7.0, src$green_dark))
  })
  new_path <- ease((q - .75) / .25)
  with_alpha(new_path, {
    path(c(.635,.675), c(.825,.825), src$amber, new_path)
    path(c(.920,cin_x+.110), c(.775,cin_y+.045), src$magenta, new_path, lwd = 1.6)
    grid.text("Stress-induced\nCIN", .855,.750, gp = src$gp_text(7.0, src$magenta_dark, "bold"))
    path(c(cin_x+.135*(1-a),.675+.135*a), c(cin_y-.045*a,.545+.035*a),
         src$magenta, new_path, lwd = 1.5)
  })
  selection <- ease((split-.5)/.5)
  with_alpha(selection, {
    path(c(.730,.675,.365,.305), c(.580,.625,.625,.580), src$green, selection)
    grid.text("Selection and retention", .5,.645, gp = src$gp_text(7.0, src$green_dark))
  })
  with_alpha(feedback, {
    path(feedback_curve[,1], feedback_curve[,2], src$green, feedback, inhibitory = TRUE)
    src$draw_arrow(.080,.685,.110,.685,col=src$muted,length_mm=1.8)
    grid.text("Activation",.122,.685,just="left",gp=src$gp_text(7.0,src$muted))
    src$draw_arrow(.080,.660,.110,.660,col=src$muted,ends="none")
    src$draw_inhibitory_bar(.110,.660,col=src$muted,length_mm=2.4)
    grid.text("Inhibition",.122,.660,just="left",gp=src$gp_text(7.0,src$muted))
  })
  with_alpha(growth, {
    path(c(.365,.331),c(.825,.825),src$blue,growth,inhibitory=TRUE)
    path(c(.240,.240),c(.580,.790),src$green,growth)
  })
  with_alpha(loss, {
    path(c(.705,.705),c(.725,.775),src$amber,loss,lwd=1.6)
    grid.text("Daughter\ncell loss",.755,.750,gp=src$gp_text(7.0,src$amber_dark,"bold"))
  })
  with_alpha(wgd, path(c(.635,.675),c(.545,.545),src$muted,wgd,lwd=1.25))

  src$draw_labeled_card(.5,.825,.27,.080,"Resource limitation","modeled through O\u2082",
    fill=src$blue_light,border=src$blue,title_size=8.6,subtitle_size=7.0,
    title_col=src$navy,subtitle_col=src$blue_dark)
  with_alpha(ease((q-.25)/.60), src$draw_labeled_card(.810,.825,.27,.095,
    if (q < 4.5) "Death hazard" else "Cell death",
    if (q < 4.5) "experienced resource stress" else "resource stress + CIN-associated loss",
    fill=src$amber_light,border=src$amber,title_size=7.5,subtitle_size=6.8,title_col=src$amber_dark))
  src$draw_labeled_card(cin_x,cin_y,.27,.085,
    if (q < .5) "CIN" else "Chromosome missegregation\nincreases",
    if (q < .5) "chromosome missegregation" else "effective per-chromosome probability",
    fill=src$magenta_light,border=src$magenta,title_size=if(q<.5) 9.5 else 7.1,
    subtitle_size=6.4,title_col=src$magenta_dark)
  if (split < .5) {
    src$draw_labeled_card(.810,.545,.27,.065,"Ploidy","mean chromosome number",
      fill=src$green_light,border=src$green,title_size=8.5,subtitle_size=7.0,
      title_col=src$green_dark,subtitle_col=src$green_dark)
  } else {
    src$draw_card(.810,.545,.27,.065,"Karyotype variation",fill=src$green_light,
      border=src$green,size=7.7,col=src$green_dark,fontface="bold")
  }
  with_alpha(ease((split-.30)/.70), {
    src$draw_labeled_card(.810-.620*split,.545,.27,.065,
      "Adapted karyotype /","ploidy composition",fill=src$green_light,border=src$green,
      title_size=7.8,subtitle_size=7.0,title_col=src$green_dark,subtitle_col=src$green_dark)
  })
  with_alpha(growth, src$draw_card(.190,.825,.27,.070,"Proliferation rate",
    fill="#F4F5F6",border=src$muted,size=7.8,col=src$ink,fontface="bold"))
  with_alpha(wgd, src$draw_labeled_card(.500,.545,.27,.065,
    "WGD generation","constant probability per division",fill="#F4F5F6",border=src$muted,
    title_size=7.7,subtitle_size=7.0))
}

draw_frame <- function(q) {
  grid.newpage()
  grid.rect(gp=gpar(fill="white",col=NA))
  # Preserve the physical geometry of Figure 2A; leave room below for one takeaway.
  pushViewport(viewport(x=.5,y=unit(-2.60475,"in"),width=unit(7.1,"in"),
    height=unit(6.35,"in"),just=c("centre","bottom"),clip="off"))
  if (q >= 6) src$draw_panel_a() else draw_partial(q)
  popViewport()
  s <- floor(q+.5)+1L
  grid.text(stages$takeaway[s],.5,unit(.245,"in"),
    gp=src$gp_text(8.6,if (s==4) src$green_dark else src$ink,if(s==4) "bold" else "plain"))
}
render_png <- function(q, file) {
  png(file,width=1920,height=1080,res=1920/7.1,type="cairo",bg="white")
  draw_frame(q)
  invisible(dev.off())
}
for (i in seq_len(nrow(stages))) render_png(stages$stage[i],file.path(stills,stages$png[i]))
cairo_pdf(file.path(out,"figure2_progressive_build.pdf"),width=7.1,height=7.1*9/16,family="Helvetica")
for (q in 0:6) draw_frame(q)
invisible(dev.off())
writeLines(c("# Speaker Notes",unlist(lapply(seq_len(nrow(stages)),function(i) {
  c("",paste0("## ",i,". ",stages$title[i]),"",stages$notes[i])
}))),file.path(out,"speaker_notes.md"))

if ("--stills" %in% args) return(invisible(out))

concat <- character()
add_clip <- function(file,duration) {
  concat <<- c(concat,paste0("file '",normalizePath(file),"'"),
               sprintf("duration %.8f",duration))
}
for (i in seq_len(nrow(stages))) {
  add_clip(file.path(stills,stages$png[i]),stages$hold[i])
  if (i < nrow(stages)) {
    n <- round(transition*fps)
    for (j in seq_len(n)) {
      file <- file.path(frames,sprintf("transition_%02d_%03d.png",i,j))
      render_png(i-1+j/n,file)
      add_clip(file,1/fps)
    }
  }
  message("Rendered stage ",i," of ",nrow(stages))
}
concat <- c(concat,paste0("file '",normalizePath(file.path(stills,tail(stages$png,1))),"'"))
concat_path <- file.path(frames,"concat.txt")
writeLines(concat,concat_path)
video <- file.path(out,"figure2_progressive_build.mp4")
code <- system2("ffmpeg",c("-hide_banner","-loglevel","warning","-y","-f","concat","-safe","0",
  "-i",shQuote(concat_path),"-an","-vf",paste0("fps=",fps),"-c:v","libx264","-threads","2",
  "-preset","medium","-crf","18","-pix_fmt","yuv420p","-movflags","+faststart",
  "-t",as.character(tail(stages$end_seconds,1)),shQuote(video)))
if (code != 0) stop("Video encoding failed")

ppt <- officer::read_pptx()
for (i in seq_len(nrow(stages))) {
  ppt <- officer::add_slide(ppt,layout="Blank",master="Office Theme")
  ppt <- officer::ph_with(ppt,officer::external_img(file.path(stills,stages$png[i]),
    width=13.333333,height=7.5,alt=stages$takeaway[i]),
    location=officer::ph_location(left=0,top=0,width=13.333333,height=7.5))
  ppt <- officer::set_notes(ppt,stages$notes[i],location=officer::notes_location_type("body"))
}
ppt_path <- file.path(out,"figure2_progressive_build.pptx")
print(ppt,target=ppt_path)
# Set the slide canvas and add click-only fades using structured OOXML APIs.
tmp <- file.path(scratch, "pptx")
dir.create(tmp)
unzip(ppt_path,exdir=tmp)
ns <- c(p="http://schemas.openxmlformats.org/presentationml/2006/main")
presentation_file <- file.path(tmp,"ppt","presentation.xml")
presentation <- xml2::read_xml(presentation_file)
size <- xml2::xml_find_first(presentation,"//p:sldSz",ns)
xml2::xml_set_attr(size,"cx","12192000")
xml2::xml_set_attr(size,"cy","6858000")
xml2::xml_set_attr(size,"type","screen16x9")
xml2::write_xml(presentation,presentation_file)
for (i in seq_len(nrow(stages))) {
  file <- file.path(tmp,"ppt","slides",paste0("slide",i,".xml"))
  slide <- xml2::read_xml(file)
  tr <- xml2::read_xml(paste0('<p:transition xmlns:p="',ns[[1]],'" spd="slow" advClick="1"><p:fade/></p:transition>'))
  xml2::xml_add_child(xml2::xml_root(slide),tr)
  xml2::write_xml(slide,file)
}
zip::zipr(ppt_path,files=list.files(tmp,all.files=TRUE,no..=TRUE),root=tmp)
unlink(tmp,recursive=TRUE)
message("Finished: ",video,"\n",ppt_path)
invisible(out)
}

if (sys.nframe() == 0L) main()
