#!/usr/bin/env Rscript

# Shared, side-effect-free HTML text escaping for report assemblers.

o2sd_html_escape_standard <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

o2sd_html_escape_full <- function(x) {
  if (is.null(x) || !length(x)) x <- ""
  x <- o2sd_html_escape_standard(x)
  gsub("'", "&#39;", x, fixed = TRUE)
}

o2sd_html_escape_minimal <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

# Shared, display-only image lightbox for generated HTML reports.  The report
# body, figures, captions, and embedded image payloads remain unchanged; the
# component is injected immediately before </body> and binds to existing <img>
# elements after the page loads.
o2sd_report_image_lightbox_version <- "1"

o2sd_report_image_lightbox_assets <- function() {
  r"---(
<!-- O2SD_REPORT_IMAGE_LIGHTBOX_START version=1 -->
<style id="o2sd-report-image-lightbox-style" data-o2sd-report-image-lightbox-version="1">
img[data-o2sd-lightbox-bound="1"]{cursor:zoom-in}
#o2sd-report-image-lightbox[hidden]{display:none!important}
#o2sd-report-image-lightbox{position:fixed;z-index:2147483000;inset:0;display:flex;flex-direction:column;background:rgba(8,10,14,.94);color:#fff;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
#o2sd-report-lightbox-toolbar{display:flex;align-items:center;gap:8px;min-height:54px;padding:8px 12px;box-sizing:border-box;border-bottom:1px solid rgba(255,255,255,.18);background:rgba(15,17,22,.96)}
#o2sd-report-lightbox-title{min-width:0;flex:1;overflow:hidden;color:#f8fafc;font:500 13px/18px -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;text-overflow:ellipsis;white-space:nowrap}
#o2sd-report-lightbox-hint{color:#cbd5e1;font:500 12px/18px -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;white-space:nowrap}
.o2sd-report-lightbox-button{height:34px;min-width:34px;padding:0 10px;border:1px solid rgba(255,255,255,.28);border-radius:8px;background:rgba(255,255,255,.08);color:#fff;font:600 12px/18px -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;cursor:pointer}
.o2sd-report-lightbox-button:hover{background:rgba(255,255,255,.16)}
.o2sd-report-lightbox-button:focus-visible{outline:2px solid #93c5fd;outline-offset:2px}
#o2sd-report-lightbox-viewport{position:relative;flex:1;min-height:0;overflow:hidden;touch-action:none;cursor:zoom-out}
#o2sd-report-lightbox-viewport.is-dragging{cursor:grabbing}
#o2sd-report-lightbox-image{position:absolute;top:50%;left:50%;display:block;max-width:none!important;max-height:none!important;width:auto!important;height:auto!important;transform-origin:center center;cursor:grab;user-select:none;-webkit-user-drag:none;will-change:transform}
#o2sd-report-lightbox-viewport.is-dragging #o2sd-report-lightbox-image{cursor:grabbing}
body.o2sd-report-lightbox-open{overflow:hidden!important}
@media(max-width:760px){#o2sd-report-lightbox-hint{display:none}#o2sd-report-lightbox-toolbar{gap:5px;padding:7px}.o2sd-report-lightbox-button{padding:0 8px}}
</style>
<div id="o2sd-report-image-lightbox" role="dialog" aria-modal="true" aria-labelledby="o2sd-report-lightbox-title" hidden>
  <div id="o2sd-report-lightbox-toolbar">
    <div id="o2sd-report-lightbox-title">Report figure</div>
    <div id="o2sd-report-lightbox-hint">Wheel: zoom · Drag: pan · Esc: close</div>
    <button class="o2sd-report-lightbox-button" id="o2sd-report-lightbox-zoom-out" type="button" aria-label="Zoom out">−</button>
    <button class="o2sd-report-lightbox-button" id="o2sd-report-lightbox-zoom-in" type="button" aria-label="Zoom in">+</button>
    <button class="o2sd-report-lightbox-button" id="o2sd-report-lightbox-fit" type="button">Fit</button>
    <button class="o2sd-report-lightbox-button" id="o2sd-report-lightbox-actual" type="button">1:1</button>
    <button class="o2sd-report-lightbox-button" id="o2sd-report-lightbox-close" type="button" aria-label="Close image viewer">Close</button>
  </div>
  <div id="o2sd-report-lightbox-viewport">
    <img id="o2sd-report-lightbox-image" alt="">
  </div>
</div>
<script id="o2sd-report-image-lightbox-script" data-o2sd-report-image-lightbox-version="1">
(function(){
  "use strict";
  if(window.__o2sdReportImageLightboxVersion==="1")return;
  window.__o2sdReportImageLightboxVersion="1";
  var lightbox=document.getElementById("o2sd-report-image-lightbox");
  var viewport=document.getElementById("o2sd-report-lightbox-viewport");
  var image=document.getElementById("o2sd-report-lightbox-image");
  var title=document.getElementById("o2sd-report-lightbox-title");
  var closeButton=document.getElementById("o2sd-report-lightbox-close");
  var scale=1,fitScale=1,panX=0,panY=0,dragging=false,pointerId=null;
  var dragStartX=0,dragStartY=0,dragPanX=0,dragPanY=0,backgroundPress=false,pointerMoved=false;
  var opener=null;
  function updateTransform(){image.style.transform="translate(calc(-50% + "+panX+"px),calc(-50% + "+panY+"px)) scale("+scale+")";}
  function imageWidth(){return image.naturalWidth||Number(image.getAttribute("width"))||1;}
  function imageHeight(){return image.naturalHeight||Number(image.getAttribute("height"))||1;}
  function reset(mode){
    if(!viewport||!image.src)return;
    var rect=viewport.getBoundingClientRect();
    var paddedWidth=Math.max(120,rect.width-32);
    var paddedHeight=Math.max(120,rect.height-32);
    fitScale=Math.min(paddedWidth/imageWidth(),paddedHeight/imageHeight(),1);
    scale=mode==="actual"?1:fitScale;
    panX=0;panY=0;updateTransform();
  }
  function zoom(factor,clientX,clientY){
    if(!viewport||!image.src)return;
    var previous=scale;
    var minimum=Math.max(fitScale*.35,.025);
    var next=Math.min(12,Math.max(minimum,previous*factor));
    if(Math.abs(next-previous)<.0001)return;
    var rect=viewport.getBoundingClientRect();
    var originX=clientX==null?rect.left+rect.width/2:clientX;
    var originY=clientY==null?rect.top+rect.height/2:clientY;
    var dx=originX-rect.left-rect.width/2-panX;
    var dy=originY-rect.top-rect.height/2-panY;
    var ratio=next/previous;
    panX+=dx*(1-ratio);panY+=dy*(1-ratio);scale=next;updateTransform();
  }
  function close(){
    if(!lightbox||lightbox.hidden)return;
    lightbox.hidden=true;document.body.classList.remove("o2sd-report-lightbox-open");
    dragging=false;backgroundPress=false;pointerMoved=false;pointerId=null;viewport.classList.remove("is-dragging");
    if(opener&&typeof opener.focus==="function")opener.focus({preventScroll:true});
  }
  function labelFor(img){
    var figure=img.closest?img.closest("figure"):null;
    var caption=figure?figure.querySelector("figcaption"):null;
    return (caption&&caption.textContent.trim())||img.getAttribute("alt")||img.getAttribute("title")||"Report figure";
  }
  function open(img){
    if(!lightbox||!viewport||!image)return;
    var src=img.currentSrc||img.getAttribute("src");if(!src)return;
    opener=img;title.textContent=labelFor(img);image.alt=title.textContent;
    image.onload=function(){reset("fit");};image.src=src;
    lightbox.hidden=false;document.body.classList.add("o2sd-report-lightbox-open");
    requestAnimationFrame(function(){reset("fit");closeButton.focus({preventScroll:true});});
  }
  function eligible(img){
    if(!img||lightbox.contains(img)||img.dataset.o2sdLightboxIgnore==="1")return false;
    if(img.getAttribute("role")==="presentation"||img.getAttribute("aria-hidden")==="true")return false;
    return Boolean(img.currentSrc||img.getAttribute("src"));
  }
  function bind(img){
    if(!eligible(img)||img.dataset.o2sdLightboxBound==="1")return;
    img.dataset.o2sdLightboxBound="1";img.setAttribute("role","button");img.setAttribute("tabindex","0");
    if(!img.getAttribute("title"))img.setAttribute("title","Click to inspect this image");
    img.addEventListener("click",function(event){event.preventDefault();event.stopPropagation();open(img);});
    img.addEventListener("keydown",function(event){if(event.key==="Enter"||event.key===" "){event.preventDefault();open(img);}});
  }
  function bindAll(root){
    if(root&&root.nodeType===1&&root.matches&&root.matches("img"))bind(root);
    var scope=root&&root.querySelectorAll?root:document;
    Array.prototype.forEach.call(scope.querySelectorAll("img"),bind);
  }
  closeButton.addEventListener("click",close);
  document.getElementById("o2sd-report-lightbox-zoom-in").addEventListener("click",function(){zoom(1.25);});
  document.getElementById("o2sd-report-lightbox-zoom-out").addEventListener("click",function(){zoom(.8);});
  document.getElementById("o2sd-report-lightbox-fit").addEventListener("click",function(){reset("fit");});
  document.getElementById("o2sd-report-lightbox-actual").addEventListener("click",function(){reset("actual");});
  viewport.addEventListener("wheel",function(event){event.preventDefault();zoom(event.deltaY<0?1.18:1/1.18,event.clientX,event.clientY);},{passive:false});
  viewport.addEventListener("pointerdown",function(event){
    if(event.button!==0)return;pointerId=event.pointerId;pointerMoved=false;dragStartX=event.clientX;dragStartY=event.clientY;
    if(event.target===image){event.preventDefault();dragging=true;backgroundPress=false;dragPanX=panX;dragPanY=panY;viewport.setPointerCapture(event.pointerId);viewport.classList.add("is-dragging");}
    else if(event.target===viewport){dragging=false;backgroundPress=true;viewport.setPointerCapture(event.pointerId);}
  });
  viewport.addEventListener("pointermove",function(event){
    if(event.pointerId!==pointerId)return;var dx=event.clientX-dragStartX,dy=event.clientY-dragStartY;
    if(Math.hypot(dx,dy)>4)pointerMoved=true;if(!dragging)return;panX=dragPanX+dx;panY=dragPanY+dy;updateTransform();
  });
  function finishPointer(event,cancelled){
    if(event.pointerId!==pointerId)return;var closeFromBackground=!cancelled&&backgroundPress&&!pointerMoved;
    dragging=false;backgroundPress=false;pointerMoved=false;pointerId=null;viewport.classList.remove("is-dragging");
    if(viewport.hasPointerCapture(event.pointerId))viewport.releasePointerCapture(event.pointerId);if(closeFromBackground)close();
  }
  viewport.addEventListener("pointerup",function(event){finishPointer(event,false);});
  viewport.addEventListener("pointercancel",function(event){finishPointer(event,true);});
  document.addEventListener("keydown",function(event){
    if(lightbox.hidden)return;
    if(event.key==="Escape"){event.preventDefault();close();}
    else if(event.key==="+"||event.key==="="){event.preventDefault();zoom(1.25);}
    else if(event.key==="-"||event.key==="_"){event.preventDefault();zoom(.8);}
    else if(event.key==="0"){event.preventDefault();reset("fit");}
    else if(event.key==="1"){event.preventDefault();reset("actual");}
  });
  window.addEventListener("resize",function(){if(!lightbox.hidden)reset("fit");},{passive:true});
  bindAll(document);
  new MutationObserver(function(records){records.forEach(function(record){Array.prototype.forEach.call(record.addedNodes,bindAll);});}).observe(document.documentElement,{childList:true,subtree:true});
})();
</script>
<!-- O2SD_REPORT_IMAGE_LIGHTBOX_END -->
)---"
}

o2sd_inject_report_image_lightbox <- function(html_path, require_images = TRUE) {
  html_path <- normalizePath(path.expand(html_path), mustWork = TRUE)
  info <- file.info(html_path)
  con <- file(html_path, open = "rb")
  on.exit(close(con), add = TRUE)
  html <- readChar(con, nchars = info$size, useBytes = TRUE)
  close(con)
  on.exit(NULL, add = FALSE)

  marker <- "id=\"o2sd-report-image-lightbox-script\""
  if (grepl(marker, html, fixed = TRUE)) {
    return(invisible(list(path = html_path, status = "already_current", changed = FALSE)))
  }
  if (isTRUE(require_images) && !grepl("<img", html, fixed = TRUE)) {
    return(invisible(list(path = html_path, status = "skipped_no_images", changed = FALSE)))
  }

  assets <- o2sd_report_image_lightbox_assets()
  body_pos <- regexpr("</body>", html, fixed = TRUE)[[1L]]
  if (body_pos < 0L) {
    updated <- paste0(html, "\n", assets)
  } else {
    updated <- paste0(
      substr(html, 1L, body_pos - 1L),
      "\n", assets, "\n",
      substr(html, body_pos, nchar(html, type = "chars"))
    )
  }

  tmp <- tempfile(pattern = paste0(".", basename(html_path), "."), tmpdir = dirname(html_path))
  tmp_con <- file(tmp, open = "wb")
  on.exit({
    try(close(tmp_con), silent = TRUE)
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)
  writeChar(updated, tmp_con, eos = NULL, useBytes = TRUE)
  close(tmp_con)
  Sys.chmod(tmp, mode = info$mode)
  if (!file.rename(tmp, html_path)) {
    if (!file.copy(tmp, html_path, overwrite = TRUE, copy.mode = TRUE)) {
      stop("Failed to replace HTML report after lightbox injection: ", html_path, call. = FALSE)
    }
    unlink(tmp)
  }
  invisible(list(path = html_path, status = "injected", changed = TRUE))
}
