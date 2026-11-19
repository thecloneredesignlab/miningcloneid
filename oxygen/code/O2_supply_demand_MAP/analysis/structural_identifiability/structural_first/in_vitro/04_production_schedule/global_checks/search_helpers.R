export_changes <- function(outputs,pairs,schedule,path) {
  all <- list()
  oxygen <- sort(unique(schedule$oxygen_pct))
  masks <- expand.grid(N=22:154,O=oxygen)
  m <- lapply(masks$N,function(n)seq.int(0L,n))
  mask_labels <- paste(rep(masks$O,lengths(m)),rep(masks$N,lengths(m)),unlist(m),sep="|")
  for (name in names(pairs)) {
    id <- as.character(pairs[[name]])
    a <- outputs[[id[1]]];b <- outputs[[id[2]]]
    append_diff <- function(kind,keys,x,y) {
      changed <- which(!is.na(x)&!is.na(y)&x!=y)
      if(length(changed)) all[[length(all)+1L]]<<-data.frame(pair=name,component=kind,key=keys[changed],
            left=as.character(x[changed]),right=as.character(y[changed]))
    }
    append_diff("endpoint_day",names(a$days),a$days,b$days)
    append_diff("parent_day",names(a$parent_days),a$parent_days,b$parent_days)
    append_diff("eligible_fallback",names(a$days),a$fallback,b$fallback)
    append_diff("pruning_keep_O_N_m",mask_labels,a$mask,b$mask)
  }
  sf_emit(do.call(rbind,all),path)
}

emit_native_pair_payload <- function(outputs,value_path,branch_path) {
  output_rows <- do.call(rbind,lapply(names(outputs),function(id) {
    r<-outputs[[id]]
    data.frame(id=id,stream=c(rep("growth",114L),rep("kary",length(r$kary))),passage=c(names(r$growth),rep(colnames(r$kary),each=133L)),
      chr=c(rep("",114L),rep(rownames(r$kary),12L)),value=c(r$growth,as.numeric(r$kary)))
  }))
  branch_rows <- do.call(rbind,lapply(names(outputs),function(id) {
    r<-outputs[[id]]
    k<-strsplit(names(r$days),"|",fixed=TRUE)
    data.frame(id=id,cohort=vapply(k,function(a)a[1L],""),segment=vapply(k,function(a)a[2L],""),
       endpoint_day=unname(r$days),parent_day=unname(r$parent_days),inoculum_fallback=r$fallback)
  }))
  sf_emit(output_rows,value_path)
  sf_emit(branch_rows,branch_path)
}
