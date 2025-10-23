##origins must be sorted according to timepoint of sample acquisition
# origins=c( "SNU-668_C_A4_seed"   , "SNU-668_P2_A18k_seed", "SNU-668_P0_A11K_seed")
# origins=c("SNU-668_C_A24_seed","SNU-668_C_A4_seed","SNU-668_G1_A4_seed","SNU-668_G1_A10_seed" )
clusterKaryotypes <- function(
    origins,
    whichP          = "GenomePerspective",
    depth           = 1,
    path2lanscape   = NULL,
    numClusters     = NULL,
    capping         = NULL,
    method          = "complete",
    chrwhole        = NULL
){
  if(is.null(chrwhole)){
    devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
    x <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz", 
               col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
    chrarms=x[ , .(length = sum(chromEnd - chromStart)),by = .(chrom, arm = substring(name, 1, 1)) ]
    chrwhole=grpstats(as.matrix(chrarms$length),chrarms$chrom, "sum")$sum
  }
  
  source("~/Repositories/ALFA-K/utils/sim_setup_functions.R")
  source("~/Repositories/ALFA-K/utils/ALFA-K.R")
  ploidy  <- 2
  min_obs <- 5
  dt      <- 5
  
  mydb  <- cloneid::connect2DB()
  X     <- list()
  
  # 1) Query subprofiles from DB and store them in X
  for(biopsy in origins){
    stmt <- paste0(
      "SELECT cloneID, size, alias, parent 
       FROM Perspective
       WHERE size=1 AND whichPerspective='", whichP, "' AND origin='", biopsy, "'"
    )
    rs <- suppressWarnings(DBI::dbSendQuery(mydb, stmt))
    sps <- DBI::fetch(rs, n = -1)
    
    # If depth > 1, also pull children of these clones
    if(depth > 1) {
      stmt <- paste0(
        "SELECT cloneID, size, alias, parent 
         FROM Perspective 
         WHERE parent IN (", paste(sps$cloneID, collapse=","), ")"
      )
      rs <- suppressWarnings(DBI::dbSendQuery(mydb, stmt))
      sps <- DBI::fetch(rs, n=-1)
    }
    
    # Get the copy-number (or other) profiles for each clone
    x <- sapply(
      sps$cloneID, 
      function(cid) cloneid::getSubProfiles(cloneID_or_sampleName = cid, whichP = whichP), 
      simplify = FALSE
    )
    X[[biopsy]] <- do.call(cbind, x)
  }
  
  ##############################################################################
  # 2) Merge across origins for karyotyping/clustering
  ##############################################################################
  #   - We get CN calls from getKaryo() 
  #   - We cluster them with findBestClustering()
  
  # cnts is a list of data frames (one per origin) with copy-number calls
  cnts_list  <- sapply(X, function(mat) getKaryo(t(mat), ploidy)$cn, simplify = FALSE)
  
  # 'sampleID' identifies the biopsy/origin for each row in combined data
  sampleID <- unlist(
    sapply(names(cnts_list), function(x) rep(x, nrow(cnts_list[[x]])))
  )
  
  # Combine all rows from all origins
  cnts_combined <- do.call(rbind, cnts_list)
  
  
  # Dendrogram parameters
  hFun <- function(x) stats::hclust(x, method = method); #ward.D2
  dFun <- function(x) chrWeightedDist(x, chrwhole=chrwhole); #function(x) stats::dist(x, method = "manhattan")
  
  # If the user specified a number of clusters, we cluster with that many groups.
  # findBestClustering() presumably returns cluster labels from 0..(numClusters-1);
  # The +1 is presumably to shift them to 1..numClusters range.
  clusters <- findBestClustering(cnts_combined, numClusters = numClusters, hFun=hFun, dFun = dFun) + 1
  
  ##############################################################################
  # 3) Plot heatmap with hierarchical clustering to extract dendrogram 
  ##############################################################################
  tmp <- substr(paste(origins, collapse = "__"), 1, 90)
  pdf(paste0(tmp, ".pdf"))
  
  # Color each sampleID differently (for RowSideColors)
  uniqueIDs <- unique(sampleID)
  colVec    <- rep("NA", length(uniqueIDs))
  names(colVec) <- uniqueIDs
  
  # Example scheme: control is gray, everything else from a Brewer palette
  idxControl <- grep("C_", names(colVec))  # or any other pattern for control
  colVec[idxControl] <- gray.colors(length(idxControl))
  
  # The rest get colored from a palette
  remaining <- setdiff(seq_along(colVec), idxControl)
  if(length(remaining) > 0) {
    colPalette <- brewer.pal(min(length(remaining), 12), "Paired")
    if(length(remaining) > length(colPalette)) {
      # Extend the palette if needed
      colPalette <- colorRampPalette(colPalette)(length(remaining))
    }
    colVec[remaining] <- colPalette[seq_along(remaining)]
  }
  
  # Draw the heatmap
  tmp=as.matrix(cnts_combined)
  if(!is.null(capping)){
    tmp[tmp>capping] = capping 
  }
  hm <- heatmap.2(
    x           = tmp,
    margins     = c(15,15),
    colRow      = clusters[rownames(cnts_combined)], 
    trace       = 'none',
    Colv        = TRUE,
    dendrogram  = "row",
    RowSideColors = colVec[sampleID],
    key.xlab    = "copy number",
    key.title   = "",
    col         = matlab::fliplr(rainbow(20))[5:12],
    hclustfun   = hFun,
    distfun     = dFun
  )
  
  legend(
    "topright",
    legend = names(colVec),
    fill   = colVec,
    cex    = 0.5
  )
  
  # Boxplot of ploidy by cluster
  ploidy_vals <- calcPloidy(cnts_combined,chrwhole)
  boxplot(
    ploidy_vals ~ factor(clusters[names(ploidy_vals)], levels = unique(clusters)),
    xlab   = "Cluster",
    ylab   = "Ploidy",
    main   = "",
    col    = unique(clusters)
  )
  
  dev.off()
  
  ##############################################################################
  # 4) Summarize and return
  ##############################################################################
  
  # Combine the results by sample
  # grpstats() presumably aggregates mean/median by sample ID
  # (This was in your original code.)
  cnts_summary <- grpstats(cnts_combined, sampleID, statscols = c("mean","median"))
  
  # Name the cluster vector by sample to keep track
  names(ploidy_vals) <- names(clusters) <- sampleID
  rownames(cnts_combined) = paste0(rownames(cnts_combined),"_",sampleID)
  
  return(list(
    clusters    = clusters,        # numeric cluster label per row in cnts_combined
    cnts        = cnts_summary,    # aggregated means/medians
    CN      = cnts_combined,      # the hierarchical clustering object
    distanceFun = dFun,            # for reference if needed
    origins     = origins,          # keep track of which origins we used
    ploidy_vals = ploidy_vals
  ))
}


# weighted Manhattan distance for an *entire* matrix
chrWeightedDist <- function(mat, chrwhole) {
  # vector of chromosome weights
  w <- chrwhole[paste0("chr", 1:22), 1]
  
  # Check if matrix has been transposed by heatmap.2
  if (ncol(mat) == length(w)) {
    # It's the original matrix (e.g., 86 samples x 22 chromosomes)
    # Apply weights to COLUMNS (MARGIN = 2)
    mat.w <- sweep(mat, 2, w, `*`)
    
  } else if (nrow(mat) == length(w)) {
    # It's the transposed matrix (e.g., 22 chromosomes x 86 samples)
    # Apply weights to ROWS (MARGIN = 1)
    mat.w <- sweep(mat, 1, w, `*`)
    
  } else {
    # Stop if dimensions don't match, as a safeguard
    stop("Matrix dimensions do not match the length of the weight vector.")
  }
  
  # Calculate the final weighted distance
  dist(mat.w, method = "manhattan") / sum(w)
}

findBestClustering<-function(allKaryo, numClusters=NULL, hFun=function(x) hclust(x, method="ward.D2"), dFun = function(x) dist(x, method="manhattan")){
  library(cluster)  
  hm=heatmap.2(allKaryo, hclustfun=hFun,distfun=dFun)
  if(!is.null(numClusters)){
    k=numClusters
  }else{
    silhouettecoefs=rep(NA,nrow(allKaryo))
    for(k in 2:(nrow(allKaryo)-1)){
      clusters=cutree(as.hclust(hm$rowDendrogram), k=k)
      sil <- summary(silhouette(clusters, dist(allKaryo)))
      silhouettecoefs[k]= sil$si.summary["Median"]
    }
    k = which.max(silhouettecoefs)
  }
  clusters=cutree(as.hclust(hm$rowDendrogram), k=k)
  return(clusters)
}

getKaryo<-function(cn,ploidy){
  ## set copy number of chromosome to copy number of largest segment for that chromosome
  segments= sapply(sapply(strsplit(colnames(cn),":"),"[[",2), function(x) strsplit(x[[1]],"-")[[1]],simplify = F)
  segments= as.data.frame(do.call(rbind,sapply(segments, as.numeric,simplify = F)))
  rownames(segments) = colnames(cn)
  colnames(segments) = c("start","end")
  segments$length=1+segments$end-segments$start
  segments$chr = as.numeric(sapply(strsplit(colnames(cn),":"),"[[",1))
  chrsegments=sapply(unique(segments$chr), function(x) segments[segments$chr==x,,drop=FALSE],simplify = F)
  chrsegments=sapply(chrsegments, function(x) x[which.max(x$length),,drop=F],simplify = F)
  chrsegments = do.call(rbind,chrsegments)
  cn=cn[,rownames(chrsegments)]
  colnames(cn)=chrsegments$chr
  
  ## all other chromosomes have copy number equal to ploidy for all cells
  otherchr = setdiff(1:22,colnames(cn))
  cn_ = matrix(ploidy,nrow(cn),length(otherchr))
  colnames(cn_)=otherchr
  cn = cbind(cn,cn_)
  gplots::heatmap.2(cn,trace='n',symbreaks = F,symkey=F)
  
  ## karyotype frequency across timepoints
  cn=round(cn)
  # cn[,apply(cn==0,2,all)]=1
  karyo=apply(cn,1,paste0,collapse=".");
  names(karyo) = rownames(cn)
  karyo_in= plyr::count(karyo)
  rownames(karyo_in)=karyo_in$x
  return(list(karyo=karyo_in[,'freq',drop=F], cn=cn ))
  
}

calcPloidy<-function(cnts, chrwhole){
  ii=paste0('chr',colnames(cnts));
  ploidy=apply(cnts,1, function(x) sum(chrwhole[ii,]*x)/sum(chrwhole[ii,]))
  return(ploidy)
}


library(data.table)
library(cluster)
library(gplots)
cloneids <- c("SUM-159_NLS_4N_A5M_K_harvest", "SUM159_NLS_4N_O2_A7K_harvest","SUM159_NLS_4N_O1_A7K_harvest","SUM-159_NLS_4N_O2_A17_seedT1", "SUM-159_NLS_4N_O1_A17_seedT1")
out=clusterKaryotypes(cloneids,whichP = "GenomePerspective", depth=1, path2lanscape=NULL, numClusters=4,method = "ward.D2")
