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


calcPloidy<-function(cnts, chrwhole){
  ii=paste0('chr',colnames(cnts));
  ploidy=apply(cnts,1, function(x) sum(chrwhole[ii,]*x)/sum(chrwhole[ii,]))
  return(ploidy)
}


library(data.table)
library(cluster)
library(RColorBrewer)
library(gplots)
source("code/Utils.R")

cloneids <- c("SUM-159_NLS_4N_A5M_K_harvest", "SUM159_NLS_4N_O2_A7K_harvest","SUM159_NLS_4N_O1_A7K_harvest","SUM-159_NLS_4N_O2_A17_seedT1", "SUM-159_NLS_4N_O1_A17_seedT1")
out=clusterKaryotypes(cloneids,whichP = "GenomePerspective", depth=1, path2lanscape=NULL, numClusters=4,method = "ward.D2")
