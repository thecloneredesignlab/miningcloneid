# Karyotypes_inVitro_cloneid.R This script is executes a karyotype
# clustering analysis on in vitro cancer cell data. Its main purpose is to 
# identify and group distinct cell subpopulations based on their chr. copy number profiles. 
# It defines custom functions to calculate ploidy and a chromosome-length-weighted Manhattan 
# distance. The script then uses these functions along with the clusterKaryotypes workflow 
# from Utils.R to analyze a predefined list of clone IDs, forcing the outcome into four 
# distinct clusters and generating visualizations like heatmaps to display the results.


#' Calculate a Weighted Manhattan Distance for a Matrix
#'
#' This function computes the Manhattan distance between rows of a matrix, but
#' weights each column's contribution by a corresponding value (e.g., chromosome
#' length). This ensures that changes in larger chromosomes contribute more to the
#' overall distance. It is designed to work correctly even if `heatmap.2`
#' internally transposes the input matrix.
#'
#' @param mat The numeric input matrix (e.g., cells x chromosomes).
#' @param chrwhole A named vector or single-column dataframe of weights, where names
#'   match the columns of `mat`.
#' @return A `dist` object representing the weighted distances between rows.
chrWeightedDist <- function(mat, chrwhole) {
  # Create a simple named vector of weights from the input dataframe.
  w <- chrwhole[paste0("chr", 1:22), 1]
  
  # `heatmap.2` can sometimes pass a transposed matrix to the distance function.
  # This block checks the dimensions to apply the weights correctly.
  
  # Case 1: The matrix has not been transposed (e.g., 86 cells x 22 chromosomes).
  if (ncol(mat) == length(w)) {
    # `sweep` applies the weights `w` to the COLUMNS (MARGIN = 2) of the matrix.
    # Each column is multiplied by its corresponding weight.
    mat.w <- sweep(mat, 2, w, `*`)
    
    # Case 2: The matrix has been transposed (e.g., 22 chromosomes x 86 cells).
  } else if (nrow(mat) == length(w)) {
    # `sweep` applies the weights `w` to the ROWS (MARGIN = 1) of the matrix.
    mat.w <- sweep(mat, 1, w, `*`)
    
  } else {
    # If dimensions don't match, something is wrong. Stop with an error.
    stop("Matrix dimensions do not match the length of the weight vector.")
  }
  
  # Calculate the standard Manhattan distance on the NOW-WEIGHTED matrix.
  # Finally, normalize the entire distance matrix by the sum of all weights.
  # This makes the distance metric independent of the absolute scale of the weights.
  dist(mat.w, method = "manhattan") / sum(w)
}

#' Calculate Ploidy based on Chromosome-Weighted Copy Numbers
#'
#' This function estimates the ploidy for each cell (row) in a copy number matrix.
#' It calculates a weighted average of the copy numbers, where each chromosome's
#' copy number is weighted by its relative size.
#'
#' @param cnts A matrix of copy number data (cells x chromosomes).
#' @param chrwhole A dataframe containing chromosome lengths.
#' @return A numeric vector of calculated ploidy values, one for each cell.
calcPloidy<-function(cnts, chrwhole){
  # Ensure chromosome names match between the copy number matrix and the weights.
  ii=paste0('chr',colnames(cnts));
  # For each cell (row), calculate the sum of (copy_number * chr_length)
  # and divide by the sum of all chromosome lengths. This is the weighted average.
  ploidy=apply(cnts,1, function(x) sum(chrwhole[ii,]*x)/sum(chrwhole[ii,]))
  return(ploidy)
}

# --- Main Script ---

# Load required libraries for data manipulation, clustering, and plotting.
library(data.table)
library(cluster)
library(RColorBrewer)
library(gplots)
# Source the utility functions defined in Utils.R, including `clusterKaryotypes`.
source("Utils.R")

# Define the set of clone IDs from the database that will be analyzed.
cloneids <- c("SUM-159_NLS_4N_A5M_K_harvest", "SUM159_NLS_4N_O2_A7K_harvest","SUM159_NLS_4N_O1_A7K_harvest","SUM-159_NLS_4N_O2_A17_seedT1", "SUM-159_NLS_4N_O1_A17_seedT1")

# Call the main analysis function with the specified clone IDs and parameters.
# - whichP: Specifies the data 'perspective' in the database.
# - depth: Controls how many related clones are fetched.
# - numClusters: Forces the clustering into 4 groups.
# - method: Uses the "ward.D2" algorithm for hierarchical clustering.
out=clusterKaryotypes(cloneids,whichP = "GenomePerspective", depth=1, path2lanscape=NULL, numClusters=4,method = "ward.D2")