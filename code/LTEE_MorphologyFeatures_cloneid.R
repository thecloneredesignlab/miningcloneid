# LTEE_MorphologyFeatures_cloneid.R
# This script is designed to process and analyze cell morphology data from a 
# Long-Term Evolution Experiment (LTEE). It begins by querying a database for 
# passaging metadata and systematically assigns lineage and experimental phase 
# labels (e.g., "early", "late") based on the sample's ID and date. The core 
# function of the script is to loop through these lineages, identify representative 
# data points at a standardized cell confluence, and then use the readCellFeatures 
# function to load the detailed morphological data for those points. Finally, it 
# saves these curated datasets as .RObj files, preparing them for downstream 
# analysis such as generating UMAP plots to visualize how cell morphology evolves over time.


library(DBI)
library(RPostgres)
library(dplyr)
library(ggplot2)
library(shiny)
###############################################################
###############################################################

# Define the base directory where cell feature CSV files are located.
baseDir='../data/S3Buckets/CellSegmentation_InVitro/output/DetectionResults'; 

# Load the helper functions from the Utils.R script.
source("Utils.R")

# --- 1. Database Query and Initial Filtering ---
# Connect to the database.
mydb = cloneid::connect2DB()
# Define and execute a SQL query to fetch passaging metadata.
q <- "SELECT id, comment, event, passaged_from_id1,cellLine, correctedCount,passage, date,lastModified,owner,areaOccupied_um2, cellSize_um2 from Passaging"
pass <- dbGetQuery(mydb,q)
# Set the row names of the dataframe to the unique passaging 'id'.
rownames(pass) <- pass$id
# Convert the date column to a proper Date object.
pass$date <- as.Date(pass$date)
# Filter out entries with comments indicating that images were not available.
pass=pass[!grepl("images unavailable", pass$comment),]
# Filter to keep only entries related to the Long-Term Evolution Experiment (LTEE).
pass=pass[grepl("LTEE", pass$comment),]


# --- 2. Main Analysis Loop ---
# Clean up and create a directory for output files.
unlink("~/Downloads/LTEEs", force = F, recursive=T)
dir.create("~/Downloads/LTEEs")
# Define a target confluence (area occupied) to standardize cell size comparisons.
the_target_area <- 1E9

# Loop through each cell line to be analyzed.
for (cellLine in c("SUM-159")){ 
  print(cellLine)
  # Define ploidy variations for the given cell line.
  ploidies = c("")
  if(cellLine=="SUM-159"){
    ploidies=c("_2N","_4N")
  }
  # Loop through each ploidy type.
  for(ploidy in ploidies){
    ## --- 2a. Lineage and Phase Assignment ---
    # Get all IDs for the current cell line and ploidy.
    origin = pass[pass$cellLine == cellLine,'id']
    origin = origin[which(grepl(ploidy,origin))]
    
    # Define a 'lineage' based on parsing the sample ID string. This is a critical
    # step that groups related experiments together.
    lineage = sapply(strsplit(origin,"_"), function(x) paste0(x[1:2],collapse = "_"))
    ## Special parsing rule for SUM-159 NLS samples.
    ii=which(grepl("_NLS_", origin))
    lineage[ii] = sapply(strsplit(origin[ii],"_"), function(x) paste0(x[1:4],collapse = "_"))
    # Further clean up the lineage names using a series of substitutions.
    lineage=gsub("SUM","",gsub("SUM-","",gsub("_A5M","_A",gsub("_A7M","_A",lineage))))
    
    # Assign the generated lineages back to the main 'pass' dataframe.
    pass$lineage=NA
    # Count the number of entries per lineage to focus on those with enough data.
    fr=plyr::count(lineage)
    fr=fr[fr$freq>4,] # Keep lineages with more than 4 data points.
    
    # For each well-represented lineage, assign an experimental phase ('early' or 'late').
    for (l in fr$x) {
      ii <- grep(l, pass$id, fixed = TRUE)         
      # Control lineages are not assigned a phase.
      if(endsWith(l,"_A") || endsWith(l,"_RPMI")){ 
        pass$lineage[ii] <- l
        next
      }
      # Calculate the temporal rank of each sample within its lineage based on date.
      p  <- rank(pass$date[ii], ties.method = "average") / length(ii)
      # Divide the experiment into phases based on this temporal rank. Here, the
      # first half of experiments are 'early' and the second half are 'late'.
      phase <- cut(p, breaks = c(0, 1/2, 1), include.lowest = TRUE, labels = c("early", "late"), right = TRUE)
      # Combine the lineage name and phase to create the final lineage label.
      pass$lineage[ii] <- paste0(l, "_", phase)
    }
    
    ## --- 2b. Modeling and Feature Extraction ---
    # Use the utility function to plot cell size vs. confluence for all defined lineages.
    res <- plot_size_vs_area_by_lineage(pass, event = "harvest")
    res$plot # Display the plot.
    rownames(res$stats)=res$stats$lineage
    
    # For each lineage, find the data point closest to the target confluence
    # and read its detailed cell feature data.
    for(l in fr$x){
      target_area <- the_target_area
      # Use a different target area for control lineages if necessary.
      if(l %in% c("159_NLS_4N_C","159_NLS_2N_C")){
        target_area <- 7E9; 
      }
      
      # For each experiment (`passaged_from_id1`), find the single measurement
      # whose `areaOccupied_um2` is closest to the `target_area`. This
      # normalizes the comparison by selecting cells at a similar confluence.
      closest_rows <- pass %>%
        filter(
          !is.na(lineage),
          abs(areaOccupied_um2 - target_area) < 0.75E9 # Pre-filter to a reasonable window.
        ) %>%
        group_by(passaged_from_id1) %>%
        # `slice_min` selects the row with the minimum value for the given expression.
        slice_min(abs(areaOccupied_um2 - target_area), with_ties = FALSE) %>%
        ungroup()
      
      # Further data cleaning on the selected rows.
      closest_rows$id=gsub("-NLS-","_NLS_",closest_rows$id)
      closest_rows$lineage=closest_rows$passaged_from_id1
      closest_rows = closest_rows[grep(l,closest_rows$id),]
      closest_rows = closest_rows["SUM-159_NLS_2N_O2_A3_seedT1"!=closest_rows$id,]
      closest_rows=closest_rows[order(closest_rows$passage),]
      
      # Read the detailed cell feature CSVs for these selected, representative wells.
      # `select_by_global_density = T` triggers the smart quartile selection.
      csv_list=readCellFeatures(closest_rows,select_by_global_density = T, bin_width = 25, base_dir=baseDir)
      names(csv_list) = closest_rows$lineage
      
      # Save the resulting list of dataframes to a file for later use.
      save('csv_list', file=paste0("~/Downloads/LTEEs/",l,ploidy,".RObj"))
    }
    
  }
}


# # ###############
# # #### UMAPs #### (This entire section is commented out in the original script)
# # The code below would be used to load the saved .RObj files and generate UMAP plots.
# 
# library(stringr)
# library(matlab)
# # List all the saved RObj files.
# f=list.files("~/Downloads/LTEEs/",pattern = "RObj", full.names = T)
# f=grep("_",f,value = T)
# dir.create("~/Downloads/LTEEs/umaps")
# # Loop through each file.
# for(x in f){
#   load(x); ## loads the 'csv_list' object.
#   if(length(csv_list)<3){
#     next; ## Skip if there are too few passages to visualize.
#   }
#   print(x)
#   print(names(csv_liFst))
#   
#   # Call the projection plotting function from Utils.R.
#   my_umap_plot <- plot_projection_by_passage(csv_list,method = "umap", n_neighbors=90, title=fileparts(x)$name)
#   
#   # Special case: for one specific dataset, launch the interactive clustering tool.
#   if(grepl("159_NLS_2N_O1_2N",x)){
#     # This would launch the Shiny app and wait for the user to finish clustering.
#     clustered_data <- interactive_cluster_projection(csv_list,method = "umap", n_neighbors=10, title=fileparts(x)$name)
#     
#     # After closing the app, process the returned data with cluster labels.
#     # Calculate the proportion of each manual cluster at each passage number.
#     distribution_data <- clustered_data %>%
#       filter(cluster != "Unassigned") %>%
#       count(cluster, passage_number, name = "n_cells") %>%
#       group_by(passage_number) %>%
#       mutate(proportion = n_cells / sum(n_cells)) %>%
#       ungroup()
# 
#     # Create a line plot showing how the cluster proportions change over time (passages).
#     final_plot <- ggplot(distribution_data, aes(x = passage_number, y = proportion, color = cluster, group = cluster)) +
#       geom_line(linewidth = 1.2) +
#       geom_point(size = 4) +
#       labs(
#         title = "Change in Manually-Defined Cluster Proportions",
#         x = "Passage Number",
#         y = "Proportion of Cells in Passage",
#         color = "Cluster"
#       ) +
#       theme_bw(base_size = 14) +
#       scale_x_continuous(breaks = unique(distribution_data$passage_number))
# 
#     print(final_plot)
#   }
#   # Save the standard UMAP plot to a PNG file.
#   ggsave(paste0("~/Downloads/LTEEs/umaps/",fileparts(x)$name,".png"), plot = my_umap_plot,width = 7.52, height = 7.52)
# }