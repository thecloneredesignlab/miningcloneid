library(DBI)
library(RPostgres)
library(dplyr)
library(ggplot2)
library(shiny)
###############################################################
###############################################################

baseDir='data/DetectionResults'; 
source("code/Utils.R")
mydb = cloneid::connect2DB()
q <- "SELECT id, comment, event, passaged_from_id1,cellLine, correctedCount,passage, date,lastModified,owner,areaOccupied_um2, cellSize_um2 from Passaging"
pass <- dbGetQuery(mydb,q)
rownames(pass) <- pass$id
pass$date <- as.Date(pass$date)
pass=pass[!grepl("images unavailable", pass$comment),]
pass=pass[grepl("LTEE", pass$comment),]


dir.create("data/LTEEs")
the_target_area <- 1E9
# the_target_area <- 3E9
for (cellLine in c("SUM-159")){ #"SNU-668","HGC-27","MDAMB231",
  print(cellLine)
  ploidies = c("")
  if(cellLine=="SUM-159"){
    ploidies=c("_2N","_4N")
  }
  for(ploidy in ploidies){
    ## compare  cell size accounting for confluence
    origin = pass[pass$cellLine == cellLine,'id']
    origin = origin[which(grepl(ploidy,origin))]
    # first assign lineage
    lineage = sapply(strsplit(origin,"_"), function(x) paste0(x[1:2],collapse = "_"))
    ## For SUM-159
    ii=which(grepl("_NLS_", origin))
    lineage[ii] = sapply(strsplit(origin[ii],"_"), function(x) paste0(x[1:4],collapse = "_"))
    lineage=gsub("SUM","",gsub("SUM-","",gsub("_A5M","_A",gsub("_A7M","_A",lineage))))
    # lineage=gsub("MB","",gsub("MB-","",gsub("MDA","",gsub("MDA-","",lineage))))
    pass$lineage=NA
    fr=plyr::count(lineage)
    fr=fr[fr$freq>4,]
    for (l in fr$x) {
      ii <- grep(l, pass$id, fixed = TRUE)         
      if(endsWith(l,"_A") || endsWith(l,"_RPMI")){ ## these are controls
        pass$lineage[ii] <- l
        next
      }
      p  <- rank(pass$date[ii], ties.method = "average") / length(ii)
      # Map to early/middle/late by tertiles
      # phase <- cut(    p, breaks = c(0, 1/3, 1/2, 2/3, 1), include.lowest = TRUE, labels = c("early","mid1","mid2", "late"),    right = TRUE  )
      phase <- cut(    p, breaks = c(0, 1/2, 1), include.lowest = TRUE, labels = c("early", "late"),    right = TRUE  )
      # if(endsWith(l,"MDA-MB-231_P")){ ## 4N O2 evolved later
      #   phase <- cut(    p, breaks = c(0, 2/3, 1), include.lowest = TRUE, labels = c("early", "late"),    right = TRUE  )
      # }
      pass$lineage[ii] <- paste0(l, "_", phase)
    }
    ## get size for a given confluence for each lineage
    res <- plot_size_vs_area_by_lineage(pass, event = "harvest")
    res$plot
    rownames(res$stats)=res$stats$lineage
    # plot_model_comparison(res$models$`159-NLS-2N_O2_early`,res$models$`159-NLS-4N_O2_early`,"2N","4N")
    target_area <- the_target_area
    
    ## plot average size over time for a given lineage
    for(l in fr$x){
      target_area <- the_target_area
      if(l %in% c("159_NLS_4N_C","159_NLS_2N_C")){
        target_area <- 7E9; ## @TODO: this should not be necessary, control not passaged at consistent densities as test
      }
      closest_rows <- pass %>%
        filter(
          !is.na(lineage),
          abs(areaOccupied_um2 - target_area) < 0.75E9        # keep only near matches
        ) %>%
        group_by(passaged_from_id1) %>%
        slice_min(abs(areaOccupied_um2 - target_area), with_ties = FALSE) %>%
        ungroup()
      closest_rows$id=gsub("-NLS-","_NLS_",closest_rows$id)
      closest_rows$lineage=closest_rows$passaged_from_id1
      closest_rows = closest_rows[grep(l,closest_rows$id),]
      closest_rows = closest_rows["SUM-159_NLS_2N_O2_A3_seedT1"!=closest_rows$id,]
      closest_rows=closest_rows[order(closest_rows$passage),]
      csv_list=readCellFeatures(closest_rows,select_by_global_density = T, bin_width = 25, base_dir=baseDir)
      names(csv_list) = closest_rows$lineage
      save('csv_list', file=paste0("data/LTEEs/",l,ploidy,".RObj"))
    }
    
  }
}


# # ###############
# # #### UMAPs ####
# library(stringr)
# library(matlab)
# f=list.files("~/Downloads/LTEEs/",pattern = "RObj", full.names = T)
# f=grep("_",f,value = T)
# dir.create("~/Downloads/LTEEs/umaps")
# for(x in f){
#   load(x); ## loads csv_list
#   if(length(csv_list)<3){
#     next; ## likely not an LTEE
#   }
#   print(x)
#   print(names(csv_list))
#   # --- 4. Example Usage ---
#   # Call the function with your list of data frames
#   # Replace 'csv_list' with the name of your actual list object
#   my_umap_plot <- plot_projection_by_passage(csv_list,method = "umap", n_neighbors=90, title=fileparts(x)$name)
#   if(grepl("159_NLS_2N_O1_2N",x)){
#     # clustering_results <- interactive_cluster_projection(csv_list,method = "umap", n_neighbors=90, title=fileparts(x)$name)
#     clustered_data <- interactive_cluster_projection(csv_list,method = "umap", n_neighbors=10, title=fileparts(x)$name)
#     distribution_data <- clustered_data %>%
#       # Ignore any points you didn't assign
#       filter(cluster != "Unassigned") %>%
#       count(cluster, passage_number, name = "n_cells") %>%
#       group_by(passage_number) %>%
#       mutate(proportion = n_cells / sum(n_cells)) %>%
#       ungroup()
# 
#     # Create the final plot
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
#   # Print the final plot
#   ggsave(paste0("~/Downloads/LTEEs/umaps/",fileparts(x)$name,".png"), plot = my_umap_plot,width = 7.52, height = 7.52)
# }