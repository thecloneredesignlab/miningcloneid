library(DBI)
library(RPostgres)
library(dplyr)
library(ggplot2)
library(shiny)
plot_size_vs_area_by_lineage <- function(pass, lineage_vec = NULL, event = "harvest", free_scales = TRUE) {
  # Subset
  keep <- pass$event == event
  if (!is.null(lineage_vec)) {
    keep <- keep & pass$lineage %in% unique(lineage_vec)
  }
  df <- pass[keep, , drop = FALSE] %>%
    filter(!is.na(cellSize_um2), !is.na(areaOccupied_um2))
  
  # Build per-lineage models: cellSize ~ areaOccupied
  split_df <- split(df, df$lineage)
  models <- lapply(split_df, function(d) {
    # Need at least 2 points and variation in predictor
    if (nrow(d) >= 2 && length(unique(d$areaOccupied_um2)) > 1) {
      lm(cellSize_um2 ~ areaOccupied_um2, data = d)
    } else {
      NULL
    }
  })
  
  # Stats for annotation (slope & R^2) for lineages with valid models
  valid_lineages <- names(models)[!vapply(models, is.null, logical(1))]
  stats_df <- do.call(
    rbind,
    lapply(valid_lineages, function(ln) {
      m <- models[[ln]]
      sm <- summary(m)
      data.frame(
        lineage = ln,
        slope = unname(coef(m)[2]),
        r2 = sm$r.squared,
        n = length(m$residuals),
        label = sprintf("slope = %.3f\nR² = %.3f", unname(coef(m)[2]), sm$r.squared),
        x = Inf,
        y = Inf,
        stringsAsFactors = FALSE
      )
    })
  )
  
  # Faceted plot
  p <- ggplot(df, aes(x = areaOccupied_um2, y = cellSize_um2)) +
    geom_point(alpha = 0.8, aes(color = lineage)) +
    geom_smooth(method = "lm", se = FALSE, aes(color = lineage)) +
    facet_wrap(~ lineage, scales = if (free_scales) "free" else "fixed") +
    geom_text(data = stats_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = 1.2) +
    labs(
      x = "Area occupied (µm²)",
      y = "Cell size (µm²)",
      color = "Lineage"
    ) +
    theme_bw()
  
  invisible(list(plot = p, models = models, stats = stats_df))
}

#' Generate a UMAP plot from a named list of data frames
#'
#' @param data_list A named list where each element is a data frame of cell metrics.
#'   The names of the list must contain a passage number formatted as 'A#'.
#' @param seed An integer for reproducibility of the UMAP algorithm.
#' @return A ggplot object representing the UMAP projection.
# Required Libraries
# Add 'dbscan' for the new clustering feature
# install.packages(c("uwot", "Rtsne", "dbscan"))
plot_projection_by_passage <- function(data_list, 
                                       method = "umap", 
                                       seed = 42, 
                                       n_neighbors = 15, 
                                       title = "Projection by Passage Number",
                                       cluster = FALSE,      # New: Toggle for clustering
                                       dbscan_eps = 0.5,
                                       dbscan_minPts = 5) {    # New: DBSCAN 'eps' parameter
  
  # --- Robust Data Preparation ---
  common_cols <- Reduce(intersect, lapply(data_list, colnames))
  cols_without_X <- grep("X", common_cols, value = TRUE, invert = TRUE)
  final_cols <- grep("Y", cols_without_X, value = TRUE, invert = TRUE)
  subsetted_list <- lapply(data_list, function(df) df[, final_cols, drop = FALSE])
  
  all_data <- bind_rows(subsetted_list, .id = "source_id") %>%
    mutate(
      passage_number = as.numeric(str_extract(source_id, "(?<=A)[0-9]+"))
    )
  
  # --- Isolate and Scale Numeric Data for Projection ---
  numeric_data <- all_data %>%
    select(where(is.numeric), -passage_number) %>%
    mutate(across(everything(), log1p)) %>% 
    scale()
  
  # --- Perform Dimensionality Reduction ---
  method <- tolower(method)
  if (method == "umap") {
    set.seed(seed)
    projection_results <- uwot::umap(numeric_data, n_neighbors = n_neighbors)
  } else if (method == "pca") {
    pca_results <- prcomp(numeric_data)
    projection_results <- pca_results$x[, 1:2]
  } else if (method == "tsne") {
    set.seed(seed)
    tsne_results <- Rtsne::Rtsne(numeric_data, dims = 2, check_duplicates = FALSE)
    projection_results <- tsne_results$Y
  } else {
    stop("Error: Method must be 'umap', 'pca', or 'tsne'.")
  }
  
  # --- Prepare Data for Plotting ---
  plot_data <- as.data.frame(projection_results)
  colnames(plot_data) <- c("Dim1", "Dim2")
  plot_data$passage_number <- all_data$passage_number
  
  # --- (Optional) Clustering and Final Output Generation ---
  if (cluster) {
    # Check if dbscan package is installed
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      stop("Package 'dbscan' is required. Please run: install.packages('dbscan')", call. = FALSE)
    }
    
    # Perform DBSCAN clustering on the UMAP coordinates.
    # DBSCAN is great for discovering arbitrarily shaped clusters.
    # The 'eps' parameter is crucial and may require tuning.
    set.seed(seed)
    db_clusters <- dbscan::dbscan(plot_data[, c("Dim1", "Dim2")], eps = dbscan_eps,minPts = dbscan_minPts)
    plot_data$cluster <- as.factor(db_clusters$cluster) # Add cluster assignments
    
    # --- Shuffle Data for Plotting ---
    set.seed(seed)
    plot_data <- plot_data[sample(nrow(plot_data)), ]
    
    # --- Visualize Clustered Results ---
    # Points are colored by cluster ID. Cluster '0' represents noise points.
    projection_plot <- ggplot(plot_data, aes(x = Dim1, y = Dim2, color = cluster)) +
      geom_point(alpha = 0.8, size = 2) +
      scale_color_brewer(palette = "Paired", name = "Cluster") +
      guides(color = guide_legend(override.aes = list(size = 5))) + # Make legend points larger
      labs(
        title = paste(toupper(method), title, "(Clustered)"),
        x = paste(toupper(method), "Dimension 1"),
        y = paste(toupper(method), "Dimension 2")
      ) +
      theme_minimal()
    
    # --- Calculate Passage Distribution within Clusters ---
    # This dataframe shows the count and proportion of cells from each passage
    # that fall into each cluster, enabling downstream analysis of cluster evolution.
    cluster_distribution <- plot_data %>%
      count(cluster, passage_number, name = "n_cells") %>%
      group_by(passage_number) %>%
      mutate(proportion = n_cells / sum(n_cells)) %>%
      ungroup() %>%
      arrange(cluster, passage_number)
    
    # Return a list containing both the plot and the distribution data
    return(list(plot = projection_plot, distribution_data = cluster_distribution))
    
  } else {
    # --- Shuffle Data for Plotting (Original Path) ---
    set.seed(seed)
    plot_data <- plot_data[sample(nrow(plot_data)), ]
    
    # --- Visualize Results (Original Path) ---
    projection_plot <- ggplot(plot_data, aes(x = Dim1, y = Dim2, color = passage_number)) +
      geom_point(alpha = 0.8, size = 2) +
      scale_color_viridis_c() +
      labs(
        title = paste(toupper(method), title),
        x = paste(toupper(method), "Dimension 1"),
        y = paste(toupper(method), "Dimension 2"),
        color = "Passage"
      ) +
      theme_minimal()
    
    return(projection_plot)
  }
}

interactive_cluster_projection <- function(data_list, 
                                           method = "umap", 
                                           seed = 42, 
                                           n_neighbors = 15,
                                           title = "Interactive Projection") {
  
  # --- [1] Data Preparation and Dimensionality Reduction (same as before) ---
  common_cols <- Reduce(intersect, lapply(data_list, colnames))
  cols_without_X <- grep("X", common_cols, value = TRUE, invert = TRUE)
  final_cols <- grep("Y", cols_without_X, value = TRUE, invert = TRUE)
  subsetted_list <- lapply(data_list, function(df) df[, final_cols, drop = FALSE])
  
  all_data <- bind_rows(subsetted_list, .id = "source_id") %>%
    mutate(
      passage_number = as.numeric(str_extract(source_id, "(?<=A)[0-9]+"))
    )
  
  numeric_data <- all_data %>%
    select(where(is.numeric), -passage_number) %>%
    mutate(across(everything(), log1p)) %>% 
    scale()
  
  method <- tolower(method)
  if (method == "umap") {
    set.seed(seed)
    projection_results <- uwot::umap(numeric_data, n_neighbors = n_neighbors)
  } else {
    # For simplicity, this example focuses on UMAP, but others could be added back.
    stop("This interactive version is currently configured for UMAP only.")
  }
  
  plot_data <- as.data.frame(projection_results)
  colnames(plot_data) <- c("Dim1", "Dim2")
  plot_data$passage_number <- all_data$passage_number
  # Initialize the cluster column
  plot_data$cluster <- "Unassigned"
  
  # --- [2] Define the Interactive Shiny Application ---
  
  # UI Definition
  ui <- fluidPage(
    titlePanel("Interactive UMAP Clustering"),
    sidebarLayout(
      sidebarPanel(
        h4("Instructions"),
        p("1. Click and drag on the plot to select points."),
        p("2. Enter a name for your selection below."),
        p("3. Click 'Assign' to save the cluster."),
        p("4. Repeat for all desired clusters."),
        hr(),
        textInput("cluster_name", "Cluster Name:", "Cluster_1"),
        actionButton("assign_btn", "Assign Selected to Cluster", icon = icon("paint-brush")),
        actionButton("reset_btn", "Reset All Assignments", icon = icon("undo")),
        hr(),
        p("Click 'Done' when finished to return the data to your R session."),
        actionButton("done_btn", "Done", icon = icon("check"), class = "btn-success")
      ),
      mainPanel(
        # The 'brush' argument is the key to enabling selection
        plotOutput("umap_plot", height = "600px", brush = brushOpts(id = "plot_brush")),
        h4("Cluster Distribution Summary"),
        DT::dataTableOutput("summary_table")
      )
    )
  )
  
  # Server Logic
  server <- function(input, output, session) {
    
    # Reactive value to store the data with cluster assignments
    rv <- reactiveValues(data = plot_data, brush_selection = NULL)
    
    # Render the main UMAP plot
    output$umap_plot <- renderPlot({
      
      # Determine points currently under the brush for highlighting
      brushed_pts <- brushedPoints(rv$data, input$plot_brush)
      
      ggplot() +
        # Plot unassigned points in grey
        geom_point(data = filter(rv$data, cluster == "Unassigned"), aes(x = Dim1, y = Dim2), color = "grey80", alpha = 0.7) +
        # Plot assigned points with color
        geom_point(data = filter(rv$data, cluster != "Unassigned"), aes(x = Dim1, y = Dim2, color = cluster), alpha = 0.9, size = 2) +
        # Highlight points currently being selected in red
        geom_point(data = brushed_pts, aes(x = Dim1, y = Dim2), color = "#d9534f", size = 2.5) +
        coord_cartesian(xlim = range(rv$data$Dim1), ylim = range(rv$data$Dim2)) +
        scale_color_brewer(palette = "Set1", name = "Cluster") +
        labs(
          title = title,
          x = "UMAP Dimension 1",
          y = "UMAP Dimension 2"
        ) +
        theme_minimal(base_size = 14)
    })
    
    # Handle the "Assign" button click
    observeEvent(input$assign_btn, {
      req(input$plot_brush, input$cluster_name)
      
      selected_rows <- brushedPoints(rv$data, input$plot_brush, allRows = TRUE)$selected_
      
      # Update the cluster label for the selected rows
      rv$data$cluster[selected_rows] <- input$cluster_name
      
      # Increment the default cluster name (e.g., from Cluster_1 to Cluster_2)
      new_num <- as.numeric(sub(".*_", "", input$cluster_name)) + 1
      if (!is.na(new_num)) {
        updateTextInput(session, "cluster_name", value = paste0("Cluster_", new_num))
      }
      
      # Clear the brush selection from the plot
      session$resetBrush("plot_brush")
    })
    
    # Handle the "Reset" button click
    observeEvent(input$reset_btn, {
      rv$data$cluster <- "Unassigned"
    })
    
    # Render the summary table
    output$summary_table <- DT::renderDataTable({
      req(rv$data)
      rv$data %>%
        filter(cluster != "Unassigned") %>%
        count(cluster, passage_number, name = "cell_count") %>%
        arrange(cluster, passage_number)
    }, options = list(pageLength = 5))
    
    # Handle the "Done" button click
    observeEvent(input$done_btn, {
      # When "Done" is clicked, stop the app and return the final data frame
      stopApp(returnValue = rv$data)
    })
  }
  
  # --- [3] Run the Application ---
  # This launches the app in a viewer or new window.
  # The function will pause here until the user clicks "Done".
  cat("Launching interactive session...\n")
  return(shiny::runApp(shinyApp(ui, server)))
}

# --- Helper Functions (to be replaced with your actual lab functions) ---

#' A placeholder function to simulate reading cell feature CSVs.
#' In a real scenario, this would read actual data from a file system.
#'
#' @param metadata_rows A dataframe of selected passaging events.
#' @param base_dir The base directory where feature files are stored.
#' @return A named list of dataframes, each containing cell features.
readCellFeatures <- function(closest_rows,
                             base_dir = "~/CellSegmentations/output/DetectionResults",
                             select_by_global_density = FALSE,
                             area_per_quartile_um2 = NULL,   # if provided, density = count / area
                             bin_width = 25,                  # binning step for mode estimation
                             tie_order = c("tl","tr","bl","br")) {
  
  # deps
  require(purrr)
  require(dplyr)
  
  dir.create("~/Downloads/tmpLTE", showWarnings = FALSE)
  
  suffixes <- c("tl","tr","bl","br")
  
  # Build file paths per id
  file_groups <- purrr::map(closest_rows$id, ~{
    paste0(file.path(base_dir, .x), "_10x_ph_", suffixes, ".csv")
  })
  names(file_groups) <- closest_rows$id
  
  # Read and combine (add suffix)
  csv_list <- purrr::map(file_groups, function(files) {
    dfs <- purrr::map(files, ~{
      if (file.exists(.x)) {
        df <- read.csv(.x, sep = "\t")
        if(ncol(df)<3){
          df <- read.csv(.x)
        }
        df$suffix <- sub(".*_ph_(.*)\\.csv$", "\\1", basename(.x))
        df
      } else {
        warning("File not found: ", .x)
        NULL
      }
    })
    dplyr::bind_rows(dfs)
  })
  
  names(csv_list) <- closest_rows$id
  
  # Early return if not selecting by global density
  if (!isTRUE(select_by_global_density)) {
    return(csv_list)
  }
  
  # ---- Compute densities per entry (per suffix) ----
  # density = count or count / area_per_quartile_um2 (if provided)
  dens_list <- purrr::map(csv_list, function(df) {
    if (is.null(df) || nrow(df) == 0) return(setNames(numeric(0), character(0)))
    counts <- table(df$suffix)
    counts <- counts[intersect(names(counts), suffixes)]
    counts <- counts[order(match(names(counts), suffixes))] # tl,tr,bl,br order
    
    if (!is.null(area_per_quartile_um2) && is.finite(area_per_quartile_um2) && area_per_quartile_um2 > 0) {
      dens <- as.numeric(counts) / area_per_quartile_um2
    } else {
      dens <- as.numeric(counts)  # simple counts
    }
    names(dens) <- names(counts)
    dens
  })
  
  # Flatten all densities for global mode estimation
  all_dens <- unlist(dens_list, use.names = FALSE)
  all_dens <- all_dens[is.finite(all_dens)]
  
  if (length(all_dens) == 0) {
    warning("No densities found to compute a global target; returning original list.")
    return(csv_list)
  }
  
  # ---- Estimate the most common density x (binned mode) ----
  # Bin the densities to stabilize the mode calculation
  # If densities are counts, bin_width works in "cells"; if normalized, in "cells/um^2"
  bins <- floor(all_dens / bin_width) * bin_width
  # mode of binned values
  bin_tab <- sort(table(bins), decreasing = TRUE)
  top_bins <- names(bin_tab)[bin_tab == max(bin_tab)]
  # If multiple top bins, pick the one whose center is closest to the median density
  bin_centers <- as.numeric(top_bins) + bin_width/2
  x <- bin_centers[which.min(abs(bin_centers - median(all_dens)))]
  
  # ---- For each entry: keep only the quartile closest to x ----
  chosen_suffix <- list()
  csv_list_filtered <- purrr::imap(csv_list, function(df, entry_name) {
    dens <- dens_list[[entry_name]]
    if (length(dens) == 0 || is.null(df) || nrow(df) == 0) return(df)
    
    # distances to global target
    dists <- abs(dens - x)
    
    # choose suffix with minimum distance; break ties by preferred order
    best <- names(dists)[which(dists == min(dists))]
    best <- best[order(match(best, tie_order))][1]
    
    chosen_suffix[[entry_name]] <<- best
    
    # subset rows
    df %>% dplyr::filter(.data$suffix == best)
  })
  
  # Attach attributes so you can inspect the selection later
  attr(csv_list_filtered, "global_target_density") <- x
  attr(csv_list_filtered, "bin_width") <- bin_width
  attr(csv_list_filtered, "area_per_quartile_um2") <- area_per_quartile_um2
  attr(csv_list_filtered, "chosen_suffix") <- chosen_suffix
  attr(csv_list_filtered, "densities") <- dens_list
  
  csv_list_filtered
}


###############################################################
###############################################################

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


clusterKaryotypes <- function(origins, whichP = "GenomePerspective", depth  = 1,  path2lanscape   = NULL, numClusters = NULL, capping = NULL, method = "complete", chrwhole = NULL){
  if(is.null(chrwhole)){
    devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
    x <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz", 
               col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
    chrarms=x[ , .(length = sum(chromEnd - chromStart)),by = .(chrom, arm = substring(name, 1, 1)) ]
    chrwhole=grpstats(as.matrix(chrarms$length),chrarms$chrom, "sum")$sum
  }
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
