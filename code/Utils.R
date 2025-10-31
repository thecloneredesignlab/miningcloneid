library(DBI)
library(RPostgres)
library(dplyr)
library(ggplot2)
library(shiny)

#' Plot Cell Size vs. Area Occupied, Faceted by Lineage
#'
#' This function generates a scatter plot of cell size against the occupied area
#' for different cell lineages. It fits a linear model for each lineage and
#' annotates the plot with the slope and R-squared value.
#'
#' @param pass A dataframe containing passaging data, which must include columns
#'   'event', 'lineage', 'cellSize_um2', and 'areaOccupied_um2'.
#' @param lineage_vec An optional character vector of lineage names to include.
#'   If NULL, all lineages are used.
#' @param event The type of event to filter for (e.g., "harvest").
#' @param free_scales A logical value passed to `facet_wrap`. If TRUE, each
#'   facet will have its own axis scales.
#' @return An invisible list containing the ggplot object ('plot'), a list of
#'   linear models ('models'), and a dataframe of model statistics ('stats').
plot_size_vs_area_by_lineage <- function(pass, lineage_vec = NULL, event = "harvest", free_scales = TRUE) {
  # --- 1. Data Subsetting ---
  # Filter the dataframe to only include the specified event type.
  keep <- pass$event == event
  # If a vector of specific lineages is provided, further subset the data.
  if (!is.null(lineage_vec)) {
    keep <- keep & pass$lineage %in% unique(lineage_vec)
  }
  # Apply the filter and remove any rows with NA values in the key metric columns.
  df <- pass[keep, , drop = FALSE] %>%
    filter(!is.na(cellSize_um2), !is.na(areaOccupied_um2))
  
  # --- 2. Per-Lineage Modeling ---
  # Split the dataframe into a list of dataframes, one for each lineage.
  split_df <- split(df, df$lineage)
  # Apply a function to each lineage's dataframe to build a linear model.
  models <- lapply(split_df, function(d) {
    # A linear model (y ~ x) requires at least 2 data points and the predictor 'x'
    # must have more than one unique value. This check prevents errors.
    if (nrow(d) >= 2 && length(unique(d$areaOccupied_um2)) > 1) {
      lm(cellSize_um2 ~ areaOccupied_um2, data = d)
    } else {
      # If the conditions aren't met, return NULL for this lineage.
      NULL
    }
  })
  
  # --- 3. Extract Model Statistics for Plot Annotation ---
  # Identify which lineages successfully produced a model.
  valid_lineages <- names(models)[!vapply(models, is.null, logical(1))]
  # For each valid model, extract the slope, R-squared, and number of points (n).
  # `do.call(rbind, ...)` efficiently combines the resulting dataframes.
  stats_df <- do.call(
    rbind,
    lapply(valid_lineages, function(ln) {
      m <- models[[ln]]
      sm <- summary(m)
      data.frame(
        lineage = ln,
        slope = unname(coef(m)[2]), # coef(m)[2] is the slope for areaOccupied_um2
        r2 = sm$r.squared,
        n = length(m$residuals),
        # Create a formatted label for the plot.
        label = sprintf("slope = %.3f\nR² = %.3f", unname(coef(m)[2]), sm$r.squared),
        # Position the label in the top-right corner of each facet.
        # 'Inf' is a ggplot trick to refer to the edge of the plot area.
        x = Inf,
        y = Inf,
        stringsAsFactors = FALSE
      )
    })
  )
  
  # --- 4. Generate the Plot ---
  p <- ggplot(df, aes(x = areaOccupied_um2, y = cellSize_um2)) +
    # Add the points, colored by lineage.
    geom_point(alpha = 0.8, aes(color = lineage)) +
    # Add the linear regression line for each lineage. `se = FALSE` hides the confidence interval.
    geom_smooth(method = "lm", se = FALSE, aes(color = lineage)) +
    # Create a separate panel for each lineage.
    facet_wrap(~ lineage, scales = if (free_scales) "free" else "fixed") +
    # Add the text annotations using the stats_df dataframe created earlier.
    # `hjust` and `vjust` are used to nudge the text away from the corner.
    geom_text(data = stats_df, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 1.05, vjust = 1.2) +
    labs(
      x = "Area occupied (µm²)",
      y = "Cell size (µm²)",
      color = "Lineage"
    ) +
    theme_bw()
  
  # Return the plot and model objects invisibly, so they don't print to the console
  # but are available for assignment (e.g., `results <- plot_size_vs_area_by_lineage(...)`).
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
                                       dbscan_eps = 0.5,     # New: DBSCAN 'eps' parameter
                                       dbscan_minPts = 5) {   # New: DBSCAN 'minPts' parameter
  
  # --- Robust Data Preparation ---
  # Find column names that are common to ALL dataframes in the list.
  common_cols <- Reduce(intersect, lapply(data_list, colnames))
  # Remove columns containing 'X' or 'Y', which are likely coordinates and not features.
  cols_without_X <- grep("X", common_cols, value = TRUE, invert = TRUE)
  final_cols <- grep("Y", cols_without_X, value = TRUE, invert = TRUE)
  # Subset each dataframe to only include these final common feature columns.
  subsetted_list <- lapply(data_list, function(df) df[, final_cols, drop = FALSE])
  
  # Combine all dataframes into a single large dataframe.
  # `.id = "source_id"` creates a new column with the name of the original list element.
  all_data <- bind_rows(subsetted_list, .id = "source_id") %>%
    mutate(
      # Extract the passage number from the source_id string (e.g., from 'A5' to 5).
      passage_number = as.numeric(str_extract(source_id, "(?<=A)[0-9]+"))
    )
  
  # --- Isolate and Scale Numeric Data for Projection ---
  # Select only the numeric feature columns for the dimensionality reduction algorithm.
  numeric_data <- all_data %>%
    select(where(is.numeric), -passage_number) %>%
    # Apply a log1p transformation (log(x+1)) to reduce the influence of extreme values.
    mutate(across(everything(), log1p)) %>% 
    # Scale the data (center to mean 0, scale to standard deviation 1). This is
    # crucial for distance-based algorithms like UMAP to work correctly.
    scale()
  
  # --- Perform Dimensionality Reduction ---
  method <- tolower(method)
  if (method == "umap") {
    set.seed(seed) # Ensure UMAP results are reproducible.
    projection_results <- uwot::umap(numeric_data, n_neighbors = n_neighbors)
  } else if (method == "pca") {
    pca_results <- prcomp(numeric_data)
    projection_results <- pca_results$x[, 1:2] # Take the first two principal components.
  } else if (method == "tsne") {
    set.seed(seed) # Ensure t-SNE results are reproducible.
    tsne_results <- Rtsne::Rtsne(numeric_data, dims = 2, check_duplicates = FALSE)
    projection_results <- tsne_results$Y
  } else {
    stop("Error: Method must be 'umap', 'pca', or 'tsne'.")
  }
  
  # --- Prepare Data for Plotting ---
  plot_data <- as.data.frame(projection_results)
  colnames(plot_data) <- c("Dim1", "Dim2")
  # Add back the passage number for coloring the points.
  plot_data$passage_number <- all_data$passage_number
  
  # --- (Optional) Clustering and Final Output Generation ---
  if (cluster) {
    # Check if the 'dbscan' package is installed before trying to use it.
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      stop("Package 'dbscan' is required. Please run: install.packages('dbscan')", call. = FALSE)
    }
    
    # Perform DBSCAN clustering on the 2D projection coordinates.
    # DBSCAN is useful here because it can find clusters of arbitrary shapes,
    # which UMAP often produces. 'eps' is the most important parameter to tune.
    set.seed(seed)
    db_clusters <- dbscan::dbscan(plot_data[, c("Dim1", "Dim2")], eps = dbscan_eps, minPts = dbscan_minPts)
    plot_data$cluster <- as.factor(db_clusters$cluster) # Add cluster assignments to the data.
    
    # --- Shuffle Data for Plotting ---
    # Shuffling prevents points from one passage or cluster from being systematically
    # plotted on top of others, which can obscure the true distribution.
    set.seed(seed)
    plot_data <- plot_data[sample(nrow(plot_data)), ]
    
    # --- Visualize Clustered Results ---
    # Points are colored by cluster ID. By DBSCAN convention, cluster '0' represents noise points.
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
    # This creates a summary table showing how many cells from each passage
    # fall into each discovered cluster. This is useful for analyzing how
    # the population structure changes over passages.
    cluster_distribution <- plot_data %>%
      count(cluster, passage_number, name = "n_cells") %>%
      group_by(passage_number) %>%
      mutate(proportion = n_cells / sum(n_cells)) %>%
      ungroup() %>%
      arrange(cluster, passage_number)
    
    # Return a list containing both the plot and the summary data.
    return(list(plot = projection_plot, distribution_data = cluster_distribution))
    
  } else {
    # --- Shuffle Data for Plotting (Original Path) ---
    set.seed(seed)
    plot_data <- plot_data[sample(nrow(plot_data)), ]
    
    # --- Visualize Results (Original Path) ---
    # If not clustering, color the points by passage number.
    projection_plot <- ggplot(plot_data, aes(x = Dim1, y = Dim2, color = passage_number)) +
      geom_point(alpha = 0.8, size = 2) +
      scale_color_viridis_c() + # Use a continuous color scale.
      labs(
        title = paste(toupper(method), title),
        x = paste(toupper(method), "Dimension 1"),
        y = paste(toupper(method), "Dimension 2"),
        color = "Passage"
      ) +
      theme_minimal()
    
    # Return only the plot object.
    return(projection_plot)
  }
}

#' Launch an Interactive Shiny App for Manual Clustering of Projection Data
#'
#' This function first computes a UMAP projection and then launches a Shiny
#' web application. The app allows the user to manually select groups of points
#' on the plot and assign them to named clusters.
#'
#' @param data_list A named list of data frames containing cell metrics.
#' @param method The dimensionality reduction method to use (currently only 'umap').
#' @param seed An integer for reproducibility.
#' @param n_neighbors UMAP parameter controlling local vs. global structure.
#' @param title A title for the plot within the Shiny app.
#' @return A dataframe containing the UMAP coordinates and the final, manually-
#'   assigned cluster labels for each cell.
interactive_cluster_projection <- function(data_list, 
                                           method = "umap", 
                                           seed = 42, 
                                           n_neighbors = 15,
                                           title = "Interactive Projection") {
  
  # --- [1] Data Preparation and Dimensionality Reduction (same as before) ---
  # This section is identical to the data preparation in `plot_projection_by_passage`.
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
    stop("This interactive version is currently configured for UMAP only.")
  }
  
  plot_data <- as.data.frame(projection_results)
  colnames(plot_data) <- c("Dim1", "Dim2")
  plot_data$passage_number <- all_data$passage_number
  # Initialize the cluster column that will be updated by the user in the app.
  plot_data$cluster <- "Unassigned"
  
  # --- [2] Define the Interactive Shiny Application ---
  
  # UI Definition: Describes the layout and appearance of the web app.
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
        # `brush = brushOpts(...)` is the key piece that enables rectangular selection on the plot.
        # The selected area information will be available in `input$plot_brush`.
        plotOutput("umap_plot", height = "600px", brush = brushOpts(id = "plot_brush")),
        h4("Cluster Distribution Summary"),
        DT::dataTableOutput("summary_table")
      )
    )
  )
  
  # Server Logic: Contains the instructions for how the app reacts to user input.
  server <- function(input, output, session) {
    
    # A `reactiveValues` object is a special list that can be updated and will
    # trigger re-rendering of outputs that depend on it.
    rv <- reactiveValues(data = plot_data, brush_selection = NULL)
    
    # This block defines how to render the main UMAP plot.
    output$umap_plot <- renderPlot({
      
      # `brushedPoints` identifies which rows in `rv$data` are inside the user's selection.
      brushed_pts <- brushedPoints(rv$data, input$plot_brush)
      
      ggplot() +
        # Plot all the points that have not yet been assigned a cluster.
        geom_point(data = filter(rv$data, cluster == "Unassigned"), aes(x = Dim1, y = Dim2), color = "grey80", alpha = 0.7) +
        # Plot the points that HAVE been assigned a cluster, with color.
        geom_point(data = filter(rv$data, cluster != "Unassigned"), aes(x = Dim1, y = Dim2, color = cluster), alpha = 0.9, size = 2) +
        # Add another layer to highlight the points currently being selected in red.
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
    
    # `observeEvent` defines what happens when an input (like a button click) changes.
    observeEvent(input$assign_btn, {
      # `req()` ensures this code only runs if a selection has been made and a name is provided.
      req(input$plot_brush, input$cluster_name)
      
      # Get the indices of the selected rows.
      selected_rows <- brushedPoints(rv$data, input$plot_brush, allRows = TRUE)$selected_
      
      # Update the 'cluster' column in our reactive data for these selected rows.
      # This change will automatically trigger `output$umap_plot` to re-render.
      rv$data$cluster[selected_rows] <- input$cluster_name
      
      # A user-friendly feature: automatically increment the cluster name for the next selection.
      new_num <- as.numeric(sub(".*_", "", input$cluster_name)) + 1
      if (!is.na(new_num)) {
        updateTextInput(session, "cluster_name", value = paste0("Cluster_", new_num))
      }
      
      # Clear the brush selection from the plot UI.
      session$resetBrush("plot_brush")
    })
    
    # Handle the "Reset" button click.
    observeEvent(input$reset_btn, {
      # Reset all assignments back to the initial "Unassigned" state.
      rv$data$cluster <- "Unassigned"
    })
    
    # Render the summary table below the plot.
    output$summary_table <- DT::renderDataTable({
      req(rv$data) # Ensure data is available.
      # Create a summary of cell counts per cluster and passage number.
      rv$data %>%
        filter(cluster != "Unassigned") %>%
        count(cluster, passage_number, name = "cell_count") %>%
        arrange(cluster, passage_number)
    }, options = list(pageLength = 5))
    
    # Handle the "Done" button click.
    observeEvent(input$done_btn, {
      # `stopApp` closes the Shiny application and returns a value.
      # Here, we return the final dataframe with all the user's cluster assignments.
      stopApp(returnValue = rv$data)
    })
  }
  
  # --- [3] Run the Application ---
  # This line launches the app. The R script will pause its execution here
  # until the user clicks the "Done" button in the app (which calls `stopApp`).
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
  
  # Load required packages. Using `require` inside a function is generally
  # discouraged in favor of `library` at the top of a script, but is seen here.
  require(purrr)
  require(dplyr)
  
  # Create a temporary directory if it doesn't exist.
  dir.create("~/Downloads/tmpLTE", showWarnings = FALSE)
  
  # Define the standard suffixes for the four image quartiles.
  suffixes <- c("tl","tr","bl","br")
  
  # --- 1. Build File Paths and Read Data ---
  # For each ID in the metadata, create the full paths to the four corresponding CSV files.
  file_groups <- purrr::map(closest_rows$id, ~{
    paste0(file.path(base_dir, .x), "_10x_ph_", suffixes, ".csv")
  })
  names(file_groups) <- closest_rows$id
  
  # Read all CSV files, combine the four quartiles for each ID, and add a 'suffix' column.
  csv_list <- purrr::map(file_groups, function(files) {
    dfs <- purrr::map(files, ~{
      if (file.exists(.x)) {
        # Attempt to read as tab-separated, then fall back to comma-separated.
        df <- read.csv(.x, sep = "\t")
        if(ncol(df)<3){
          df <- read.csv(.x)
        }
        # Extract the suffix (e.g., 'tl') from the filename and add it as a column.
        df$suffix <- sub(".*_ph_(.*)\\.csv$", "\\1", basename(.x))
        df
      } else {
        warning("File not found: ", .x)
        NULL
      }
    })
    # Combine the list of dataframes (tl, tr, bl, br) into a single dataframe for this ID.
    dplyr::bind_rows(dfs)
  })
  
  # Name the final list of dataframes by their corresponding IDs.
  names(csv_list) <- closest_rows$id
  
  # If we are not performing the density-based selection, return the combined data now.
  if (!isTRUE(select_by_global_density)) {
    return(csv_list)
  }
  
  # ---- 2. Compute Densities per Quartile for Each Entry ----
  # Density can be either the raw cell count or count normalized by area.
  dens_list <- purrr::map(csv_list, function(df) {
    if (is.null(df) || nrow(df) == 0) return(setNames(numeric(0), character(0)))
    counts <- table(df$suffix)
    # Ensure a consistent order for the quartiles.
    counts <- counts[intersect(names(counts), suffixes)]
    counts <- counts[order(match(names(counts), suffixes))] 
    
    # If an area is provided, calculate density as cells/area; otherwise, use raw counts.
    if (!is.null(area_per_quartile_um2) && is.finite(area_per_quartile_um2) && area_per_quartile_um2 > 0) {
      dens <- as.numeric(counts) / area_per_quartile_um2
    } else {
      dens <- as.numeric(counts) 
    }
    names(dens) <- names(counts)
    dens
  })
  
  # --- 3. Estimate the Global Target Density (Mode) ---
  # Pool all calculated densities from all quartiles of all samples.
  all_dens <- unlist(dens_list, use.names = FALSE)
  all_dens <- all_dens[is.finite(all_dens)]
  
  if (length(all_dens) == 0) {
    warning("No densities found to compute a global target; returning original list.")
    return(csv_list)
  }
  
  # To find the most common density, we first bin the values. This makes the
  # mode calculation more stable and less sensitive to minor fluctuations.
  bins <- floor(all_dens / bin_width) * bin_width
  bin_tab <- sort(table(bins), decreasing = TRUE)
  # It's possible for multiple bins to have the same max count.
  top_bins <- names(bin_tab)[bin_tab == max(bin_tab)]
  # To break this tie, we choose the bin whose center is closest to the overall median density.
  bin_centers <- as.numeric(top_bins) + bin_width/2
  x <- bin_centers[which.min(abs(bin_centers - median(all_dens)))]
  
  # ---- 4. Filter Each Entry to Keep Only the Best Quartile ----
  chosen_suffix <- list()
  csv_list_filtered <- purrr::imap(csv_list, function(df, entry_name) {
    dens <- dens_list[[entry_name]]
    if (length(dens) == 0 || is.null(df) || nrow(df) == 0) return(df)
    
    # Calculate the absolute difference between each quartile's density and the global target `x`.
    dists <- abs(dens - x)
    
    # Find the suffix(es) with the minimum distance.
    best <- names(dists)[which(dists == min(dists))]
    # If there's a tie, use the predefined `tie_order` to pick one consistently.
    best <- best[order(match(best, tie_order))][1]
    
    # Store which suffix was chosen for this entry.
    chosen_suffix[[entry_name]] <<- best
    
    # Filter the original dataframe to keep only the rows from the chosen suffix.
    df %>% dplyr::filter(.data$suffix == best)
  })
  
  # Attach metadata about the selection process as attributes to the returned object.
  # This is useful for downstream inspection and debugging.
  attr(csv_list_filtered, "global_target_density") <- x
  attr(csv_list_filtered, "bin_width") <- bin_width
  attr(csv_list_filtered, "area_per_quartile_um2") <- area_per_quartile_um2
  attr(csv_list_filtered, "chosen_suffix") <- chosen_suffix
  attr(csv_list_filtered, "densities") <- dens_list
  
  csv_list_filtered
}


###############################################################
###############################################################

#' Find the Best Clustering for Karyotype Data
#'
#' This function performs hierarchical clustering and determines the optimal
#' number of clusters using the silhouette method, unless a specific number
#' of clusters is provided.
#'
#' @param allKaryo A numeric matrix of karyotype data (cells x chromosomes).
#' @param numClusters Optional: an integer specifying the number of clusters. If NULL,
#'   the optimal number is determined automatically.
#' @param hFun The hierarchical clustering function to use.
#' @param dFun The distance function to use.
#' @return A numeric vector of cluster assignments for each row in `allKaryo`.
findBestClustering<-function(allKaryo, numClusters=NULL, hFun=function(x) hclust(x, method="ward.D2"), dFun = function(x) dist(x, method="manhattan")){
  library(cluster)  
  # First, perform hierarchical clustering to get the dendrogram.
  # The `heatmap.2` function is used here primarily for its clustering capabilities.
  hm=heatmap.2(allKaryo, hclustfun=hFun,distfun=dFun)
  
  # If the user has specified the number of clusters, use that value.
  if(!is.null(numClusters)){
    k=numClusters
  }else{
    # Otherwise, find the optimal 'k' by maximizing the median silhouette coefficient.
    silhouettecoefs=rep(NA,nrow(allKaryo))
    # Iterate through a range of possible cluster numbers (from 2 to n-1).
    for(k in 2:(nrow(allKaryo)-1)){
      # For each 'k', cut the dendrogram to get cluster assignments.
      clusters=cutree(as.hclust(hm$rowDendrogram), k=k)
      # Calculate the silhouette information for this clustering.
      sil <- summary(silhouette(clusters, dist(allKaryo)))
      # Store the median silhouette score. A higher score indicates better clustering.
      silhouettecoefs[k]= sil$si.summary["Median"]
    }
    # Choose the 'k' that resulted in the highest median silhouette score.
    k = which.max(silhouettecoefs)
  }
  # Cut the dendrogram using the chosen (or determined) 'k'.
  clusters=cutree(as.hclust(hm$rowDendrogram), k=k)
  return(clusters)
}

#' Process Raw Copy Number Data into Karyotypes
#'
#' This function simplifies segmented copy number data into whole-chromosome
#' copy numbers and generates a summary of unique karyotypes.
#'
#' @param cn A matrix of copy number data (cells x segments). Column names are
#'   expected in the format "chr:start-end".
#' @param ploidy The baseline ploidy to assume for chromosomes not present in the data.
#' @return A list containing `karyo` (a frequency table of unique karyotypes) and
#'   `cn` (the simplified cell x chromosome copy number matrix).
getKaryo<-function(cn,ploidy){
  ## For each chromosome, find its largest segment and use the copy number of
  ## that segment as the representative copy number for the entire chromosome.
  
  # Parse the "chr:start-end" column names into a structured dataframe.
  segments= sapply(sapply(strsplit(colnames(cn),":"),"[[",2), function(x) strsplit(x[[1]],"-")[[1]],simplify = F)
  segments= as.data.frame(do.call(rbind,sapply(segments, as.numeric,simplify = F)))
  rownames(segments) = colnames(cn)
  colnames(segments) = c("start","end")
  segments$length=1+segments$end-segments$start
  segments$chr = as.numeric(sapply(strsplit(colnames(cn),":"),"[[",1))
  
  # For each chromosome, identify the segment with the maximum length.
  chrsegments=sapply(unique(segments$chr), function(x) segments[segments$chr==x,,drop=FALSE],simplify = F)
  chrsegments=sapply(chrsegments, function(x) x[which.max(x$length),,drop=F],simplify = F)
  chrsegments = do.call(rbind,chrsegments)
  
  # Subset the original copy number matrix to only include these representative segments.
  cn=cn[,rownames(chrsegments)]
  # Rename the columns to be just the chromosome number.
  colnames(cn)=chrsegments$chr
  
  ## For any autosomes (1-22) not present in the data, assume their copy number
  ## is equal to the baseline ploidy for all cells.
  otherchr = setdiff(1:22,colnames(cn))
  cn_ = matrix(ploidy,nrow(cn),length(otherchr))
  colnames(cn_)=otherchr
  cn = cbind(cn,cn_)
  gplots::heatmap.2(cn,trace='n',symbreaks = F,symkey=F)
  
  ## Create a unique string representation for each cell's karyotype.
  cn=round(cn)
  # Collapse the copy numbers for each cell (row) into a single string (e.g., "2.2.3.2...").
  karyo=apply(cn,1,paste0,collapse=".");
  names(karyo) = rownames(cn)
  # Count the occurrences of each unique karyotype string.
  karyo_in= plyr::count(karyo)
  rownames(karyo_in)=karyo_in$x
  return(list(karyo=karyo_in[,'freq',drop=F], cn=cn ))
  
}

#' Cluster Karyotypes from Multiple Origins
#'
#' A high-level wrapper function that queries copy number profiles from a
#' database, processes them into whole-chromosome karyotypes, performs
#' hierarchical clustering, and generates visualizations.
#'
#' @param origins A vector of sample/biopsy identifiers to query from the database.
#' @param whichP The perspective to use in the database query.
#' @param depth The query depth for fetching related clones.
#' @param path2lanscape (Not used in this version) Path to landscape data.
#' @param numClusters The desired number of clusters.
#' @param capping An optional upper limit to cap copy number values for visualization.
#' @param method The agglomeration method for hierarchical clustering (e.g., "complete", "ward.D2").
#' @param chrwhole A dataframe with chromosome length information, used for weighted distance.
#' @return A list containing cluster assignments, summary stats, and raw data.
clusterKaryotypes <- function(origins, whichP = "GenomePerspective", depth  = 1,  path2lanscape   = NULL, numClusters = NULL, capping = NULL, method = "complete", chrwhole = NULL){
  # If chromosome length data is not provided, download it.
  if(is.null(chrwhole)){
    devtools::source_url("https://github.com/noemiandor/Utils/blob/master/grpstats.R?raw=TRUE")
    x <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz", 
               col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
    chrarms=x[ , .(length = sum(chromEnd - chromStart)),by = .(chrom, arm = substring(name, 1, 1)) ]
    chrwhole=grpstats(as.matrix(chrarms$length),chrarms$chrom, "sum")$sum
  }
  # --- Set default parameters ---
  ploidy  <- 2
  min_obs <- 5
  dt      <- 5
  
  # --- 1) Query subprofiles from DB and store them in list 'X' ---
  mydb  <- cloneid::connect2DB()
  X     <- list()
  for(biopsy in origins){
    # Query for parent clones.
    stmt <- paste0(
      "SELECT cloneID, size, alias, parent 
       FROM Perspective
       WHERE size=1 AND whichPerspective='", whichP, "' AND origin='", biopsy, "'"
    )
    rs <- suppressWarnings(DBI::dbSendQuery(mydb, stmt))
    sps <- DBI::fetch(rs, n = -1)
    
    # If depth > 1, also query for children of the parent clones.
    if(depth > 1) {
      stmt <- paste0(
        "SELECT cloneID, size, alias, parent 
         FROM Perspective 
         WHERE parent IN (", paste(sps$cloneID, collapse=","), ")"
      )
      rs <- suppressWarnings(DBI::dbSendQuery(mydb, stmt))
      sps <- DBI::fetch(rs, n=-1)
    }
    
    # For each cloneID, fetch its detailed copy number profile.
    x <- sapply(
      sps$cloneID, 
      function(cid) cloneid::getSubProfiles(cloneID_or_sampleName = cid, whichP = whichP), 
      simplify = FALSE
    )
    # Combine the profiles for this biopsy into a single matrix.
    X[[biopsy]] <- do.call(cbind, x)
  }
  
  ##############################################################################
  # 2) Merge across origins for karyotyping/clustering
  ##############################################################################
  
  # Process the raw profiles from each origin into whole-chromosome copy numbers.
  cnts_list  <- sapply(X, function(mat) getKaryo(t(mat), ploidy)$cn, simplify = FALSE)
  
  # Create a vector that tracks the origin (biopsy) of each cell.
  sampleID <- unlist(
    sapply(names(cnts_list), function(x) rep(x, nrow(cnts_list[[x]])))
  )
  
  # Combine the copy number matrices from all origins into one large matrix.
  cnts_combined <- do.call(rbind, cnts_list)
  
  
  # --- Define clustering parameters ---
  # Hierarchical clustering function
  hFun <- function(x) stats::hclust(x, method = method);
  # Custom distance function that weights by chromosome length.
  dFun <- function(x) chrWeightedDist(x, chrwhole=chrwhole);
  
  # Perform the clustering. `findBestClustering` will either use the user-supplied
  # `numClusters` or determine the optimal number itself.
  # The +1 is likely to shift cluster labels from a 0-based to a 1-based index.
  clusters <- findBestClustering(cnts_combined, numClusters = numClusters, hFun=hFun, dFun = dFun) + 1
  
  ##############################################################################
  # 3) Plot heatmap with hierarchical clustering to extract dendrogram 
  ##############################################################################
  # Create a safe filename from the origin names.
  tmp <- substr(paste(origins, collapse = "__"), 1, 90)
  pdf(paste0(tmp, ".pdf"))
  
  # Create a color mapping for the `RowSideColors` annotation, where each
  # unique sampleID (origin) gets a distinct color.
  uniqueIDs <- unique(sampleID)
  colVec    <- rep("NA", length(uniqueIDs))
  names(colVec) <- uniqueIDs
  
  # Assign specific colors (e.g., grey for controls).
  idxControl <- grep("C_", names(colVec))
  colVec[idxControl] <- gray.colors(length(idxControl))
  
  # Assign a color palette to the remaining samples.
  remaining <- setdiff(seq_along(colVec), idxControl)
  if(length(remaining) > 0) {
    colPalette <- brewer.pal(min(length(remaining), 12), "Paired")
    if(length(remaining) > length(colPalette)) {
      # If there are more samples than colors, extend the palette.
      colPalette <- colorRampPalette(colPalette)(length(remaining))
    }
    colVec[remaining] <- colPalette[seq_along(remaining)]
  }
  
  # Apply capping to the data matrix for better color scaling in the heatmap.
  tmp=as.matrix(cnts_combined)
  if(!is.null(capping)){
    tmp[tmp>capping] = capping 
  }
  
  # Generate the heatmap.
  hm <- heatmap.2(
    x           = tmp,
    margins     = c(15,15),
    # Color the rows by their assigned cluster number.
    colRow      = clusters[rownames(cnts_combined)], 
    trace       = 'none',
    Colv        = TRUE, # Cluster columns
    dendrogram  = "row",
    # Add a sidebar of colors indicating the sample origin for each row.
    RowSideColors = colVec[sampleID],
    key.xlab    = "copy number",
    key.title   = "",
    col         = matlab::fliplr(rainbow(20))[5:12],
    hclustfun   = hFun,
    distfun     = dFun
  )
  
  # Add a legend for the RowSideColors.
  legend(
    "topright",
    legend = names(colVec),
    fill   = colVec,
    cex    = 0.5
  )
  
  # --- Generate a boxplot of ploidy values per cluster ---
  # First, calculate the ploidy for each cell.
  ploidy_vals <- calcPloidy(cnts_combined,chrwhole)
  # Then, create the boxplot.
  boxplot(
    ploidy_vals ~ factor(clusters[names(ploidy_vals)], levels = unique(clusters)),
    xlab   = "Cluster",
    ylab   = "Ploidy",
    main   = "",
    col    = unique(clusters)
  )
  
  dev.off()
  
  ##############################################################################
  # 4) Summarize and return results
  ##############################################################################
  
  # Calculate summary statistics (mean, median) for each sample origin.
  cnts_summary <- grpstats(cnts_combined, sampleID, statscols = c("mean","median"))
  
  # Attach sample IDs as names to the results vectors for easy tracking.
  names(ploidy_vals) <- names(clusters) <- sampleID
  rownames(cnts_combined) = paste0(rownames(cnts_combined),"_",sampleID)
  
  # Return a structured list containing all the key outputs.
  return(list(
    clusters    = clusters,        # Cluster assignment for each cell
    cnts        = cnts_summary,    # Aggregated stats per origin
    CN      = cnts_combined,      # The full, merged copy number matrix
    distanceFun = dFun,            # The distance function used
    origins     = origins,         # The input origins
    ploidy_vals = ploidy_vals      # The calculated ploidy for each cell
  ))
}