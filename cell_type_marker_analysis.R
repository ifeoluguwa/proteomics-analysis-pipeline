#!/usr/bin/env Rscript
#
# Complete Cell Type Marker Validation Analysis
# This script combines all functions and analysis steps into a single file
#

# Load required packages
required_packages <- c("writexl", "pheatmap", "RColorBrewer", "tools")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Please run your package installation script first, then try again.\n")
    cat("Missing package:", pkg, "\n")
    stop("Required package missing.")
  }
}

#------------------------------------------------------------------------------
# FUNCTION DEFINITIONS
#------------------------------------------------------------------------------

#' Find cell type markers in datasets
#' 
#' This function searches for the specified markers in each dataset using flexible pattern matching
#' 
#' @param micro_input Data frame or matrix containing microglia proteomics data
#' @param neuro_input Data frame or matrix containing neurons proteomics data
#' @param microglia_markers Vector of microglia marker names to search for
#' @param neuron_markers Vector of neuron marker names to search for
#' @param verbose Logical; if TRUE, prints detailed matching information
#' @return List with found microglia and neuron markers
find_cell_type_markers <- function(micro_input, neuro_input, 
                                   microglia_markers, neuron_markers,
                                   verbose = TRUE) {
  microglia_markers_found <- c()
  neuron_markers_found <- c()
  
  # Convert matrices to data frames if needed for consistent handling
  if (!is.data.frame(micro_input)) {
    micro_input <- as.data.frame(micro_input)
  }
  if (!is.data.frame(neuro_input)) {
    neuro_input <- as.data.frame(neuro_input)
  }
  
  # First, verify that protein names match between datasets
  microglia_proteins <- rownames(micro_input)
  neuron_proteins <- rownames(neuro_input)
  
  # Check name consistency
  protein_overlap <- intersect(microglia_proteins, neuron_proteins)
  proteins_only_in_microglia <- setdiff(microglia_proteins, neuron_proteins)
  proteins_only_in_neurons <- setdiff(neuron_proteins, microglia_proteins)
  
  if (verbose) {
    cat("Proteins in both datasets:", length(protein_overlap), "\n")
    cat("Proteins only in microglia:", length(proteins_only_in_microglia), "\n")
    cat("Proteins only in neurons:", length(proteins_only_in_neurons), "\n\n")
  }
  
  # Search for microglia markers
  if (verbose) {
    cat("Searching for microglia markers...\n")
  }
  
  for (marker in microglia_markers) {
    # Try exact match first
    if (marker %in% microglia_proteins) {
      microglia_markers_found <- c(microglia_markers_found, marker)
      if (verbose) cat("  Found exact match for:", marker, "\n")
    } else {
      # More flexible pattern matching
      marker_dash_to_under <- gsub("-", "_", marker)
      marker_under_to_dash <- gsub("_", "-", marker)
      
      patterns <- c(
        paste0("^", marker, "$"),  # Exact match (case-insensitive)
        paste0("^", marker),       # Starts with marker
        paste0(marker, "$"),       # Ends with marker
        marker,                    # Contains marker
        paste0("^", marker_dash_to_under, "$"),
        paste0("^", marker_dash_to_under),
        paste0(marker_dash_to_under, "$"),
        marker_dash_to_under,
        paste0("^", marker_under_to_dash, "$"),
        paste0("^", marker_under_to_dash),
        paste0(marker_under_to_dash, "$"),
        marker_under_to_dash
      )
      
      found_match <- FALSE
      for (pattern in patterns) {
        if (found_match) break
        
        matches <- grep(pattern, microglia_proteins, ignore.case = TRUE, value = TRUE)
        
        if (length(matches) > 0) {
          # Prioritize shortest matches as they're likely more specific
          matches <- matches[order(nchar(matches))]
          microglia_markers_found <- c(microglia_markers_found, matches[1])
          if (verbose) cat("  Found fuzzy match for:", marker, "->", matches[1], "\n")
          found_match <- TRUE
        }
      }
      
      if (!found_match && verbose) {
        cat("  No match found for microglia marker:", marker, "\n")
      }
    }
  }
  
  # Search for neuron markers
  if (verbose) {
    cat("\nSearching for neuron markers...\n")
  }
  
  for (marker in neuron_markers) {
    # Try exact match first
    if (marker %in% neuron_proteins) {
      neuron_markers_found <- c(neuron_markers_found, marker)
      if (verbose) cat("  Found exact match for:", marker, "\n")
    } else {
      # More flexible pattern matching
      marker_dash_to_under <- gsub("-", "_", marker)
      marker_under_to_dash <- gsub("_", "-", marker)
      
      patterns <- c(
        paste0("^", marker, "$"),  # Exact match (case-insensitive)
        paste0("^", marker),       # Starts with marker
        paste0(marker, "$"),       # Ends with marker
        marker,                    # Contains marker
        paste0("^", marker_dash_to_under, "$"),
        paste0("^", marker_dash_to_under),
        paste0(marker_dash_to_under, "$"),
        marker_dash_to_under,
        paste0("^", marker_under_to_dash, "$"),
        paste0("^", marker_under_to_dash),
        paste0(marker_under_to_dash, "$"),
        marker_under_to_dash
      )
      
      found_match <- FALSE
      for (pattern in patterns) {
        if (found_match) break
        
        matches <- grep(pattern, neuron_proteins, ignore.case = TRUE, value = TRUE)
        
        if (length(matches) > 0) {
          # Prioritize shortest matches as they're likely more specific
          matches <- matches[order(nchar(matches))]
          neuron_markers_found <- c(neuron_markers_found, matches[1])
          if (verbose) cat("  Found fuzzy match for:", marker, "->", matches[1], "\n")
          found_match <- TRUE
        }
      }
      
      if (!found_match && verbose) {
        cat("  No match found for neuron marker:", marker, "\n")
      }
    }
  }
  
  # Return unique matches to avoid duplicates
  return(list(
    microglia = unique(microglia_markers_found),
    neuron = unique(neuron_markers_found)
  ))
}

#' Create mean expression table for cell type markers
#' 
#' @param micro_input Matrix or data frame with microglia protein expression data (log2-normalized)
#' @param neuro_input Matrix or data frame with neuron protein expression data (log2-normalized)
#' @param microglia_markers_present Vector of microglia markers found in the data
#' @param neuron_markers_present Vector of neuron markers found in the data
#' @param microglia_cols Vector of column names to use from microglia data
#' @param neuron_cols Vector of column names to use from neuron data
#' @return Data frame with mean expression values for each marker in each cell type
create_mean_expression_table <- function(micro_input, neuro_input,
                                         microglia_markers_present, neuron_markers_present,
                                         microglia_cols, neuron_cols) {
  
  # Convert to data frame if needed
  if (!is.data.frame(micro_input)) {
    micro_input <- as.data.frame(micro_input)
  }
  if (!is.data.frame(neuro_input)) {
    neuro_input <- as.data.frame(neuro_input)
  }
  
  # Initialize results data frame
  all_markers <- c(microglia_markers_present, neuron_markers_present)
  mean_data <- data.frame(
    Marker = all_markers,
    Microglia = NA,
    Neuron = NA,
    Marker_Type = NA,
    stringsAsFactors = FALSE
  )
  rownames(mean_data) <- all_markers
  
  # Calculate mean expression for microglia markers
  for (marker in microglia_markers_present) {
    if (marker %in% rownames(micro_input)) {
      mean_data[marker, "Microglia"] <- mean(as.numeric(micro_input[marker, microglia_cols]), na.rm = TRUE)
      mean_data[marker, "Marker_Type"] <- "Microglia Marker"
    }
    
    # Also calculate expression in neuron dataset if available
    if (marker %in% rownames(neuro_input)) {
      mean_data[marker, "Neuron"] <- mean(as.numeric(neuro_input[marker, neuron_cols]), na.rm = TRUE)
    }
  }
  
  # Calculate mean expression for neuron markers
  for (marker in neuron_markers_present) {
    if (marker %in% rownames(neuro_input)) {
      mean_data[marker, "Neuron"] <- mean(as.numeric(neuro_input[marker, neuron_cols]), na.rm = TRUE)
      mean_data[marker, "Marker_Type"] <- "Neuron Marker"
    }
    
    # Also calculate expression in microglia dataset if available
    if (marker %in% rownames(micro_input)) {
      mean_data[marker, "Microglia"] <- mean(as.numeric(micro_input[marker, microglia_cols]), na.rm = TRUE)
    }
  }
  
  # For any missing values, replace with a low expression value
  min_value <- min(c(mean_data$Microglia, mean_data$Neuron), na.rm = TRUE)
  mean_data$Microglia[is.na(mean_data$Microglia)] <- min_value - 0.5
  mean_data$Neuron[is.na(mean_data$Neuron)] <- min_value - 0.5
  
  # Remove rows where both values are NA or Marker_Type is NA
  mean_data <- mean_data[!is.na(mean_data$Marker_Type), ]
  
  # Calculate differential metrics properly for log2-normalized data
  # Since the input data is already log2-transformed:
  # - The difference of log2 values is the log2 fold change
  # - The actual fold change is 2 raised to the power of this difference
  mean_data$Log2FC <- mean_data$Neuron - mean_data$Microglia
  mean_data$Fold_Change <- 2^mean_data$Log2FC
  
  # Clip extreme log2FC values for better visualization
  mean_data$Log2FC_Clipped <- pmin(pmax(mean_data$Log2FC, -3), 3)
  
  return(mean_data)
}

#' Create fold change heatmap to visualize cell type specificity
#' 
#' @param mean_data Data frame with mean expression values from create_mean_expression_table
#' @param save_file Optional filename to save the heatmap as a PNG file
#' @return pheatmap object
create_fold_change_heatmap <- function(mean_data, save_file = NULL) {
  # Load required package if not already loaded
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("pheatmap")
    library(pheatmap)
  }
  
  # Create fold change matrix
  fold_change_matrix <- as.matrix(mean_data$Log2FC_Clipped)
  rownames(fold_change_matrix) <- rownames(mean_data)
  colnames(fold_change_matrix) <- "Neuron/Microglia Log2FC"
  
  # Determine marker order
  microglia_indices <- which(mean_data$Marker_Type == "Microglia Marker")
  neuron_indices <- which(mean_data$Marker_Type == "Neuron Marker")
  marker_order <- c(rownames(mean_data)[microglia_indices], rownames(mean_data)[neuron_indices])
  
  # Reorder fold change matrix
  fold_change_matrix <- fold_change_matrix[marker_order, , drop = FALSE]
  
  # Create blue-white-red gradient centered at 0
  fold_change_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  breaks_fc <- seq(-3, 3, length.out = 101)
  
  # Create row annotation for markers
  row_annotation <- data.frame(
    Marker_Type = mean_data[marker_order, "Marker_Type"]
  )
  rownames(row_annotation) <- marker_order
  
  # Define annotation colors
  ann_colors <- list(
    Marker_Type = c("Microglia Marker" = "#FF9999", "Neuron Marker" = "#99CCFF")
  )
  
  # Create the fold change heatmap
  p_fc <- pheatmap::pheatmap(
    fold_change_matrix,
    color = fold_change_colors,
    breaks = breaks_fc,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    annotation_colors = ann_colors,
    fontsize_row = 10,
    fontsize_col = 12,
    cellwidth = 60,
    cellheight = 10,
    main = "Cell Type Specificity (Log2 Fold Change: Neuron/Microglia)",
    border_color = NA,
    annotation_names_row = FALSE,
    gaps_row = length(microglia_indices),
    treeheight_row = 0,
    treeheight_col = 0,
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black",
    fontsize_number = 7,
    angle_col = 0
  )
  
  # Save the plot if a filename is provided
  if (!is.null(save_file)) {
    png(save_file, width = 8, height = 10, units = "in", res = 300)
    print(p_fc)
    dev.off()
    cat("Fold change heatmap saved to:", save_file, "\n")
  }
  
  return(p_fc)
}

#' Create and display expression heatmap for cell type validation
#' 
#' @param mean_data Data frame with mean expression values from create_mean_expression_table
#' @param title_suffix String to append to heatmap title
#' @param use_scaled Logical; if TRUE, uses z-scaled data for the heatmap
#' @param save_file Optional filename to save the heatmap as a PNG file
#' @return pheatmap object
create_expression_heatmap <- function(mean_data, title_suffix = "", use_scaled = FALSE, save_file = NULL) {
  # Load required package if not already loaded
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("pheatmap")
    library(pheatmap)
  }
  
  # Prepare data matrix for heatmap
  heatmap_data <- mean_data[, c("Microglia", "Neuron")]
  rownames(heatmap_data) <- rownames(mean_data)
  
  if (use_scaled) {
    # Scale data (emphasizes differences between cell types)
    display_data <- t(scale(t(heatmap_data)))
    color_breaks <- seq(-2, 2, length.out = 100)
    heatmap_title <- paste0("Cell Type Validation ", title_suffix, " (Z-scaled)")
  } else {
    # Use a simple min-max normalization for more intuitive display
    min_val <- min(heatmap_data, na.rm = TRUE)
    max_val <- max(heatmap_data, na.rm = TRUE)
    
    # Normalize to range [-1, 1] for the color scale
    display_data <- 2 * (heatmap_data - min_val) / (max_val - min_val) - 1
    
    # Create custom breaks to ensure more color variation
    color_breaks <- seq(-1, 1, length.out = 100)
    heatmap_title <- paste0("Cell Type Validation ", title_suffix)
  }
  
  # Reorder rows to group markers by cell type
  microglia_indices <- which(mean_data$Marker_Type == "Microglia Marker")
  neuron_indices <- which(mean_data$Marker_Type == "Neuron Marker")
  marker_order <- c(rownames(mean_data)[microglia_indices], rownames(mean_data)[neuron_indices])
  
  # Make sure we only include markers that are in the data
  display_data <- display_data[marker_order, ]
  
  # Create row annotation to identify which markers are for which cell type
  row_annotation <- data.frame(
    Marker_Type = mean_data[marker_order, "Marker_Type"]
  )
  rownames(row_annotation) <- marker_order
  
  # Define annotation colors for the sidebar
  ann_colors <- list(
    Marker_Type = c("Microglia Marker" = "#FF9999", "Neuron Marker" = "#99CCFF")
  )
  
  # Create the heatmap
  p <- pheatmap::pheatmap(
    display_data,
    color = colorRampPalette(c("navy", "royalblue", "lightblue", "white", 
                               "mistyrose", "salmon", "firebrick3"))(100),
    breaks = color_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    annotation_colors = ann_colors,
    fontsize_row = 10,
    fontsize_col = 12,
    cellwidth = 60,
    cellheight = 10,
    main = heatmap_title,
    border_color = NA,
    annotation_names_row = FALSE,
    gaps_row = length(microglia_indices),
    treeheight_row = 0,
    treeheight_col = 0,
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black",
    fontsize_number = 7
  )
  
  # Save the plot if a filename is provided
  if (!is.null(save_file)) {
    png(save_file, width = 8, height = 10, units = "in", res = 300)
    print(p)
    dev.off()
    cat("Heatmap saved to:", save_file, "\n")
  }
  
  return(p)
}


#------------------------------------------------------------------------------
# MAIN ANALYSIS CODE
#------------------------------------------------------------------------------

# Set output directory
output_dir <- "cell_type_validation_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#------------------------------------------------------------------------------
# Part 1: Load and verify data
#------------------------------------------------------------------------------
cat("Loading data and verifying protein name consistency...\n")

# Convert matrices to data frames with new variable names to avoid conflicts
microglia_df_markers <- as.data.frame(microglia_combat_corrected)
neurons_df_markers <- as.data.frame(neurons_std_sva_corrected)

# Verify data dimensions and structure
cat("Microglia dataset dimensions:", nrow(microglia_df_markers), "proteins x", ncol(microglia_df_markers), "samples\n")
cat("Neurons dataset dimensions:", nrow(neurons_df_markers), "proteins x", ncol(neurons_df_markers), "samples\n")

# Get the intersection of proteins for cross-cell type comparison
proteins_in_both <- intersect(rownames(microglia_df_markers), rownames(neurons_df_markers))
cat("Proteins present in both datasets:", length(proteins_in_both), "\n")

# Display some sample protein names to verify naming conventions
cat("\nSample protein names in microglia dataset (first 5):\n")
print(head(rownames(microglia_df_markers), 5))
cat("\nSample protein names in neurons dataset (first 5):\n")
print(head(rownames(neurons_df_markers), 5))

#------------------------------------------------------------------------------
# Part 2: Define cell-type marker lists
#------------------------------------------------------------------------------
cat("\nDefining cell-type marker lists...\n")

# Define the marker genes for each cell type
microglia_markers <- c(
  # Common names
  "IBA-1", "TMEM119", "P2Y12", "CD68", "CD11b", "CD14", "CD45", "CD80", "CD115", "CX3CR1", 
  "FCER1G", "FCRLS", "Sirpa", "Siglec", "Glut5", "HexB", "Vimentin", "Ferritin", "Sall1", 
  "CD16", "CD32", "CD40", "CD86", "F4/80", "MHC II", "CD163", "CD206", "iNOS", "ARG1", 
  "Ym1", "FIZZ1", "MMP9", "MMP12", "NIBAN1", "PYCARD",
  # UniProt style IDs
  "AIF1_HUMAN", "TMEM119_HUMAN", "P2RY12_HUMAN", "CD68_HUMAN", "ITGAM_HUMAN", "CD14_HUMAN",
  "PTPRC_HUMAN", "CD80_HUMAN", "CSF1R_HUMAN", "CX3CR1_HUMAN", "FCER1G_HUMAN", "HEXB_HUMAN",
  "VIM_HUMAN", "FTH1_HUMAN", "SALL1_HUMAN", "FCGR3A_HUMAN", "FCGR2A_HUMAN", "CD40_HUMAN",
  "CD86_HUMAN", "EMR1_HUMAN", "CD163_HUMAN", "MRC1_HUMAN", "NOS2_HUMAN", "ARG1_HUMAN",
  "MMP9_HUMAN", "MMP12_HUMAN", "PU1_HUMAN", "SPI1_HUMAN", "NIBA1_HUMAN", "TNR5_HUMAN", "ASC_HUMAN"
)

neuron_markers <- c(
  # Common names
  "MAP2", "NCAM1", "TUBB3", "vGlut", "CamKII", "DBH", "EBF3", "GRIA1", "GRIA2", "GRIA3", 
  "GRIA4", "GRIN1", "GRIN2A", "GRIN2B", "DLG4", "SYN1", "BCL11B", "TBR1", "CUX1", "POU3F2", 
  "SATB2", "BRN2", "SLC17A6", "SLC17A7", "SLC32A1", "TENM2", "TTR", "VAT1L", "NEFL",
  "NEUROD6", "GAD1", "GAD2", "CNR1", "MAPT", "SYN2", "GPHN", "GABRA2", "GABRA5", "GABRB1", 
  # UniProt style IDs
  "MAP2_HUMAN", "NCAM1_HUMAN", "TBB3_HUMAN", "TUBB3_HUMAN", "CAMK2_HUMAN", "DBH_HUMAN", 
  "EBF3_HUMAN", "GRIA1_HUMAN", "GRIA2_HUMAN", "GRIA3_HUMAN", "GRIA4_HUMAN", "NMDZ1_HUMAN", 
  "NMDE1_HUMAN", "NMDE2_HUMAN", "DLG4_HUMAN", "SYN1_HUMAN", "BCL11B_HUMAN", "TBR1_HUMAN", 
  "CUX1_HUMAN", "PO3F2_HUMAN", "SATB2_HUMAN", "BRN2_HUMAN", "VGLU2_HUMAN", "VGLU1_HUMAN", 
  "VIAAT_HUMAN", "TEN2_HUMAN", "TTHY_HUMAN", "VAT1L_HUMAN", "NFL_HUMAN", "NDF6_HUMAN", 
  "DCE1_HUMAN", "DCE2_HUMAN", "CNR1_HUMAN", "TAU_HUMAN", "SYN_HUMAN", "GEPH_HUMAN", 
  "GBRA2_HUMAN", "GBRA5_HUMAN", "GBRB1_HUMAN"
)

cat("Defined", length(microglia_markers), "potential microglia markers\n")
cat("Defined", length(neuron_markers), "potential neuron markers\n")

#------------------------------------------------------------------------------
# Part 3: Identify markers in the datasets
#------------------------------------------------------------------------------
cat("\nSearching for markers in datasets...\n")

# Select vehicle columns for each dataset (assuming vehicle columns end with _V)
microglia_vehicle_columns <- grep("_V[0-9]?$|_V$", colnames(microglia_df_markers), value = TRUE)
neuron_vehicle_columns <- grep("_V[0-9]?$|_V$", colnames(neurons_df_markers), value = TRUE)

cat("Selected", length(microglia_vehicle_columns), "vehicle columns for microglia\n")
cat("Selected", length(neuron_vehicle_columns), "vehicle columns for neurons\n")

# Find markers in the datasets
marker_results <- find_cell_type_markers(
  microglia_df_markers, neurons_df_markers,
  microglia_markers, neuron_markers,
  verbose = TRUE
)

# Extract found markers
microglia_markers_found <- marker_results$microglia
neuron_markers_found <- marker_results$neuron

cat("\nSummary: Found", length(microglia_markers_found), "microglia markers and", 
    length(neuron_markers_found), "neuron markers\n")

#------------------------------------------------------------------------------
# Part 4: Create expression table and visualize results
#------------------------------------------------------------------------------
cat("\nCalculating mean expression values for each marker...\n")

# Calculate mean expression for each marker in each cell type
mean_expression_data <- create_mean_expression_table(
  microglia_df_markers, neurons_df_markers,
  microglia_markers_found, neuron_markers_found,
  microglia_vehicle_columns, neuron_vehicle_columns
)

# Save expression data to Excel
expression_data_file <- file.path(output_dir, "cell_type_marker_expression.xlsx")
writexl::write_xlsx(mean_expression_data, expression_data_file)
cat("Expression data saved to:", expression_data_file, "\n")

# Create expression heatmap
cat("\nCreating expression heatmap...\n")
expression_heatmap_file <- file.path(output_dir, "cell_type_marker_expression_heatmap.png")
p_expression <- create_expression_heatmap(
  mean_expression_data,
  title_suffix = "- Cell-Specific Markers",
  use_scaled = FALSE,
  save_file = expression_heatmap_file
)

# Create fold change heatmap
cat("\nCreating fold change heatmap...\n")
fold_change_heatmap_file <- file.path(output_dir, "cell_type_marker_fold_change_heatmap.png")
p_fold_change <- create_fold_change_heatmap(
  mean_expression_data,
  save_file = fold_change_heatmap_file
)

#------------------------------------------------------------------------------
# Part 5: Create intersection analysis (Proteins in both datasets)
#------------------------------------------------------------------------------
cat("\nPerforming intersection analysis...\n")

# Filter markers to only those present in both datasets
microglia_markers_intersection <- intersect(microglia_markers_found, proteins_in_both)
neuron_markers_intersection <- intersect(neuron_markers_found, proteins_in_both)

cat("Microglia markers in intersection:", length(microglia_markers_intersection), "\n")
cat("Neuron markers in intersection:", length(neuron_markers_intersection), "\n")

# Calculate expression for intersection markers
mean_expression_intersection <- create_mean_expression_table(
  microglia_df_markers[proteins_in_both, ], neurons_df_markers[proteins_in_both, ],
  microglia_markers_intersection, neuron_markers_intersection,
  microglia_vehicle_columns, neuron_vehicle_columns
)

# Save intersection data
intersection_data_file <- file.path(output_dir, "cell_type_marker_intersection.xlsx")
writexl::write_xlsx(mean_expression_intersection, intersection_data_file)
cat("Intersection data saved to:", intersection_data_file, "\n")

# Create intersection heatmaps if enough markers are found
if (nrow(mean_expression_intersection) > 3) {
  cat("\nCreating intersection heatmaps...\n")
  
  # Expression heatmap for intersection markers
  intersection_heatmap_file <- file.path(output_dir, "cell_type_marker_intersection_heatmap.png")
  p_intersection <- create_expression_heatmap(
    mean_expression_intersection,
    title_suffix = "- Markers in Both Datasets",
    use_scaled = FALSE,
    save_file = intersection_heatmap_file
  )
  
  # Fold change heatmap for intersection markers
  intersection_fc_heatmap_file <- file.path(output_dir, "cell_type_marker_intersection_fc_heatmap.png")
  p_intersection_fc <- create_fold_change_heatmap(
    mean_expression_intersection,
    save_file = intersection_fc_heatmap_file
  )
}

#------------------------------------------------------------------------------
# Part 6: Save a marker list summary for reference
#------------------------------------------------------------------------------
cat("\nGenerating marker summary report...\n")

# Create a summary data frame of all found markers
marker_summary <- data.frame(
  Marker = c(microglia_markers_found, neuron_markers_found),
  Cell_Type = c(rep("Microglia", length(microglia_markers_found)), 
                rep("Neuron", length(neuron_markers_found))),
  In_Both_Datasets = c(microglia_markers_found %in% proteins_in_both,
                       neuron_markers_found %in% proteins_in_both),
  stringsAsFactors = FALSE
)

# Save the summary
summary_file <- file.path(output_dir, "cell_type_marker_summary.xlsx")
writexl::write_xlsx(marker_summary, summary_file)
cat("Marker summary saved to:", summary_file, "\n")

# Create a markdown report with textual summary
report_file <- file.path(output_dir, "cell_type_marker_report.md")
cat("# Cell Type Marker Validation Report\n\n", file = report_file)
cat("## Dataset Summary\n\n", file = report_file, append = TRUE)
cat(paste("- Microglia dataset: ", nrow(microglia_df_markers), " proteins, ", 
          ncol(microglia_df_markers), " samples\n", sep = ""), 
    file = report_file, append = TRUE)
cat(paste("- Neuron dataset: ", nrow(neurons_df_markers), " proteins, ", 
          ncol(neurons_df_markers), " samples\n", sep = ""), 
    file = report_file, append = TRUE)
cat(paste("- Proteins in both datasets: ", length(proteins_in_both), "\n", sep = ""), 
    file = report_file, append = TRUE)

cat("\n## Marker Summary\n\n", file = report_file, append = TRUE)
cat(paste("- Microglia markers found: ", length(microglia_markers_found), 
          " out of ", length(microglia_markers), " tested\n", sep = ""), 
    file = report_file, append = TRUE)
cat(paste("- Neuron markers found: ", length(neuron_markers_found), 
          " out of ", length(neuron_markers), " tested\n", sep = ""), 
    file = report_file, append = TRUE)
cat(paste("- Microglia markers in both datasets: ", length(microglia_markers_intersection), "\n", sep = ""), 
    file = report_file, append = TRUE)
cat(paste("- Neuron markers in both datasets: ", length(neuron_markers_intersection), "\n", sep = ""), 
    file = report_file, append = TRUE)

cat("\n## Top Markers by Cell Type Specificity\n\n", file = report_file, append = TRUE)

# Add top 5 most specific markers for each cell type
top_microglia <- mean_expression_data[mean_expression_data$Marker_Type == "Microglia Marker", ]
top_microglia <- top_microglia[order(top_microglia$Log2FC), ]
top_neuron <- mean_expression_data[mean_expression_data$Marker_Type == "Neuron Marker", ]
top_neuron <- top_neuron[order(top_neuron$Log2FC, decreasing = TRUE), ]

cat("### Most Specific Microglia Markers\n\n", file = report_file, append = TRUE)
cat("| Marker | Microglia Expression | Neuron Expression | Log2 Fold Change |\n", 
    file = report_file, append = TRUE)
cat("|--------|---------------------|-------------------|------------------|\n", 
    file = report_file, append = TRUE)
for (i in 1:min(5, nrow(top_microglia))) {
  cat(paste("| ", rownames(top_microglia)[i], " | ", 
            round(top_microglia$Microglia[i], 2), " | ", 
            round(top_microglia$Neuron[i], 2), " | ", 
            round(top_microglia$Log2FC[i], 2), " |\n", sep = ""), 
      file = report_file, append = TRUE)
}

cat("\n### Most Specific Neuron Markers\n\n", file = report_file, append = TRUE)
cat("| Marker | Microglia Expression | Neuron Expression | Log2 Fold Change |\n", 
    file = report_file, append = TRUE)
cat("|--------|---------------------|-------------------|------------------|\n", 
    file = report_file, append = TRUE)
for (i in 1:min(5, nrow(top_neuron))) {
  cat(paste("| ", rownames(top_neuron)[i], " | ", 
            round(top_neuron$Microglia[i], 2), " | ", 
            round(top_neuron$Neuron[i], 2), " | ", 
            round(top_neuron$Log2FC[i], 2), " |\n", sep = ""), 
      file = report_file, append = TRUE)
}

# Add execution information to report
cat("\n## Analysis Information\n\n", file = report_file, append = TRUE)
cat(paste("- Analysis date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = ""), 
    file = report_file, append = TRUE)
cat(paste("- Vehicle columns used in microglia: ", paste(microglia_vehicle_columns, collapse = ", "), "\n", sep = ""), 
    file = report_file, append = TRUE)
cat(paste("- Vehicle columns used in neurons: ", paste(neuron_vehicle_columns, collapse = ", "), "\n", sep = ""), 
    file = report_file, append = TRUE)

cat("\nAnalysis complete. Results saved to:", output_dir, "\n")
