#!/usr/bin/env Rscript
#===============================================================================
# Proteomics Data Batch Correction and Merging Script (Post-Duplicate Resolution)
# 
# This script performs:
# 1. Loading of filtered datasets
# 2. Duplicate resolution for protein entries
# 3. Proper batch identification for all neuron datasets
# 4. Merging of neuron datasets along common proteins with tracking
# 5. Batch correction of microglia and neuron data using:
#    - ComBat
#    - Standard SVA+Limma
#    - Manual SVA
# 6. QC plots and protein retention validation
#===============================================================================

#-------------------------------------------------------------------------------
# 1. Load required libraries
#-------------------------------------------------------------------------------
# Load required packages 
# source("install_packages_updated.R") 

library(sva)        # For batch correction with SVA
library(limma)      # For batch correction with removeBatchEffect
library(ggplot2)    # For plotting
library(writexl)    # For Excel output
library(dplyr)      # For data manipulation
library(stringr)    # For string manipulation
library(factoextra) # For PCA visualization
library(gridExtra)  # For arranging multiple plots

#-------------------------------------------------------------------------------
# 2. Load filtered datasets 
#-------------------------------------------------------------------------------
# Set timestamp for the output files to prevent overwriting
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_suffix <- "post_dupl_res" # Add this suffix to all output files

# Load the filtered datasets
load("QR2_contaminant_filtered_data_no_neuron_20250412_231227.RData")
# This should load microglia_data_filtered, neuron_data_batch1_filtered, and neuron_data_batch2_filtered

# Function to resolve duplicates by taking the mean
resolve_duplicates <- function(dataset) {
  proteins <- rownames(dataset)
  unique_proteins <- unique(proteins)
  result <- matrix(NA, nrow=length(unique_proteins), ncol=ncol(dataset))
  rownames(result) <- unique_proteins
  colnames(result) <- colnames(dataset)
  
  for(p in unique_proteins) {
    rows <- which(proteins == p)
    if(length(rows) == 1) {
      result[p,] <- dataset[rows,]
    } else {
      result[p,] <- colMeans(dataset[rows,,drop=FALSE])
    }
  }
  return(result)
}

# Apply the function to resolve duplicates
microglia_clean <- resolve_duplicates(microglia_data_filtered)
neuron1_clean <- resolve_duplicates(neuron_data_batch1_filtered)
neuron2_clean <- resolve_duplicates(neuron_data_batch2_filtered)

# Check protein counts before and after resolution
cat("\n===========================================================\n")
cat("DUPLICATE RESOLUTION VALIDATION\n")
cat("===========================================================\n")
cat("Microglia - Original unique proteins:", length(unique(rownames(microglia_data_filtered))), 
    ", After resolution:", nrow(microglia_clean), "\n")
cat("Neuron1 - Original unique proteins:", length(unique(rownames(neuron_data_batch1_filtered))), 
    ", After resolution:", nrow(neuron1_clean), "\n")
cat("Neuron2 - Original unique proteins:", length(unique(rownames(neuron_data_batch2_filtered))), 
    ", After resolution:", nrow(neuron2_clean), "\n")

# Confirm no duplicates remain
cat("\nDuplicates remaining - Microglia:", sum(duplicated(rownames(microglia_clean))), "\n")
cat("Duplicates remaining - Neuron1:", sum(duplicated(rownames(neuron1_clean))), "\n")
cat("Duplicates remaining - Neuron2:", sum(duplicated(rownames(neuron2_clean))), "\n")
cat("===========================================================\n\n")

# Create directories for output files if they don't exist
dir.create("resultssss", showWarnings = FALSE)
dir.create("resultssss/data", showWarnings = FALSE)
dir.create("resultssss/plots", showWarnings = FALSE)
dir.create("resultssss/tracking", showWarnings = FALSE)

#-------------------------------------------------------------------------------
# 3. Define the manual SVA function (MISSING FROM ORIGINAL)
#-------------------------------------------------------------------------------
# Create a function for manual SVA
manual_sva <- function(dat, mod, mod0, n.sv = NULL) {
  n <- dim(dat)[1]
  m <- dim(dat)[2]
  
  # Calculate residuals
  H <- mod %*% solve(t(mod) %*% mod) %*% t(mod)
  res <- dat - H %*% dat
  
  # Calculate residuals for null model
  H0 <- mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0)
  res0 <- dat - H0 %*% dat
  
  # Perform SVD on residuals
  svd_res <- svd(res0)
  
  # Determine number of significant SVs if not provided
  if (is.null(n.sv)) {
    n.sv <- 2  # Default to 2 surrogate variables
    print(paste("Using", n.sv, "surrogate variables for manual SVA"))
  }
  
  # Extract surrogate variables
  sv <- svd_res$u[, 1:n.sv, drop = FALSE]
  
  # Calculate adjustment
  X <- cbind(mod, sv)
  Hat <- X %*% solve(t(X) %*% X) %*% t(X)
  adj <- Hat %*% dat - H %*% dat
  
  # Return results
  return(list(sv = sv, n.sv = n.sv, adjustment = adj))
}

#-------------------------------------------------------------------------------
# 4. Automatic batch identification for neuron datasets
#-------------------------------------------------------------------------------
# AUTOMATIC BATCH ASSIGNMENT (like microglia)
# Extract batch information automatically from sample names

cat("\n--- AUTOMATIC BATCH ASSIGNMENT ---\n")

# For neuron1 dataset (assuming all are batch 1)
neuron1_batches <- rep(1, ncol(neuron1_clean))

# For neuron2 dataset - extract batch number from sample names
# Pattern: Neuron_X_CellLine_Treatment where X is the batch number
neuron2_batches <- as.numeric(gsub("Neuron_([0-9])_.*", "\\1", colnames(neuron2_clean)))

# Combine all neuron batch information
all_neuron_batches <- c(neuron1_batches, neuron2_batches)

# Verify extraction worked
cat("Neuron1 batch assignment:\n")
print(table(neuron1_batches))
cat("Neuron2 batch assignment:\n")
print(table(neuron2_batches))
cat("Combined neuron batch assignment:\n")
print(table(all_neuron_batches))

all_neuron_samples <- c(colnames(neuron1_clean), colnames(neuron2_clean))

# Create a batch assignment summary
batch_summary <- data.frame(
  Sample = all_neuron_samples,
  Dataset = c(rep("Neuron1", ncol(neuron1_clean)), 
              rep("Neuron2", ncol(neuron2_clean))),
  Batch = all_neuron_batches,
  stringsAsFactors = FALSE
)

# Save batch assignment for review
write_xlsx(batch_summary, 
           file.path("resultssss", "tracking", paste0("neuron_batch_assignment_", output_suffix, "_", timestamp, ".xlsx")))

# Print batch distribution
cat("\nBatch distribution:\n")
print(table(batch_summary$Dataset, batch_summary$Batch))

cat("===========================================================\n\n")

#-------------------------------------------------------------------------------
# 5. Data inspection and tracking preparation
#-------------------------------------------------------------------------------
# Function to track protein counts and save to Excel
track_proteins <- function(data_list, step_name, output_file) {
  tracking_df <- data.frame(
    Step = character(),
    Dataset = character(),
    ProteinCount = integer(),
    stringsAsFactors = FALSE
  )
  
  for (dataset_name in names(data_list)) {
    dataset <- data_list[[dataset_name]]
    protein_count <- nrow(dataset)
    tracking_df <- rbind(tracking_df, data.frame(
      Step = step_name,
      Dataset = dataset_name,
      ProteinCount = protein_count,
      stringsAsFactors = FALSE
    ))
  }
  
  # Append to existing tracking file if it exists, otherwise create new
  if (file.exists(output_file)) {
    existing_df <- read.csv(output_file)
    tracking_df <- rbind(existing_df, tracking_df)
  }
  
  write.csv(tracking_df, output_file, row.names = FALSE)
  return(tracking_df)
}

# Initial tracking of protein counts
initial_data <- list(
  "microglia" = microglia_clean,
  "neuron_batch1" = neuron1_clean,
  "neuron_batch2" = neuron2_clean
)
track_proteins(initial_data, "Post-Duplicate Resolution", 
               file.path("resultssss", "tracking", paste0("protein_tracking_", output_suffix, "_", timestamp, ".csv")))

#-------------------------------------------------------------------------------
# 6. Merge neuron datasets along common proteins with tracking
#-------------------------------------------------------------------------------
# Find common proteins between neuron datasets
common_proteins_neurons <- intersect(rownames(neuron1_clean), rownames(neuron2_clean))

cat("Protein overlap analysis:\n")
cat("Neuron1 proteins:", nrow(neuron1_clean), "\n")
cat("Neuron2 proteins:", nrow(neuron2_clean), "\n")
cat("Common proteins:", length(common_proteins_neurons), "\n")

# Track proteins dropped in neuron merging
dropped_proteins_batch1 <- setdiff(rownames(neuron1_clean), common_proteins_neurons)
dropped_proteins_batch2 <- setdiff(rownames(neuron2_clean), common_proteins_neurons)

cat("Proteins dropped from Neuron1:", length(dropped_proteins_batch1), "\n")
cat("Proteins dropped from Neuron2:", length(dropped_proteins_batch2), "\n\n")

# Save lists of dropped proteins
write.csv(data.frame(Protein = dropped_proteins_batch1), 
          file.path("resultssss", "tracking", paste0("neuron_batch1_dropped_proteins_", output_suffix, "_", timestamp, ".csv")), 
          row.names = FALSE)
write.csv(data.frame(Protein = dropped_proteins_batch2), 
          file.path("resultssss", "tracking", paste0("neuron_batch2_dropped_proteins_", output_suffix, "_", timestamp, ".csv")), 
          row.names = FALSE)

# Subset to common proteins for merging
neuron_batch1_common <- neuron1_clean[common_proteins_neurons, ]
neuron_batch2_common <- neuron2_clean[common_proteins_neurons, ]

# Merge the neuron datasets
neurons_combined <- cbind(neuron_batch1_common, neuron_batch2_common)

# Extract experimental factors for all samples
extract_factors <- function(sample_names) {
  # Adjust these patterns based on your naming convention
  cell_line <- gsub(".*_(Ctrl|sAD|Bio|UK)_.*", "\\1", sample_names)
  treatment <- gsub(".*_(V|I)[0-9]*$", "\\1", sample_names)
  
  # Handle cases where pattern doesn't match
  cell_line[!grepl("_(Ctrl|sAD|Bio|UK)_", sample_names)] <- "Unknown"
  treatment[!grepl("_(V|I)[0-9]*$", sample_names)] <- "Unknown"
  
  return(list(cell_line = cell_line, treatment = treatment))
}

# Extract factors for combined neuron data
neuron_factors <- extract_factors(colnames(neurons_combined))
neuron_combined_cell_line <- neuron_factors$cell_line
neuron_combined_treatment <- neuron_factors$treatment

#-------------------------------------------------------------------------------
# 7. Prepare batch and condition information for microglia
#-------------------------------------------------------------------------------
# Extract batch information and treatment from column names for microglia
microglia_batch <- as.numeric(gsub(".*_([0-9])_.*", "\\1", colnames(microglia_clean)))
microglia_factors <- extract_factors(colnames(microglia_clean))
microglia_cell_line <- microglia_factors$cell_line
microglia_treatment <- microglia_factors$treatment

# Check for confounding in microglia
cat("Microglia batch vs cell line table:\n")
print(table(microglia_batch, microglia_cell_line))
cat("Microglia batch vs treatment table:\n")
print(table(microglia_batch, microglia_treatment))

# Check for confounding in neurons with updated batch assignments
cat("Neuron batch vs cell line table:\n")
print(table(all_neuron_batches, neuron_combined_cell_line))
cat("Neuron batch vs treatment table:\n")
print(table(all_neuron_batches, neuron_combined_treatment))

# Create and save sample information
microglia_sample_info <- data.frame(
  Sample = colnames(microglia_clean),
  Batch = microglia_batch,
  CellLine = microglia_cell_line,
  Treatment = microglia_treatment
)

neuron_sample_info <- data.frame(
  Sample = colnames(neurons_combined),
  Batch = all_neuron_batches,
  CellLine = neuron_combined_cell_line,
  Treatment = neuron_combined_treatment
)

write_xlsx(list("Microglia" = microglia_sample_info, 
                "Neurons" = neuron_sample_info),
           file.path("resultssss", "tracking", paste0("sample_info_", output_suffix, "_", timestamp, ".xlsx")))

#-------------------------------------------------------------------------------
# 8. Apply ComBat batch correction to microglia data
#-------------------------------------------------------------------------------
# Convert to matrix for ComBat
microglia_matrix <- as.matrix(microglia_clean)

# Apply ComBat
print("Applying ComBat to microglia data...")
microglia_combat_corrected <- ComBat(
  dat = microglia_matrix,
  batch = microglia_batch,
  ref.batch = 2,  # Use batch 2 as reference (assuming it's most balanced as per your comment)
  par.prior = TRUE,
  prior.plots = FALSE
)

# Track protein counts after ComBat
combat_data <- list("microglia_combat_corrected" = microglia_combat_corrected)
track_proteins(combat_data, "After ComBat Correction (Post-Duplicate Resolution)", 
               file.path("resultssss", "tracking", paste0("protein_tracking_", output_suffix, "_", timestamp, ".csv")))

# Verify that no proteins were dropped
identical_rownames_combat <- identical(sort(rownames(microglia_clean)), sort(rownames(microglia_combat_corrected)))
print(paste("ComBat correction retained all proteins for microglia:", identical_rownames_combat))

# Save ComBat-corrected data
microglia_combat_for_excel <- as.data.frame(microglia_combat_corrected)
microglia_combat_for_excel <- cbind(Protein = rownames(microglia_combat_for_excel), microglia_combat_for_excel)
write_xlsx(microglia_combat_for_excel, 
           file.path("resultssss", "data", paste0("microglia_combat_corrected_", output_suffix, "_", timestamp, ".xlsx")))

#-------------------------------------------------------------------------------
# 9. Apply Standard SVA + Limma batch correction to microglia data
#-------------------------------------------------------------------------------
# Prepare model matrices for SVA
print("Applying Standard SVA+Limma to microglia data...")
microglia_cell_line_factor <- as.factor(microglia_cell_line)
microglia_treatment_factor <- as.factor(microglia_treatment)
microglia_batch_factor <- as.factor(microglia_batch)

# Create model matrices
mod <- model.matrix(~ 0 + microglia_cell_line_factor + microglia_treatment_factor)
mod_with_batch <- model.matrix(~ 0 + microglia_cell_line_factor + microglia_treatment_factor + microglia_batch_factor)

# Estimate the number of surrogate variables
n_sv <- num.sv(microglia_matrix, mod, method = "be")
print(paste("Estimated number of surrogate variables for microglia:", n_sv))

# If n_sv is 0, set it to a small number
if (n_sv == 0) {
  n_sv <- 2
  print("Setting n_sv to 2 since estimated value was 0")
}

# Create the SVA model
svobj <- try(sva(microglia_matrix, mod, mod_with_batch, n.sv = n_sv))

# If SVA fails, try a simplified approach
if (inherits(svobj, "try-error")) {
  print("Standard SVA failed for microglia, trying alternative approach with fixed number of SVs")
  n_sv <- 2  # Fixed number
  
  # Create a simplified design
  mod_simple <- model.matrix(~ microglia_cell_line_factor + microglia_treatment_factor)
  mod0_simple <- model.matrix(~ 1, data = data.frame(x = rep(1, length(microglia_cell_line_factor))))
  
  # Use svaseq instead which can be more stable
  svobj <- svaseq(microglia_matrix, mod_simple, mod0_simple, n.sv = n_sv)
}

# Get surrogate variables
sv <- svobj$sv

# Apply limma's removeBatchEffect
microglia_std_sva_corrected <- removeBatchEffect(
  microglia_matrix,
  covariates = sv,
  design = mod
)

# Track protein counts after Standard SVA+Limma
std_sva_data <- list("microglia_std_sva_corrected" = microglia_std_sva_corrected)
track_proteins(std_sva_data, "After Standard SVA+Limma Correction (Post-Duplicate Resolution)", 
               file.path("resultssss", "tracking", paste0("protein_tracking_", output_suffix, "_", timestamp, ".csv")))

# Verify that no proteins were dropped
identical_rownames_std_sva <- identical(sort(rownames(microglia_clean)), sort(rownames(microglia_std_sva_corrected)))
print(paste("Standard SVA+Limma correction retained all proteins for microglia:", identical_rownames_std_sva))

# Save Standard SVA+Limma-corrected data
microglia_std_sva_for_excel <- as.data.frame(microglia_std_sva_corrected)
microglia_std_sva_for_excel <- cbind(Protein = rownames(microglia_std_sva_for_excel), microglia_std_sva_for_excel)
write_xlsx(microglia_std_sva_for_excel, 
           file.path("resultssss", "data", paste0("microglia_std_sva_corrected_", output_suffix, "_", timestamp, ".xlsx")))

#-------------------------------------------------------------------------------
# 10. Apply Manual SVA batch correction to microglia data
#-------------------------------------------------------------------------------
# Prepare for manual SVA on microglia
print("Applying Manual SVA to microglia data...")

# Transpose the data for Manual SVA
microglia_matrix_t <- t(microglia_matrix)

# Create design matrices for microglia
microglia_design_df <- data.frame(
  cell_line = microglia_cell_line,
  treatment = microglia_treatment
)

microglia_mod <- model.matrix(~ cell_line + treatment, data = microglia_design_df)
microglia_mod0 <- matrix(1, nrow = nrow(microglia_matrix_t), ncol = 1)

# Apply manual SVA to microglia
microglia_manual_svs <- manual_sva(microglia_matrix_t, microglia_mod, microglia_mod0, n.sv = 2)

# Apply correction
microglia_manual_corrected_t <- microglia_matrix_t - microglia_manual_svs$adjustment
microglia_manual_corrected <- t(microglia_manual_corrected_t)

# Track protein counts after Manual SVA
microglia_manual_sva_data <- list("microglia_manual_corrected" = microglia_manual_corrected)
track_proteins(microglia_manual_sva_data, "After Manual SVA Correction (Post-Duplicate Resolution)", 
               file.path("resultssss", "tracking", paste0("protein_tracking_", output_suffix, "_", timestamp, ".csv")))

# Verify that no proteins were dropped
identical_rownames_microglia_manual <- identical(sort(rownames(microglia_clean)), 
                                                 sort(rownames(microglia_manual_corrected)))
print(paste("Manual SVA correction retained all proteins for microglia:", identical_rownames_microglia_manual))

# Save Manual SVA-corrected microglia data
microglia_manual_for_excel <- as.data.frame(microglia_manual_corrected)
microglia_manual_for_excel <- cbind(Protein = rownames(microglia_manual_for_excel), microglia_manual_for_excel)
write_xlsx(microglia_manual_for_excel, 
           file.path("resultssss", "data", paste0("microglia_manual_sva_corrected_", output_suffix, "_", timestamp, ".xlsx")))

#-------------------------------------------------------------------------------
# 11. Apply batch correction methods to neuron data
#-------------------------------------------------------------------------------
# Prepare neuron data matrix
neurons_combined_matrix <- as.matrix(neurons_combined)

#-------------------------------------------------------------------------------
# 11a. Apply Manual SVA batch correction to neuron data
#-------------------------------------------------------------------------------
# Prepare for manual SVA on neurons
print("Applying Manual SVA to neuron data...")

# Transpose the data for Manual SVA
neurons_combined_t <- t(neurons_combined_matrix)

# Create design matrices for neurons
neuron_design_df <- data.frame(
  cell_line = neuron_combined_cell_line,
  treatment = neuron_combined_treatment
)

neuron_mod <- model.matrix(~ cell_line + treatment, data = neuron_design_df)
neuron_mod0 <- matrix(1, nrow = nrow(neurons_combined_t), ncol = 1)

# Apply manual SVA to neurons
neuron_manual_svs <- manual_sva(neurons_combined_t, neuron_mod, neuron_mod0, n.sv = 2)

# Apply correction
neurons_manual_corrected_t <- neurons_combined_t - neuron_manual_svs$adjustment
neurons_manual_corrected <- t(neurons_manual_corrected_t)

# Track protein counts after Manual SVA
manual_sva_data <- list("neurons_manual_corrected" = neurons_manual_corrected)
track_proteins(manual_sva_data, "After Manual SVA Correction (Post-Duplicate Resolution)", 
               file.path("resultssss", "tracking", paste0("protein_tracking_", output_suffix, "_", timestamp, ".csv")))

# Verify that no proteins were dropped
identical_rownames_manual_sva <- identical(sort(rownames(neurons_combined)), sort(rownames(neurons_manual_corrected)))
print(paste("Manual SVA correction retained all proteins for neurons:", identical_rownames_manual_sva))

# Save Manual SVA-corrected neuron data
neurons_manual_for_excel <- as.data.frame(neurons_manual_corrected)
neurons_manual_for_excel <- cbind(Protein = rownames(neurons_manual_for_excel), neurons_manual_for_excel)
write_xlsx(neurons_manual_for_excel, 
           file.path("resultssss", "data", paste0("neurons_manual_sva_corrected_", output_suffix, "_", timestamp, ".xlsx")))

#-------------------------------------------------------------------------------
# 11b. Apply ComBat to neuron data
#-------------------------------------------------------------------------------
cat("Applying ComBat to neuron data with multiple batches...\n")
neurons_combat_corrected <- try(
  ComBat(
    dat = neurons_combined_matrix,
    batch = all_neuron_batches,  # Use the updated batch assignments
    ref.batch = 1,
    par.prior = FALSE,
    prior.plots = FALSE
  )
)

if (inherits(neurons_combat_corrected, "try-error")) {
  cat("ComBat failed for neurons. This might be due to:\n")
  cat("1. Insufficient samples per batch\n")
  cat("2. Perfect confounding between batch and biological factors\n")
  cat("3. Incorrect batch assignments\n")
  cat("Please check your batch assignments and sample distribution.\n\n")
} else {
  # Track protein counts after ComBat
  neurons_combat_data <- list("neurons_combat_corrected" = neurons_combat_corrected)
  track_proteins(neurons_combat_data, "After ComBat Correction (Post-Duplicate Resolution)", 
                 file.path("resultssss", "tracking", paste0("protein_tracking_", output_suffix, "_", timestamp, ".csv")))
  
  # Verify that no proteins were dropped
  identical_rownames_neurons_combat <- identical(sort(rownames(neurons_combined)), 
                                                 sort(rownames(neurons_combat_corrected)))
  print(paste("ComBat correction retained all proteins for neurons:", identical_rownames_neurons_combat))
  
  # Save ComBat-corrected neuron data
  neurons_combat_for_excel <- as.data.frame(neurons_combat_corrected)
  neurons_combat_for_excel <- cbind(Protein = rownames(neurons_combat_for_excel), neurons_combat_for_excel)
  write_xlsx(neurons_combat_for_excel, 
             file.path("resultssss", "data", paste0("neurons_combat_corrected_", output_suffix, "_", timestamp, ".xlsx")))
}

#-------------------------------------------------------------------------------
# 11c. Apply Standard SVA + Limma to neuron data
#-------------------------------------------------------------------------------
print("Applying Standard SVA+Limma to neuron data...")

# Create model matrices for Standard SVA
neuron_cell_line_factor <- as.factor(neuron_combined_cell_line)
neuron_treatment_factor <- as.factor(neuron_combined_treatment)
neuron_batch_factor <- as.factor(all_neuron_batches)  # FIXED: use correct variable name

# Create model matrices
neuron_mod_std <- model.matrix(~ 0 + neuron_cell_line_factor + neuron_treatment_factor)
neuron_mod_with_batch <- model.matrix(~ 0 + neuron_cell_line_factor + neuron_treatment_factor + neuron_batch_factor)

# Estimate the number of surrogate variables
neuron_n_sv <- num.sv(neurons_combined_matrix, neuron_mod_std, method = "be")
print(paste("Estimated number of surrogate variables for neurons:", neuron_n_sv))

# If n_sv is 0, set it to a small number
if (neuron_n_sv == 0) {
  neuron_n_sv <- 2
  print("Setting neuron_n_sv to 2 since estimated value was 0")
}

# Create the SVA model
neuron_svobj <- try(sva(neurons_combined_matrix, neuron_mod_std, neuron_mod_with_batch, n.sv = neuron_n_sv))

# If SVA fails, try a simplified approach
if (inherits(neuron_svobj, "try-error")) {
  print("Standard SVA failed for neurons, trying alternative approach with fixed number of SVs")
  neuron_n_sv <- 2  # Fixed number
  
  # Create a simplified design
  neuron_mod_simple <- model.matrix(~ neuron_cell_line_factor + neuron_treatment_factor)
  neuron_mod0_simple <- model.matrix(~ 1, data = data.frame(x = rep(1, length(neuron_cell_line_factor))))
  
  # Use svaseq instead which can be more stable
  neuron_svobj <- svaseq(neurons_combined_matrix, neuron_mod_simple, neuron_mod0_simple, n.sv = neuron_n_sv)
}

# Get surrogate variables
neuron_sv <- neuron_svobj$sv

# Apply limma's removeBatchEffect using the SVA variables
neurons_std_sva_corrected <- removeBatchEffect(
  neurons_combined_matrix,
  covariates = neuron_sv,
  design = neuron_mod_std
)

# Track protein counts after Standard SVA+Limma
neurons_std_sva_data <- list("neurons_std_sva_corrected" = neurons_std_sva_corrected)
track_proteins(neurons_std_sva_data, "After Standard SVA+Limma Correction (Post-Duplicate Resolution)", 
               file.path("resultssss", "tracking", paste0("protein_tracking_", output_suffix, "_", timestamp, ".csv")))

# Verify that no proteins were dropped
identical_rownames_neurons_std_sva <- identical(sort(rownames(neurons_combined)), 
                                                sort(rownames(neurons_std_sva_corrected)))
print(paste("Standard SVA+Limma correction retained all proteins for neurons:", identical_rownames_neurons_std_sva))

# Save Standard SVA+Limma-corrected neuron data
neurons_std_sva_for_excel <- as.data.frame(neurons_std_sva_corrected)
neurons_std_sva_for_excel <- cbind(Protein = rownames(neurons_std_sva_for_excel), neurons_std_sva_for_excel)
write_xlsx(neurons_std_sva_for_excel, 
           file.path("resultssss", "data", paste0("neurons_std_sva_corrected_", output_suffix, "_", timestamp, ".xlsx")))

#-------------------------------------------------------------------------------
# 12. Perform PCA and create plots for microglia data
#-------------------------------------------------------------------------------
# PCA for raw microglia data
microglia_raw_pca <- prcomp(t(microglia_clean), scale. = TRUE)
microglia_raw_pca_data <- as.data.frame(microglia_raw_pca$x)
microglia_raw_pca_data$Batch <- as.factor(microglia_batch)
microglia_raw_pca_data$CellLine <- as.factor(microglia_cell_line)
microglia_raw_pca_data$Treatment <- as.factor(microglia_treatment)

# PCA for ComBat-corrected microglia data
microglia_combat_pca <- prcomp(t(microglia_combat_corrected), scale. = TRUE)
microglia_combat_pca_data <- as.data.frame(microglia_combat_pca$x)
microglia_combat_pca_data$Batch <- as.factor(microglia_batch)
microglia_combat_pca_data$CellLine <- as.factor(microglia_cell_line)
microglia_combat_pca_data$Treatment <- as.factor(microglia_treatment)

# PCA for Standard SVA+Limma-corrected microglia data
microglia_std_sva_pca <- prcomp(t(microglia_std_sva_corrected), scale. = TRUE)
microglia_std_sva_pca_data <- as.data.frame(microglia_std_sva_pca$x)
microglia_std_sva_pca_data$Batch <- as.factor(microglia_batch)
microglia_std_sva_pca_data$CellLine <- as.factor(microglia_cell_line)
microglia_std_sva_pca_data$Treatment <- as.factor(microglia_treatment)

# PCA for Manual SVA-corrected microglia data
microglia_manual_sva_pca <- prcomp(t(microglia_manual_corrected), scale. = TRUE)
microglia_manual_sva_pca_data <- as.data.frame(microglia_manual_sva_pca$x)
microglia_manual_sva_pca_data$Batch <- as.factor(microglia_batch)
microglia_manual_sva_pca_data$CellLine <- as.factor(microglia_cell_line)
microglia_manual_sva_pca_data$Treatment <- as.factor(microglia_treatment)

# Calculate variance explained
microglia_raw_var <- summary(microglia_raw_pca)$importance[2, 1:5] * 100
microglia_combat_var <- summary(microglia_combat_pca)$importance[2, 1:5] * 100
microglia_std_sva_var <- summary(microglia_std_sva_pca)$importance[2, 1:5] * 100
microglia_manual_sva_var <- summary(microglia_manual_sva_pca)$importance[2, 1:5] * 100

# Create a data frame for variance comparison
variance_comparison_microglia <- data.frame(
  PC = paste0("PC", 1:5),
  Raw = microglia_raw_var,
  ComBat = microglia_combat_var,
  Standard_SVA_Limma = microglia_std_sva_var,
  Manual_SVA = microglia_manual_sva_var
)

# Save variance comparison
write_xlsx(variance_comparison_microglia, 
           file.path("resultssss", "data", paste0("microglia_variance_comparison_", output_suffix, "_", timestamp, ".xlsx")))

# Create PCA plots for microglia
# Batch effect visualization
p1_raw <- ggplot(microglia_raw_pca_data, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Raw Microglia Data (Post-Duplicate Resolution)", color = "Batch") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p1_combat <- ggplot(microglia_combat_pca_data, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "ComBat-Corrected Microglia (Post-Duplicate Resolution)", color = "Batch") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p1_std_sva <- ggplot(microglia_std_sva_pca_data, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Standard SVA+Limma-Corrected Microglia (Post-Duplicate Resolution)", color = "Batch") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p1_manual_sva <- ggplot(microglia_manual_sva_pca_data, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Manual SVA-Corrected Microglia (Post-Duplicate Resolution)", color = "Batch") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Combine plots (arrange in a 2x2 grid)
batch_comparison <- grid.arrange(p1_raw, p1_combat, p1_std_sva, p1_manual_sva, ncol = 2)

# Save the comparison plot
ggsave(
  file.path("resultssss", "plots", paste0("microglia_batch_comparison_", output_suffix, "_", timestamp, ".png")),
  batch_comparison,
  width = 12,
  height = 10,
  dpi = 300
)

# Cell line visualization for all corrected methods
p2_combat <- ggplot(microglia_combat_pca_data, aes(x = PC1, y = PC2, color = CellLine)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "ComBat-Corrected Microglia (Post-Duplicate Resolution)", color = "Cell Line") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p2_std_sva <- ggplot(microglia_std_sva_pca_data, aes(x = PC1, y = PC2, color = CellLine)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Standard SVA+Limma-Corrected Microglia (Post-Duplicate Resolution)", color = "Cell Line") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p2_manual_sva <- ggplot(microglia_manual_sva_pca_data, aes(x = PC1, y = PC2, color = CellLine)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Manual SVA-Corrected Microglia (Post-Duplicate Resolution)", color = "Cell Line") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Combine cell line plots
cellline_comparison <- grid.arrange(p2_combat, p2_std_sva, p2_manual_sva, ncol = 3)

# Save the cell line comparison plot
ggsave(
  file.path("resultssss", "plots", paste0("microglia_cellline_comparison_", output_suffix, "_", timestamp, ".png")),
  cellline_comparison,
  width = 15,
  height = 5,
  dpi = 300
)

# Treatment visualization for all corrected methods
p3_combat <- ggplot(microglia_combat_pca_data, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "ComBat-Corrected Microglia (Post-Duplicate Resolution)", color = "Treatment") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p3_std_sva <- ggplot(microglia_std_sva_pca_data, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Standard SVA+Limma-Corrected Microglia (Post-Duplicate Resolution)", color = "Treatment") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p3_manual_sva <- ggplot(microglia_manual_sva_pca_data, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Manual SVA-Corrected Microglia (Post-Duplicate Resolution)", color = "Treatment") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# Combine treatment plots
treatment_comparison <- grid.arrange(p3_combat, p3_std_sva, p3_manual_sva, ncol = 3)

# Save the treatment comparison plot
ggsave(
  file.path("resultssss", "plots", paste0("microglia_treatment_comparison_", output_suffix, "_", timestamp, ".png")),
  treatment_comparison,
  width = 15,
  height = 5,
  dpi = 300
)

#-------------------------------------------------------------------------------
# 13. Function for PCA and plots for neuron data
#-------------------------------------------------------------------------------
# Function to handle PCA for neuron data, allowing for different correction methods
neuron_pca_plots <- function(raw_data, corrected_data, correction_method, 
                             batch, cell_line, treatment, timestamp, output_suffix) {
  # PCA for raw data
  raw_pca <- prcomp(t(raw_data), scale. = TRUE)
  raw_pca_data <- as.data.frame(raw_pca$x)
  raw_pca_data$Batch <- as.factor(batch)
  raw_pca_data$CellLine <- as.factor(cell_line)
  raw_pca_data$Treatment <- as.factor(treatment)
  
  # PCA for corrected data
  corrected_pca <- prcomp(t(corrected_data), scale. = TRUE)
  corrected_pca_data <- as.data.frame(corrected_pca$x)
  corrected_pca_data$Batch <- as.factor(batch)
  corrected_pca_data$CellLine <- as.factor(cell_line)
  corrected_pca_data$Treatment <- as.factor(treatment)
  
  # Calculate variance explained
  raw_var <- summary(raw_pca)$importance[2, 1:5] * 100
  corrected_var <- summary(corrected_pca)$importance[2, 1:5] * 100
  
  # Create comparison dataframe
  variance_df <- data.frame(
    PC = paste0("PC", 1:5),
    Raw = raw_var,
    Corrected = corrected_var
  )
  colnames(variance_df)[3] <- correction_method
  
  # Save variance data
  write_xlsx(variance_df, 
             file.path("resultssss", "data", 
                       paste0("neurons_variance_", correction_method, "_", output_suffix, "_", timestamp, ".xlsx")))
  
  # Create PCA plots
  # Batch effect visualization
  p1_raw <- ggplot(raw_pca_data, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = "Raw Neuron Data (Post-Duplicate Resolution)", color = "Batch") +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  
  p1_corrected <- ggplot(corrected_pca_data, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = paste0(correction_method, "-Corrected Neurons (Post-Duplicate Resolution)"), color = "Batch") +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  
  # Combine plots
  batch_comparison <- grid.arrange(p1_raw, p1_corrected, ncol = 2)
  
  # Save the comparison plot
  ggsave(
    file.path("resultssss", "plots", paste0("neurons_batch_comparison_", correction_method, "_", output_suffix, "_", timestamp, ".png")),
    batch_comparison,
    width = 10,
    height = 5,
    dpi = 300
  )
  
  # Cell line and treatment plots only for corrected data
  p2_corrected <- ggplot(corrected_pca_data, aes(x = PC1, y = PC2, color = CellLine)) +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = paste0(correction_method, "-Corrected Neurons (Post-Duplicate Resolution)"), color = "Cell Line") +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  
  p3_corrected <- ggplot(corrected_pca_data, aes(x = PC1, y = PC2, color = Treatment)) +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = paste0(correction_method, "-Corrected Neurons (Post-Duplicate Resolution)"), color = "Treatment") +
    theme(plot.title = element_text(hjust = 0.5, size = 12))
  
  # Combine plots
  factor_plots <- grid.arrange(p2_corrected, p3_corrected, ncol = 2)
  
  # Save the plots
  ggsave(
    file.path("resultssss", "plots", paste0("neurons_factors_", correction_method, "_", output_suffix, "_", timestamp, ".png")),
    factor_plots,
    width = 10,
    height = 5,
    dpi = 300
  )
  
  # Return the variance dataframe for comparison
  return(variance_df)
}

# Generate PCA plots for Manual SVA-corrected neurons
manual_sva_variance <- neuron_pca_plots(
  raw_data = neurons_combined,
  corrected_data = neurons_manual_corrected,
  correction_method = "Manual_SVA",
  batch = all_neuron_batches,
  cell_line = neuron_combined_cell_line,
  treatment = neuron_combined_treatment,
  timestamp = timestamp,
  output_suffix = output_suffix
)

# Generate PCA plots for ComBat-corrected neurons if available
if (exists("neurons_combat_corrected") && !inherits(neurons_combat_corrected, "try-error")) {
  combat_variance <- neuron_pca_plots(
    raw_data = neurons_combined,
    corrected_data = neurons_combat_corrected,
    correction_method = "ComBat",
    batch = all_neuron_batches,
    cell_line = neuron_combined_cell_line,
    treatment = neuron_combined_treatment,
    timestamp = timestamp,
    output_suffix = output_suffix
  )
}

# Generate PCA plots for Standard SVA+Limma-corrected neurons if available
if (exists("neurons_std_sva_corrected") && !inherits(neurons_std_sva_corrected, "try-error")) {
  std_sva_variance <- neuron_pca_plots(
    raw_data = neurons_combined,
    corrected_data = neurons_std_sva_corrected,
    correction_method = "Standard_SVA_Limma",
    batch = all_neuron_batches,
    cell_line = neuron_combined_cell_line,
    treatment = neuron_combined_treatment,
    timestamp = timestamp,
    output_suffix = output_suffix
  )
}

#-------------------------------------------------------------------------------
# 14. Compare protein retention across all methods
#-------------------------------------------------------------------------------
# Create a summary of protein counts at each step
tracking_df <- read.csv(file.path("resultssss", "tracking", paste0("protein_tracking_", output_suffix, "_", timestamp, ".csv")))

# Create a data frame for protein count verification
protein_verification <- data.frame(
  Dataset = character(),
  Original_Count = integer(),
  Corrected_Count = integer(),
  All_Retained = logical(),
  stringsAsFactors = FALSE
)

# Protein verification for microglia
protein_verification <- rbind(protein_verification, data.frame(
  Dataset = "Microglia (ComBat) - Post-Duplicate Resolution",
  Original_Count = nrow(microglia_clean),
  Corrected_Count = nrow(microglia_combat_corrected),
  All_Retained = identical(sort(rownames(microglia_clean)), sort(rownames(microglia_combat_corrected))),
  stringsAsFactors = FALSE
))

protein_verification <- rbind(protein_verification, data.frame(
  Dataset = "Microglia (Standard SVA+Limma) - Post-Duplicate Resolution",
  Original_Count = nrow(microglia_clean),
  Corrected_Count = nrow(microglia_std_sva_corrected),
  All_Retained = identical(sort(rownames(microglia_clean)), sort(rownames(microglia_std_sva_corrected))),
  stringsAsFactors = FALSE
))

protein_verification <- rbind(protein_verification, data.frame(
  Dataset = "Microglia (Manual SVA) - Post-Duplicate Resolution",
  Original_Count = nrow(microglia_clean),
  Corrected_Count = nrow(microglia_manual_corrected),
  All_Retained = identical(sort(rownames(microglia_clean)), sort(rownames(microglia_manual_corrected))),
  stringsAsFactors = FALSE
))

# Verify neuron proteins
protein_verification <- rbind(protein_verification, data.frame(
  Dataset = "Neurons (Manual SVA) - Post-Duplicate Resolution",
  Original_Count = nrow(neurons_combined),
  Corrected_Count = nrow(neurons_manual_corrected),
  All_Retained = identical(sort(rownames(neurons_combined)), sort(rownames(neurons_manual_corrected))),
  stringsAsFactors = FALSE
))

if (exists("neurons_combat_corrected") && !inherits(neurons_combat_corrected, "try-error")) {
  protein_verification <- rbind(protein_verification, data.frame(
    Dataset = "Neurons (ComBat) - Post-Duplicate Resolution",
    Original_Count = nrow(neurons_combined),
    Corrected_Count = nrow(neurons_combat_corrected),
    All_Retained = identical(sort(rownames(neurons_combined)), sort(rownames(neurons_combat_corrected))),
    stringsAsFactors = FALSE
  ))
}

if (exists("neurons_std_sva_corrected") && !inherits(neurons_std_sva_corrected, "try-error")) {
  protein_verification <- rbind(protein_verification, data.frame(
    Dataset = "Neurons (Standard SVA+Limma) - Post-Duplicate Resolution",
    Original_Count = nrow(neurons_combined),
    Corrected_Count = nrow(neurons_std_sva_corrected),
    All_Retained = identical(sort(rownames(neurons_combined)), sort(rownames(neurons_std_sva_corrected))),
    stringsAsFactors = FALSE
  ))
}

# Save the verification summary
write_xlsx(protein_verification, 
           file.path("resultssss", "tracking", paste0("protein_retention_verification_", output_suffix, "_", timestamp, ".xlsx")))

print("Protein retention verification:")
print(protein_verification)

#-------------------------------------------------------------------------------
# 15. Compare variance explained by different methods
#-------------------------------------------------------------------------------
# Create a combined variance comparison dataframe for neurons
if (exists("combat_variance") && exists("std_sva_variance")) {
  # If all three neuron methods are available
  neuron_variance_comparison <- data.frame(
    PC = paste0("PC", 1:5),
    Raw = manual_sva_variance$Raw,  # Raw is the same in all dataframes
    Manual_SVA = manual_sva_variance$Manual_SVA,
    ComBat = combat_variance$ComBat,
    Standard_SVA_Limma = std_sva_variance$Standard_SVA_Limma
  )
} else if (exists("combat_variance")) {
  # If only Combat and Manual SVA are available
  neuron_variance_comparison <- data.frame(
    PC = paste0("PC", 1:5),
    Raw = manual_sva_variance$Raw,
    Manual_SVA = manual_sva_variance$Manual_SVA,
    ComBat = combat_variance$ComBat
  )
} else if (exists("std_sva_variance")) {
  # If only Standard SVA+Limma and Manual SVA are available
  neuron_variance_comparison <- data.frame(
    PC = paste0("PC", 1:5),
    Raw = manual_sva_variance$Raw,
    Manual_SVA = manual_sva_variance$Manual_SVA,
    Standard_SVA_Limma = std_sva_variance$Standard_SVA_Limma
  )
} else {
  # If only Manual SVA is available
  neuron_variance_comparison <- manual_sva_variance
}

# Save the combined variance comparison
write_xlsx(neuron_variance_comparison, 
           file.path("resultssss", "data", paste0("neurons_variance_comparison_all_methods_", output_suffix, "_", timestamp, ".xlsx")))

# Save a complete summary of all results
summary_list <- list(
  "Microglia Variance Comparison" = variance_comparison_microglia,
  "Neuron Variance Comparison" = neuron_variance_comparison,
  "Protein Retention Verification" = protein_verification,
  "Protein Tracking" = tracking_df
)

write_xlsx(summary_list, 
           file.path("resultssss", paste0("batch_correction_summary_", output_suffix, "_", timestamp, ".xlsx")))

#-------------------------------------------------------------------------------
# 16. Create a combined batch effect visualization for all methods
#-------------------------------------------------------------------------------
# Function to create a comprehensive plot showing batch effects before/after correction
create_batch_effect_summary <- function(cell_type, raw_data, corrected_data_list, 
                                        batch_vector, method_names, timestamp, output_suffix) {
  # PCA for raw data
  raw_pca <- prcomp(t(raw_data), scale. = TRUE)
  raw_pca_data <- as.data.frame(raw_pca$x)
  raw_pca_data$Batch <- as.factor(batch_vector)
  raw_pca_data$Method <- "Raw"
  
  # Create a combined dataframe for visualization
  combined_pca <- raw_pca_data[, c("PC1", "PC2", "Batch", "Method")]
  
  # Add PCA data for each correction method
  for (i in 1:length(corrected_data_list)) {
    method_name <- method_names[i]
    corrected_data <- corrected_data_list[[i]]
    
    if (!is.null(corrected_data)) {
      corrected_pca <- prcomp(t(corrected_data), scale. = TRUE)
      corrected_pca_data <- as.data.frame(corrected_pca$x)
      corrected_pca_data$Batch <- as.factor(batch_vector)
      corrected_pca_data$Method <- method_name
      
      # Combine with existing data
      combined_pca <- rbind(combined_pca, 
                            corrected_pca_data[, c("PC1", "PC2", "Batch", "Method")])
    }
  }
  
  # Create a faceted plot showing all methods
  summary_plot <- ggplot(combined_pca, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    facet_wrap(~ Method, scales = "free") +
    theme_bw() +
    labs(title = paste0(cell_type, " Batch Effect Removal Comparison (Post-Duplicate Resolution)"),
         color = "Batch") +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          strip.background = element_rect(fill = "lightgray"),
          strip.text = element_text(size = 12, face = "bold"))
  
  # Save the summary plot
  ggsave(
    file.path("resultssss", "plots", paste0(tolower(cell_type), "_batch_effect_summary_", output_suffix, "_", timestamp, ".png")),
    summary_plot,
    width = 15,
    height = 10,
    dpi = 300
  )
  
  return(summary_plot)
}

# Create summary plots if all methods worked successfully
if (all(exists("microglia_combat_corrected"), 
        exists("microglia_std_sva_corrected"), 
        exists("microglia_manual_corrected"))) {
  
  microglia_methods <- c("ComBat", "Standard SVA+Limma", "Manual SVA")
  microglia_corrected_list <- list(
    microglia_combat_corrected,
    microglia_std_sva_corrected,
    microglia_manual_corrected
  )
  
  microglia_summary_plot <- create_batch_effect_summary(
    cell_type = "Microglia",
    raw_data = microglia_clean,
    corrected_data_list = microglia_corrected_list,
    batch_vector = microglia_batch,
    method_names = microglia_methods,
    timestamp = timestamp,
    output_suffix = output_suffix
  )
}

# Create a combined summary for neurons
neuron_corrected_list <- list()
neuron_methods <- c()

if (exists("neurons_manual_corrected")) {
  neuron_corrected_list <- c(neuron_corrected_list, list(neurons_manual_corrected))
  neuron_methods <- c(neuron_methods, "Manual SVA")
}

if (exists("neurons_combat_corrected") && !inherits(neurons_combat_corrected, "try-error")) {
  neuron_corrected_list <- c(neuron_corrected_list, list(neurons_combat_corrected))
  neuron_methods <- c(neuron_methods, "ComBat")
}

if (exists("neurons_std_sva_corrected") && !inherits(neurons_std_sva_corrected, "try-error")) {
  neuron_corrected_list <- c(neuron_corrected_list, list(neurons_std_sva_corrected))
  neuron_methods <- c(neuron_methods, "Standard SVA+Limma")
}

if (length(neuron_corrected_list) > 0) {
  neuron_summary_plot <- create_batch_effect_summary(
    cell_type = "Neurons",
    raw_data = neurons_combined,
    corrected_data_list = neuron_corrected_list,
    batch_vector = all_neuron_batches,
    method_names = neuron_methods,
    timestamp = timestamp,
    output_suffix = output_suffix
  )
}

#-------------------------------------------------------------------------------
# 17. Create basic raw neuron plot for validation
#-------------------------------------------------------------------------------
# Create a basic PCA plot to validate the raw neuron data shows batch effects
neuron_pca_raw <- prcomp(t(neurons_combined_matrix), scale. = TRUE)
neuron_pca_data <- as.data.frame(neuron_pca_raw$x)
neuron_pca_data$Batch <- as.factor(all_neuron_batches)
neuron_pca_data$CellLine <- as.factor(neuron_combined_cell_line)
neuron_pca_data$Treatment <- as.factor(neuron_combined_treatment)

# Plot to check batch effects before correction
p_neuron_raw <- ggplot(neuron_pca_data, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  labs(title = "Raw Neuron Data - Batch Effects") +
  theme_bw()

ggsave(file.path("resultssss", "plots", paste0("neurons_raw_batch_effects_", output_suffix, "_", timestamp, ".png")),
       p_neuron_raw, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------
# 18. Save individual protein lists for reference
#-------------------------------------------------------------------------------
# Save lists of proteins in each dataset for reference
write.csv(data.frame(Protein = rownames(microglia_clean)), 
          file.path("resultssss", "tracking", paste0("microglia_proteins_", output_suffix, "_", timestamp, ".csv")), 
          row.names = FALSE)
write.csv(data.frame(Protein = rownames(neuron1_clean)), 
          file.path("resultssss", "tracking", paste0("neuron_batch1_proteins_", output_suffix, "_", timestamp, ".csv")), 
          row.names = FALSE)
write.csv(data.frame(Protein = rownames(neuron2_clean)), 
          file.path("resultssss", "tracking", paste0("neuron_batch2_proteins_", output_suffix, "_", timestamp, ".csv")), 
          row.names = FALSE)
write.csv(data.frame(Protein = rownames(neurons_combined)), 
          file.path("resultssss", "tracking", paste0("neurons_combined_proteins_", output_suffix, "_", timestamp, ".csv")), 
          row.names = FALSE)

#-------------------------------------------------------------------------------
# 19. Create method comparison summary table
#-------------------------------------------------------------------------------
# Create a comprehensive summary table
method_summary <- data.frame(
  Method = character(),
  Cell_Type = character(),
  Success = character(),
  Proteins_Retained = character(),
  Notes = character(),
  stringsAsFactors = FALSE
)

# Add microglia results
method_summary <- rbind(method_summary, data.frame(
  Method = "ComBat",
  Cell_Type = "Microglia",
  Success = "Yes",
  Proteins_Retained = paste0(nrow(microglia_combat_corrected), "/", nrow(microglia_clean)),
  Notes = "Successful batch correction",
  stringsAsFactors = FALSE
))

method_summary <- rbind(method_summary, data.frame(
  Method = "Standard SVA+Limma",
  Cell_Type = "Microglia",
  Success = "Yes",
  Proteins_Retained = paste0(nrow(microglia_std_sva_corrected), "/", nrow(microglia_clean)),
  Notes = "Successful batch correction",
  stringsAsFactors = FALSE
))

method_summary <- rbind(method_summary, data.frame(
  Method = "Manual SVA",
  Cell_Type = "Microglia",
  Success = "Yes",
  Proteins_Retained = paste0(nrow(microglia_manual_corrected), "/", nrow(microglia_clean)),
  Notes = "Successful batch correction",
  stringsAsFactors = FALSE
))

# Add neuron results
method_summary <- rbind(method_summary, data.frame(
  Method = "Manual SVA",
  Cell_Type = "Neurons",
  Success = "Yes",
  Proteins_Retained = paste0(nrow(neurons_manual_corrected), "/", nrow(neurons_combined)),
  Notes = "Successful batch correction",
  stringsAsFactors = FALSE
))

if (exists("neurons_combat_corrected") && !inherits(neurons_combat_corrected, "try-error")) {
  method_summary <- rbind(method_summary, data.frame(
    Method = "ComBat",
    Cell_Type = "Neurons",
    Success = "Yes",
    Proteins_Retained = paste0(nrow(neurons_combat_corrected), "/", nrow(neurons_combined)),
    Notes = "Successful batch correction",
    stringsAsFactors = FALSE
  ))
} else {
  method_summary <- rbind(method_summary, data.frame(
    Method = "ComBat",
    Cell_Type = "Neurons",
    Success = "No",
    Proteins_Retained = "N/A",
    Notes = "Failed - likely due to batch/biological confounding",
    stringsAsFactors = FALSE
  ))
}

if (exists("neurons_std_sva_corrected") && !inherits(neurons_std_sva_corrected, "try-error")) {
  method_summary <- rbind(method_summary, data.frame(
    Method = "Standard SVA+Limma",
    Cell_Type = "Neurons",
    Success = "Yes",
    Proteins_Retained = paste0(nrow(neurons_std_sva_corrected), "/", nrow(neurons_combined)),
    Notes = "Successful batch correction",
    stringsAsFactors = FALSE
  ))
} else {
  method_summary <- rbind(method_summary, data.frame(
    Method = "Standard SVA+Limma",
    Cell_Type = "Neurons",
    Success = "No",
    Proteins_Retained = "N/A",
    Notes = "Failed during SVA estimation",
    stringsAsFactors = FALSE
  ))
}

# Save method summary
write_xlsx(method_summary, 
           file.path("resultssss", "tracking", paste0("method_comparison_summary_", output_suffix, "_", timestamp, ".xlsx")))

#-------------------------------------------------------------------------------
# 20. Print final summary message
#-------------------------------------------------------------------------------
cat("\n=======================================================================\n")
cat("POST-DUPLICATE RESOLUTION BATCH CORRECTION ANALYSIS COMPLETED\n")
cat("=======================================================================\n")
cat("Results are saved in the 'resultssss' directory with timestamp:", timestamp, "\n\n")

cat("Data files are in 'resultssss/data/'\n")
cat("Plots are in 'resultssss/plots/'\n")
cat("Protein tracking information is in 'resultssss/tracking/'\n\n")

cat("Summary file: batch_correction_summary_", output_suffix, "_", timestamp, ".xlsx\n")
cat("Method comparison: method_comparison_summary_", output_suffix, "_", timestamp, ".xlsx\n")
cat("=======================================================================\n")

# Print method success summary
cat("METHOD SUCCESS SUMMARY:\n")
print(method_summary)
cat("\n")

cat("PROTEIN COUNTS:\n")
cat("Microglia original proteins:", nrow(microglia_clean), "\n")
cat("Neuron1 original proteins:", nrow(neuron1_clean), "\n")
cat("Neuron2 original proteins:", nrow(neuron2_clean), "\n")
cat("Neurons combined (common proteins):", nrow(neurons_combined), "\n")
cat("\n")

cat("BATCH INFORMATION:\n")
cat("Microglia batches:", paste(sort(unique(microglia_batch)), collapse = ", "), "\n")
cat("Neuron batches:", paste(sort(unique(all_neuron_batches)), collapse = ", "), "\n")
cat("\n")

cat("=======================================================================\n")
cat("Methods tested for microglia: ComBat, Standard SVA+Limma, Manual SVA\n")
cat("Methods tested for neurons: Manual SVA, ComBat, Standard SVA+Limma\n")
cat("All analyses performed after protein duplicate resolution\n")
cat("=======================================================================\n")

# Create a final validation check
cat("FINAL VALIDATION:\n")
cat("All output directories exist:", 
    all(dir.exists(c("resultssss", "resultssss/data", "resultssss/plots", "resultssss/tracking"))), "\n")

# List key output files
key_files <- c(
  paste0("batch_correction_summary_", output_suffix, "_", timestamp, ".xlsx"),
  paste0("method_comparison_summary_", output_suffix, "_", timestamp, ".xlsx")
)

for (file in key_files) {
  cat("Key file exists -", file, ":", file.exists(file.path("resultssss", file)), "\n")
}

cat("=======================================================================\n")
cat("SCRIPT EXECUTION COMPLETED SUCCESSFULLY!\n")
cat("=======================================================================\n")


# End of script
