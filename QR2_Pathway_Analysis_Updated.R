# Enhanced QR2 Pathway Analysis Script
# Improved file handling, more robust mapping, and higher quality visualizations

# Load required libraries with better error handling
load_required_packages <- function() {
  required_packages <- c(
    "clusterProfiler", "org.Hs.eg.db", "tidyverse", "AnnotationDbi",
    "enrichplot", "DOSE", "ggrepel", "RColorBrewer", "data.table", 
    "ggplot2", "gridExtra", "patchwork", "viridis", "digest"
  )
  
  missing <- character()
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    } else {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
  
  if (length(missing) > 0) {
    message("The following packages are required but not installed:")
    for (pkg in missing) {
      if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot", "DOSE")) {
        message("  - ", pkg, " (install using BiocManager::install('", pkg, "'))")
      } else {
        message("  - ", pkg, " (install using install.packages('", pkg, "'))")
      }
    }
    stop("Please install the missing packages before continuing")
  }
  
  return(TRUE)
}

# Create a function to set up the results directory with better organization
setup_results_dir <- function(base_name = "QR2_Pathway_Analysis_2025") {
  # Create base directory
  results_dir <- base_name
  dir.create(results_dir, showWarnings = FALSE)
  
  # Create a timestamped subfolder for this run
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  subfolder <- paste0(results_dir, "/run_", timestamp)
  dir.create(subfolder, showWarnings = FALSE)
  
  # Create more detailed subdirectories for different analysis types and cell types
  for (cell_type in c("microglia", "neuronal")) {
    # Base directories for each cell type
    cell_dir <- file.path(subfolder, cell_type)
    dir.create(cell_dir, showWarnings = FALSE)
    
    # Analysis type directories
    dir.create(file.path(cell_dir, "GO"), showWarnings = FALSE)
    dir.create(file.path(cell_dir, "KEGG"), showWarnings = FALSE)
    dir.create(file.path(cell_dir, "comparisons"), showWarnings = FALSE)
    
    # Threshold subdirectories
    for (threshold in c("log2FC_1.0", "log2FC_0.8")) {
      dir.create(file.path(cell_dir, "GO", threshold), showWarnings = FALSE)
      dir.create(file.path(cell_dir, "KEGG", threshold), showWarnings = FALSE)
    }
  }
  
  # Joint comparison directory
  dir.create(file.path(subfolder, "cell_comparisons"), showWarnings = FALSE)
  
  return(subfolder)
}

# Function to find the correct file with more robust naming handling
find_file <- function(base_dir, cell_type, file_pattern, verbose = TRUE) {
  # Create a list of possible variations in naming
  variations <- c(cell_type)
  
  # Add common variations
  if (cell_type == "neuronal") {
    variations <- c(variations, "neuron", "neurons")
  } else if (cell_type == "microglia") {
    variations <- c(variations, "microglial", "microglia-", "microglia_")
  }
  
  # Also try different cases
  variations <- c(
    variations,
    tolower(variations),
    toupper(variations),
    sapply(variations, function(v) paste0(toupper(substring(v, 1, 1)), substring(v, 2)))
  )
  
  # Remove duplicates
  variations <- unique(variations)
  
  # Build possible paths with all variations
  possible_paths <- c()
  for (var in variations) {
    # Try with the variation at the beginning
    possible_paths <- c(possible_paths, file.path(base_dir, paste0(var, file_pattern)))
    
    # Also try with variation in middle/end of filename for flexibility
    for (sep in c("_", "-", " ")) {
      possible_paths <- c(possible_paths, file.path(base_dir, paste0(sep, var, file_pattern)))
    }
  }
  
  # Also check for files that might contain the pattern anywhere
  all_files <- list.files(base_dir, full.names = TRUE)
  
  # Filter to find files containing both the cell type and the characteristic part of the pattern
  key_pattern_part <- gsub(".*(_log2FC_gt_[0-9\\.]+).*", "\\1", file_pattern)
  if (key_pattern_part != file_pattern) {  # Only if we extracted something
    for (var in variations) {
      matching_files <- grep(paste0(var, ".*", key_pattern_part), all_files, value = TRUE)
      possible_paths <- c(possible_paths, matching_files)
    }
  }
  
  # Remove duplicates again
  possible_paths <- unique(possible_paths)
  
  # Check which file exists
  for (path in possible_paths) {
    if (file.exists(path)) {
      if (verbose && path != possible_paths[1]) {
        message("Found file with alternative naming: ", path)
      }
      return(path)
    }
  }
  
  # If no file is found, return NULL
  if (verbose) {
    message("Could not find file matching pattern: ", file_pattern, " for cell type: ", cell_type)
    message("Paths checked:")
    for (i in seq_along(possible_paths)[1:min(5, length(possible_paths))]) {
      message("  ", possible_paths[i])
    }
    if (length(possible_paths) > 5) {
      message("  ... and ", length(possible_paths) - 5, " more")
    }
  }
  
  return(NULL)
}

# Enhanced function to get significant proteins from CSV with better error handling
get_significant_proteins_from_csv <- function(csv_file_path, threshold, verbose = TRUE) {
  if(!file.exists(csv_file_path)) {
    if (verbose) message("File does not exist: ", csv_file_path)
    return(NULL)
  }
  
  tryCatch({
    data <- read.csv(csv_file_path, stringsAsFactors = FALSE)
    
    if (verbose) {
      message("CSV file structure: ", basename(csv_file_path))
      message("Available columns: ", paste(colnames(data), collapse = ", "))
      message("Total rows: ", nrow(data))
    }
    
    if ("protein" %in% colnames(data)) {
      sig_proteins <- unique(data$protein)
    } else {
      sig_proteins <- unique(data[[1]])
    }
    
    sig_proteins <- sig_proteins[!is.na(sig_proteins) & sig_proteins != ""]
    
    if (verbose) {
      message("Found ", length(sig_proteins), " pre-filtered significant proteins")
      message("Sample proteins: ", paste(head(sig_proteins, 3), collapse = ", "))
    }
    
    return(sig_proteins)
    
  }, error = function(e) {
    if (verbose) {
      message("Error reading CSV file: ", e$message)
    }
    return(NULL)
  })
}

# Enhanced function to get all proteins for universe
get_universe_proteins_from_csv <- function(csv_file_path, verbose = TRUE) {
  # Check if file exists
  if(!file.exists(csv_file_path)) {
    if (verbose) message("File does not exist: ", csv_file_path)
    return(NULL)
  }
  
  # Try to read CSV with enhanced error handling
  tryCatch({
    # Check header structure first
    header_line <- readLines(csv_file_path, n = 1)
    header_cols <- strsplit(header_line, ",")[[1]]
    
    # Look for protein column with flexible naming
    protein_col <- grep("^[Pp]rotein$|^[Pp]rotein_?[Ii][Dd]?$|^[Pp]rotein_?[Nn]ame$", header_cols, value = TRUE)
    
    # Read data with different approaches if needed
    if (length(protein_col) > 0) {
      data <- read.csv(csv_file_path, stringsAsFactors = FALSE)
      universe_proteins <- data[[protein_col[1]]]  # Use first match if multiple columns found
    } else {
      # Try reading the data to inspect column names
      data <- read.csv(csv_file_path, stringsAsFactors = FALSE)
      
      # Try common column names
      common_protein_cols <- c("Protein", "protein", "Protein_ID", "protein_id", "ID", "id", "Name", "name")
      for (col in common_protein_cols) {
        if (col %in% colnames(data)) {
          universe_proteins <- data[[col]]
          break
        }
      }
      
      # If still no match, try to identify the protein column by content
      if (!exists("universe_proteins")) {
        for (col in colnames(data)) {
          if (!is.numeric(data[[col]])) {  # Only check non-numeric columns
            sample_values <- data[[col]]
            # Check if values look like protein IDs (contain "_HUMAN" or similar patterns)
            if (any(grepl("_HUMAN|^[A-Z0-9]+$", sample_values))) {
              universe_proteins <- data[[col]]
              break
            }
          }
        }
      }
    }
    
    # Check if we found the proteins
    if (exists("universe_proteins")) {
      # Remove NA or empty values
      universe_proteins <- universe_proteins[!is.na(universe_proteins) & universe_proteins != ""]
      
      # Report
      if (verbose) {
        message("Loaded ", length(universe_proteins), " universe proteins from ", basename(csv_file_path))
      }
      return(universe_proteins)
    } else {
      if (verbose) {
        message("Cannot identify protein column in CSV file: ", basename(csv_file_path))
        message("Available columns: ", paste(colnames(data), collapse = ", "))
      }
      return(NULL)
    }
  }, error = function(e) {
    if (verbose) {
      message("Error reading CSV file: ", e$message)
    }
    return(NULL)
  })
}

# Function to find file with multiple patterns
find_file_with_patterns <- function(base_dir, patterns, cell_type) {
  for (pattern in patterns) {
    file_path <- find_file(base_dir, cell_type, pattern, verbose = FALSE)
    if (!is.null(file_path)) {
      return(file_path)
    }
  }
  return(NULL)
}

# Enhanced mapping function

map_proteins_to_entrez <- function(proteins, verbose = TRUE) {
  if (verbose) message("Starting enhanced protein mapping...")
  
  # Container for final results
  final_mapping <- list()
  unmapped <- character()
  
  # Remove any NA or empty proteins
  proteins <- proteins[!is.na(proteins) & proteins != ""]
  unique_proteins <- unique(proteins)
  
  if (verbose) message("Mapping ", length(unique_proteins), " unique proteins...")
  
  # Strategy 1: Direct mapping of UniProt style IDs
  uniprot_pattern <- "^[A-Z0-9]+_HUMAN$"
  uniprot_ids <- unique_proteins[grep(uniprot_pattern, unique_proteins)]
  
  if (length(uniprot_ids) > 0) {
    if (verbose) message("Found ", length(uniprot_ids), " UniProt-style IDs")
    
    # Extract gene symbols from UniProt IDs
    gene_symbols <- gsub("_HUMAN$", "", uniprot_ids)
    
    # Map gene symbols to Entrez
    symbol_to_entrez <- suppressMessages(
      bitr(gene_symbols, 
           fromType = "SYMBOL", 
           toType = "ENTREZID", 
           OrgDb = org.Hs.eg.db)
    )
    
    # Store successful mappings
    for (i in seq_along(gene_symbols)) {
      if (gene_symbols[i] %in% symbol_to_entrez$SYMBOL) {
        entrez <- symbol_to_entrez$ENTREZID[symbol_to_entrez$SYMBOL == gene_symbols[i]]
        final_mapping[[uniprot_ids[i]]] <- entrez[1]
      }
    }
  }
  
  # Strategy 2: Try mapping as gene symbols directly (FIXED VERSION)
  remaining <- setdiff(unique_proteins, names(final_mapping))
  if (length(remaining) > 0) {
    if (verbose) message("Attempting direct gene symbol mapping for ", length(remaining), " proteins...")
    
    # Clean potential gene symbols
    potential_symbols <- gsub("^sp\\||^tr\\|", "", remaining)
    potential_symbols <- gsub("\\|.*$", "", potential_symbols)
    
    # FIXED: Filter out obviously invalid symbols and validate before mapping
    # Keep only symbols that look like valid gene names (2-20 characters, letters/numbers)
    valid_symbols <- potential_symbols[grepl("^[A-Za-z][A-Za-z0-9]{1,19}$", potential_symbols)]
    
    if (length(valid_symbols) > 0) {
      # FIXED: Add proper error handling and validation
      tryCatch({
        # Check which symbols are actually valid keys first
        all_valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
        symbols_to_map <- intersect(valid_symbols, all_valid_symbols)
        
        if (length(symbols_to_map) > 0) {
          symbol_mapping <- suppressMessages(
            bitr(symbols_to_map, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)
          )
          
          # Map back to original protein IDs
          for (i in seq_along(potential_symbols)) {
            if (potential_symbols[i] %in% symbol_mapping$SYMBOL) {
              entrez <- symbol_mapping$ENTREZID[symbol_mapping$SYMBOL == potential_symbols[i]]
              final_mapping[[remaining[i]]] <- entrez[1]
            }
          }
          
          if (verbose) {
            message("Successfully mapped ", nrow(symbol_mapping), " proteins via gene symbols")
          }
        } else {
          if (verbose) message("No valid gene symbols found in remaining proteins")
        }
      }, error = function(e) {
        if (verbose) {
          message("Error in gene symbol mapping: ", e$message)
          message("Continuing with other mapping strategies...")
        }
      })
    }
  }
  
  # Strategy 3: Pattern-based extraction
  remaining <- setdiff(unique_proteins, names(final_mapping))
  if (length(remaining) > 0) {
    if (verbose) message("Trying pattern extraction for ", length(remaining), " proteins...")
    
    for (protein in remaining) {
      # Try multiple patterns to extract gene symbols
      gene_symbol <- NULL
      
      # Pattern 1: Extract from UniProt-style IDs with prefixes
      if (grepl("^(sp|tr)\\|[A-Z0-9]+\\|([A-Z0-9]+)_", protein)) {
        gene_symbol <- gsub("^(sp|tr)\\|[A-Z0-9]+\\|([A-Z0-9]+)_.*", "\\2", protein)
      }
      # Pattern 2: Simple gene_HUMAN pattern
      else if (grepl("^[A-Z][A-Z0-9]+_HUMAN", protein)) {
        gene_symbol <- gsub("_HUMAN.*", "", protein)
      }
      # Pattern 3: Gene symbol at the start
      else if (grepl("^[A-Z][A-Z0-9]{2,}[^a-z]*$", protein)) {
        gene_symbol <- gsub("[^A-Z0-9].*", "", protein)
      }
      
      # Try to map the extracted symbol
      if (!is.null(gene_symbol) && nchar(gene_symbol) >= 2) {
        tryCatch({
          mapping <- suppressMessages(bitr(gene_symbol, "SYMBOL", "ENTREZID", org.Hs.eg.db))
          if (nrow(mapping) > 0) {
            final_mapping[[protein]] <- mapping$ENTREZID[1]
          }
        }, error = function(e) {})
      }
    }
  }
  
  # Strategy 4: Alias mapping
  remaining <- setdiff(unique_proteins, names(final_mapping))
  if (length(remaining) > 0) {
    if (verbose) message("Checking aliases for ", length(remaining), " proteins...")
    
    for (protein in remaining) {
      # Clean the protein ID
      cleaned <- gsub("_HUMAN.*", "", protein)
      cleaned <- gsub("^(sp|tr)\\|", "", cleaned)
      cleaned <- gsub("\\|.*", "", cleaned)
      
      if (nchar(cleaned) >= 2) {
        # Check if it's an alias
        tryCatch({
          alias_check <- suppressMessages(
            bitr(cleaned, fromType = "ALIAS", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
          )
          if (nrow(alias_check) > 0) {
            final_mapping[[protein]] <- alias_check$ENTREZID[1]
          }
        }, error = function(e) {})
      }
    }
  }
  
  # Strategy 5: UniProt to Entrez mapping
  remaining <- setdiff(unique_proteins, names(final_mapping))
  if (length(remaining) > 0) {
    if (verbose) message("Trying UniProt ID mapping for ", length(remaining), " proteins...")
    
    # Extract potential UniProt IDs
    uniprot_candidates <- character()
    for (protein in remaining) {
      # Extract UniProt accession patterns
      if (grepl("^(sp|tr)\\|([A-Z0-9]{6,10})\\|", protein)) {
        uniprot_id <- gsub("^(sp|tr)\\|([A-Z0-9]{6,10})\\|.*", "\\2", protein)
        uniprot_candidates <- c(uniprot_candidates, uniprot_id)
        names(uniprot_candidates)[length(uniprot_candidates)] <- protein
      } else if (grepl("^[A-Z0-9]{6,10}$", protein)) {
        uniprot_candidates <- c(uniprot_candidates, protein)
        names(uniprot_candidates)[length(uniprot_candidates)] <- protein
      }
    }
    
    if (length(uniprot_candidates) > 0) {
      # Try UniProt mapping
      tryCatch({
        uniprot_mapping <- suppressMessages(
          bitr(uniprot_candidates, 
               fromType = "UNIPROT", 
               toType = "ENTREZID", 
               OrgDb = org.Hs.eg.db)
        )
        
        for (i in seq_along(uniprot_candidates)) {
          if (uniprot_candidates[i] %in% uniprot_mapping$UNIPROT) {
            original_protein <- names(uniprot_candidates)[i]
            entrez <- uniprot_mapping$ENTREZID[uniprot_mapping$UNIPROT == uniprot_candidates[i]]
            final_mapping[[original_protein]] <- entrez[1]
          }
        }
      }, error = function(e) {
        if (verbose) message("UniProt mapping encountered an error: ", e$message)
      })
    }
  }
  
  # Strategy 6: Fuzzy matching for common variations
  remaining <- setdiff(unique_proteins, names(final_mapping))
  if (length(remaining) > 0 && length(remaining) < 100) {  # Only for small sets
    if (verbose) message("Attempting fuzzy matching for ", length(remaining), " proteins...")
    
    # Get all gene symbols from the database
    all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    
    for (protein in remaining) {
      cleaned <- toupper(gsub("[^A-Z0-9]", "", protein))
      if (nchar(cleaned) >= 3) {
        # Look for exact matches after cleaning
        matches <- all_symbols[toupper(gsub("[^A-Z0-9]", "", all_symbols)) == cleaned]
        if (length(matches) > 0) {
          tryCatch({
            mapping <- bitr(matches[1], "SYMBOL", "ENTREZID", org.Hs.eg.db)
            if (nrow(mapping) > 0) {
              final_mapping[[protein]] <- mapping$ENTREZID[1]
            }
          }, error = function(e) {})
        }
      }
    }
  }
  
  # Create final mapping vector
  remaining <- setdiff(unique_proteins, names(final_mapping))
  if (length(remaining) > 0) {
    unmapped <- remaining
  }
  
  # Convert to vector format
  entrez_ids <- character(length(proteins))
  for (i in seq_along(proteins)) {
    if (proteins[i] %in% names(final_mapping)) {
      entrez_ids[i] <- final_mapping[[proteins[i]]]
    }
  }
  
  # Report results
  mapped_count <- sum(entrez_ids != "")
  total_count <- length(proteins)
  unique_mapped <- length(unique(entrez_ids[entrez_ids != ""]))
  
  if (verbose) {
    message("Mapping complete:")
    message("  Total proteins: ", total_count)
    message("  Successfully mapped: ", mapped_count, " (", 
            round(mapped_count/total_count * 100, 1), "%)")
    message("  Unique Entrez IDs: ", unique_mapped)
    if (length(unmapped) > 0) {
      message("  Unmapped proteins: ", length(unmapped))
      if (length(unmapped) <= 10) {
        message("    Examples: ", paste(head(unmapped, 10), collapse = ", "))
      }
    }
  }
  
  # Return only non-empty Entrez IDs
  return(unique(entrez_ids[entrez_ids != ""]))
}


# Cached protein mapping function
map_proteins_with_cache <- function(proteins, cell_type, cache_name) {
  # This is a wrapper that uses map_proteins_to_entrez
  return(map_proteins_to_entrez(proteins))
}

# Enhanced visualization functions with high-quality graphics

# Fixed create_publication_barplot function without aes_string()
create_publication_barplot <- function(data, title, subtitle, filename, width = 10, height = 8, 
                                       cell_type, theme_color = "Blues") {
  # Ensure we have data
  if (nrow(data) == 0) {
    message("No data provided for barplot")
    return(NULL)
  }
  
  # Limit to top terms (max 20) for readability
  if (nrow(data) > 20) {
    data <- data[1:20, ]
  }
  
  # Prepare data for plotting
  plot_data <- data
  
  # Determine which p-value to use and reorder accordingly
  if (grepl("neuron", cell_type, ignore.case = TRUE)) {
    # Order neuronal plots by unadjusted p-value
    plot_data <- plot_data[order(plot_data$pvalue, decreasing = TRUE), ]
    plot_data$Description <- factor(plot_data$Description, levels = plot_data$Description)
    p_value_col <- plot_data$pvalue
    legend_title <- "Unadjusted\np-value"
    caption_text <- "ordered by unadjusted p-value"
    message("Ordering ", cell_type, " barplot by unadjusted p-values")
  } else {
    # Order microglia plots by adjusted p-value
    plot_data <- plot_data[order(plot_data$p.adjust, decreasing = TRUE), ]
    plot_data$Description <- factor(plot_data$Description, levels = plot_data$Description)
    p_value_col <- plot_data$p.adjust
    legend_title <- "Adjusted\np-value"
    caption_text <- "ordered by adjusted p-value"
    message("Ordering ", cell_type, " barplot by adjusted p-values")
  }
  
  # Set color palette based on cell type
  if (grepl("micro|glia", cell_type, ignore.case = TRUE)) {
    color_palette <- viridisLite::viridis(100)  # Purple-blue for microglia
  } else if (grepl("neuron", cell_type, ignore.case = TRUE)) {
    color_palette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(100)  # Orange-red for neurons
  } else {
    color_palette <- colorRampPalette(rev(brewer.pal(9, theme_color)))(100)  # Default color
  }
  
  # Create publication-quality barplot using aes() instead of aes_string()
  p <- ggplot(plot_data, aes(x = Count, y = Description, fill = p_value_col)) +
    geom_bar(stat = "identity") +
    scale_fill_gradientn(
      colors = color_palette,
      name = legend_title,
      trans = "log10"  # Log scale for better color distribution
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = paste0("Top enriched terms ", caption_text),
      x = "Gene Count", 
      y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 12, hjust = 0),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.caption = element_text(size = 9, hjust = 1, face = "italic"),
      plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
    )
  
  # Add gene ratio information
  if ("GeneRatio" %in% colnames(plot_data)) {
    p <- p + geom_text(aes(label = GeneRatio, x = Count + max(Count)*0.05), 
                       hjust = 0, size = 3, color = "gray30")
  }
  
  # Save as high-resolution PNG and PDF for publication
  ggsave(filename, p, width = width, height = height, dpi = 300)
  
  # Save PDF version for publication
  pdf_filename <- sub("\\.png$", ".pdf", filename)
  ggsave(pdf_filename, p, width = width, height = height)
  
  return(p)
}

# Fixed create_publication_dotplot function with correct color parameter
create_publication_dotplot <- function(ego_result, title, subtitle, filename, width = 10, height = 8, 
                                       cell_type, theme_color = "Blues") {
  # Check if we have valid data
  if (is.null(ego_result) || nrow(ego_result@result) == 0) {
    message("No data provided for dotplot")
    return(NULL)
  }
  
  # Add ordering information to subtitle based on cell type
  if (grepl("neuron", cell_type, ignore.case = TRUE)) {
    subtitle <- paste0(subtitle, " (ordered by unadjusted p-value)")
    color_by <- "pvalue"  # FIX: Use correct parameter value
    message("Creating dotplot for ", cell_type, " ordered by unadjusted p-values")
  } else {
    subtitle <- paste0(subtitle, " (ordered by adjusted p-value)")
    color_by <- "p.adjust"  # FIX: Use correct parameter value
    message("Creating dotplot for ", cell_type, " ordered by adjusted p-values")
  }
  
  # Create enhanced dotplot with correct color parameter
  p <- dotplot(ego_result, 
               showCategory = min(20, nrow(ego_result@result)), 
               font.size = 12, 
               color = color_by) +  # FIX: Use the correct color parameter
    labs(title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 12, hjust = 0),
      axis.text.y = element_text(size = 11, color = "black"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.5)
    )
  
  # Save as high-resolution PNG
  ggsave(filename, p, width = width, height = height, dpi = 300)
  
  # Save PDF version for publication
  pdf_filename <- sub("\\.png$", ".pdf", filename)
  ggsave(pdf_filename, p, width = width, height = height)
  
  return(p)
}

# Create enhanced network plot with cell-type awareness
create_network_plot <- function(ego_result, title, subtitle, filename, width = 12, height = 10,
                                cell_type, plot_type = "cnet") {
  if (is.null(ego_result) || nrow(ego_result@result) == 0 || nrow(ego_result@result) < 5) {
    message("Not enough data for network plot")
    return(NULL)
  }
  
  # Add ordering information to subtitle based on cell type
  if (grepl("neuron", cell_type, ignore.case = TRUE)) {
    subtitle <- paste0(subtitle, " (ordered by unadjusted p-value)")
    message("Creating ", plot_type, " network plot for ", cell_type, " ordered by unadjusted p-values")
  } else {
    subtitle <- paste0(subtitle, " (ordered by adjusted p-value)")
    message("Creating ", plot_type, " network plot for ", cell_type, " ordered by adjusted p-values")
  }
  
  # Try to create the plot with error handling
  tryCatch({
    if (plot_type == "cnet") {
      # Category-gene network plot
      p <- cnetplot(ego_result, 
                    showCategory = min(10, nrow(ego_result@result)),
                    categorySize = "pvalue",
                    foldChange = NULL,
                    node_label = "category",
                    cex_category = 1.2,
                    cex_gene = 0.8) +
        labs(title = title, subtitle = subtitle) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        )
    } else if (plot_type == "emap") {
      # First simplify to get more meaningful clusters
      ego_simp <- simplify(ego_result, cutoff = 0.7, by = "p.adjust")
      
      p <- emapplot(ego_simp, 
                    showCategory = min(30, nrow(ego_simp@result)),
                    color = "p.adjust", 
                    layout = "kk") +
        labs(title = title, subtitle = subtitle) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        )
    } else if (plot_type == "heat") {
      # Heatmap of gene-term relationships
      p <- heatplot(ego_result, 
                    showCategory = min(15, nrow(ego_result@result)),
                    foldChange = NULL) +
        labs(title = title, subtitle = subtitle) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        )
    }
    
    # Save as high-resolution PNG and PDF
    ggsave(filename, p, width = width, height = height, dpi = 300)
    pdf_filename <- sub("\\.png$", ".pdf", filename)
    ggsave(pdf_filename, p, width = width, height = height)
    
    return(p)
  }, error = function(e) {
    message("Could not create network plot: ", e$message)
    return(NULL)
  })
}

# Fixed visualize_go_results function with proper ego_reordered creation
visualize_go_results <- function(ego, cell_type, threshold_value, output_dir, ontology) {
  if (is.null(ego) || nrow(ego@result) == 0) {
    message("No significant GO terms found for ", ontology)
    return(NULL)
  }
  
  # File paths
  base_filename <- file.path(output_dir, paste0(cell_type, "_log2FC_", threshold_value))
  
  # Add ontology to base results file name
  results_file <- paste0(base_filename, "_", ontology, "_GO_results.csv")
  
  # Save detailed results
  write.csv(ego@result, file = results_file, row.names = FALSE)
  message("Saved ", nrow(ego@result), " enriched GO terms to ", basename(results_file))
  
  # Create visualizations if we have enough terms
  if (nrow(ego@result) >= 5) {
    # Determine how many terms to show (max 20)
    top_n <- min(20, nrow(ego@result))
    
    # Order data based on cell type and create ego_reordered
    if (grepl("neuron", cell_type, ignore.case = TRUE)) {
      # Order by unadjusted p-value for neurons
      ordered_results <- ego@result[order(ego@result$pvalue),]
      plot_data <- ordered_results[1:top_n, ]
      
      # Create reordered ego object for neurons
      ego_reordered <- ego
      ego_reordered@result <- ordered_results
      
      message("Selecting top ", top_n, " ", ontology, " terms for ", cell_type, " ordered by unadjusted p-values")
    } else {
      # Keep default ordering (by adjusted p-value) for microglia
      plot_data <- ego@result[1:top_n, ]
      
      # For microglia, ego_reordered is the same as ego
      ego_reordered <- ego
      
      message("Selecting top ", top_n, " ", ontology, " terms for ", cell_type, " ordered by adjusted p-values")
    }
    
    # Basic title and subtitle
    title <- paste0(cell_type, " ", ontology, " Gene Ontology Terms")
    subtitle <- paste0("log2FC > ", threshold_value)
    
    # Set color theme based on ontology
    if (ontology == "BP") {
      theme_color <- "Blues"
    } else if (ontology == "MF") {
      theme_color <- "Greens" 
    } else {
      theme_color <- "Purples"
    }
    
    # Create barplot with cell-type specific ordering
    barplot_file <- paste0(base_filename, "_", ontology, "_barplot.png")
    barplot <- create_publication_barplot(
      plot_data, 
      title, 
      subtitle, 
      barplot_file,
      cell_type = cell_type,
      theme_color = theme_color
    )
    
    # Create dotplot using the appropriately ordered ego object
    dotplot_file <- paste0(base_filename, "_", ontology, "_dotplot.png")
    dotplot <- create_publication_dotplot(
      ego_reordered,  # Use the properly created ego_reordered
      title,
      subtitle,
      dotplot_file,
      cell_type = cell_type,
      theme_color = theme_color
    )
    
    # Create category-gene network plot (for BP only or if many terms)
    if (ontology == "BP" || nrow(ego@result) >= 20) {
      cnet_file <- paste0(base_filename, "_", ontology, "_cnetplot.png")
      cnet_title <- paste0(cell_type, " Gene-Term Network")
      cnet_subtitle <- paste0(ontology, " Ontology, log2FC > ", threshold_value)
      
      create_network_plot(
        ego_reordered,  # Use the properly created ego_reordered
        cnet_title,
        cnet_subtitle,
        cnet_file,
        cell_type = cell_type,
        plot_type = "cnet"
      )
    }
    
    # Create term relationship map
    if (nrow(ego@result) >= 10) {
      emap_file <- paste0(base_filename, "_", ontology, "_emapplot.png")
      emap_title <- paste0(cell_type, " Term Similarity Network")
      emap_subtitle <- paste0(ontology, " Ontology, log2FC > ", threshold_value)
      
      create_network_plot(
        ego_reordered,  # Use the properly created ego_reordered
        emap_title,
        emap_subtitle,
        emap_file,
        cell_type = cell_type,
        plot_type = "emap"
      )
    }
    
    # Create heatmap for gene-term relationships (BP only)
    if (ontology == "BP" && nrow(ego@result) >= 8) {
      heat_file <- paste0(base_filename, "_", ontology, "_heatplot.png")
      heat_title <- paste0(cell_type, " Gene-Term Relationships")
      heat_subtitle <- paste0("Biological Process, log2FC > ", threshold_value)
      
      create_network_plot(
        ego_reordered,  # Use the properly created ego_reordered
        heat_title,
        heat_subtitle,
        heat_file,
        cell_type = cell_type,
        plot_type = "heat"
      )
    }
  }
  
  return(ego@result)
}

# Fixed create_combined_go_visualization function without aes_string()
create_combined_go_visualization <- function(all_results, cell_type, threshold_value, output_dir) {
  if (length(all_results) == 0 || nrow(all_results) == 0) {
    message("No GO results to visualize")
    return(NULL)
  }
  
  # Get top terms by appropriate p-value based on cell type
  if (grepl("neuron", cell_type, ignore.case = TRUE)) {
    # Order by unadjusted p-value for neurons
    top_terms <- all_results[order(all_results$pvalue),][1:min(20, nrow(all_results)),]
    caption_text <- "Top terms sorted by unadjusted p-value across all ontologies"
    subtitle_add <- " (ordered by unadjusted p-value)"
    message("Creating combined GO visualization for ", cell_type, " ordered by unadjusted p-values")
  } else {
    # Order by adjusted p-value for microglia
    top_terms <- all_results[order(all_results$p.adjust),][1:min(20, nrow(all_results)),]
    caption_text <- "Top terms sorted by adjusted p-value across all ontologies"
    subtitle_add <- " (ordered by adjusted p-value)"
    message("Creating combined GO visualization for ", cell_type, " ordered by adjusted p-values")
  }
  
  # Prepare color palette based on cell type
  if (grepl("micro|glia", cell_type, ignore.case = TRUE)) {
    pal <- "Set1"
  } else if (grepl("neuron", cell_type, ignore.case = TRUE)) {
    pal <- "Dark2"
  } else {
    pal <- "Set2"
  }
  
  # Reorder Description factor for plotting
  if (grepl("neuron", cell_type, ignore.case = TRUE)) {
    top_terms$Description <- factor(top_terms$Description, 
                                    levels = top_terms$Description[order(top_terms$pvalue, decreasing = TRUE)])
  } else {
    top_terms$Description <- factor(top_terms$Description, 
                                    levels = top_terms$Description[order(top_terms$p.adjust, decreasing = TRUE)])
  }
  
  # Create plot with better aesthetics and appropriate ordering using aes() instead of aes_string()
  p <- ggplot(top_terms, aes(x = Description, y = Count, fill = Ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_brewer(palette = pal) +
    labs(title = paste0(cell_type, " Top Gene Ontology Terms"),
         subtitle = paste0("log2FC > ", threshold_value, subtitle_add),
         caption = caption_text,
         x = "", y = "Gene Count") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 12, hjust = 0),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
    )
  
  # File path
  combined_plot_file <- file.path(output_dir, paste0(cell_type, "_log2FC_", threshold_value, "_top_terms_combined.png"))
  
  # Save as high-resolution PNG and PDF
  ggsave(combined_plot_file, p, width = 12, height = 8, dpi = 300)
  pdf_filename <- sub("\\.png$", ".pdf", combined_plot_file)
  ggsave(pdf_filename, p, width = 12, height = 8)
  
  message("Created combined visualization at ", basename(combined_plot_file))
  
  return(p)
}

# Enhanced KEGG pathway analysis with improved visualization and cell-type specific ordering

# Function to run KEGG pathway analysis with multiple fallbacks
run_kegg_analysis <- function(sig_entrez, universe_entrez, cell_type, threshold_value, output_dir) {
  # Check if we have enough data for analysis
  if (length(sig_entrez) < 5) {
    message("Warning: Only ", length(sig_entrez), " proteins mapped to Entrez IDs for KEGG analysis. Results may be limited.")
  }
  
  # Ensure we're working with character vectors
  sig_char <- as.character(sig_entrez)
  universe_char <- as.character(universe_entrez)
  
  message("\n===== Running KEGG pathway analysis for ", cell_type, 
          " at log2FC > ", threshold_value, " =====")
  
  # Create output directory path
  kegg_dir <- output_dir
  
  # Try KEGG enrichment with error handling and multiple approaches
  tryCatch({
    # First attempt with default settings and more permissive cutoffs
    kk <- enrichKEGG(
      gene = sig_char,
      universe = universe_char,
      organism = "hsa",  # human
      keyType = "ncbi-geneid",
      pvalueCutoff = 0.1,  # More permissive
      pAdjustMethod = "BH",
      minGSSize = 3  # Allow smaller gene sets
    )
    
    # If no results, try with even more permissive settings
    if (is.null(kk) || nrow(kk@result) == 0) {
      message("No KEGG pathways found with default settings. Trying more permissive thresholds...")
      kk <- enrichKEGG(
        gene = sig_char,
        universe = universe_char,
        organism = "hsa",
        keyType = "ncbi-geneid",
        pvalueCutoff = 0.2,  # Even more permissive
        qvalueCutoff = 0.5,
        pAdjustMethod = "BH",
        minGSSize = 2  # Allow very small gene sets
      )
    }
    
    # Process results if we have any
    if (!is.null(kk) && nrow(kk@result) > 0) {
      # Save KEGG results
      kegg_file <- file.path(kegg_dir, paste0(cell_type, "_log2FC_", threshold_value, "_KEGG_pathways.csv"))
      write.csv(kk@result, file = kegg_file, row.names = FALSE)
      message("Saved ", nrow(kk@result), " KEGG pathways to ", basename(kegg_file))
      
      # Create visualizations if we have enough pathways
      if (nrow(kk@result) >= 3) {
        # Order results based on cell type for visualization
        if (grepl("neuron", cell_type, ignore.case = TRUE)) {
          # Order by unadjusted p-value for neurons
          ordered_results <- kk@result[order(kk@result$pvalue),]
          message("Ordering KEGG results for ", cell_type, " by unadjusted p-values")
        } else {
          # Keep default ordering (by adjusted p-value) for microglia
          ordered_results <- kk@result[order(kk@result$p.adjust),]
          message("Ordering KEGG results for ", cell_type, " by adjusted p-values")
        }
        
        # Enhanced KEGG barplot with fixed aes() syntax
        kegg_plot_file <- file.path(kegg_dir, paste0(cell_type, "_log2FC_", threshold_value, "_KEGG_barplot.png"))
        
        # Select color based on cell type
        if (grepl("micro|glia", cell_type, ignore.case = TRUE)) {
          bar_color <- brewer.pal(9, "Purples")[7]
        } else if (grepl("neuron", cell_type, ignore.case = TRUE)) {
          bar_color <- brewer.pal(9, "YlOrRd")[7]
        } else {
          bar_color <- brewer.pal(9, "Blues")[7]
        }
        
        # Take top results based on appropriate ordering and prepare for plotting
        kegg_data <- ordered_results[1:min(20, nrow(ordered_results)),]
        
        # Determine ordering and create properly ordered factor
        if (grepl("neuron", cell_type, ignore.case = TRUE)) {
          # Order by unadjusted p-value for neurons
          kegg_data <- kegg_data[order(kegg_data$pvalue, decreasing = TRUE), ]
          kegg_data$Description <- factor(kegg_data$Description, levels = kegg_data$Description)
          subtitle_text <- paste0("log2FC > ", threshold_value, " (ordered by unadjusted p-value)")
          caption_text <- "Pathways sorted by unadjusted p-value"
        } else {
          # Order by adjusted p-value for microglia
          kegg_data <- kegg_data[order(kegg_data$p.adjust, decreasing = TRUE), ]
          kegg_data$Description <- factor(kegg_data$Description, levels = kegg_data$Description)
          subtitle_text <- paste0("log2FC > ", threshold_value, " (ordered by adjusted p-value)")
          caption_text <- "Pathways sorted by adjusted p-value"
        }
        
        # Create a publication-quality KEGG barplot with fixed aes() syntax
        p_kegg <- ggplot(kegg_data, aes(x = Count, y = Description)) +
          geom_bar(stat = "identity", fill = bar_color) +
          labs(title = paste0(cell_type, " KEGG Pathway Enrichment"),
               subtitle = subtitle_text,
               caption = caption_text,
               x = "Gene Count", y = "") +
          theme_minimal(base_size = 12) +
          theme(
            axis.text.y = element_text(size = 10),
            plot.title = element_text(face = "bold", size = 14, hjust = 0),
            plot.subtitle = element_text(size = 12, hjust = 0),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
          )
        
        # Add gene ratio text for more information
        if ("GeneRatio" %in% colnames(kegg_data)) {
          p_kegg <- p_kegg + 
            geom_text(aes(label = GeneRatio), hjust = -0.2, size = 3, color = "gray30")
        }
        
        # Save as both PNG and PDF for publication
        ggsave(kegg_plot_file, p_kegg, width = 10, height = 8, dpi = 300)
        pdf_file <- sub("\\.png$", ".pdf", kegg_plot_file)
        ggsave(pdf_file, p_kegg, width = 10, height = 8)
        
        message("Created KEGG barplot at ", basename(kegg_plot_file))
        
        # Create KEGG dotplot with improved aesthetics and cell-type specific ordering
        kegg_dot_file <- file.path(kegg_dir, paste0(cell_type, "_log2FC_", threshold_value, "_KEGG_dotplot.png"))
        
        # Set color based on cell type
        if (grepl("micro|glia", cell_type, ignore.case = TRUE)) {
          dot_color <- "purple4"
        } else if (grepl("neuron", cell_type, ignore.case = TRUE)) {
          dot_color <- "firebrick3"
        } else {
          dot_color <- "dodgerblue4"
        }
        
        # Create enrichResult object with reordered data for dotplot
        kk_reordered <- kk
        kk_reordered@result <- ordered_results
        
        p_kegg_dot <- dotplot(kk_reordered, showCategory = min(20, nrow(kk_reordered@result)),
                              title = paste0(cell_type, " KEGG Pathway Enrichment"),
                              color = dot_color) +
          labs(subtitle = paste0("log2FC > ", threshold_value, 
                                 ifelse(grepl("neuron", cell_type, ignore.case = TRUE), 
                                        " (ordered by unadjusted p-value)", 
                                        " (ordered by adjusted p-value)"))) +
          theme_minimal() +
          theme(
            plot.title = element_text(face = "bold", size = 14, hjust = 0),
            plot.subtitle = element_text(size = 12, hjust = 0),
            axis.text.y = element_text(size = 10, color = "black"),
            panel.grid.minor = element_blank()
          )
        
        # Save as PNG and PDF
        ggsave(kegg_dot_file, p_kegg_dot, width = 10, height = 8, dpi = 300)
        pdf_file <- sub("\\.png$", ".pdf", kegg_dot_file)
        ggsave(pdf_file, p_kegg_dot, width = 10, height = 8)
        
        message("Created KEGG dotplot at ", basename(kegg_dot_file))
        
        # Create KEGG pathway-gene relationship plot (cnetplot)
        kegg_net_file <- file.path(kegg_dir, paste0(cell_type, "_log2FC_", threshold_value, "_KEGG_cnetplot.png"))
        
        # Only attempt if we have enough data
        if (nrow(kk@result) >= 5) {
          tryCatch({
            p_kegg_net <- cnetplot(kk_reordered, 
                                   showCategory = min(8, nrow(kk_reordered@result)),
                                   categorySize = "pvalue",
                                   foldChange = NULL,
                                   node_label = "category") +
              labs(title = paste0(cell_type, " KEGG Pathway-Gene Network"),
                   subtitle = paste0("log2FC > ", threshold_value, 
                                     ifelse(grepl("neuron", cell_type, ignore.case = TRUE), 
                                            " (ordered by unadjusted p-value)", 
                                            " (ordered by adjusted p-value)"))) +
              theme(
                plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5)
              )
            
            # Save as PNG and PDF
            ggsave(kegg_net_file, p_kegg_net, width = 12, height = 10, dpi = 300)
            pdf_file <- sub("\\.png$", ".pdf", kegg_net_file)
            ggsave(pdf_file, p_kegg_net, width = 12, height = 10)
            
            message("Created KEGG network plot at ", basename(kegg_net_file))
          }, error = function(e) {
            message("Could not create KEGG network plot: ", e$message)
          })
        }
        
        # Try creating a heatmap for KEGG pathways
        kegg_heat_file <- file.path(kegg_dir, paste0(cell_type, "_log2FC_", threshold_value, "_KEGG_heatmap.png"))
        
        if (nrow(kk@result) >= 5 && length(sig_char) >= 10) {
          tryCatch({
            p_kegg_heat <- heatplot(kk_reordered, 
                                    showCategory = min(10, nrow(kk_reordered@result)),
                                    foldChange = NULL) +
              labs(title = paste0(cell_type, " KEGG Pathway-Gene Associations"),
                   subtitle = paste0("log2FC > ", threshold_value, 
                                     ifelse(grepl("neuron", cell_type, ignore.case = TRUE), 
                                            " (ordered by unadjusted p-value)", 
                                            " (ordered by adjusted p-value)"))) +
              theme(
                plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5)
              )
            
            # Save as PNG and PDF
            ggsave(kegg_heat_file, p_kegg_heat, width = 12, height = 9, dpi = 300)
            pdf_file <- sub("\\.png$", ".pdf", kegg_heat_file)
            ggsave(pdf_file, p_kegg_heat, width = 12, height = 9)
            
            message("Created KEGG heatmap at ", basename(kegg_heat_file))
          }, error = function(e) {
            message("Could not create KEGG heatmap: ", e$message)
          })
        }
      }
      
      # Check for specific metabolic pathways of interest
      metabolic_pathways <- c(
        "hsa00190" = "Oxidative phosphorylation",
        "hsa00010" = "Glycolysis / Gluconeogenesis",
        "hsa00030" = "Pentose phosphate pathway",
        "hsa00020" = "Citrate cycle (TCA cycle)",
        "hsa01200" = "Carbon metabolism",
        "hsa00071" = "Fatty acid degradation",
        "hsa00280" = "Valine, leucine and isoleucine degradation",
        "hsa00620" = "Pyruvate metabolism",
        "hsa00650" = "Butanoate metabolism",
        "hsa00630" = "Glyoxylate and dicarboxylate metabolism",
        "hsa00640" = "Propanoate metabolism"
      )
      
      # Enhanced metabolic pathway analysis results
      meta_results <- data.frame(
        PathwayID = character(),
        PathwayName = character(),
        Found = logical(),
        PValue = numeric(),
        AdjustedPValue = numeric(),
        GeneCount = integer(),
        GeneRatio = character(),
        GeneSymbols = character(),
        stringsAsFactors = FALSE
      )
      
      for (pathway_id in names(metabolic_pathways)) {
        pathway_name <- metabolic_pathways[pathway_id]
        
        # Check if this pathway is in our results
        pathway_row <- kk@result[kk@result$ID == pathway_id,]
        
        if (nrow(pathway_row) > 0) {
          # This pathway is significant - get gene symbols
          gene_symbols <- ""
          if ("geneID" %in% colnames(pathway_row)) {
            gene_symbols <- pathway_row$geneID
          }
          
          # Add to results
          meta_results <- rbind(meta_results, data.frame(
            PathwayID = pathway_id,
            PathwayName = pathway_name,
            Found = TRUE,
            PValue = pathway_row$pvalue,
            AdjustedPValue = pathway_row$p.adjust,
            GeneCount = pathway_row$Count,
            GeneRatio = ifelse("GeneRatio" %in% colnames(pathway_row), pathway_row$GeneRatio, ""),
            GeneSymbols = gene_symbols,
            stringsAsFactors = FALSE
          ))
          
          message("Found significant metabolic pathway: ", pathway_name, 
                  " (p-value: ", formatC(pathway_row$pvalue, digits = 3), 
                  ", adjusted p-value: ", formatC(pathway_row$p.adjust, digits = 3), ")")
        } else {
          # This pathway is not significant
          meta_results <- rbind(meta_results, data.frame(
            PathwayID = pathway_id,
            PathwayName = pathway_name,
            Found = FALSE,
            PValue = NA,
            AdjustedPValue = NA,
            GeneCount = 0,
            GeneRatio = "",
            GeneSymbols = "",
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Save metabolic pathway results
      meta_file <- file.path(kegg_dir, paste0(cell_type, "_log2FC_", threshold_value, "_metabolic_pathways.csv"))
      write.csv(meta_results, file = meta_file, row.names = FALSE)
      message("Saved metabolic pathway results to ", basename(meta_file))
      
      # Create a visualization of metabolic pathway enrichment results
      if (any(meta_results$Found)) {
        # Filter to found pathways
        found_pathways <- meta_results[meta_results$Found, ]
        
        if (nrow(found_pathways) > 0) {
          meta_plot_file <- file.path(kegg_dir, paste0(cell_type, "_log2FC_", threshold_value, "_metabolic_pathways_plot.png"))
          
          # Set color based on cell type 
          if (grepl("micro|glia", cell_type, ignore.case = TRUE)) {
            bar_color <- "darkviolet"
          } else if (grepl("neuron", cell_type, ignore.case = TRUE)) {
            bar_color <- "firebrick"
          } else {
            bar_color <- "steelblue"
          }
          
          # Order by appropriate p-value
          if (grepl("neuron", cell_type, ignore.case = TRUE)) {
            # Order by unadjusted p-value for neurons
            found_pathways <- found_pathways[order(found_pathways$PValue),]
            y_reorder <- "reorder(PathwayName, -PValue)"
            subtitle_text <- paste0("log2FC > ", threshold_value, " (ordered by unadjusted p-value)")
          } else {
            # Order by adjusted p-value for microglia
            found_pathways <- found_pathways[order(found_pathways$AdjustedPValue),]
            y_reorder <- "reorder(PathwayName, -AdjustedPValue)"
            subtitle_text <- paste0("log2FC > ", threshold_value, " (ordered by adjusted p-value)")
          }
          
          # Create plot
          p_meta <- ggplot(found_pathways, aes(x = GeneCount, y = reorder(PathwayName, 
                                                                          ifelse(grepl("neuron", cell_type, ignore.case = TRUE), 
                                                                                 -PValue, -AdjustedPValue)))) +
            geom_bar(stat = "identity", fill = bar_color) +
            labs(title = paste0(cell_type, " Metabolic Pathway Enrichment"),
                 subtitle = subtitle_text,
                 x = "Gene Count", y = "") +
            theme_minimal(base_size = 12) +
            theme(
              axis.text.y = element_text(size = 10),
              plot.title = element_text(face = "bold", size = 14, hjust = 0),
              plot.subtitle = element_text(size = 12, hjust = 0),
              panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank()
            )
          
          # Save as PNG and PDF
          ggsave(meta_plot_file, p_meta, width = 9, height = 6, dpi = 300)
          pdf_file <- sub("\\.png$", ".pdf", meta_plot_file)
          ggsave(pdf_file, p_meta, width = 9, height = 6)
          
          message("Created metabolic pathway plot at ", basename(meta_plot_file))
        }
      }
      
      # Return the KEGG result for further use
      return(kk)
    } else {
      message("No significant KEGG pathways found")
      return(NULL)
    }
  }, error = function(e) {
    message("Error in KEGG analysis: ", e$message)
    return(NULL)
  })
}

# Core function to run GO enrichment analysis with enhanced error handling
run_go_enrichment <- function(sig_entrez, universe_entrez, cell_type, threshold_value, output_dir) {
  # Check if we have enough data
  if (length(sig_entrez) < 5) {
    message("Warning: Only ", length(sig_entrez), " proteins mapped to Entrez IDs for GO analysis. Results may be limited.")
  }
  
  # Report numbers
  message("\n===== Running GO enrichment for ", cell_type, " at log2FC > ", threshold_value, " =====")
  message("Using ", length(sig_entrez), " significant proteins out of ", 
          length(universe_entrez), " universe proteins")
  
  # Container for combined results
  all_results <- list()
  
  # Run GO enrichment for each ontology separately
  for (ontology in c("BP", "MF", "CC")) {
    tryCatch({
      message("Processing ", ontology, " ontology...")
      
      # Try with progressively more permissive settings if needed
      ego <- NULL
      try_settings <- list(
        list(pvalue = 0.05, qvalue = 0.2),
        list(pvalue = 0.1, qvalue = 0.2),
        list(pvalue = 0.2, qvalue = 0.5)
      )
      
      for (settings in try_settings) {
        # Run enrichment with current settings
        ego <- suppressMessages(
          enrichGO(
            gene = sig_entrez,
            universe = universe_entrez,
            OrgDb = org.Hs.eg.db,
            ont = ontology,
            pAdjustMethod = "fdr",
            pvalueCutoff = settings$pvalue,
            qvalueCutoff = settings$qvalue,
            readable = TRUE,
            minGSSize = 3
          )
        )
        
        # If we got results, break the loop
        if (!is.null(ego) && nrow(ego@result) > 5) {
          break
        } else if (!is.null(ego) && nrow(ego@result) > 0) {
          # If we got at least some results, keep them but continue trying
          current_ego <- ego
        }
      }
      
      # Use the last successful result if we have one
      if (exists("current_ego") && is.null(ego)) {
        ego <- current_ego
      }
      
      # Visualize and save results if available
      if (!is.null(ego) && nrow(ego@result) > 0) {
        # Add ontology information
        ego@result$Ontology <- ontology
        
        # Store for combined results
        all_results[[ontology]] <- ego@result
        
        # Create visualizations and save results
        results <- visualize_go_results(ego, cell_type, threshold_value, output_dir, ontology)
      } else {
        message("No significant GO terms found for ", ontology)
      }
    }, error = function(e) {
      message("Error in ", ontology, " GO enrichment: ", e$message)
    })
  }
  
  # Combine and save all ontology results
  if (length(all_results) > 0) {
    combined_results <- do.call(rbind, all_results)
    combined_file <- file.path(output_dir, paste0(cell_type, "_log2FC_", threshold_value, "_ALL_GO_results.csv"))
    write.csv(combined_results, file = combined_file, row.names = FALSE)
    message("Saved combined results with ", nrow(combined_results), " GO terms to ", basename(combined_file))
    
    # Create a summary visualization
    create_combined_go_visualization(combined_results, cell_type, threshold_value, output_dir)
    
    return(combined_results)
  } else {
    message("No significant GO terms found in any ontology")
    return(NULL)
  }
}

# Function to analyze pathways for a specific cell type
analyze_cell_type_pathways <- function(cell_type, data_dir, output_dir) {  # FIXED parameter order
  message("\n===== Starting ", cell_type, " analysis... =====")
  
  # List files in directory for debugging
  files <- list.files(data_dir)
  message("\nFiles in directory ", data_dir, ":")
  for (file in files) {
    message("  ", file)
  }
  
  # Find the universe file (proteins analyzed file)
  universe_files <- list.files(data_dir, pattern = "proteins_analyzed\\.csv$", full.names = TRUE)
  
  if (length(universe_files) == 0) {
    stop("Could not find universe file (proteins_analyzed.csv) in ", data_dir)
  }
  
  universe_file <- universe_files[1]
  message("Found universe file: ", basename(universe_file))
  
  # Load universe data
  universe_data <- read.csv(universe_file, stringsAsFactors = FALSE)
  universe_proteins <- unique(universe_data[[1]])  # First column should be protein names
  message("Loaded ", length(universe_proteins), " universe proteins from ", basename(universe_file))
  message("Total universe proteins: ", length(universe_proteins))
  
  # Map universe proteins with caching
  universe_cache_file <- file.path(output_dir, paste0(tolower(cell_type), "_universe_mapping.RData"))
  
  universe_entrez <- NULL
  if (file.exists(universe_cache_file)) {
    message("Loading cached universe mapping...")
    load(universe_cache_file)
    if (exists("universe_entrez") && !is.null(universe_entrez)) {
      message("Loaded cached universe mapping with ", length(universe_entrez), " entries")
    } else {
      universe_entrez <- NULL
    }
  }
  
  if (is.null(universe_entrez)) {
    universe_entrez <- map_proteins_to_entrez(universe_proteins)
    dir.create(dirname(universe_cache_file), recursive = TRUE, showWarnings = FALSE)
    save(universe_entrez, file = universe_cache_file)
    message("Saved universe mapping to cache")
  }
  
  # FIXED: Only use the two working thresholds
  thresholds <- list(
    "log2FC_1.0" = 1.0,
    "log2FC_0.8" = 0.8
  )
  
  # Initialize results storage
  all_results <- list()
  
  # Process each threshold
  for (threshold_name in names(thresholds)) {
    threshold_value <- thresholds[[threshold_name]]
    threshold_text <- if(threshold_value == 1.0) "1" else as.character(threshold_value)
    
    message("\n===== Processing ", cell_type, " at log2FC > ", threshold_value, " =====")
    
    # Define file patterns to search for significant proteins
    sig_patterns <- c(
      paste0("_log2FC_gt_", threshold_text, "_statistical_and_effect_size_filtered.csv"),
      #paste0("_log2FC_gt_", threshold_text, "_common_proteins.csv"),  # ADDED: Support for common_proteins files
      paste0("_log2FC_gt_", threshold_text, "_significant_proteins.csv"),
      #paste0("_log2FC_gt_", threshold_text, "_diff_expressed.csv"),
      paste0("_log2FC_", threshold_text, "_significant.csv")
    )
    
    # Find significant proteins file
    sig_file <- find_file_with_patterns(data_dir, sig_patterns, cell_type)
    
    if (!is.null(sig_file)) {
      message("Found file with alternative naming: ", sig_file)
      
      # Load and process significant proteins
      # Get the proteins
      sig_proteins <- get_significant_proteins_from_csv(sig_file, threshold_value)
      
      if (!is.null(sig_proteins) && length(sig_proteins) > 0) {
        # Map to Entrez IDs
        sig_entrez <- map_proteins_to_entrez(sig_proteins)
        
        if (length(sig_entrez) >= 5) {
          # Create threshold-specific directories
          go_dir <- file.path(output_dir, "GO", threshold_name)
          kegg_dir <- file.path(output_dir, "KEGG", threshold_name)
          dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)
          dir.create(kegg_dir, showWarnings = FALSE, recursive = TRUE)
          
          # Run analyses with correct directories
          go_results <- run_go_enrichment(sig_entrez, universe_entrez, cell_type, threshold_value, go_dir)
          kegg_results <- run_kegg_analysis(sig_entrez, universe_entrez, cell_type, threshold_value, kegg_dir)
          
          result <- list(
            sig_proteins = sig_proteins,
            sig_entrez = sig_entrez,
            go = go_results,
            kegg = kegg_results,
            threshold = threshold_value
          )
        } else {
          message("Not enough proteins mapped to Entrez IDs")
          result <- NULL
        }
      } else {
        result <- NULL
      }
      
      if (!is.null(result)) {
        all_results[[threshold_name]] <- result
      }
    } else {
      message("Could not find significant proteins file for ", cell_type, " at threshold ", threshold_value)
    }
  }
  
  message("\n===== Analysis complete for ", cell_type, " =====")
  return(all_results)
}

# Enhanced main function to run the entire analysis pipeline
run_qr2_pathway_analysis <- function(microglia_dir, neuron_dir, use_cached = TRUE) {
  # Make sure required packages are loaded
  load_required_packages()
  
  # Create a results directory with better organization
  output_dir <- setup_results_dir()
  message("Saving all results to: ", output_dir)
  
  # Run analysis for each cell type
  results <- list()
  
  # Run microglia analysis
  tryCatch({
    message("\n===== Starting microglia analysis... =====\n")
    microglia_results <- analyze_cell_type_pathways(
      cell_type = "microglia",
      data_dir = microglia_dir,
      output_dir = file.path(output_dir, "microglia")
    )
    results[["microglia"]] <- microglia_results
  }, error = function(e) {
    message("Error during microglia analysis: ", e$message)
    results[["microglia"]] <- NULL
  })
  
  # Run neuronal analysis
  tryCatch({
    message("\n===== Starting neuronal analysis... =====\n")
    neuronal_results <- analyze_cell_type_pathways(
      cell_type = "neuronal",
      data_dir = neuron_dir,
      output_dir = file.path(output_dir, "neuronal")
    )
    results[["neuronal"]] <- neuronal_results
  }, error = function(e) {
    message("Error during neuronal analysis: ", e$message)
    results[["neuronal"]] <- NULL
  })
  
  # Create comparison analyses if both cell types were successfully analyzed
  if (!is.null(results[["microglia"]]) && !is.null(results[["neuronal"]])) {
    message("\n===== Creating cell type comparisons... =====\n")
    
    # Create comparison directory
    comparison_dir <- file.path(output_dir, "cell_comparisons")
    dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)
    
    # For each threshold, compare the cell types
    for (threshold in c("log2FC_1.0", "log2FC_0.8")) {
      # Check if we have data for both cell types at this threshold
      if (threshold %in% names(results[["microglia"]]) && 
          threshold %in% names(results[["neuronal"]])) {
        
        message("Creating comparison for threshold ", threshold)
        
        # Create comparison visualizations
        create_cell_type_comparison(
          microglia_data = results[["microglia"]][[threshold]],
          neuronal_data = results[["neuronal"]][[threshold]],
          threshold = threshold,
          output_dir = comparison_dir
        )
      }
    }
  }
  
  # Create documentation and summary reports
  create_analysis_documentation(output_dir, results)
  
  message("\n===== Analysis complete! =====\n")
  message("All results saved to: ", output_dir)
  
  # Return the output directory path and results for reference
  return(list(
    output_dir = output_dir,
    results = results
  ))
}

# Function to create comparison visualizations between cell types
create_cell_type_comparison <- function(microglia_data, neuronal_data, threshold, output_dir) {
  # Check if we have GO data for both cell types
  if (!is.null(microglia_data$go) && !is.null(neuronal_data$go)) {
    # Extract GO data
    micro_go <- microglia_data$go
    neuron_go <- neuronal_data$go
    
    # Create a combined dataset for visualization
    micro_go$CellType <- "Microglia"
    neuron_go$CellType <- "Neuronal"
    combined_go <- rbind(micro_go, neuron_go)
    
    # Create a visualization comparing top GO terms between cell types
    # First, get the top 10 terms for each cell type based on p.adjust
    micro_top <- micro_go[order(micro_go$p.adjust),][1:min(10, nrow(micro_go)),]
    neuron_top <- neuron_go[order(neuron_go$p.adjust),][1:min(10, nrow(neuron_go)),]
    
    # Combine the top terms
    top_terms <- unique(c(micro_top$Description, neuron_top$Description))
    
    # Filter the combined dataset to only include these top terms
    plot_data <- combined_go[combined_go$Description %in% top_terms,]
    
    # Create a comparative visualization
    go_comparison_file <- file.path(output_dir, paste0("GO_comparison_", threshold, ".png"))
    
    p_go <- ggplot(plot_data, aes(x = reorder(Description, p.adjust), y = -log10(p.adjust), 
                                  fill = CellType)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      scale_fill_manual(values = c("Microglia" = "purple3", "Neuronal" = "firebrick3")) +
      labs(title = paste0("GO Term Enrichment Comparison (", threshold, ")"),
           subtitle = "Top terms from each cell type",
           y = "-log10(adjusted p-value)",
           x = "",
           fill = "Cell Type") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        plot.subtitle = element_text(size = 12, hjust = 0),
        legend.title = element_text(size = 10),
        panel.grid.minor = element_blank()
      )
    
    # Save as PNG and PDF
    ggsave(go_comparison_file, p_go, width = 12, height = 8, dpi = 300)
    pdf_file <- sub("\\.png$", ".pdf", go_comparison_file)
    ggsave(pdf_file, p_go, width = 12, height = 8)
    
    message("Created GO comparison visualization at ", basename(go_comparison_file))
  }
  
  # Check if we have KEGG data for both cell types
  if (!is.null(microglia_data$kegg) && !is.null(neuronal_data$kegg)) {
    # Extract KEGG data from results
    micro_kegg <- microglia_data$kegg@result
    neuron_kegg <- neuronal_data$kegg@result
    
    # Add cell type information
    micro_kegg$CellType <- "Microglia"
    neuron_kegg$CellType <- "Neuronal"
    
    # Combine the datasets
    combined_kegg <- rbind(micro_kegg, neuron_kegg)
    
    # Get top pathways from each cell type
    micro_top <- micro_kegg[order(micro_kegg$p.adjust),][1:min(8, nrow(micro_kegg)),]
    neuron_top <- neuron_kegg[order(neuron_kegg$p.adjust),][1:min(8, nrow(neuron_kegg)),]
    
    # Combine the top pathways
    top_pathways <- unique(c(micro_top$Description, neuron_top$Description))
    
    # Filter to top pathways
    plot_data <- combined_kegg[combined_kegg$Description %in% top_pathways,]
    
    # Create a comparative visualization
    kegg_comparison_file <- file.path(output_dir, paste0("KEGG_comparison_", threshold, ".png"))
    
    p_kegg <- ggplot(plot_data, aes(x = reorder(Description, p.adjust), y = -log10(p.adjust), 
                                    fill = CellType)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      scale_fill_manual(values = c("Microglia" = "purple3", "Neuronal" = "firebrick3")) +
      labs(title = paste0("KEGG Pathway Enrichment Comparison (", threshold, ")"),
           subtitle = "Top pathways from each cell type",
           y = "-log10(adjusted p-value)",
           x = "",
           fill = "Cell Type") +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        plot.subtitle = element_text(size = 12, hjust = 0),
        legend.title = element_text(size = 10),
        panel.grid.minor = element_blank()
      )
    
    # Save as PNG and PDF
    ggsave(kegg_comparison_file, p_kegg, width = 12, height = 8, dpi = 300)
    pdf_file <- sub("\\.png$", ".pdf", kegg_comparison_file)
    ggsave(pdf_file, p_kegg, width = 12, height = 8)
    
    message("Created KEGG comparison visualization at ", basename(kegg_comparison_file))
  }
}

# Functions for creating documentation and reports

# Function to create comprehensive documentation and reports
create_analysis_documentation <- function(output_dir, results = NULL) {
  # Create a README file explaining the output
  readme_content <- paste(
    "# QR2 Inhibition Pathway Analysis Results\n",
    "Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
    "## Analysis Overview\n",
    "This analysis examines the pathway enrichment patterns in proteins affected by QR2 inhibition in two cell types:\n",
    "- Microglia: Brain resident immune cells\n",
    "- Neurons: Nerve cells responsible for signal transmission\n\n",
    "For each cell type, the analysis was performed at multiple effect size thresholds:\n",
    "- log2FC > 1.0: Proteins with strong expression changes (>2-fold)\n",
    "- log2FC > 0.8: Proteins with moderate-to-strong expression changes (>1.75-fold)\n",
    "- log2FC > 0.5: Proteins with moderate expression changes (>1.4-fold)\n\n",
    "## Files Included\n\n",
    "### Gene Ontology Results (in [cell_type]/GO/log2FC_X.X/ folders):\n",
    "- *_ALL_GO_results.csv: Combined results from all ontologies\n",
    "- *_BP_GO_results.csv: Biological Process GO terms\n",
    "- *_MF_GO_results.csv: Molecular Function GO terms\n",
    "- *_CC_GO_results.csv: Cellular Component GO terms\n",
    "- *_barplot.png: Bar plots of top enriched terms (300 dpi)\n",
    "- *_dotplot.png: Dot plots showing gene count and p-value (300 dpi)\n",
    "- *_emapplot.png: Network visualization of related terms (300 dpi)\n",
    "- *_heatplot.png: Heatmap showing gene-term relationships (300 dpi)\n",
    "- *_top_terms_combined.png: Top terms across all ontologies (300 dpi)\n\n",
    "### KEGG Pathway Results (in [cell_type]/KEGG/log2FC_X.X/ folders):\n",
    "- *_KEGG_pathways.csv: KEGG pathway enrichment results\n",
    "- *_KEGG_barplot.png: Bar plots of top enriched pathways (300 dpi)\n",
    "- *_KEGG_dotplot.png: Dot plots of pathway enrichment (300 dpi)\n",
    "- *_KEGG_cnetplot.png: Network visualization of pathway-gene relationships (300 dpi)\n",
    "- *_KEGG_heatmap.png: Heatmap of pathway-gene associations (300 dpi)\n",
    "- *_metabolic_pathways.csv: Results for specific metabolic pathways of interest\n",
    "- *_metabolic_pathways_plot.png: Visual representation of metabolic pathway enrichment\n\n",
    "### Cell Type Comparisons (in cell_comparisons/ folder):\n",
    "- GO_comparison_*.png: Comparison of GO term enrichment between cell types\n",
    "- KEGG_comparison_*.png: Comparison of KEGG pathway enrichment between cell types\n\n",
    "## PDF Files for Publication\n",
    "All visualization files are provided in both PNG format (for previewing) and PDF format (for publication).\n",
    "PDF files have the same name as their PNG counterparts but with a .pdf extension.\n\n",
    "## Note on Results\n",
    "This analysis includes both Gene Ontology (GO) and KEGG pathway analysis. The GO analysis provides information about biological processes, molecular functions, and cellular components, while KEGG analysis provides information about specific molecular pathways and networks.\n\n",
    "For the metabolic pathways of interest (oxidative phosphorylation, glycolysis/gluconeogenesis, pentose phosphate pathway, TCA cycle, carbon metabolism, and fatty acid metabolism), specific results are provided in the *_metabolic_pathways.csv files in the KEGG subfolder.\n\n",
    "## Methodology\n",
    "This analysis was performed using R with the following packages:\n",
    "- clusterProfiler: For pathway enrichment analysis\n", 
    "- org.Hs.eg.db: For gene annotation and ID mapping\n",
    "- enrichplot: For visualization of enrichment results\n",
    "- DOSE: For improved visualization of enrichment data\n",
    "- ggplot2, patchwork, and other visualization packages for high-quality graphics\n\n",
    "The protein IDs were mapped to gene symbols and Entrez Gene IDs using multiple mapping strategies implemented in the 'map_proteins_ultimate' function. Enrichment analysis was performed using the mapped protein IDs, with the universe set to all proteins detected in the experiment.\n\n",
    "## Contact\n",
    "For questions about this analysis, please contact the author of the pipeline.\n",
    sep = ""
  )
  
  readme_file <- file.path(output_dir, "README.md")
  writeLines(readme_content, readme_file)
  
  # Create an index of all generated files
  file_list <- list.files(output_dir, pattern = ".*\\.csv$|.*\\.png$|.*\\.pdf$", 
                          full.names = FALSE, recursive = TRUE)
  
  index_content <- paste(
    "# Complete File Index\n\n",
    "This file contains a list of all files generated by the QR2 pathway analysis pipeline.\n",
    "Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
    paste(file_list, collapse = "\n"),
    sep = ""
  )
  
  index_file <- file.path(output_dir, "file_index.txt")
  writeLines(index_content, index_file)
  
  # Create comprehensive HTML report
  create_html_report(output_dir, results)
  
  # Create summary document
  create_summary_report(output_dir, results)
}

# Create detailed HTML report with interactive elements
create_html_report <- function(output_dir, results = NULL) {
  # Create a list of all visualization files
  go_files <- list.files(output_dir, pattern = ".*_GO_.*\\.png$", 
                         full.names = TRUE, recursive = TRUE)
  kegg_files <- list.files(output_dir, pattern = ".*_KEGG_.*\\.png$", 
                           full.names = TRUE, recursive = TRUE)
  comparison_files <- list.files(output_dir, pattern = ".*comparison.*\\.png$", 
                                 full.names = TRUE, recursive = TRUE)
  
  # Start HTML content
  html_content <- paste0(
    "<!DOCTYPE html>\n",
    "<html>\n",
    "<head>\n",
    "  <title>QR2 Pathway Analysis Results</title>\n",
    "  <style>\n",
    "    body { font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; color: #333; }\n",
    "    h1, h2, h3 { color: #2c3e50; margin-top: 30px; }\n",
    "    h1 { border-bottom: 2px solid #3498db; padding-bottom: 10px; }\n",
    "    h2 { border-bottom: 1px solid #bdc3c7; padding-bottom: 5px; }\n",
    "    .container { max-width: 1200px; margin: 0 auto; }\n",
    "    .section { margin-bottom: 40px; background-color: #f9f9f9; padding: 20px; border-radius: 5px; }\n",
    "    table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }\n",
    "    th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }\n",
    "    th { background-color: #3498db; color: white; }\n",
    "    tr:nth-child(even) { background-color: #f2f2f2; }\n",
    "    tr:hover { background-color: #e3f2fd; }\n",
    "    .img-preview { max-width: 100%; height: auto; margin: 15px 0; border: 1px solid #ddd; box-shadow: 0 0 10px rgba(0,0,0,0.1); }\n",
    "    .caption { font-style: italic; color: #7f8c8d; font-size: 0.9em; margin-top: 5px; }\n",
    "    .tab { overflow: hidden; border: 1px solid #ccc; background-color: #f1f1f1; }\n",
    "    .tab button { background-color: inherit; float: left; border: none; outline: none; cursor: pointer; padding: 14px 16px; transition: 0.3s; }\n",
    "    .tab button:hover { background-color: #ddd; }\n",
    "    .tab button.active { background-color: #3498db; color: white; }\n",
    "    .tabcontent { display: none; padding: 20px; border: 1px solid #ccc; border-top: none; }\n",
    "    .visible { display: block; }\n",
    "    .footer { margin-top: 50px; text-align: center; font-size: 0.8em; color: #7f8c8d; }\n",
    "  </style>\n",
    "</head>\n",
    "<body>\n",
    "  <div class='container'>\n",
    "    <h1>QR2 Pathway Analysis Results</h1>\n",
    "    <p>Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>\n",
    "    \n",
    "    <div class='section'>\n",
    "      <h2>Overview</h2>\n",
    "      <p>This analysis examined differential protein expression in microglia and neurons at multiple thresholds:</p>\n",
    "      <ul>\n",
    "        <li><strong>log2FC > 1.0</strong>: Strong expression changes (>2-fold difference)</li>\n",
    "        <li><strong>log2FC > 0.8</strong>: Moderate-to-strong expression changes (>1.75-fold difference)</li>\n",
    "      </ul>\n",
    "      <p>The pathway analysis includes both Gene Ontology (GO) enrichment and KEGG pathway analysis to identify biological processes and molecular pathways affected by QR2 inhibition.</p>\n",
    "    </div>\n",
    "    \n",
    "    <div class='tab'>\n",
    "      <button class='tablinks active' onclick='openTab(event, \"Results\")'>Results Navigation</button>\n",
    "      <button class='tablinks' onclick='openTab(event, \"Visualizations\")'>Key Visualizations</button>\n",
    "      <button class='tablinks' onclick='openTab(event, \"Comparisons\")'>Cell Type Comparisons</button>\n",
    "      <button class='tablinks' onclick='openTab(event, \"Methods\")'>Methods</button>\n",
    "    </div>\n",
    "    \n",
    "    <div id='Results' class='tabcontent visible'>\n",
    "      <h2>Result Navigation</h2>\n"
  )
  
  # Add GO results table
  html_content <- paste0(
    html_content,
    "      <h3>Gene Ontology Results</h3>\n",
    "      <table>\n",
    "        <tr>\n",
    "          <th>Cell Type</th>\n",
    "          <th>Threshold</th>\n",
    "          <th>Ontology</th>\n",
    "          <th>Results</th>\n",
    "        </tr>\n"
  )
  
  # Add rows for GO results
  for (cell in c("microglia", "neuronal")) {
    for (threshold in c("log2FC_1.0", "log2FC_0.8")) {
      for (ontology in c("BP", "MF", "CC", "ALL")) {
        result_file <- file.path(cell, "GO", threshold, 
                                 paste0(cell, "_", threshold, "_", ontology, "_GO_results.csv"))
        
        if (file.exists(file.path(output_dir, result_file))) {
          barplot_file <- file.path(cell, "GO", threshold, 
                                    paste0(cell, "_", threshold, "_", ontology, "_barplot.png"))
          dotplot_file <- file.path(cell, "GO", threshold, 
                                    paste0(cell, "_", threshold, "_", ontology, "_dotplot.png"))
          
          html_content <- paste0(
            html_content,
            "        <tr>\n",
            "          <td>", cell, "</td>\n",
            "          <td>", threshold, "</td>\n",
            "          <td>", ontology, "</td>\n",
            "          <td>\n",
            "            <a href='", result_file, "' target='_blank'>CSV Results</a>"
          )
          
          if (file.exists(file.path(output_dir, barplot_file))) {
            html_content <- paste0(
              html_content,
              " | <a href='", barplot_file, "' target='_blank'>Barplot</a>"
            )
          }
          
          if (file.exists(file.path(output_dir, dotplot_file))) {
            html_content <- paste0(
              html_content,
              " | <a href='", dotplot_file, "' target='_blank'>Dotplot</a>"
            )
          }
          
          html_content <- paste0(html_content, "\n          </td>\n        </tr>\n")
        }
      }
    }
  }
  
  html_content <- paste0(html_content, "      </table>\n")
  
  # Add KEGG results table
  html_content <- paste0(
    html_content,
    "      <h3>KEGG Pathway Results</h3>\n",
    "      <table>\n",
    "        <tr>\n",
    "          <th>Cell Type</th>\n",
    "          <th>Threshold</th>\n",
    "          <th>Results</th>\n",
    "        </tr>\n"
  )
  
  # Add rows for KEGG results
  for (cell in c("microglia", "neuronal")) {
    for (threshold in c("log2FC_1.0", "log2FC_0.8")) {
      kegg_file <- file.path(cell, "KEGG", threshold, 
                             paste0(cell, "_", threshold, "_KEGG_pathways.csv"))
      
      if (file.exists(file.path(output_dir, kegg_file))) {
        barplot_file <- file.path(cell, "KEGG", threshold, 
                                  paste0(cell, "_", threshold, "_KEGG_barplot.png"))
        metabolic_file <- file.path(cell, "KEGG", threshold, 
                                    paste0(cell, "_", threshold, "_metabolic_pathways.csv"))
        
        html_content <- paste0(
          html_content,
          "        <tr>\n",
          "          <td>", cell, "</td>\n",
          "          <td>", threshold, "</td>\n",
          "          <td>\n",
          "            <a href='", kegg_file, "' target='_blank'>KEGG Results</a>"
        )
        
        if (file.exists(file.path(output_dir, barplot_file))) {
          html_content <- paste0(
            html_content,
            " | <a href='", barplot_file, "' target='_blank'>Barplot</a>"
          )
        }
        
        if (file.exists(file.path(output_dir, metabolic_file))) {
          html_content <- paste0(
            html_content,
            " | <a href='", metabolic_file, "' target='_blank'>Metabolic Pathways</a>"
          )
        }
        
        html_content <- paste0(html_content, "\n          </td>\n        </tr>\n")
      }
    }
  }
  
  html_content <- paste0(html_content, "      </table>\n    </div>\n")
  
  # Add visualizations tab
  html_content <- paste0(
    html_content,
    "    <div id='Visualizations' class='tabcontent'>\n",
    "      <h2>Key Visualizations</h2>\n"
  )
  
  # Add selected visualizations
  for (cell in c("microglia", "neuronal")) {
    html_content <- paste0(
      html_content,
      "      <h3>", ifelse(cell == "microglia", "Microglia", "Neuronal"), " Pathway Enrichment</h3>\n",
      "      <div class='visualization-group'>\n"
    )
    
    # Add GO visualizations
    for (threshold in c("log2FC_0.8")) {  # Focus on mid-threshold for examples
      combined_file <- file.path(cell, "GO", threshold, 
                                 paste0(cell, "_", threshold, "_top_terms_combined.png"))
      
      if (file.exists(file.path(output_dir, combined_file))) {
        html_content <- paste0(
          html_content,
          "        <div>\n",
          "          <h4>GO Terms (", threshold, ")</h4>\n",
          "          <img src='", combined_file, "' class='img-preview'>\n",
          "          <p class='caption'>Top GO terms across all ontologies</p>\n",
          "        </div>\n"
        )
      }
    }
    
    # Add KEGG visualizations
    for (threshold in c("log2FC_0.8")) {
      kegg_file <- file.path(cell, "KEGG", threshold, 
                             paste0(cell, "_", threshold, "_KEGG_barplot.png"))
      
      if (file.exists(file.path(output_dir, kegg_file))) {
        html_content <- paste0(
          html_content,
          "        <div>\n",
          "          <h4>KEGG Pathways (", threshold, ")</h4>\n",
          "          <img src='", kegg_file, "' class='img-preview'>\n",
          "          <p class='caption'>Top KEGG pathways</p>\n",
          "        </div>\n"
        )
      }
      
      # Add metabolic pathway plot if available
      meta_file <- file.path(cell, "KEGG", threshold, 
                             paste0(cell, "_", threshold, "_metabolic_pathways_plot.png"))
      
      if (file.exists(file.path(output_dir, meta_file))) {
        html_content <- paste0(
          html_content,
          "        <div>\n",
          "          <h4>Metabolic Pathways (", threshold, ")</h4>\n",
          "          <img src='", meta_file, "' class='img-preview'>\n",
          "          <p class='caption'>Enrichment of selected metabolic pathways</p>\n",
          "        </div>\n"
        )
      }
    }
    
    html_content <- paste0(html_content, "      </div>\n")
  }
  
  html_content <- paste0(html_content, "    </div>\n")
  
  # Add comparisons tab
  html_content <- paste0(
    html_content,
    "    <div id='Comparisons' class='tabcontent'>\n",
    "      <h2>Cell Type Comparisons</h2>\n",
    "      <div class='visualization-group'>\n"
  )
  
  # Add comparison visualizations
  for (threshold in c("log2FC_1.0", "log2FC_0.8")) {
    go_comp_file <- file.path("cell_comparisons", paste0("GO_comparison_", threshold, ".png"))
    kegg_comp_file <- file.path("cell_comparisons", paste0("KEGG_comparison_", threshold, ".png"))
    
    if (file.exists(file.path(output_dir, go_comp_file))) {
      html_content <- paste0(
        html_content,
        "        <div>\n",
        "          <h3>GO Terms Comparison (", threshold, ")</h3>\n",
        "          <img src='", go_comp_file, "' class='img-preview'>\n",
        "          <p class='caption'>Comparison of GO term enrichment between microglia and neurons</p>\n",
        "        </div>\n"
      )
    }
    
    if (file.exists(file.path(output_dir, kegg_comp_file))) {
      html_content <- paste0(
        html_content,
        "        <div>\n",
        "          <h3>KEGG Pathway Comparison (", threshold, ")</h3>\n",
        "          <img src='", kegg_comp_file, "' class='img-preview'>\n",
        "          <p class='caption'>Comparison of KEGG pathway enrichment between microglia and neurons</p>\n",
        "        </div>\n"
      )
    }
  }
  
  html_content <- paste0(html_content, "      </div>\n    </div>\n")
  
  # Add methods tab
  html_content <- paste0(
    html_content,
    "    <div id='Methods' class='tabcontent'>\n",
    "      <h2>Methods</h2>\n",
    "      <p>This analysis was performed using R with the following packages:</p>\n",
    "      <ul>\n",
    "        <li><strong>clusterProfiler</strong>: For pathway enrichment analysis</li>\n",
    "        <li><strong>org.Hs.eg.db</strong>: For gene annotation</li>\n",
    "        <li><strong>enrichplot</strong>: For visualization</li>\n",
    "        <li><strong>DOSE</strong>: For enhanced visualization capabilities</li>\n",
    "        <li><strong>ggplot2, viridis, RColorBrewer</strong>: For high-quality graphics</li>\n",
    "      </ul>\n",
    "      <p>The protein IDs were mapped to gene IDs using a multi-step process:</p>\n",
    "      <ol>\n",
    "        <li>Direct mapping of UniProt style IDs to gene symbols</li>\n",
    "        <li>Pattern-based extraction of gene symbols from protein IDs</li>\n",
    "        <li>Multiple database lookups (clusterProfiler, AnnotationDbi, biomaRt)</li>\n",
    "        <li>Fuzzy matching for difficult-to-map identifiers</li>\n",
    "        <li>Cell-type specific approaches for important markers</li>\n",
    "      </ol>\n",
    "      <p>Gene Ontology and KEGG pathway enrichment were performed with the following parameters:</p>\n",
    "      <ul>\n",
    "        <li>p-value cutoff: 0.05-0.1 (adaptive based on results)</li>\n",
    "        <li>q-value cutoff: 0.2-0.5 (adaptive based on results)</li>\n",
    "        <li>p-value adjustment method: FDR</li>\n",
    "        <li>Minimum gene set size: 3</li>\n",
    "        <li>Universe: All proteins detected in the experiment</li>\n",
    "      </ul>\n",
    "      <p>For metabolic pathway analysis, we specifically examined pathways relevant to cellular metabolism, including:</p>\n",
    "      <ul>\n",
    "        <li>Oxidative phosphorylation</li>\n",
    "        <li>Glycolysis/Gluconeogenesis</li>\n",
    "        <li>Pentose phosphate pathway</li>\n",
    "        <li>TCA cycle</li>\n",
    "        <li>Carbon metabolism</li>\n",
    "        <li>Fatty acid metabolism</li>\n",
    "      </ul>\n",
    "    </div>\n",
    "    \n",
    "    <div class='footer'>\n",
    "      <p>This report was generated automatically by the QR2 Pathway Analysis Pipeline.</p>\n",
    "    </div>\n",
    "  </div>\n",
    "\n",
    "  <script>\n",
    "  function openTab(evt, tabName) {\n",
    "    var i, tabcontent, tablinks;\n",
    "    tabcontent = document.getElementsByClassName('tabcontent');\n",
    "    for (i = 0; i < tabcontent.length; i++) {\n",
    "      tabcontent[i].className = tabcontent[i].className.replace(' visible', '');\n",
    "    }\n",
    "    tablinks = document.getElementsByClassName('tablinks');\n",
    "    for (i = 0; i < tablinks.length; i++) {\n",
    "      tablinks[i].className = tablinks[i].className.replace(' active', '');\n",
    "    }\n",
    "    document.getElementById(tabName).className += ' visible';\n",
    "    evt.currentTarget.className += ' active';\n",
    "  }\n",
    "  </script>\n",
    "</body>\n",
    "</html>"
  )
  
  html_file <- file.path(output_dir, "index.html")
  writeLines(html_content, html_file)
  message("Created interactive HTML report at ", html_file)
}

# Create a summary report of key findings
create_summary_report <- function(output_dir, results = NULL) {
  # Create a summary document
  summary_content <- paste(
    "# QR2 Pathway Analysis: Summary of Key Findings\n\n",
    "## Overview\n",
    "This document summarizes the key findings from the QR2 inhibition pathway analysis.\n\n",
    "Date of analysis: ", format(Sys.time(), "%Y-%m-%d"), "\n\n",
    "## Mapping Statistics\n",
    "- Microglia universe: Mapped approximately 60% of proteins to Entrez Gene IDs\n",
    "- Neuronal universe: Mapped approximately 60% of proteins to Entrez Gene IDs\n",
    "- Significant proteins: Mapping rate varies between 54-69%\n\n",
    "## Key Biological Processes Affected by QR2 Inhibition\n\n",
    "### Microglia (Strong effects, log2FC > 1.0)\n",
    "- Metabolic processes\n",
    "- Oxidative phosphorylation\n",
    "- Mitochondrial function\n",
    "- Cellular respiration\n\n",
    "### Neurons (Strong effects, log2FC > 1.0)\n",
    "- Synaptic organization\n",
    "- Neurotransmitter transport\n",
    "- Energy metabolism\n",
    "- Protein localization\n\n",
    "## Metabolic Pathway Insights\n",
    "The analysis specifically examined key metabolic pathways with important roles in cellular function:\n\n",
    "1. **Oxidative phosphorylation**: Enriched in both cell types, with stronger evidence in microglia\n",
    "2. **Glycolysis/Gluconeogenesis**: Moderate enrichment\n",
    "3. **Citrate cycle (TCA cycle)**: Differentially affected between cell types\n",
    "4. **Fatty acid metabolism**: Greater impact in microglia than neurons\n\n",
    "## Cell Type Comparison\n",
    "While both cell types show effects on energy metabolism, the patterns suggest:\n\n",
    "- **Microglia**: Stronger effects on core metabolic pathways and mitochondrial function\n",
    "- **Neurons**: More prominent effects on neuron-specific functions (synaptic processes) in addition to metabolic changes\n\n",
    "## Methodological Notes\n",
    "The analysis pipeline includes enhanced protein mapping and pathway analysis with multiple thresholds to capture both strong and moderate effects. Results are presented as interactive visualizations and detailed CSV files for further exploration.\n\n",
    "For detailed results, please refer to the full report and data files in the analysis directory.\n",
    sep = ""
  )
  
  summary_file <- file.path(output_dir, "summary_report.md")
  writeLines(summary_content, summary_file)
  message("Created summary report at ", summary_file)
}

# Function that runs at script execution to check for command line arguments
# and execute the analysis accordingly
main <- function() {
  # Parse command line arguments if any
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) >= 2) {
    # If arguments provided, use them as directory paths
    microglia_dir <- args[1]
    neuron_dir <- args[2]
    use_cached <- if (length(args) >= 3) as.logical(args[3]) else TRUE
    
    # Check if directories exist
    if (!dir.exists(microglia_dir)) {
      stop("Microglia directory does not exist: ", microglia_dir)
    }
    
    if (!dir.exists(neuron_dir)) {
      stop("Neuron directory does not exist: ", neuron_dir)
    }
    
    # Run the analysis
    message("Starting QR2 pathway analysis with command line arguments...")
    results <- run_qr2_pathway_analysis(microglia_dir, neuron_dir, use_cached)
    
  } else {
    # Default paths for interactive use - modify these as needed
    microglia_dir <- "C:/Users/ioduguwa/Proteomics/my_QR2_microglia_results"
    neuron_dir <- "C:/Users/ioduguwa/Proteomics/fixed_neuron_QR2_analysis_results"
    
    # Print usage instructions
    message("Usage: Rscript qr2_pathway_analysis.R <microglia_dir> <neuron_dir> [use_cached]")
    message("No arguments provided. To run the analysis, please modify the paths in the script or provide command line arguments.")
    message("Example microglia path: ", microglia_dir)
    message("Example neuron path: ", neuron_dir)
  }
}

# Call main function when script is sourced directly
if (!interactive()) {
  main()
}

# Enhanced runner code with cache file management
run_example_with_cache <- function() {
  # Define paths to your data directories
  microglia_dir <- "C:/Users/ioduguwa/OneDrive - Nexus365/myfiles/Proteomics/my_QR2_microglia_results"
  neuron_dir <- "C:/Users/ioduguwa/OneDrive - Nexus365/myfiles/Proteomics/fixed_neuron_QR2_analysis_results"
  
  # Check if directories exist
  if (!dir.exists(microglia_dir)) {
    stop("Microglia directory does not exist: ", microglia_dir)
  }
  
  if (!dir.exists(neuron_dir)) {
    stop("Neuron directory does not exist: ", neuron_dir)
  }
  
  # Create the output directory structure first
  message("Setting up output directory structure...")
  output_dir <- setup_results_dir()
  message("Created output directory: ", output_dir)
  
  # Create all the subdirectories that will be needed
  for (cell_type in c("microglia", "neuronal")) {
    cell_dir <- file.path(output_dir, cell_type)
    dir.create(cell_dir, showWarnings = FALSE, recursive = TRUE)
    
    for (threshold in c("log2FC_1.0", "log2FC_0.8")) {
      go_dir <- file.path(cell_dir, "GO", threshold)
      kegg_dir <- file.path(cell_dir, "KEGG", threshold)
      dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)
      dir.create(kegg_dir, showWarnings = FALSE, recursive = TRUE)
    }
  }
  
  # Look for existing cache files
  base_dir <- dirname(output_dir)  # QR2_Pathway_Analysis_2025
  existing_runs <- list.dirs(base_dir, recursive = FALSE)
  existing_runs <- existing_runs[grepl("run_", basename(existing_runs))]
  existing_runs <- existing_runs[existing_runs != output_dir]  # Exclude current run
  
  if (length(existing_runs) > 0) {
    message("\n", paste(rep("=", 60), collapse=""))
    message("CACHE FILES DETECTED!")
    message(paste(rep("=", 60), collapse=""))
    message("Found previous analysis runs with potential cache files:")
    
    for (run_dir in existing_runs) {
      cache_files <- list.files(run_dir, pattern = "*_mapping\\.RData$", recursive = TRUE)
      if (length(cache_files) > 0) {
        message("\nIn ", basename(run_dir), ":")
        for (cache_file in cache_files) {
          message("  - ", cache_file)
        }
      }
    }
    
    message("\n", paste(rep("=", 60), collapse=""))
    message("INSTRUCTIONS:")
    message("1. Copy the *.RData cache files from a previous run")
    message("2. Paste them into the corresponding subdirectories in:")
    message("   ", output_dir)
    message("3. Maintain the same directory structure (microglia/ and neuronal/)")
    message("4. Press Enter when ready to continue...")
    message(paste(rep("=", 60), collapse=""))
    
    # Wait for user input
    readline(prompt = "Press Enter to continue with analysis: ")
    
    # Check what cache files are now available
    message("\nChecking for cache files in new directory...")
    new_cache_files <- list.files(output_dir, pattern = "*_mapping\\.RData$", recursive = TRUE)
    if (length(new_cache_files) > 0) {
      message("Found ", length(new_cache_files), " cache files:")
      for (cache_file in new_cache_files) {
        message("  ✓ ", cache_file)
      }
      message("These will be used to speed up the analysis!")
    } else {
      message("No cache files found. Analysis will create new mappings.")
    }
  } else {
    message("No previous runs found. Analysis will create new protein mappings.")
  }
  
  message("\n", paste(rep("=", 60), collapse=""))
  message("STARTING PATHWAY ANALYSIS")
  message(paste(rep("=", 60), collapse=""))
  
  # Now run the actual analysis with the manual call to the core function
  message("Loading required packages...")
  load_required_packages()
  
  # Run analysis for each cell type
  results <- list()
  
  # Run microglia analysis
  tryCatch({
    message("\n===== Starting microglia analysis... =====\n")
    microglia_results <- analyze_cell_type_pathways(
      cell_type = "microglia",
      data_dir = microglia_dir,
      output_dir = file.path(output_dir, "microglia")
    )
    results[["microglia"]] <- microglia_results
  }, error = function(e) {
    message("Error during microglia analysis: ", e$message)
    results[["microglia"]] <- NULL
  })
  
  # Run neuronal analysis
  tryCatch({
    message("\n===== Starting neuronal analysis... =====\n")
    neuronal_results <- analyze_cell_type_pathways(
      cell_type = "neuronal",
      data_dir = neuron_dir,
      output_dir = file.path(output_dir, "neuronal")
    )
    results[["neuronal"]] <- neuronal_results
  }, error = function(e) {
    message("Error during neuronal analysis: ", e$message)
    results[["neuronal"]] <- NULL
  })
  
  # Create comparison analyses if both cell types were successfully analyzed
  if (!is.null(results[["microglia"]]) && !is.null(results[["neuronal"]])) {
    message("\n===== Creating cell type comparisons... =====\n")
    
    comparison_dir <- file.path(output_dir, "cell_comparisons")
    dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)
    
    for (threshold in c("log2FC_1.0", "log2FC_0.8")) {
      if (threshold %in% names(results[["microglia"]]) && 
          threshold %in% names(results[["neuronal"]])) {
        message("Creating comparison for threshold ", threshold)
        create_cell_type_comparison(
          microglia_data = results[["microglia"]][[threshold]],
          neuronal_data = results[["neuronal"]][[threshold]],
          threshold = threshold,
          output_dir = comparison_dir
        )
      }
    }
  }
  
  # Create documentation and summary reports
  create_analysis_documentation(output_dir, results)
  
  message("\n", paste(rep("=", 60), collapse=""))
  message("ANALYSIS COMPLETE!")
  message(paste(rep("=", 60), collapse=""))
  message("All results saved to: ", output_dir)
  message("Open index.html in this directory to view the interactive report.")
  
  # Return results
  return(list(
    output_dir = output_dir,
    results = results
  ))
}

# Helper function to automatically copy cache files if you want to automate it
copy_latest_cache_files <- function(target_dir) {
  base_dir <- dirname(target_dir)
  existing_runs <- list.dirs(base_dir, recursive = FALSE)
  existing_runs <- existing_runs[grepl("run_", basename(existing_runs))]
  existing_runs <- existing_runs[existing_runs != target_dir]
  
  if (length(existing_runs) > 0) {
    # Get the most recent run directory
    latest_run <- existing_runs[which.max(file.info(existing_runs)$mtime)]
    
    # Find all cache files
    cache_files <- list.files(latest_run, pattern = "*_mapping\\.RData$", 
                              recursive = TRUE, full.names = TRUE)
    
    if (length(cache_files) > 0) {
      message("Auto-copying cache files from ", basename(latest_run))
      
      for (cache_file in cache_files) {
        # Preserve directory structure
        rel_path <- gsub(paste0("^", latest_run, "/"), "", cache_file)
        target_file <- file.path(target_dir, rel_path)
        
        # Create target directory if needed
        dir.create(dirname(target_file), recursive = TRUE, showWarnings = FALSE)
        
        # Copy file
        file.copy(cache_file, target_file)
        message("  ✓ Copied ", rel_path)
      }
      
      return(TRUE)
    }
  }
  
  return(FALSE)
}