###############################################################################
#
# Title: Complex Heatmap of methylation data of several clocks with PAM optimization
# Date: 26/07/2024
# Author: Catia Antunes
#
# Description:
# Complex heatmap of methylation data of top CpG sites identified by Clocks
# 1, 2, 3 and overlaps further grouped by tissue and including PAM clustering 
# optimization..
#
# Based on the code from:
# https://github.com/kevinblighe/E-MTAB-6141?tab=readme-ov-file
#
###############################################################################



#### Load Libraries 
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)
library(circlize)
library(cluster)
library(magick)
library(grid)



# #### Load Data

# Get the path to the directory containing the script
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the script's directory
setwd(script_dir)

# List all CSV files in the "Clocks Methylation Data" folder
files <- list.files(path = "../Clocks Methylation Data",
                    pattern = "Clock[1-3]_methylation_data_complete.csv|Clock2_3_overlap_methylation_data_complete.csv",
                    # pattern = "Clock[2-3]_methylation_data_complete.csv",
                    full.names = TRUE)

# Read metadata
metadata <- read.csv("../complete_metadata_GSE223748.csv", header = TRUE)

# Define the list of tissues to include
tissues_to_include <- c('Blood', 'Skin', 'Liver', 'Ear', 'Cortex', 'Muscle', 'Heart',
                        'Cerebellum', 'Striatum', 'Kidney', 'Lung', 'Brain', 'Spleen', 'Fibroblast')

# Iterate over each CSV file and create a heatmap for each CSV file
for (file in files) {
  
  # Load Methylation Data
  clock_data <- read.table(file, header = TRUE, sep = ",")
  
  # Filter the metadata to include only the specified tissues
  filtered_metadata <- metadata[metadata$Tissue %in% tissues_to_include, ]
  
  # Subset the clock data to include only the columns corresponding to the filtered metadata
  sample_ids <- filtered_metadata$Sample  # Extracting sample IDs from the metadata
  sample_ids <- as.character(sample_ids)  # Ensuring sample IDs are character type
  
  # Subset clock_data based on sample IDs
  # Including the first column which is CpG IDs
  filtered_clock_data <- clock_data[, c('CpG', sample_ids), drop = FALSE]
  
  # Extract relevant data and convert to matrices
  mat <- as.matrix(filtered_clock_data[, -1])
  rownames(mat) <- clock_data[, 1]  # Ensure CpG IDs are row names
  filtered_metadata <- as.data.frame(filtered_metadata)
  
  # Preprocess Data (without Scaling)
  heat <- mat
  
  # Clustering Analysis (PAM)
  silhouette_scores <- sapply(2:10, function(k) {
    pam_result <- pam(heat, k = k)
    if (is.null(pam_result$silinfo$avg.width)) {
      return(NA)
    } else {
      return(mean(pam_result$silinfo$avg.width))
    }
  })
  
  optimal_k <- which.max(silhouette_scores) + 1  # Adding 1 because index starts from 2
  pamClusters <- pam(heat, k = optimal_k)
  pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)
  pamClusters$clustering <- factor(pamClusters$clustering,
                                   levels = paste0('Cluster ', 1:optimal_k))
  
  # Set Color Scheme of Methylation Values
  myCol <- colorRampPalette(c('#1368aa', 'white', '#ef3c2d'))(100)
  myBreaks <- seq(0, 1, length.out = 100)
  
  # Create Annotations
  ann <- data.frame(
    Tissue = filtered_metadata$Tissue,
    Sex = filtered_metadata$Sex,
    Age.Group = filtered_metadata$Age.Group,
    stringsAsFactors = FALSE
  )
  
  colours <- list(
    Tissue = c(
      "Blood" = "#9e0142",
      "Skin" =  "#d53e4f",
      "Liver" = "#f46d43",
      "Ear" =  "#fdae61",
      "Cortex" = "#fee08b",
      "Muscle" =  "#e6f598",
      "Heart" = "#abdda4",
      "Cerebellum" = "#66c2a5",
      "Striatum" = "#3288bd",
      "Kidney" = "#5e4fa2",
      "Lung" = "#32285f",
      "Brain" = "#542788",
      "Spleen" = "#8073ac",
      "Fibroblast" = "#b2abd2"),
    
    Sex = c("Male" = "#22223b",
            "Female" = "#bc4749",
            "Missing" = "#f8f3eb"),
    
    Age.Group = c("0-30"  =   "#c9f2c7",
                  "30-60" =  "#96be8c" ,
                  "60-90" =  "#629460",
                  "90+" = "#243119",
                  "Missing" = "#f8f3eb")
  )
  
  colAnn <- HeatmapAnnotation(
    df = ann,
    which = 'col',
    na_col = 'white',
    col = colours,
    annotation_name_gp = gpar(fontsize = 11, fontface = 'bold', family = 'serif'),
    annotation_height = 5,
    annotation_width = unit(1, 'cm'), 
    gap = unit(1, 'mm'), 
    annotation_legend_param = list(
      Tissue = list(
        nrow = 6,
        title = 'Tissue',
        title_position = 'topleft',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold', family = 'serif'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain', family = 'serif')),
      Sex = list(
        nrow = 2,
        title = 'Sex',
        title_position = 'topleft',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold', family = 'serif'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain', family = 'serif')),
      Age.Group = list(
        nrow = 4,
        title = 'Age\nGroup',
        title_position = 'topleft',
        legend_direction = 'vertical',
        title_gp = gpar(fontsize = 12, fontface = 'bold', family = 'serif'),
        labels_gp = gpar(fontsize = 10, fontface = 'plain', family = 'serif'))
    )
  )
  
  # Create Saving Instructions
  heatmap_dir <- file.path(script_dir, "Heatmaps")
  if (!file.exists(heatmap_dir)) 
    dir.create(heatmap_dir)
  
  # Group Between Groups (Tissues)
  tissue_factor <- factor(filtered_metadata$Tissue)
  tissue_numeric <- as.numeric(tissue_factor)
  dend1 = cluster_between_groups(heat, tissue_numeric)
  
  # Create Heatmap Object
  hmap <- Heatmap(heat,
                  split = pamClusters$clustering,
                  cluster_row_slices = FALSE,
                  col = colorRamp2(myBreaks, myCol),
                  use_raster = TRUE,  # Set use_raster explicitly
                  heatmap_legend_param = list(
                    title = 'CpGs\nBeta-value',
                    color_bar = 'continuous',
                    legend_direction = 'vertical',
                    legend_width = unit(7.0, 'cm'),
                    legend_height = unit(4.0, 'cm'),
                    title_position = 'topleft',
                    title_gp=gpar(fontsize = 10, fontface = 'bold', family = 'serif'),
                    labels_gp=gpar(fontsize = 10, fontface = 'plain', family = 'serif'),
                    padding = unit(c(2, 2, 2, 2), 'cm')),
                  cluster_rows = TRUE,
                  show_row_dend = TRUE,
                  row_title = 'PAM Clustered CpGs',
                  row_title_side = 'left',
                  row_title_gp = gpar(fontsize = 12,  fontface = 'bold', family = 'serif'),
                  row_title_rot = 90,
                  show_row_names = FALSE,
                  row_names_gp = gpar(fontsize = 10, fontface = 'plain', family = 'serif'),
                  row_names_side = 'left',
                  row_dend_width = unit(25,'mm'),
                  cluster_columns = dend1,
                  show_column_dend = FALSE,
                  column_title = '',
                  column_title_side = 'bottom',
                  column_title_gp = gpar(fontsize = 12, fontface = 'bold', family = 'serif'),
                  column_title_rot = 0,
                  show_column_names = FALSE,
                  column_names_gp = gpar(fontsize = 10, fontface = 'bold', family = 'serif'),
                  column_names_max_height = unit(10, 'cm'),
                  column_dend_height = unit(25,'mm'),
                  clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                  clustering_method_columns = 'ward.D2',
                  top_annotation = colAnn
  )
  
  # Save the heatmap as a high-resolution PDF
  pdf_filename <- file.path(heatmap_dir, sub(".csv$", "_Tissues_Grouped_PAM.pdf", basename(file)))
  
  # Create directories if they don't exist
  dir.create(dirname(pdf_filename), recursive = TRUE, showWarnings = FALSE)
  ht_opt$message = FALSE  # Suppress messages
  
  pdf(pdf_filename, width = 12, height = 8)
  draw(hmap)
  dev.off()
}