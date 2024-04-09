
###############################################################################
#
# Title: Complex Heatmap of methylation data of Clock1 
# Date: 28/03/2024
# Author: Catia Antunes
#
# Description:
# Complex heatmap of methylation data of top 100 CpG sites identified by Clock1
# further grouped by tissue
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



#### Load Data

# Get the path to the directory containing the script
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the script's directory
setwd(script_dir)

# Load data & Sneak a peek
metadata <- read.table("human_metadata_184211.csv", header = TRUE, sep = ",")
metadata[1:5,1:5]
dim(metadata)

clock1_top100_methylation <- read.table("./Clocks Methylation Data/Clock1_methylation_data_top100.csv", header = TRUE, sep = ",")
clock1_top100_methylation[1:5,1:5]
dim(clock1_top100_methylation)

# human_methylation_data <- read.table("GSE184211_datBetaNormalized.csv", header = TRUE, sep = ",")
# human_methylation_data[1:5,1:5]
# dim(human_methylation_data)



#### Read metadata & Load Methylation Data

# Extract relevant data and convert to matrixes
mat <- as.matrix(clock1_top100_methylation[, -1])
# mat <- as.matrix(human_methylation_data[, -1])
metadata <- as.data.frame(metadata)

# print the unique names of the tissues
unique(metadata$Tissue)

# Check if number of rows (samples) in the metadata is the same as the cols (samples) in the methylation data
all(rownames(metadata) == colnames(mat)) # prints false because of the index col

# print number of rows and columns
dim(metadata)
dim(mat)
mat[1:5,1:5]



#### Preprocess Data (Scaling)

# Scale the data
heat <- t(scale(t(mat)))



#### Clusetirng Analysis (PAM)

# Perform partitioning around medoids (PAM) to identify clusters in the data
pamClusters <- cluster::pam(heat, k = 4) # pre-select k = 4 centers
pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)

# fix order of the clusters to have 1 to 4, top to bottom
pamClusters$clustering <- factor(pamClusters$clustering,
                                 levels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'))

# print the number of CpGs in each cluster
table(pamClusters$clustering)



####  Set Color Scheme of Methylation values

# Set color scheme
myCol <- colorRampPalette(c('#1368aa', 'white', '#ef3c2d'))(100)
myBreaks <- seq(0, 1, length.out = 100)  # Adjust the breaks to range from 0 to



#### Create annotations

# create a df with the annotations (Tissues, Sex, Age Group)
ann <- data.frame(
  Tissue = metadata$Tissue,
  Sex = metadata$Sex,
  Age.Group = metadata$Age.Group,
  stringsAsFactors = FALSE
)

# Set colors for the annotations
colours <- list(
  Tissue = c(
    "Adipose" = "#9e0142",
    "Lung" =  "#d53e4f",
    "Blood" = "#f46d43",
    "BoneMarrow" =  "#fdae61",
    "Heart" = "#fee08b",
    "Kidney" =  "#e6f598",
    "Liver" = "#abdda4",
    "Spleen" = "#66c2a5",
    "LymphNode" = "#3288bd",
    "Muscle" = "#5e4fa2",
    "Pituitary" = "#32285f"),

  Sex = c("Male" = "#22223b",
          "Female" = "#bc4749"),

  Age.Group = c("0-30"  =   "#c9f2c7",
                "30-60" =  "#96be8c" ,
                "60-90" =  "#629460",
                "90+" = "#243119")
)

# create the ComplexHeatmap annotation object
colAnn <- HeatmapAnnotation( # annotation object for the columns
  df = ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = colours,
  annotation_name_gp = gpar(fontsize = 11, fontface = 'bold', family = 'serif'),
  annotation_height = 5, # height of the annotation in the heatmap
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Tissue = list(
      nrow = 6, # number of rows across which the legend will be arranged
      title = 'Tissue', # title of the legend
      title_position = 'topleft', # position of the title, options are
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



#### Create saving instructions

# Create a directory for saving heatmaps if it doesn't exist
heatmap_dir <- file.path("Heatmaps")
if (!file.exists(heatmap_dir))
  dir.create(heatmap_dir)

# create the filename for the heatmap
filename <- file.path(heatmap_dir, sub(".csv$", "_Heatmap.png", basename(file))) # Dynamically set the filename and save in the heatmap directory

# Instructions for the png file
png(filename, width=10,height=7,units="in",res=1200)



#### Groupe between groups (Tissues)

# Generate a factor vector based on metadata$Tissue
tissue_factor <- factor(metadata$Tissue)

# Assign numeric values to each level
tissue_numeric <- as.numeric(tissue_factor)

#Create a grouping instance
dend1 = cluster_between_groups(heat, tissue_numeric)



#### Create heatmap object

hmap <- Heatmap(heat,

                # split the CpGs / rows according to the PAM clusters
                split = pamClusters$clustering,
                cluster_row_slices = FALSE,


                # Set colour scheme of the range of values
                col = colorRamp2(myBreaks, myCol),

                # parameters for the colour-bar that represents gradient of methylation?
                heatmap_legend_param = list(
                  title = 'CpGs\nBeta-value',
                  color_bar = 'continuous',
                  legend_direction = 'vertical', # set the direction of the legend
                  legend_width = unit(7.0, 'cm'), # set the width of the legend
                  legend_height = unit(4.0, 'cm'), # set the height of the legend
                  title_position = 'topleft', # set the position of the title
                  title_gp=gpar(fontsize = 10, fontface = 'bold', family = 'serif'),
                  labels_gp=gpar(fontsize = 10, fontface = 'plain', family = 'serif'),
                  padding = unit(c(2, 2, 2, 2), 'cm')), # (top, right, bottom, left)


                # row (CpG) parameters
                cluster_rows = TRUE, # No clustering of rows
                show_row_dend = TRUE, # No row dendrogram

                row_title = 'PAM Clustered CpGs',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 12,  fontface = 'bold', family = 'serif'),
                row_title_rot = 90,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 10, fontface = 'plain', family = 'serif'),
                row_names_side = 'left',
                row_dend_width = unit(25,'mm'),


                # column (sample) parameters
                cluster_columns = dend1, # Cluster columns (cluster samples)
                show_column_dend = FALSE, # No column dendrogram
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 12, fontface = 'bold', family = 'serif'),
                column_title_rot = 0,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 10, fontface = 'bold', family = 'serif'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'),

                # # cluster methods for rows and columns
                clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                clustering_method_columns = 'ward.D2',
                # clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                # clustering_method_rows = 'ward.D2',

                # specify top and bottom annotations
                top_annotation = colAnn,
                # bottom_annotation = boxplotCol

)



#### Create heatmap object
draw(
  hmap,
  heatmap_legend_side = 'right',
  annotation_legend_side = 'left',
  row_sub_title_side = 'left',
  legend_gap = unit(0.5, "cm"),
  padding = unit(c(2, 1, 1, 1), "cm"),
  column_title="Methylation Profile of Top 100 CpG's Identified by Clock 1",
  column_title_gp=grid::gpar(fontsize=13, fontface="bold", family="serif"),
)



### Save the heatmap
dev.off()