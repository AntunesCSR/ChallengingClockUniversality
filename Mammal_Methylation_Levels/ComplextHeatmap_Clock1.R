
###############################################################################
#
# Title: Complex Heatmap of methylation data of Clock1 
# Date: 28/03/2024
# Author: Catia Antunes
#
# Description:
# Complex heatmap of methylation data of top 100 CpG sites identified by Clock1
#
###############################################################################



library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)
        
# https://github.com/kevinblighe/E-MTAB-6141?tab=readme-ov-file

# Load data & Sneak a peek
metadata <- read.table("human_metadata_184211.csv", header = TRUE, sep = ",")
metadata[1:5,1:5]
dim(metadata)

clock1_top100_methylation <- read.table("Clock1_methylation_data.csv", header = TRUE, sep = ",")
clock1_top100_methylation[1:5,1:5]
dim(clock1_top100_methylation)


# Extract relevant data and convert to matrixes
mat <- as.matrix(clock1_top100_methylation[, -1])
metadata <- as.data.frame(metadata)

# print the unique names of the tissues
unique(metadata$Tissue)


# Check if number of rows (samples) in the metadata is the same as the cols (samples) in the methylation data
all(rownames(metadata) == colnames(mat))



# Scale the data 
heat <- t(scale(t(mat)))



# Set color scheme
myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)



# Create annotation for the heatmap

# Tissues
tissues <- metadata$Tissue
pick.col <- brewer.pal(9, 'Greens')
col.tissues <- colorRampPalette(pick.col)(length(unique(tissues)))
unique(col.tissues)

# Gender
sex <- metadata$Sex
pick.col <- brewer.pal(9, 'Purples')
col.sex <- colorRampPalette(pick.col)(length(unique(sex)))
unique(col.sex)

# Age group
age.group <- metadata$Age.Group
pick.col <- brewer.pal(9, 'Oranges')
col.age.group <- colorRampPalette(pick.col)(length(unique(age.group)))
unique(col.age.group)


# Create annotation for the heatmap

ann <- data.frame(
  Tissue = metadata$Tissue,
  Sex = metadata$Sex,
  Age.Group = metadata$Age.Group,
  stringsAsFactors = FALSE
)

colours <- list(
  Tissue = c(
    "Adipose" = "#FFC300",
    "Lung" =  "#900C3F",
    "Blood" = "#C70039",
    "BoneMarrow" =  "#FF5733",
    "Heart" = "#581845",
    "Kidney" =  "#DAF7A6",
    "Liver" = "#36454F",
    "Spleen" = "#A7C7E7",
    "LymphNode" = "#6F4E37",
    "Muscle" = "#8A9A5B",
    "Pituitary" = "#FF5795"),
  
  Sex = c("Male" = "#6082B6",
          "Female" = "#800020"),
            
  Age.Group = c("30-60" =  "#f9e79f" ,
                "60-90" =  "#abebc6",
                "0-30"  =   "#aed6f1",
                "90+" = "#d2b4de")
)

# create the ComplexHeatmap annotation object
colAnn <- HeatmapAnnotation(
  df = ann,
  which = 'col', # 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
  col = colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Tissue = list(
      nrow = 10, # number of rows across which the legend will be arranged
      title = 'Tissue', # title of the legend
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    Sex = list(
      nrow = 2,
      title = 'Sex',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    Age.Group = list(
      nrow = 4,
      title = 'Age Group',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold'))
  )
)


# Perform partitioning around medoids (PAM) to identify clusters in the data

pamClusters <- cluster::pam(heat, k = 4) # pre-select k = 4 centers
pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)


# fix order of the clusters to have 1 to 4, top to bottom
pamClusters$clustering <- factor(pamClusters$clustering,
                                 levels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'))

# Create heatmap object

hmap <- Heatmap(heat,
                
                # split the CpGs / rows according to the PAM clusters
                split = pamClusters$clustering,
                cluster_row_slices = FALSE, 
                
              
                name = 'CpG\nBeta-\nvalue',
                
                col = colorRamp2(myBreaks, myCol), # colorRamp2(myBreaks, myCol)
                
                # parameters for the colour-bar that represents gradient of expression/methylation?
                heatmap_legend_param = list(
                  color_bar = 'continuous',
                  legend_direction = 'vertical',
                  legend_width = unit(8, 'cm'),
                  legend_height = unit(5.0, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 12, fontface = 'bold'),
                  labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                
                # row (CpG) parameters
                cluster_rows = TRUE,
                show_row_dend = TRUE,
                
                #row_title = 'Statistically significant CpGs',
                row_title_side = 'left',
                row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                row_title_rot = 90,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                row_names_side = 'left',
                row_dend_width = unit(25,'mm'),
                
                # column (sample) parameters
                cluster_columns = TRUE,
                show_column_dend = TRUE,
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                column_title_rot = 0,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'),
                
                # cluster methods for rows and columns
                clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                clustering_method_columns = 'ward.D2',
                clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                clustering_method_rows = 'ward.D2',
                
                # specify top and bottom annotations
                top_annotation = colAnn,
                # bottom_annotation = boxplotCol
                )

# Draw the heatmap                
draw(hmap,
     heatmap_legend_side = 'left',
     annotation_legend_side = 'right',
     row_sub_title_side = 'left')



# Save the heatmap
save_heatmap(hmap, "heatmap_methylation_data_Clock1_top100.png")

# save the figure
dev.copy(png, "heatmap_methylation_data_Clock1_top100.png", width = 10, height = 8, units = "in", res = 300)






















