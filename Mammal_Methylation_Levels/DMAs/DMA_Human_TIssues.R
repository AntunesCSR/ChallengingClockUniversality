# Load necessary libraries
library(limma) # for differential methylation analysis
library(ComplexHeatmap)
library(colorRamp2)
library(cluster)
library(circlize)


### Load the methylation data and metadata

# Load methylation data
methylation_data <- read.csv("../GSE184211_datBetaNormalized.csv", header = TRUE, row.names = 1) # Assuming the first column contains row names
head(methylation_data)
dim(methylation_data)

# Load metadata
metadata <- read.csv("../human_metadata_184211.csv", header = TRUE)
head(metadata)
dim(metadata)


# Extract tissue types from metadata
tissue_types <- metadata$Tissue
head(tissue_types)
unique(tissue_types)
length(unique(tissue_types)) # Should be 11


# Design matrix
design <- model.matrix(~0 + factor(tissue_types)) # Create a design matrix with tissue types as factors
colnames(design) <- levels(factor(tissue_types))
head(design)


# Fit the linear model
fit <- lmFit(methylation_data, design)

# Perform empirical Bayes moderation
fit <- eBayes(fit)

# Get the differential methylation results
results <- topTable(fit, coef = colnames(design), adjust.method = "BH", sort.by = "none")
results

# Extract the significant CpG sites
significant_cpgs <- results[results$adj.P.Val < 0.05, ]
significant_cpgs
dim(significant_cpgs)




# 
# ###################### DMA with ANOVA
# 
# # Define the model formula
# formula <- as.formula(paste("methylation_data ~ tissue_types"))
# 
# # Fit ANOVA model
# fit <- lmFit(methylation_data, formula)
# 
# # Perform empirical Bayes moderation of standard errors
# fit <- eBayes(fit)
# 
# # Extract differential methylation results
# results <- topTable(fit, coef = "tissue_types", number = Inf)
# 
# # Adjust p-values for multiple testing
# results$p.adjust <- p.adjust(results$P.Value, method = "fdr")
# 
# # Filter significant results
# significant_results <- subset(results, p.adjust < 0.05)
# 
# # Output significant results
# write.csv(significant_results, file = "significant_results_anova.csv")
# 
# 
# 
# 
# 
# 
# ######################  DMA with Kruskall-Wallis test
# 
# # Kruskal-Wallis test for each CpG site
# kw_results <- apply(methylation_data, 1, function(x) kruskal.test(x ~ tissue_types)$p.value)
# 
# # Adjust p-values for multiple testing using the false discovery rate (FDR) method
# kw_results_adj <- p.adjust(kw_results, method = "fdr")
# 
# # Filter significant results
# significant_kw_results <- which(kw_results_adj < 0.05)
# 
# # Output significant results
# write.csv(significant_kw_results, file = "significant_kw_results.csv")
# 
# # Print some of the significant results
# head(significant_kw_results)
# 
# # Count the number of significant CpG sites
# num_significant_kw <- length(significant_kw_results)
# 
# # Print the number of significant CpG sites
# print(num_significant_kw)



############# Heatmap of Diffrerentially Methylated CpG Sites (TISSUE GROUPING)

#### Read metadata & Load Methylation Data

# Extract relevant data and convert to matrixes
mat <- as.matrix(methylation_data[rownames(significant_cpgs), ])
metadata <- as.data.frame(metadata)


#### Preprocess Data (Scaling)

# Scale the data
heat <- t(scale(t(mat)))
dim(heat)


#### PAM Clustering Analysis to determine optimal number of clusters

## Silhouette Plot to determine optimal number of clusters

silhouette_scores <- sapply(2:9, function(k) { # modify the range of k as needed (k must be between 1 and n-1, here n =10 [number of datapoints])
  pam_result <- pam(heat, k = k)
  if (is.null(pam_result$silinfo$avg.width)) {
    return(NA)
  } else {
    return(mean(pam_result$silinfo$avg.width))
  }
})

# Plot silhouette scores
plot(2:9, silhouette_scores, type = "b", pch = 19, xlab = "Number of Clusters (k)", 
     ylab = "Average Silhouette Width", main = "Silhouette Analysis")

# Choose the value of k that maximizes the average silhouette width
optimal_k <- which.max(silhouette_scores) + 1  # Adding 1 because index starts from 2

# Perform clustering analysis with the optimal number of clusters
pamClusters <- pam(heat, k = optimal_k)
# Continue with your analysis using the optimal number of cluster


## Permorm PAM Clustering according to optimal number of clusters

# Perform partitioning around medoids (PAM) to identify clusters in the data
pamClusters <- cluster::pam(heat, optimal_k)  # Pre-select k = 4 centers
pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)

# Rename clusters to 'Cluster 1', 'Cluster 2', etc.
pamClusters$clustering <- factor(pamClusters$clustering,
                                 levels = paste0('Cluster ', 1:optimal_k))

# Print the number of CpGs in each cluster
table(pamClusters$clustering)


####  Set Color Scheme of Methylation values

# Set color scheme
myCol <- colorRampPalette(c('#1368aa', 'white', '#ef3c2d'))(100)
myBreaks <- seq(0, 1, length.out = 100)  # Adjust the breaks to range from 0 to



#### Create annotations

# create a df with the annotations (Tissues, Sex, Age Group)
ann <- data.frame(
  Tissue = metadata$Tissue,
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
    "Pituitary" = "#32285f")
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
      labels_gp = gpar(fontsize = 10, fontface = 'plain', family = 'serif'))
  )
)



#### Create saving instructions

# Create a directory for saving heatmaps if it doesn't exist
heatmap_dir <- file.path("Heatmaps")
if (!file.exists(heatmap_dir))
  dir.create(heatmap_dir)

# create the filename for the heatmap
filename <- file.path(heatmap_dir, "GSE184211_DMC_Heatmap.png")


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
                
                # # split the CpGs / rows according to the PAM clusters
                # split the CpGs / rows according to the PAM clusters
                split = pamClusters$clustering,
                cluster_row_slices = TRUE,
                
                
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
                show_row_names = TRUE,
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
  column_title="Methylation Profile of Significant CpG's from DMAs of Human Samples (GSE184211)",
  column_title_gp=grid::gpar(fontsize=13, fontface="bold", family="serif"),
)



### Save the heatmap
dev.off()











################## Heatmap of Diffrerentially Methylated CpG Sites (NO GROUPING)


# Extract relevant data and convert to matrices
mat <- as.matrix(methylation_data[rownames(significant_cpgs), ])

metadata <- as.data.frame(metadata)

# print the unique names of the tissues
unique(metadata$Tissue)


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


#### Clusetirng Analysis (PAM)

# # Perform partitioning around medoids (PAM) to identify clusters in the data
# pamClusters <- cluster::pam(heat, k = 4) # pre-select k = 4 centers
# pamClusters$clustering <- paste0('Cluster ', pamClusters$clustering)
# 
# # fix order of the clusters to have 1 to 4, top to bottom
# pamClusters$clustering <- factor(pamClusters$clustering,
#                                  levels = c('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'))
# 
# # print the number of CpGs in each cluster
# table(pamClusters$clustering)

####  Set Color Scheme of Methylation values

# Set color scheme
myCol <- colorRampPalette(c('#1368aa', 'white', '#ef3c2d'))(100)
myBreaks <- seq(0, 1, length.out = 100)  # Adjust the breaks to range from 0 to



#### Create annotations

# create a df with the annotations (Tissues, Sex, Age Group)
ann <- data.frame(
  Tissue = metadata$Tissue,
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
    "Pituitary" = "#32285f")
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
      labels_gp = gpar(fontsize = 10, fontface = 'plain', family = 'serif'))
  )
)



#### Create saving instructions

# # Create a directory for saving heatmaps if it doesn't exist
# heatmap_dir <- file.path("Heatmaps")
# if (!file.exists(heatmap_dir)) 
#   dir.create(heatmap_dir)
# 
# # create the filename for the heatmap
# filename <- file.path(heatmap_dir, sub(".csv$", "_Heatmap.png", basename(file))) # Dynamically set the filename and save in the heatmap directory
# 
# # Instructions for the png file
# png(filename, width=10,height=7,units="in",res=1200)

png("Test.png" , width=10,height=7,units="in",res=1200)


#### Create heatmap object

hmap <- Heatmap(heat,
                
                # split the CpGs / rows according to the PAM clusters
                # split = pamClusters$clustering,
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
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 10, fontface = 'plain', family = 'serif'),
                row_names_side = 'left',
                row_dend_width = unit(25,'mm'),
                
                
                # column (sample) parameters
                cluster_columns = TRUE, # Cluster columns (cluster samples)
                show_column_dend = FALSE, # No column dendrogram
                column_title = '',
                column_title_side = 'bottom',
                column_title_gp = gpar(fontsize = 12, fontface = 'bold', family = 'serif'),
                column_title_rot = 0,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 10, fontface = 'bold', family = 'serif'), 
                column_names_max_height = unit(10, 'cm'),
                column_dend_height = unit(25,'mm'), 
                
                
                
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
  column_title="Methylation Profile of Significant CpG's from Differential Methylation Analysis",
  column_title_gp=grid::gpar(fontsize=13, fontface="bold", family="serif"),
)



### Save the heatmap
dev.off()
