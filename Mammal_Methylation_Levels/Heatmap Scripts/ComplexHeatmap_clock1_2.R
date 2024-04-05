# Load required libraries
library(ComplexHeatmap)

# Load data
metadata <- read.table("human_metadata_184211.csv", header = TRUE, sep = ",")
clock1_top100_methylation <- read.table("Clock1_methylation_data.csv", header = TRUE, sep = ",")

# Extract relevant data
clock1_top100_methylation_data <- as.matrix(clock1_top100_methylation[, -1])
rownames(clock1_top100_methylation_data) <- clock1_top100_methylation[, 1]

# Aggregate samples by tissue
tissue_data <- split(metadata$Sample, metadata$Tissue)
tissue_names <- names(tissue_data)

# Create top annotation aggregated by tissue
top_annotation_list <- lapply(tissue_names, function(tissue) {
  samples <- tissue_data[[tissue]]
  anno_block(gp = gpar(fill = rep(1:5, length.out = length(samples))),
             labels = rep(LETTERS[1:5], length.out = length(samples)))
})

# Create heatmap with aggregated top annotation
Heatmap(clock1_top100_methylation_data, name = "clock1_top100_methylation_data",
        top_annotation = HeatmapAnnotation(tissue = anno_blocks(top_annotation_list, align_to = tissue_samples),
        column_title = "Samples"))
