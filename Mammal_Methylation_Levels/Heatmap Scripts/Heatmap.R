# Heatmap of methylation data of Clock1 

### Load libraries
library(gplots)
library(RColorBrewer)

### Load data

# Metada of 661 human samples (contains Sample, GEO_ID, Tissue, Age, Sex, Age Group)
metadata <- read.table("human_metadata_184211.csv", header = TRUE, sep = ",")

# Full human methylation data (contains beta vals for aprox 37,000 CpGs of 661 human samples)
#methylation <- read.table("GSE184211_datbetaNormalized.csv", header = TRUE, sep = ",") 

# Methylation data of the to 100 CpG sites identified by Clock 1 (Contains beta vals for 100 CpGs of 661 human samples)
clock1_top100_methylation <- read.table("Clock1_methylation_data.csv", header = True, sep = ",") # Top 100 CpG sites identified by clock 1


### Extract relevant data

# Extract the ....
#methylation_data <- as.matrix(methylation[, -1])  # Assuming first column is row names
#rownames(methylation_data) <- methylation[, 1]    # Assign row names
clock1_top100_methylation_data <- as.matrix(clock1_top100_methylation[, -1]) # 
rownames(clock1_top100_methylation_data) <- clock1_top100_methylation[, 1]
metadata_data <- metadata        

# Create heatmap
heatmap_data <- heatmap.2(methylation_data,
                          scale = "row",         # Scale rows
                          Colv = NA,             # Don't cluster columns
                          dendrogram = "row",    # Add row dendrogram
                          trace = "none",        # No trace lines
                          margins = c(12, 9),    # Adjust margins
                          col = colorRampPalette(c("blue", "white", "red"))(100),  # Color scheme
                          key = TRUE,            # Show color key
                          keysize = 1.5,         # Adjust size of color key
                          symkey = FALSE,        # Don't show symmetric key
                          density.info = "none", # Don't show density plot
                          cexRow = 0.5,          # Adjust row label size
                          cexCol = 0.7,          # Adjust column label size
                          main = "Methylation Heatmap")  # Title

# Add metadata as color side bar
heatmap_data$colSideColors <- as.vector(metadata_data$Color)

# Save heatmap and name it according to the input data
#save_heatmap(heatmap_data, "heatmap_methylation_data.png")


png("heatmap_full_human_dataset.png", width = 10, height = 8, units = "in", res = 300)
heatmap_data
dev.off()

