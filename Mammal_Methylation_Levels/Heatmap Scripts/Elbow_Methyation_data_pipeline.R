# Set the working directory to "Mammal_Methylation_Levels"
# setwd("Mammal_Methylation_Levels")

# List all CSV files in the "Clocks Methylation Data" folder
files <- list.files(path = "Clocks Methylation Data", pattern = "*.csv", full.names = TRUE)

# Vector to store WCSS values
wcss_values <- numeric()

# Maximum number of clusters to consider
max_clusters <- 10

# Iterate over each CSV file
for (file in files) {
  # Read metadata
  metadata <- read.table("human_metadata_184211.csv", header = TRUE, sep = ",")
  
  # Read methylation data
  clock_data <- read.table(file, header = TRUE, sep = ",")
  
  # Extract relevant data and convert to matrices
  mat <- as.matrix(clock_data[, -1])
  metadata <- as.data.frame(metadata)
  
  # Scale the data
  heat <- t(scale(t(mat)))
  
  # Calculate WCSS for different number of clusters
  wcss <- numeric(max_clusters)
  for (k in 1:max_clusters) {
    pamClusters <- cluster::pam(heat, k = k)
    wcss[k] <- sum(pamClusters$medoids.dist)
  }
  
  # Append WCSS values to the vector
  wcss_values <- cbind(wcss_values, wcss)
}

# Plot the elbow curve
plot(1:max_clusters, colMeans(wcss_values), type = "b", pch = 19, xlab = "Number of Clusters", ylab = "Within-Cluster Sum of Squares (WCSS)")
