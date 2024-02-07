#### LOAD THE DATA ####

# List of file names
file_names <- c("top_25_cpgs_clock_1.csv", "top_25_cpgs_clock_2.csv", "top_25_cpgs_clock_3.csv", "top_25_cpgs_overlap_1_2_3.csv", "top_25_cpgs_overlap_2_3.csv")

# List to store data frames
cpgs_data_list <- list()

# Loop through each file and read into a data frame
for (file_name in file_names) {
  # Create the file path
  csv_file_path <- file.path("FunctionalAnalysis", file_name)
  
  # Read the CSV file into a data frame and add it to the list
  cpgs_data_list[[file_name]] <- read.csv(csv_file_path)
}

# Access the data frames from the list
clock_1_df <- cpgs_data_list[["top_25_cpgs_clock_1.csv"]]
clock_2_df <- cpgs_data_list[["top_25_cpgs_clock_2.csv"]]
clock_3_df <- cpgs_data_list[["top_25_cpgs_clock_3.csv"]]
overlap_1_2_3_df <- cpgs_data_list[["top_25_cpgs_overlap_1_2_3.csv"]]
overlap_2_3_df <- cpgs_data_list[["top_25_cpgs_overlap_2_3.csv"]]



#### GENE SET ENRICHMENT ANALYSIS ####

# Load the required libraries
library(clusterProfiler)

# Create a list of data frames for the gene sets
gene_set_data_list <- list(
  "clock_1" = clock_1_df,
  "clock_2" = clock_2_df,
  "clock_3" = clock_3_df,
  "overlap_1_2_3" = overlap_1_2_3_df,
  "overlap_2_3" = overlap_2_3_df
)

# Create a list to store the results of the gene set enrichment analysis
gsea_results_list <- list()

# Loop through each data frame and perform gene set enrichment analysis
for (gene_set_name in names(gene_set_data_list)) {
  # Get the data frame for the gene set
  gene_set_df <- gene_set_data_list[[gene_set_name]]
  
  # Extract the gene symbols from the data frame
  gene_symbols <- gene_set_df$Gene # symbol
  gene_entrez <- gene_set_df$ENTREZID
  
  # Perform gene set enrichment analysis using the clusterProfiler package
  gsea_results_KEGG <- enrichKEGG(gene = gene_symbols, organism = "hsa", keyType = "symbol", pvalueCutoff = 0.05)
  gsea_results_DAVID <- enrichDAVID(gene = gene_entrez, pvalueCutoff = 0.05)
  

  
  
  # Add the results to the list
  gsea_results_list[[gene_set_name]] <- gsea_results
}

# Access the results of the gene set enrichment analysis
clock_1_gsea_results <- gsea_results_list[["clock_1"]]
clock_2_gsea_results <- gsea_results_list[["clock_2"]]
clock_3_gsea_results <- gsea_results_list[["clock_3"]]
overlap_1_2_3_gsea_results <- gsea_results_list[["overlap_1_2_3"]]
overlap_2_3_gsea_results <- gsea_results_list[["overlap_2_3"]]

# Print the results
print(clock_1_gsea_results)
print(clock_2_gsea_results)
print(clock_3_gsea_results)
print(overlap_1_2_3_gsea_results)
print(overlap_2_3_gsea_results)

# Save the results to CSV files
write.csv(clock_1_gsea_results, file = "FunctionalAnalysis/clock_1_gsea_GSEA_results.csv")
write.csv(clock_2_gsea_results, file = "FunctionalAnalysis/clock_2_gsea_GSEA_results.csv")
write.csv(clock_3_gsea_results, file = "FunctionalAnalysis/clock_3_gsea_GSEA_results.csv")
write.csv(overlap_1_2_3_gsea_results, file = "FunctionalAnalysis/overlap_1_2_3_gsea_GSEA_results.csv")
write.csv(overlap_2_3_gsea_results, file = "FunctionalAnalysis/overlap_2_3_gsea_GSEA_results.csv")

#### NOTE TO SELF.... error with KEGG ### 


