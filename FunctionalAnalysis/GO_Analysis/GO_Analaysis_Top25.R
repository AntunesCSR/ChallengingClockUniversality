#### LOAD LIBRARIES ####
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

#### LOAD THE DATA ####

# List of file names
file_names <- c("top_25_cpgs_clock_1.csv", 
                "top_25_cpgs_clock_2.csv", 
                "top_25_cpgs_clock_3.csv", 
                "top_25_cpgs_overlap_1_2_3.csv", 
                "top_25_cpgs_overlap_2_3.csv")

# List to store data frames
cpgs_data_list <- list()

# Loop through each file and read into a data frame
for (file_name in file_names) {
  # Create the file path
  csv_file_path <- file.path("../../FunctionalAnalysis", file_name)
  
  # Read the CSV file into a data frame and add it to the list
  cpgs_data_list[[file_name]] <- read.csv(csv_file_path)
}

# Access the data frames from the list
clock_1_df <- cpgs_data_list[["top_25_cpgs_clock_1.csv"]]
clock_2_df <- cpgs_data_list[["top_25_cpgs_clock_2.csv"]]
clock_3_df <- cpgs_data_list[["top_25_cpgs_clock_3.csv"]]
overlap_1_2_3_df <- cpgs_data_list[["top_25_cpgs_overlap_1_2_3.csv"]]
overlap_2_3_df <- cpgs_data_list[["top_25_cpgs_overlap_2_3.csv"]]



#### GENE ONTOLOGY ENRICHMENT ANALYSIS ####


# Create a list of data frames for the gene sets
gene_set_data_list <- list(
  "clock_1" = clock_1_df,
  "clock_2" = clock_2_df,
  "clock_3" = clock_3_df,
  "overlap_1_2_3" = overlap_1_2_3_df,
  "overlap_2_3" = overlap_2_3_df
)

# Create a list to store the results of the gene set enrichment analysis
GOA_results_list <- list()

# Loop through each data frame and perform gene set enrichment analysis
for (gene_set_name in names(gene_set_data_list)) {
  # Get the data frame for the gene set
  gene_set_df <- gene_set_data_list[[gene_set_name]]
  
  # Extract the gene symbols from the data frame
  gene_symbols <- gene_set_df$Gene # symbol
  gene_entrez <- gene_set_df$ENTREZID
  
  # Perform gene set enrichment analysis using the DOSE package 
  GOA_results <- enrichGO(gene = gene_symbols, 
                              OrgDb = org.Hs.eg.db,
                              keyType = "SYMBOL",
                              ont = "BP", # GO is structured into three main ontologies: Biological Process (BP), Molecular Function (MF), and Cellular Component (CC).
                              pAdjustMethod = "BH", # Benjamini-Hochberg (BH) coerrection
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05) # use 0.2 or 0.05 depending on the number of genes
  
  
  # Add the results to the list
  GOA_results_list[[gene_set_name]] <- GOA_results
}

# Access the results of the gene set enrichment analysis
clock_1_GOA_results <- GOA_results_list[["clock_1"]]
clock_2_GOA_results <- GOA_results_list[["clock_2"]]
clock_3_GOA_results <- GOA_results_list[["clock_3"]]
overlap_1_2_3_GOA_results <- GOA_results_list[["overlap_1_2_3"]]
overlap_2_3_GOA_results <- GOA_results_list[["overlap_2_3"]]

# Print the results
print(clock_1_GOA_results) # 0 enriched terms found
print(clock_2_GOA_results)
print(clock_3_GOA_results)
print(overlap_1_2_3_GOA_results)
print(overlap_2_3_GOA_results)


# Save the results to CSV files
# write.csv(clock_1_GOA_results, file = "FunctionalAnalysis/GO_Analysis/clock_1_GOA_25_results.csv") # 0 enriched terms found
write.csv(clock_2_GOA_results, file = "FunctionalAnalysis/GO_Analysis/clock_2_GOA_25_results.csv")
write.csv(clock_3_GOA_results, file = "FunctionalAnalysis/GO_Analysis/clock_3_GOA_25_results.csv")
write.csv(overlap_1_2_3_GOA_results, file = "FunctionalAnalysis/GO_Analysis/overlap_1_2_3_GOA_25_results.csv")
write.csv(overlap_2_3_GOA_results, file = "FunctionalAnalysis/GO_Analysis/overlap_2_3_GOA_25_results.csv")


# Visualize the results
# dotplot(clock_1_GOA_results, showCategory = 10, title = "Top 10 GO terms for clock 1 (25 Genes)") # 0 enriched terms found
dotplot(clock_2_GOA_results, showCategory = 10, title = "Top 10 GO terms for clock 2 (25 Genes)")
dotplot(clock_3_GOA_results, showCategory = 10, title = "Top 10 GO terms for clock 3 (25 Genes)")
dotplot(overlap_1_2_3_GOA_results, showCategory = 10, title = "Top 10 GO terms for overlap of clocks 1, 2 and 3 (25 Genes)")
dotplot(overlap_2_3_GOA_results, showCategory = 10, title = "Top 10 GO terms for overlap 2 and 3 (25 Genes)")

# Need to change wd to GO annot

# Export the plots as PNG files
# ggsave("clock_1_GOA_25_results.png", plot = dotplot(clock_1_GOA_results, showCategory = 10, title = "Top 10 GO terms for clock 1 (25 Genes)")) # 0 enriched terms found
ggsave("clock_2_GOA_25_results.png", plot = dotplot(clock_2_GOA_results, showCategory = 10, title = "Top 10 GO terms for clock 2 (25 Genes)"))
ggsave("clock_3_GOA_25_results.png", plot = dotplot(clock_3_GOA_results, showCategory = 10, title = "Top 10 GO terms for clock 3 (25 Genes)"))
ggsave("overlap_1_2_3_GOA_25_results.png", plot = dotplot(overlap_1_2_3_GOA_results, showCategory = 10, title = "Top 10 GO terms for overlap of clocks 1, 2 and 3 (25 Genes)"))
ggsave("overlap_2_3_GOA_25_results.png", plot = dotplot(overlap_2_3_GOA_results, showCategory = 10, title = "Top 10 GO terms for overlap 2 and 3 (25 Genes)"))













