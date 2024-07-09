### Load necessary packages

library(clusterProfiler)
library(ggplot2)
library(here)



### Define Functions

# Function to perform GO enrichment analysis and save results
perform_GO_analysis <- function(gene_set_df, gene_set_name, ont_type, output_folder) {
  gene_symbols <- gene_set_df$Gene # assuming the column is named 'Gene'
  
  # Perform GO enrichment analysis
  GOA_results <- enrichGO(
    gene = gene_symbols,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = ont_type,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  # Set the working directory to the output folder
  setwd(output_folder)
  
  # Save results to CSV file
  write.csv(GOA_results, file = paste0(gene_set_name, "_GOA_", ont_type, "_results.csv"))
  
  # Visualize with dotplot (if GOA_results is not empty)
  if (nrow(GOA_results) > 0) {
    dotplot(GOA_results, showCategory = 10, title = paste("Top 10 GO", ont_type, "terms for", gene_set_name))
    # Save dotplot visualization to PNG file
    ggsave(paste0(gene_set_name, "_GOA_", ont_type, "_results_dotplot.png"), 
           plot = last_plot())
  }
  
  # Visualize with barplot (if GOA_results is not empty)
  if (nrow(GOA_results) > 0) {
    barplot(GOA_results, showCategory = 10, title = paste("Top 10 GO", ont_type, "terms for", gene_set_name))
    # Save barplot visualization to PNG file
    ggsave(paste0(gene_set_name, "_GOA_", ont_type, "_results_barplot.png"), 
           plot = last_plot())
  }
  
  # Return the GO enrichment results
  return(GOA_results)
}



### Access the data frames from the list
clock_1_df <- cpgs_data_list[["top_100_cpgs_clock_1.csv"]]
clock_2_df <- cpgs_data_list[["top_100_cpgs_clock_2.csv"]]
clock_3_df <- cpgs_data_list[["top_100_cpgs_clock_3.csv"]]
overlap_1_2_3_df <- cpgs_data_list[["top_100_cpgs_overlap_1_2_3.csv"]]
overlap_2_3_df <- cpgs_data_list[["top_100_cpgs_overlap_2_3.csv"]]



### Create a list of data frames for the gene sets
gene_set_data_list <- list(
  "Clock 1" = clock_1_df,
  "Clock 2" = clock_2_df,
  "Clock 3" = clock_3_df,
  "Overlap of Clocks 1, 2 & 3" = overlap_1_2_3_df,
  "overlap of Clocks 2 & 3" = overlap_2_3_df
)



### Set the output folder
output_folder <- here("FunctionalAnalysis/Functional_Enrichment/GO_Analysis")



### Create a list to store the results of the gene set enrichment analysis
GOA_results_list <- list()



### Loop through each data frame and perform gene set enrichment analysis for BP, MF, and CC
for (gene_set_name in names(gene_set_data_list)) {
  gene_set_df <- gene_set_data_list[[gene_set_name]]
  
  # Perform GO enrichment analysis for Biological Process (BP)
  BP_results <- perform_GO_analysis(gene_set_df, gene_set_name, "BP", output_folder)
  GOA_results_list[[paste(gene_set_name, "BP", sep = "_")]] <- BP_results
  
  # Perform GO enrichment analysis for Molecular Function (MF)
  MF_results <- perform_GO_analysis(gene_set_df, gene_set_name, "MF", output_folder)
  GOA_results_list[[paste(gene_set_name, "MF", sep = "_")]] <- MF_results
  
  # Perform GO enrichment analysis for Cellular Component (CC)
  CC_results <- perform_GO_analysis(gene_set_df, gene_set_name, "CC", output_folder)
  GOA_results_list[[paste(gene_set_name, "CC", sep = "_")]] <- CC_results
}

















