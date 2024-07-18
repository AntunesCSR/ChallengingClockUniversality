######################################################################################################
###################### Gene Ontology (GO) Enrichment Analysis for Relevant CpGs ######################
###################### identified by the clustering analysis of the clock CpGs  ######################
######################################################################################################

### Load necessary packages
library(clusterProfiler)
library(ggplot2)
library(here)
library(org.Hs.eg.db)

### Define Functions ### 
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
  
  # Create a subfolder for the specific gene set
  gene_set_folder <- file.path(output_folder, gene_set_name)
  if (!dir.exists(gene_set_folder)) {
    dir.create(gene_set_folder, recursive = TRUE)
  }
  
  # Save results to CSV file
  write.csv(GOA_results, file = file.path(gene_set_folder, paste0(gene_set_name, "_GOA_", ont_type, "_results.csv")))
  
  # Visualize with dotplot (if GOA_results is not empty)
  if (nrow(GOA_results) > 0) {
    dotplot(GOA_results, showCategory = 10, title = paste("Top 10 GO", ont_type, "terms for", gene_set_name))
    # Save dotplot visualization to PNG file
    ggsave(file.path(gene_set_folder, paste0(gene_set_name, "_GOA_", ont_type, "_results_dotplot.png")), 
           plot = last_plot())
  }
  
  # Visualize with barplot (if GOA_results is not empty)
  if (nrow(GOA_results) > 0) {
    barplot(GOA_results, showCategory = 10, title = paste("Top 10 GO", ont_type, "terms for", gene_set_name))
    # Save barplot visualization to PNG file
    ggsave(file.path(gene_set_folder, paste0(gene_set_name, "_GOA_", ont_type, "_results_barplot.png")), 
           plot = last_plot())
  }
  
  # Return the GO enrichment results
  return(GOA_results)
}

#### LOAD THE DATA ####
# List of filenames
file_names <- c(
  "all_clocks_high_methylated_cpgs_annotated_relevant_CpGs.csv",
  "all_clocks_low_methylated_cpgs_annotated_relevant_CpGs.csv",
  "all_clocks_variable_methylated_cpgs_annotated_relevant_CpGs.csv",
  "anova_all_clocks_significant_cpgs_annotated_relevant_CpGs.csv"
) 

# Loop through the filenames and read the data into a list
cpgs_data_list <- list() 

for (file_name in file_names) {
  file_path <- here("FunctionalAnalysis/Annotations/Annotations Results/Cross-Tissue Methylation", file_name)
  cpgs_data_list[[file_name]] <- read.csv(file_path, stringsAsFactors = FALSE)
}

# Create a list of data frames for the genes seets
gene_set_data_list <- list(
  "high_methylated_cpgs" = cpgs_data_list[[file_names[1]]],
  "low_methylated_cpgs" = cpgs_data_list[[file_names[2]]],
  "variable_methylated_cpgs" = cpgs_data_list[[file_names[3]]],
  "sig_variable_cpgs" = cpgs_data_list[[file_names[4]]]
)

### Set the output folder
output_folder <- here("FunctionalAnalysis/Functional_Enrichment/Go_Analysis_CrossTissue/Results/")

### Create a list to store the results of the gene set enrichment analysis
GOA_results_list <- list()

### Loop through each data frame and perform gene set enrichment analysis for BP, MF, and CC
for (gene_set_name in names(gene_set_data_list)) {
  gene_set_df <- gene_set_data_list[[gene_set_name]]
  
  # Modify the col gene to Gene for each of the data frames
  colnames(gene_set_df)[colnames(gene_set_df) == "gene"] <- "Gene"
  
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













