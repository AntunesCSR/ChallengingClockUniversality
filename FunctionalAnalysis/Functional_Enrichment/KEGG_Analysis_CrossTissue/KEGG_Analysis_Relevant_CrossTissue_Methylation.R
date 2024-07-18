### KEGG Analysis for relevant CpGs 

# AUthor: Dr. Catia Antunes
# Date: 15/07/2024
# Description: This script performs KEGG analysis relevant CpGs 


### Load necessary packages
library(clusterProfiler)
library(ggplot2)
library(here)
library(org.Hs.eg.db)

### Define Functions ###
perform_KEGG_analysis <- function(gene_set_df, gene_set_name, output_folder) {
  if (nrow(gene_set_df) == 0) {
    warning(paste("The gene set for", gene_set_name, "is empty. Skipping KEGG analysis."))
    return(NULL)
  }
  
  entrez_ids <- gene_set_df$ENTREZID
  entrez_ids <- unique(entrez_ids[!is.na(entrez_ids) & entrez_ids != ""])
  
  print(paste("Performing KEGG analysis for", gene_set_name))
  print(head(entrez_ids))
  
  if (length(entrez_ids) == 0) {
    warning(paste("No valid ENTREZ IDs for", gene_set_name, ". Skipping KEGG analysis."))
    return(NULL)
  }
  
  KEGG_results <- enrichKEGG(
    gene = entrez_ids,
    organism = 'hsa',
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  if (is.null(KEGG_results) || nrow(KEGG_results) == 0) {
    warning(paste("No KEGG pathways found for", gene_set_name))
    return(NULL)
  }
  
  gene_set_folder <- file.path(output_folder, gene_set_name)
  if (!dir.exists(gene_set_folder)) {
    dir.create(gene_set_folder, recursive = TRUE)
  }
  
  write.csv(KEGG_results, file = file.path(gene_set_folder, paste0(gene_set_name, "_KEGG_results.csv")))
  
  if (nrow(KEGG_results) > 0) {
    dotplot(KEGG_results, showCategory = 10, title = paste("Top 10 KEGG pathways for", gene_set_name))
    ggsave(file.path(gene_set_folder, paste0(gene_set_name, "_KEGG_results_dotplot.png")), plot = last_plot())
  }
  
  if (nrow(KEGG_results) > 0) {
    barplot(KEGG_results, showCategory = 10, title = paste("Top 10 KEGG pathways for", gene_set_name))
    ggsave(file.path(gene_set_folder, paste0(gene_set_name, "_KEGG_results_barplot.png")), plot = last_plot())
  }
  
  return(KEGG_results)
}




#### LOAD THE DATA ####

# List of filenames
file_names <- c(
  "all_clocks_high_methylated_cpgs_annotated_relevant_CpGs.csv",
  "all_clocks_low_methylated_cpgs_annotated_relevant_CpGs.csv",
  "all_clocks_variable_methylated_cpgs_annotated_relevant_CpGs.csv",
  "anova_all_clocks_significant_cpgs_annotated_relevant_CpGs.csv"
) 

# Load data into a list taking into account the location of the data
cpgs_data_list <- list()
for (file_name in file_names) {
  file_path <- here("FunctionalAnalysis/Annotations/Annotations Results/Cross-Tissue Methylation", file_name)
  cpgs_data_list[[file_name]] <- read.csv(file_path, stringsAsFactors = FALSE)
}

# Create a list of data frames for the genes of each clock dataset
gene_set_data_list <- list(
  "high_methylated_cpgs" = cpgs_data_list[[file_names[1]]],
  "low_methylated_cpgs" = cpgs_data_list[[file_names[2]]],
  "variable_methylated_cpgs" = cpgs_data_list[[file_names[3]]],
  "sig_variable_cpgs" = cpgs_data_list[[file_names[4]]]
)

# Select the output folder
output_folder <- here("FunctionalAnalysis/Functional_Enrichment/KEGG_Analysis_CrossTissue")
KEGG_results_list <- list()

# Run the KEGG analysis for each gene set
for (gene_set_name in names(gene_set_data_list)) {
  gene_set_df <- gene_set_data_list[[gene_set_name]]
  
  KEGG_results <- perform_KEGG_analysis(gene_set_df, gene_set_name, output_folder)
  KEGG_results_list[[gene_set_name]] <- KEGG_results
}

