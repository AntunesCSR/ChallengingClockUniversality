######################################################################################################
###################### Gene Ontology (GO) Enrichment Analysis for Relevant CpGs ######################
###################### identified by the clustering analysis of the clock CpGs  ######################
######################################################################################################

### Load necessary packages
library(clusterProfiler)
library(ggplot2)
library(here)
library(org.Hs.eg.db)


# 
# ### Define Functions ### 
# # Function that allow multiple GO enrichment analysis to be performed. 
# # The function takes the gene set data frame, gene set name, ontology type, and output folder as input parameters.
# 
# # Function to perform GO enrichment analysis and save results
# perform_GO_analysis <- function(gene_set_df, gene_set_name, ont_type, output_folder) {
#   gene_symbols <- gene_set_df$Gene # assuming the column is named 'Gene'
#   
#   # Perform GO enrichment analysis
#   GOA_results <- enrichGO(
#     gene = gene_symbols,
#     OrgDb = org.Hs.eg.db,
#     keyType = "SYMBOL",
#     ont = ont_type,
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
#   )
#   
#   # Set the working directory to the output folder
#   setwd(output_folder)
#   
#   # Save results to CSV file
#   write.csv(GOA_results, file = paste0(gene_set_name, "_GOA_", ont_type, "_results.csv"))
#   
#   # Visualize with dotplot (if GOA_results is not empty)
#   if (nrow(GOA_results) > 0) {
#     dotplot(GOA_results, showCategory = 10, title = paste("Top 10 GO", ont_type, "terms for", gene_set_name))
#     # Save dotplot visualization to PNG file
#     ggsave(paste0(gene_set_name, "_GOA_", ont_type, "_results_dotplot.png"), 
#            plot = last_plot())
#   }
#   
#   # Visualize with barplot (if GOA_results is not empty)
#   if (nrow(GOA_results) > 0) {
#     barplot(GOA_results, showCategory = 10, title = paste("Top 10 GO", ont_type, "terms for", gene_set_name))
#     # Save barplot visualization to PNG file
#     ggsave(paste0(gene_set_name, "_GOA_", ont_type, "_results_barplot.png"), 
#            plot = last_plot())
#   }
#   
#   # Return the GO enrichment results
#   return(GOA_results)
# }
# 
# 
# #### LOAD THE DATA ####
# # Load the data frames containing the relevant CpGs and associated annotations identified by the clustering analysis of the clock CpGs
# # In this case, load the files produced by the Annotation step. The data frames should contain the relevant CpGs and associated genes.
# # Files should look like this clockX_annotated_relevant_cpgs.csv.
# 
# # List of filenames (when many files are required. Add the file names to the list as required)
# file_names <- c(
#   "clock1_annotated_relevant_CpGs.csv",
#   "clock2_annotated_relevant_CpGs.csv",
#   "clock3_annotated_relevant_CpGs.csv",
#   "clock2_3_annotated_relevant_CpGs.csv"
# ) 
# 
# # Loop through the filenames and read the data into a list
# cpgs_data_list <- list() 
# 
# for (file_name in file_names) {
#   # create file path
#   file_path <- here("FunctionalAnalysis/Annotations/Annotations Results", file_name)
#   # read the data into a data frame
#   cpgs_data_list[[file_name]] <- read.csv(file_path, stringsAsFactors = FALSE) # StringsAsFactors = FALSE to avoid conversion of character to factor
# }
# 
# cpgs_data_list
# 
# # Create a list of data frames for the gene sets
# gene_set_data_list <- list(
#   "Clock 1" = cpgs_data_list[[1]],
#   "Clock 2" = cpgs_data_list[[2]],
#   "Clock 3" = cpgs_data_list[[3]],
#   "overlap of Clocks 2 & 3" = cpgs_data_list[[4]]
# )
# 
# # gene_set_data_list[["Clock 3"]]
# 
# ### Set the output folder
# output_folder <- here("FunctionalAnalysis/Functional_Enrichment/GO_Analysis_Relevant_CpGs/Results/")
# 
# 
# ### Create a list to store the results of the gene set enrichment analysis
# GOA_results_list <- list()
# 
# 
# ### Loop through each data frame and perform gene set enrichment analysis for BP, MF, and CC
# for (dataframe in names(gene_set_data_list)) {
#   gene_set_df <- gene_set_data_list[[dataframe]]
# 
#   # Modify the col gene to gene for each of the data frames
#   colnames(gene_set_df)[colnames(gene_set_df) == "gene"] <- "Gene"
#   
#   # Perform GO enrichment analysis for Biological Process (BP)
#   BP_results <- perform_GO_analysis(gene_set_df, gene_set_name, "BP", output_folder)
#   GOA_results_list[[paste(gene_set_name, "BP", sep = "_")]] <- BP_results
#   
#   # Perform GO enrichment analysis for Molecular Function (MF)
#   MF_results <- perform_GO_analysis(gene_set_df, gene_set_name, "MF", output_folder)
#   GOA_results_list[[paste(gene_set_name, "MF", sep = "_")]] <- MF_results
#   
#   # Perform GO enrichment analysis for Cellular Component (CC)
#   CC_results <- perform_GO_analysis(gene_set_df, gene_set_name, "CC", output_folder)
#   GOA_results_list[[paste(gene_set_name, "CC", sep = "_")]] <- CC_results
# }


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

#### LOAD THE DATA ####
# List of filenames
file_names <- c(
  "clock1_annotated_relevant_CpGs.csv",
  "clock2_annotated_relevant_CpGs.csv",
  "clock3_annotated_relevant_CpGs.csv",
  "clock2_3_annotated_relevant_CpGs.csv"
) 

# Loop through the filenames and read the data into a list
cpgs_data_list <- list() 

for (file_name in file_names) {
  file_path <- here("FunctionalAnalysis/Annotations/Annotations Results", file_name)
  cpgs_data_list[[file_name]] <- read.csv(file_path, stringsAsFactors = FALSE)
}

# Create a list of data frames for the gene sets
gene_set_data_list <- list(
  "Clock_1" = cpgs_data_list[[1]],
  "Clock_2" = cpgs_data_list[[2]],
  "Clock_3" = cpgs_data_list[[3]],
  "Clock_2_3_overlap" = cpgs_data_list[[4]]
)

### Set the output folder
output_folder <- here("FunctionalAnalysis/Functional_Enrichment/GO_Analysis_Relevant_CpGs/Results/")

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
















