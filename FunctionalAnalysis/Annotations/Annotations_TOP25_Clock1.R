# Author: Catia Antunes
# Date: 2024-02-29
# Description: Compiles biological data for the top 25 genes of all clocks,
# for further network analysis. The biological data is retrieved from BioGrid (PPI),
# DAVID (functional enrichment), Ensembl (gene annotations), ENTREZ (gene annotations),
# KEGG (pathway data), STRING (PPI), and UniProt (gene annotations).


#### IMPORT LIBRARIES ####
library(biomaRt)
library(STRINGdb)
library(dplyr)
library(purrr)
library(stringr)
library(httr)
library(org.Hs.eg.db)
library(KEGGREST)
library(KEGGgraph)
library(AnnotationHub)
library(AnnotationDbi)
library(jsonlite)
library(clusterProfiler)  # For pathway enrichment analysis
library(igraph)  # For network analysis

#### LOAD THE DATA ####

# Load top 25 genes from the clock1, clock2, clock3, overlap1_2_3 and overlap _2_3
top25_clk1 <- read.table("./FunctionalAnalysis/top_25_cpgs_clock_1.csv", header=TRUE, sep=",", quote="\"")
top25_clk1_s <- top25_clk1$Gene
top25_clk1_e <- top25_clk1$ENTREZID

top25_clk2 <- read.table("./FunctionalAnalysis/top_25_cpgs_clock_2.csv", header=TRUE, sep=",", quote="\"")
top25_clk2_s <- top25_clk2$Gene
top25_clk2_e <- top25_clk2$ENTREZID

top25_clk3 <- read.table("./FunctionalAnalysis/top_25_cpgs_clock_3.csv", header=TRUE, sep=",", quote="\"")
top25_clk3_s <- top25_clk3$Gene
top25_clk3_e <- top25_clk3$ENTREZID

top25_ov1_2_3 <- read.table("./FunctionalAnalysis/top_25_cpgs_overlap_1_2_3.csv", header=TRUE, sep=",", quote="\"")
top25_ov1_2_3_s <- top25_ov1_2_3$Gene
top25_ov1_2_3_e <- top25_ov1_2_3$ENTREZID

top25_ov2_3 <- read.table("./FunctionalAnalysis/top_25_cpgs_overlap_2_3.csv", header=TRUE, sep=",", quote="\"")
top25_ov2_3_s <- top25_ov2_3$Gene
top25_ov2_3_e <- top25_ov2_3$ENTREZID # Empty. WHY??????? T^T



#### FUNCTIONS ####


# Function to retrieve protein-protein interaction data from BioGRID
retrieve_ppi_data_biogrid <- function(gene_list) {
  api_key <- "e250241e4295ebc3a4aa90057f51cc3d"
  url <- paste0("https://webservice.thebiogrid.org/interactions/?accesskey=", api_key,
                "&geneList=", paste(gene_list, collapse = "|"), # List of genes to retrieve interactions for
                "&searchNames=TRUE", # Search for gene names in addition to gene 
                "&searchSynonyms=TRUE", # Search for gene synonyms in addition to gene IDs
                "&includeInteractors=TRUE", # Include interactors in the response
                "&includeInteractorsInteractions=FALSE", # TO NOT Include interactors interactions in the response
                "&includeHeader=TRUE", # Include the header in the response
                "&format=json", # Return the data in JSON format. standard option .tab2
                "&max=10" # Maximum number of interactions to retrieve
  )
  response <- httr::GET(url)
  biogrid_ppi_data <- httr::content(response, "parsed") 
  
  # # Extracting only the required fields and storing in a list
  # extracted_data <- lapply(biogrid_ppi_data, function(x) {
  #   list(
  #     BIOGRID_INTERACTION_ID = x$BIOGRID_INTERACTION_ID,
  #     ENTREZ_GENE_A = x$ENTREZ_GENE_A,
  #     ENTREZ_GENE_B = x$ENTREZ_GENE_B
  #   )
  # })
  # 
  # return(extracted_data)
  
  # Extracting only the required fields and storing in an edge list
  edge_list <- lapply(biogrid_ppi_data, function(x) {
    c(source = x$ENTREZ_GENE_A, target = x$ENTREZ_GENE_B)
  })
  
  return(edge_list)

}


# Function to retrieve gene annotations from NCBI Entrez
retrieve_annotations_entrez <- function(gene_list) {
  entrez <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  entrez_annotations_data <- getBM(attributes = c("entrezgene_id", "external_gene_name", "description"), # external_gene_name correct
                       filters = "entrezgene_id", # To use the EntrezID, otherwise "external_gene_name" for symbol
                       values = gene_list,
                       mart = entrez)
  return(entrez_annotations_data)
}


# Function to retrieve pathway data from KEGG (REQUIRES NCBI GENE ID's)
retrieve_pathway_data_kegg <- function(gene_list) {
  
  # Convert NCBI Gene IDs to KEGG IDs
  kegg_ids <- lapply(gene_list, function(gene_id) {
    conv_result <- tryCatch(keggConv("genes", paste0("ncbi-geneid:", gene_id)), error = function(e) NA)
    return(conv_result)
  })
  
  # Print the KEGG IDs
  print(kegg_ids)
  
  # Filter out NAs
  kegg_ids <- kegg_ids[!sapply(kegg_ids, function(x) identical(x, NA))]
  
  # Retrieve pathway data using the converted KEGG IDs
  if (length(kegg_ids) > 0) {
    kegg_pathway_data <- keggLink("pathway", kegg_ids)
    return(kegg_pathway_data)
  } else {
    warning("No valid KEGG IDs found.")
    return(NULL)
  }
}




#### RETRIEVE ANNOTATIONS - NO LOOOP - JSUT TESTING ####

# Define the gene lists
gene_lists_s <- list('ANK1', 'RPS6KA3', 'RPS6KA3', 'RPS6KA3', 'RPS6KA3')
gene_lists_e <- list('268', '6195', '6195', '6195', '6195') 

# 1.1 Retrieve PPI data from BioGrid
ppi_data <- lapply(gene_lists_s, retrieve_ppi_data_biogrid_edge_list)
print(ppi_data)


# 1.2. Retrieve Gene Annotations from ENTREZ - NCBI GENE ID's
annotations_data <- lapply(gene_lists_e, retrieve_annotations_entrez)
print(annotations_data)


# 1.3. Retrieve pathway from KEGG - NEEDS NCBI GENE ID's
pathway_data <- lapply(gene_lists_e, retrieve_pathway_data_kegg)
print(pathway_data)


#### LOOP ###

# Define the gene lists
gene_lists_s <- list('ANK1', 'RPS6KA3', 'RPS6KA3', 'RPS6KA3', 'RPS6KA3')
gene_lists_e <- list('268', '6195', '6195', '6195', '6195')

# Initialize empty lists to store results for each gene list
ppi_data_list <- list()
annotations_data_list <- list()
pathway_data_list <- list()

# Iterate over each gene list
for (i in seq_along(gene_lists_e)) {
  # 1.1 Retrieve PPI data from BioGrid
  ppi_data <- retrieve_ppi_data_biogrid_edge_list(gene_lists_s[[i]])
  ppi_data_list[[i]] <- ppi_data
  
  # 1.2 Retrieve Gene Annotations from ENTREZ - NCBI GENE ID's
  annotations_data <- retrieve_annotations_entrez(gene_lists_e[[i]])
  annotations_data_list[[i]] <- annotations_data
  
  # 1.3 Retrieve pathway from KEGG - NEEDS NCBI GENE ID's
  pathway_data <- retrieve_pathway_data_kegg(gene_lists_e[[i]])
  pathway_data_list[[i]] <- pathway_data
}

# Print the results for each gene list
print(ppi_data_list)
print(annotations_data_list)
print(pathway_data_list)





