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
library(jsonlite)

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


# Function to retrieve protein-protein interaction data from BioGRI
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
                       filters = "external_gene_name", # To use the EntrezID, otherwise "external_gene_name" for symbol
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


# Function to retrieve protein-protein interaction data from STRING          ############# FIX THIS ONE
retrieve_ppi_data_string <- function(gene_list) {
  string_db <- STRINGdb$new(version = "11", score_threshold = 700)
  string_ppi_data <- string_db$get_interactions(gene_list)
  return(string_ppi_data)
}


############# FIX FROM HERE
# Function to clean and combine gene annotations
clean_and_combine_annotations <- function(source_list) {
  # Merge data frames based on the common gene identifier
  combined_annotations <- Reduce(function(df1, df2) {
    merge(df1, df2, by = "common_gene_identifier", all = TRUE)   #change the common_gene_identifier to the common column name
    
  }, annotation_list)
  
  # Combine columns and handle missing values
  combined_annotations <- combined_annotations %>%
    mutate(ensembl_gene_id = coalesce(ensembl_gene_id.x, ensembl_gene_id.y),
           description = coalesce(description.x, description.y)) %>%
    select(common_gene_identifier, ensembl_gene_id, description)
  
  return(combined_annotations)
}



#### RETRIEVE ANNOTATIONS ####

# Define the gene lists
gene_lists_s <- list('ANK1')
gene_lists_e <- list('268') 

# 1.1 Retrieve Gene Annotations from BioGrid
gene_annotations_biogrid <- lapply(gene_lists_s, retrieve_ppi_data_biogrid)
gene_ppi_biogrid_edge_list <- lapply(gene_lists_s, retrieve_ppi_data_biogrid_edge_list)
print(gene_annotations_biogrid)
print(gene_ppi_biogrid_edge_list)


# # 1.2. Retrieve Gene Annotations from DAVID
# gene_annotations_david <- lapply(gene_lists, retrieve_functional_classification_david)
# print(gene_annotations_david)


# 1.3. Retrieve Gene Annotations from Ensembl
gene_annotations_ensembl <- lapply(gene_lists, retrieve_annotations_ensembl)
print(gene_annotations_ensembl)

# 1.4. Retrieve Gene Annotations from ENTREZ
gene_annotations_entrez <- lapply(gene_lists, retrieve_annotations_entrez)
print(gene_annotations_entrez)


# 1.5. Retrieve Gene Annotations from KEGG
gene_pathway_kegg <- lapply(gene_lists_e, retrieve_pathway_data_kegg)
print(gene_pathway_kegg)


# 1.6. Retrieve Gene Annotations from STRING
gene_annotations_string <- lapply(gene_lists, retrieve_ppi_data_string)  
print(gene_annotations_string)   

# 1.7. Retrieve Gene Annotations from UniProt
gene_annotations_uniprot <- lapply(gene_lists, retrieve_annotations_uniprot)
print(gene_annotations_uniprot)


# 2. Define the sources. Here we can remove any unwanted source from the list so it doesn't appear
# in the final combined annotations.
source_list <- list(biogrid_ppi_data, # BioGrid PPI data
                    david_enrichment_data, # DAVID functional enrichment data
                    ensembl_annotations_data, # ENSEMBL gene annotations
                    entrez_annotations_data, # NCBI ENTREZ gene annotations
                    kegg_pathway_data, # KEGG pathway data
                    string_ppi_data, # STRING PPI data
                    uniprot_annotations_data # UniProt gene annotations
)



# 3. Clean and Combine Gene Annotations
final_combined_annotations <- clean_and_combine_annotations(source_list)


# 4. Print the combined annotations
print("Final Combined Gene Annotations:")
print(final_combined_annotations)


# 5. Save the combined annotations to a file
write.csv(final_combined_annotations, "combined_annotations.csv", row.names = FALSE)



# #### RETRIEVE ANNOTATIONS - LOOP####
# 
# # Define the gene lists for each clock separately
# gene_lists <- list(top25_clk1_s) #, top25_clk2_s, top25_clk3_s, top25_ov1_2_3_s, top25_ov2_3_s)
# 
# # Loop over each clock separately
# for (i in seq_along(gene_lists)) {
#   
#   # Retrieve Gene Annotations from BioGrid 
#   gene_annotations_biogrid <- retrieve_ppi_data_biogrid(gene_lists[[i]])
#   print(gene_annotations_biogrid)
#   
#   # Retrieve Gene Annotations from DAVID
#   gene_annotations_david <- retrieve_functional_enrichment_david(gene_lists[[i]])
#   print(gene_annotations_david)
#   
#   # Retrieve Gene Annotations from Ensembl 
#   gene_annotations_ensembl <- retrieve_annotations_ensembl(gene_lists[[i]])
#   print(gene_annotations_ensembl)
#   
#   # Retrieve Gene Annotations from ENTREZ
#   gene_annotations_entrez <- retrieve_annotations_entrez(gene_lists[[i]])
#   print(gene_annotations_entrez)
#   
#   # Retrieve Gene Annotations from KEGG
#   gene_annotations_kegg <- retrieve_pathway_data_kegg(gene_lists[[i]])
#   print(gene_annotations_kegg)
#   
#   # Retrieve Gene Annotations from STRING
#   gene_annotations_string <- retrieve_ppi_data_string(gene_lists[[i]])
#   print(gene_annotations_string)
#   
#   # Retrieve Gene Annotations from UniProt
#   gene_annotations_uniprot <- retrieve_annotations_uniprot(gene_lists[[i]])
#   print(gene_annotations_uniprot)
#   
#   # Define the sources for each clock
#   source_list <- list(
#     biogrid_ppi_data = gene_annotations_biogrid,
#     david_enrichment_data = gene_annotations_david,
#     ensembl_annotations_data = gene_annotations_ensembl,
#     entrez_annotations_data = gene_annotations_entrez,
#     kegg_pathway_data = gene_annotations_kegg,
#     string_ppi_data = gene_annotations_string,
#     uniprot_annotations_data = gene_annotations_uniprot
#   )
#   
#   # Clean and Combine Gene Annotations
#   final_combined_annotations <- clean_and_combine_annotations(source_list)
#   
#   # Print the combined annotations
#   print("Final Combined Gene Annotations:")
#   print(final_combined_annotations)
#   
#   # Save the combined annotations to a file with a name reflecting the clock being processed
#   write.csv(final_combined_annotations, paste0("combined_annotations_clock", i, ".csv"), row.names = FALSE)
# }
# 
