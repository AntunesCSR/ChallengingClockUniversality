#### LOAD LIBRARIES ####
library(igraph)


#### LOAD THE DATA ####

# Load genes from csv file for network analysis, select column Gene
genes <- read.csv("./FunctionalAnalysis/top_25_cpgs_clock_1.csv", header = TRUE, sep = ",")
genes_symbol <- genes$Gene
genes_entrezid <- genes$ENTREZID

# Create a fully connected graph
gene_symbol_graph <- graph.full(n = length(genes_symbol), directed = FALSE)
genes_entrezid_graph <- graph.full(n = length(genes_entrezid), directed = FALSE)

# Plot the graph
plot(gene_symbol_graph, layout = layout.circle(gene_symbol_graph), vertex.label = genes_symbol, main = "Gene Network (Symbols)")
plot(genes_entrezid_graph, layout = layout.circle(genes_entrezid_graph), vertex.label = genes_entrezid, main = "Gene Network (EntrezID)")

# Save the graph
png("gene_symbol_graph.png")
png("genes_entrezid_graph.png")