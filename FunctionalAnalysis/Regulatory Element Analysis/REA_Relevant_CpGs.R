# Load required libraries
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

# Set up AnnotationHub
ah <- AnnotationHub()

# Function to create GRanges object from CpG data
create_cpg_granges <- function(cpg_data) {
  GRanges(
    seqnames = cpg_data$chr,
    ranges = IRanges(start = cpg_data$start, end = cpg_data$end),
    mcols = cpg_data[, !(colnames(cpg_data) %in% c("chr", "start", "end"))]
  )
}

# Load your CpG data
# Replace this with your actual data loading code
cpg_data <- read.csv("your_cpg_data.csv")
cpg_granges <- create_cpg_granges(cpg_data)

# Load regulatory element data from AnnotationHub
# Example: ENCODE transcription factor binding sites
encode_tfbs <- ah[["AH33540"]]

# Overlap CpGs with regulatory elements
overlaps <- findOverlaps(cpg_granges, encode_tfbs)

# Summarize overlaps
overlap_summary <- table(encode_tfbs$name[subjectHits(overlaps)])
overlap_summary <- sort(overlap_summary, decreasing = TRUE)

# Plot top overlapping transcription factors
top_n <- 20
plot_data <- data.frame(
  TF = names(overlap_summary)[1:top_n],
  Count = as.numeric(overlap_summary)[1:top_n]
)

ggplot(plot_data, aes(x = reorder(TF, Count), y = Count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top Overlapping Transcription Factors",
       x = "Transcription Factor",
       y = "Number of CpG Overlaps") +
  theme_minimal()

# Function to get nearby genes
get_nearby_genes <- function(gr, distance = 5000) {
  genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  nearby <- distanceToNearest(gr, genes, ignore.strand = TRUE)
  nearby_filtered <- nearby[mcols(nearby)$distance <= distance]
  unique(genes$gene_id[subjectHits(nearby_filtered)])
}

# Get nearby genes for CpGs overlapping regulatory elements
overlapping_cpgs <- cpg_granges[queryHits(overlaps)]
nearby_genes <- get_nearby_genes(overlapping_cpgs)

# Perform GO enrichment analysis on nearby genes
go_enrichment <- enrichGO(gene = nearby_genes,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

# Plot GO enrichment results
dotplot(go_enrichment, showCategory = 20) +
  ggtitle("GO Enrichment of Genes Near CpGs Overlapping Regulatory Elements")

# Save results
write.csv(as.data.frame(overlap_summary), "tf_overlap_summary.csv")
write.csv(as.data.frame(go_enrichment), "go_enrichment_results.csv")
