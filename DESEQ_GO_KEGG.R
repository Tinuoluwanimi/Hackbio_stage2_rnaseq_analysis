# DE_analysis_plant.R
# Differential expression + GO/KEGG enrichment for Arabidopsis
# Assumes:
#  - counts.txt (featureCounts output) in working dir
#  - metadata.csv (sample,condition) in working dir
#  - working dir set to C:/Users/User/Desktop/planta

# Set working directory
setwd("C:/Users/User/Desktop/planta")   

# Load libraries
  library(DESeq2)
  library(pheatmap)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.At.tair.db)   
  library(enrichplot)

#  Read counts + metadata
# Read featureCounts output 
fc <- read.delim("counts.txt", comment.char = "#", header = TRUE, stringsAsFactors = FALSE)

# Read metadata (CSV); should have columns: sample,condition
meta <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# Check matching
if(!all(c("sample","condition") %in% colnames(meta))){
  stop("metadata.csv must have columns named 'sample' and 'condition'.")
}
# Extract only count columns and clean colnames
# keep Geneid + all sample columns (7:ncol)
counts_raw <- fc[, c(1, 7:ncol(fc))]
rownames(counts_raw) <- counts_raw$Geneid
counts_raw <- counts_raw[, -1, drop = FALSE]   # now pure numeric counts but as characters

# Clean column names to match metadata sample IDs
# Example original names: plants.alignments.IGV.SRR12808497_Aligned.sortedByCoord.out.bam
colnames(counts_raw) <- colnames(counts_raw) %>%
  gsub("plants\\.alignments\\.IGV\\.", "", .) %>%
  gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", .)

# ensure metadata sample names line up
meta$sample <- as.character(meta$sample)

# Check matching
if(!all(colnames(counts_raw) %in% meta$sample)){
  cat("WARNING: Some count columns do not match metadata sample IDs.\n")
  cat("Count columns:\n"); print(colnames(counts_raw))
  cat("Metadata samples:\n"); print(meta$sample)
  stop("Fix sample names so counts and metadata match.")
}

# Reorder metadata rows to match count columns
meta <- meta[match(colnames(counts_raw), meta$sample), , drop = FALSE]

#  Convert to numeric matrix and clean NAs / zeros
counts_mat <- apply(counts_raw, 2, as.integer)
rownames(counts_mat) <- rownames(counts_raw)

# Replace any NA with 0 
if(any(is.na(counts_mat))){
  warning("NA values found in count matrix — replacing NAs with 0.")
  counts_mat[is.na(counts_mat)] <- 0
}

# Remove genes with zero counts across all samples
counts_mat <- counts_mat[rowSums(counts_mat) > 0, , drop = FALSE]

# Final checks
if(any(is.na(counts_mat))) stop("NA remain in counts matrix after cleaning.")
message("Count matrix dims (genes x samples): ", paste(dim(counts_mat), collapse = " x "))

# Build DESeq2 object and run DE

# make condition factor with control first
meta$condition <- factor(meta$condition, levels = c("control", "treatment"))

dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = meta,
                              design    = ~ condition)

dds <- DESeq(dds)
res <- results(dds)

# Save raw results
write.csv(as.data.frame(res), "Final_DESeq2_Results.csv", row.names = TRUE)

# Identify DEGs (up/down)
upregulated   <- subset(res, padj < 0.05 & log2FoldChange > 1)
downregulated <- subset(res, padj < 0.05 & log2FoldChange < -1)

write.csv(as.data.frame(upregulated), "Upregulated_Genes.csv", row.names = TRUE)
write.csv(as.data.frame(downregulated), "Downregulated_Genes.csv", row.names = TRUE)

#  Volcano plot (save PNG)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

png("VolcanoPlot.png", width = 1400, height = 1000, res = 150)
plot(res_df$log2FoldChange, -log10(res_df$padj),
     pch = 19, cex = 0.5, col = "grey",
     xlab = "Log2 Fold Change", ylab = "-log10 Adjusted P-value",
     main = "Volcano Plot: Treatment vs Control")
abline(v = c(-1, 1), col = "darkblue", lty = 2)
abline(h = -log10(0.05), col = "darkred", lty = 2)

# highlight
with(subset(res_df, padj < 0.05 & log2FoldChange > 1),
     points(log2FoldChange, -log10(padj), col = "salmon", pch = 19, cex = 0.7))
with(subset(res_df, padj < 0.05 & log2FoldChange < -1),
     points(log2FoldChange, -log10(padj), col = "skyblue", pch = 19, cex = 0.7))

legend("topright", legend = c("Upregulated", "Downregulated"),
       col = c("salmon", "skyblue"), pch = 19)
dev.off()

# Heatmap of top 50 DEGs 
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# get significant genes, order by padj, select top 50
sig_degs <- res[!is.na(res$padj) & res$padj < 0.05, , drop = FALSE]
if(nrow(sig_degs) == 0) {
  warning("No significant DEGs (padj < 0.05). Using top 50 by smallest padj nonetheless.")
  order_idx <- order(res$padj, na.last = NA)
  top50_genes <- rownames(res)[head(order_idx, 50)]
} else {
  top50_genes <- rownames(sig_degs)[order(sig_degs$padj)][1:min(50, nrow(sig_degs))]
}
vsd_top50 <- vsd_mat[top50_genes, , drop = FALSE]

# Prepare annotation (samples)
annotation_col <- data.frame(Condition = meta$condition)
rownames(annotation_col) <- meta$sample

png("Heatmap_Top50.png", width = 1400, height = 1000, res = 150)
pheatmap(vsd_top50,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         color = colorRampPalette(c("blue", "white", "red"))(60),
         main = "Top 50 DEGs Heatmap")
dev.off()

# Extract DEGs and run GO/KEGG enrichment

# Use all significant DEGs (padj < 0.05) for enrichment
sig_degs_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  filter(!is.na(padj) & padj < 0.05)

# Save top 100 by abs(log2FC) for reporting
top100_degs <- sig_degs_df %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(100)
write.csv(top100_degs, "Top100_DEGs.csv", row.names = FALSE)

# Map TAIR IDs to ENTREZ IDs
gene_list <- sig_degs_df$gene
gene2entrez <- mapIds(org.At.tair.db, keys = gene_list,
                      column = "ENTREZID", keytype = "TAIR", multiVals = "first")
entrez_ids <- na.omit(gene2entrez)

if(length(entrez_ids) == 0){
  warning("No ENTREZ IDs mapped; skipping GO/KEGG enrichment.")
} else {
  # -----------------------------
  # GO Enrichment (BP, MF, CC)
  # -----------------------------
  ego_BP <- enrichGO(gene = entrez_ids, OrgDb = org.At.tair.db,
                     keyType = "ENTREZID", ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  ego_MF <- enrichGO(gene = entrez_ids, OrgDb = org.At.tair.db,
                     keyType = "ENTREZID", ont = "MF",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  ego_CC <- enrichGO(gene = entrez_ids, OrgDb = org.At.tair.db,
                     keyType = "ENTREZID", ont = "CC",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

  # Save results 
  write.csv(as.data.frame(ego_BP), "GO_BP_results.csv", row.names = FALSE)
  write.csv(as.data.frame(ego_MF), "GO_MF_results.csv", row.names = FALSE)
  write.csv(as.data.frame(ego_CC), "GO_CC_results.csv", row.names = FALSE)

  # Dotplots
  if(nrow(as.data.frame(ego_BP)) > 0){
    png("GO_BP_Dotplot.png", width = 1400, height = 1000, res = 150)
    print(dotplot(ego_BP, showCategory = 15) + ggtitle("GO: Biological Process"))
    dev.off()
  }

  if(nrow(as.data.frame(ego_MF)) > 0){
    png("GO_MF_Dotplot.png", width = 1400, height = 1000, res = 150)
    print(dotplot(ego_MF, showCategory = 15) + ggtitle("GO: Molecular Function"))
    dev.off()
  }

  if(nrow(as.data.frame(ego_CC)) > 0){
    png("GO_CC_Dotplot.png", width = 1400, height = 1000, res = 150)
    print(dotplot(ego_CC, showCategory = 15) + ggtitle("GO: Cellular Component"))
    dev.off()
  }

  #  KEGG Enrichment for Top 5 pathways using TAIR IDs
ekegg <- tryCatch({
  enrichKEGG(gene         = sig_degs_df$gene,   # use TAIR IDs directly
             organism     = "ath",
             keyType      = "kegg",             # for KEGG, use TAIR
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             qvalueCutoff = 0.05)
}, error = function(e) {
  warning("KEGG enrichment failed: ", conditionMessage(e))
  return(NULL)
})

if(!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0){
  kegg_df <- as.data.frame(ekegg)
  write.csv(kegg_df, "KEGG_Enrichment.csv", row.names = FALSE)

  top5_kegg <- kegg_df %>% arrange(p.adjust) %>% head(5)
  write.csv(top5_kegg, "Top5_KEGG_Pathways.csv", row.names = FALSE)

  png("KEGG_Dotplot.png", 1400, 1000, res = 150)
  print(dotplot(ekegg, showCategory = 5) + ggtitle("Top 5 KEGG Pathways"))
  dev.off()
} else {
  warning("No significant KEGG pathways found.")
}

