#!/usr/bin/env Rscript

# =========================
# LIBRARIES
# =========================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# =========================
# ARGUMENTS
# =========================
args <- commandArgs(trailingOnly = TRUE)
deseq_file <- args[1]

# =========================
# LOAD DATA
# =========================
res <- read_tsv(deseq_file)

# Zmień nazwę kolumny z genami, jeśli trzeba
gene_col <- if ("gene" %in% colnames(res)) "gene" else "gene_id"

res <- res %>%
  filter(!is.na(padj), !is.na(log2FoldChange))

# =========================
# ID CONVERSION (ENSEMBL → ENTREZ)
# =========================
gene_map <- bitr(
  res[[gene_col]],
  fromType = "REFSEQ",
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Hs.eg.db
)


res <- res %>%
  inner_join(gene_map, by = c("gene" = "REFSEQ"))


# =========================
# DIFFERENTIALLY EXPRESSED GENES
# =========================
deg <- res %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

deg_ids <- deg$ENTREZID

# =========================
# GO ENRICHMENT (BP)
# =========================
ego <- enrichGO(
  gene          = deg_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

write_tsv(as.data.frame(ego), "GO_BP_enrichment.tsv")

png("GO_BP_dotplot.png", width = 2000, height = 1600, res = 300)
dotplot(ego, showCategory = 20) + ggtitle("GO Biological Process")
dev.off()

# =========================
# KEGG ENRICHMENT
# =========================
ekk <- enrichKEGG(
  gene         = deg_ids,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

write_tsv(as.data.frame(ekk), "KEGG_enrichment.tsv")

png("KEGG_dotplot.png", width = 2000, height = 1600, res = 300)
dotplot(ekk, showCategory = 20) + ggtitle("KEGG Pathways")
dev.off()

# -------------------------------
# REMOVE DUPLICATE ENTREZ IDS
# -------------------------------
res_unique <- res %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1) %>%
  ungroup()

gene_list <- res_unique$log2FoldChange
names(gene_list) <- res_unique$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_go <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  keyType      = "ENTREZID",
  pvalueCutoff = 0.05
)

write_tsv(as.data.frame(gsea_go), "GSEA_GO_BP.tsv")

png("GSEA_GO_BP.png", width = 2000, height = 1600, res = 300)
ridgeplot(gsea_go) + ggtitle("GSEA GO BP")
dev.off()

# =========================
# DONE
# =========================
message("Functional enrichment analysis finished successfully.")
