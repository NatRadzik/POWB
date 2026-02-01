#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

#parsowanie argumentów
args <- commandArgs(trailingOnly = TRUE)
deseq_file <- args[1]
pval_threshold  <- as.numeric(args[2])
logfc_threshold <- as.numeric(args[3])

#wczytaj wyniki
res <- read_tsv(deseq_file)

# Nazwy kolumn się upweniam
stopifnot(all(c("padj", "log2FoldChange", "gene") %in% colnames(res)))


#ID CONVERSION
# Mamy Refseq Id ("NM_033031") checmy mieć nazwy "ENTREZID"
gene_map <- bitr(
  unique(res$gene),
  fromType = "REFSEQ",
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Hs.eg.db
)

# JOin z mapowaniem
res <- res %>%
  inner_join(gene_map, by = c("gene" = "REFSEQ"))

message("Liczba genów po mapowaniu: ", nrow(res))


#Filtrowanie DEG
deg <- res %>%
  filter(
    padj < pval_threshold,
    abs(log2FoldChange) > logfc_threshold
  )

deg_ids <- unique(deg$ENTREZID)
message("Liczba DEG po mapowaniu: ", length(deg_ids))

# GO ENRICHMENT (BP)
if (length(deg_ids) >= 20) {

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

# KEGG ENRICHMENT
ekk <- enrichKEGG(
  gene         = deg_ids,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

write_tsv(as.data.frame(ekk), "KEGG_enrichment.tsv")

png("KEGG_dotplot.png", width = 2000, height = 1600, res = 300)
dotplot(ekk, showCategory = 20) + ggtitle("KEGG Pathways")
dev.off()
}


#GSEA

# Entrez usuń duplikaty
res_unique <- res %>%
  arrange(desc(abs(log2FoldChange))) %>%
  distinct(ENTREZID, .keep_all = TRUE)


#Przygotowanie danych
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
