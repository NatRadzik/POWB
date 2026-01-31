#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

counts_file <- args[1]
meta_file   <- args[2]
pval_threshold  <- as.numeric(args[3])
logfc_threshold <- as.numeric(args[4])


suppressMessages({
  library(DESeq2)
  library(ggplot2)
})

#Wczytanie i przygotowanie countów
counts <- read.delim(
  counts_file,
  comment.char = "#",
  stringsAsFactors = FALSE
)

#GeneID jako rownames
rownames(counts) <- counts$Geneid

#Usunięcie kolumn technicznych
counts <- counts[, -(1:6)]

#Uproszczenie nazw kolumn (BAM → sample ID)
colnames(counts) <- gsub(
  "Aligned.sortedByCoord.out.bam",
  "",
  colnames(counts)
)
colnames(counts) <- gsub("\\.$", "", colnames(counts))

#Wczytanie metadata - co jest kontrolą a co treated plus uporządkowanie
meta <- read.delim(meta_file, header = TRUE, row.names = 1, sep = "", stringsAsFactors = FALSE)
meta <- meta[colnames(counts), , drop = FALSE]

#Sprawdzenie zgodności
stopifnot(all(colnames(counts) == rownames(meta)))

#DESeq2 
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

#Zapis tabeli wyników
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)


#Volcano plot
res_df$significant <- with(res_df, ifelse(!is.na(padj) & padj <= pval_threshold & abs(log2FoldChange) >= logfc_threshold, "yes", "no"))


#Write table with all
write.table(
  res_df,
  "deseq2_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#Filter down genes
res_df_significant <- subset(res_df, significant=="yes")  


#Write table with significant
write.table(
  res_df_significant,
  "deseq2_results_significant.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


png("volcano_plot.png", width = 800, height = 600)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("no" = "black", "yes" = "magenta")) +
  geom_vline(xintercept = c(-logfc_threshold , logfc_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Differential gene expression (DESeq2)", color = "Significant")

dev.off()


