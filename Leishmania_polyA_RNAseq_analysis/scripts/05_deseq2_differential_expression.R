library(DESeq2)
library(tidyverse)

counts_file <- "data/counts/htseq_counts_matrix.tsv"
metadata_file <- "data/metadata/sample_metadata.tsv"

counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_file, row.names = 1, check.names = FALSE)

# Ensure sample order matches
counts <- counts[, rownames(metadata)]

# Keep only mRNAs with at least one read in all libraries
keep <- rowSums(counts >= 1) == ncol(counts)
counts_filtered <- counts[keep, ]

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "treatment", "control"))

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    fold_change = 2^log2FoldChange,
    significant = ifelse(
      abs(log2FoldChange) >= 1 & padj < 0.01,
      "significant",
      "not_significant"
    )
  )

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

write.table(
  res_df,
  file = "results/deseq2/DESeq2_all_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sig_df <- res_df %>%
  filter(significant == "significant")

write.table(
  sig_df,
  file = "results/deseq2/DESeq2_significant_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
