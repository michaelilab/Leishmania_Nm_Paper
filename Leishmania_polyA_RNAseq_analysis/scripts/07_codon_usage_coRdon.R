library(coRdon)
library(tidyverse)
library(Biostrings)

cds_file <- "data/reference/GenesByGeneModelChars_LM_11July2022.fasta"
de_file <- "results/deseq2/DESeq2_all_results.tsv"

cds <- readSet(cds_file)

de <- read.delim(de_file)

codon_table <- codonTable(cds)

scuo_values <- SCUO(codon_table)

scuo_df <- data.frame(
  gene_id = names(scuo_values),
  SCUO = as.numeric(scuo_values)
)

merged <- de %>%
  left_join(scuo_df, by = "gene_id") %>%
  mutate(
    translation_group = case_when(
      padj < 0.01 & log2FoldChange >= 1 ~ "Upregulated",
      padj < 0.01 & log2FoldChange <= -1 ~ "Downregulated",
      TRUE ~ "Not_significant"
    )
  )

wilcox_result <- wilcox.test(
  SCUO ~ translation_group,
  data = merged %>%
    filter(translation_group %in% c("Upregulated", "Downregulated"))
)

dir.create("results/codon_usage", recursive = TRUE, showWarnings = FALSE)

write.table(
  merged,
  file = "results/codon_usage/SCUO_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    test = "Wilcoxon rank-sum test",
    p_value = wilcox_result$p.value
  ),
  file = "results/codon_usage/codon_usage_wilcoxon_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
