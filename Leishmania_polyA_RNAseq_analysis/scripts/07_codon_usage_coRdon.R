#!/usr/bin/env Rscript

library(coRdon)
library(tidyverse)
library(Biostrings)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggforce)

# =========================
# Input files
# =========================

cds_file <- "data/reference/GenesByGeneModelChars_LM_11July2022.fasta"
de_file  <- "results/deseq2/DESeq2_all_results.tsv"

# =========================
# Output directory
# =========================

outdir <- "results/codon_usage"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Read input data
# =========================

cds <- readSet(cds_file)
de  <- read.delim(de_file)

# Make sure DESeq2 table has gene_id column
# If your gene IDs are rownames instead, uncomment this:
# de <- de %>% rownames_to_column("gene_id")

# =========================
# Codon usage table
# =========================

codon_table <- codonTable(cds)

# =========================
# SCUO analysis
# =========================

scuo_values <- SCUO(codon_table)

scuo_df <- data.frame(
  gene_id = names(scuo_values),
  SCUO = as.numeric(scuo_values)
)

merged <- de %>%
  left_join(scuo_df, by = "gene_id") %>%
  mutate(
    translation_group = case_when(
      padj < 0.01 & log2FoldChange >= 1  ~ "Upregulated",
      padj < 0.01 & log2FoldChange <= -1 ~ "Downregulated",
      TRUE ~ "Not_significant"
    )
  )

write.table(
  merged,
  file = file.path(outdir, "SCUO_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Wilcoxon test for SCUO: Upregulated vs Downregulated
wilcox_data <- merged %>%
  filter(translation_group %in% c("Upregulated", "Downregulated")) %>%
  filter(!is.na(SCUO))

if (length(unique(wilcox_data$translation_group)) == 2) {
  wilcox_result <- wilcox.test(
    SCUO ~ translation_group,
    data = wilcox_data
  )

  wilcox_out <- data.frame(
    test = "Wilcoxon rank-sum test",
    comparison = "SCUO: Upregulated vs Downregulated",
    p_value = wilcox_result$p.value
  )
} else {
  wilcox_out <- data.frame(
    test = "Wilcoxon rank-sum test",
    comparison = "SCUO: Upregulated vs Downregulated",
    p_value = NA
  )
}

write.table(
  wilcox_out,
  file = file.path(outdir, "codon_usage_wilcoxon_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================
# Codon frequency plot
# =========================

# Convert codon table to data frame
cu_table <- as.data.frame(codon_table)

# Add gene ID column
cu_table$gene_id <- rownames(cu_table)

# Convert from wide codon-frequency table to long format
# Assumes codon_table columns are codons and values are frequencies/counts.
# If coRdon gives extra non-codon columns, this keeps only standard codon-like columns.
codon_cols <- grep("^[ACGTU]{3}$", colnames(cu_table), value = TRUE)

codon_long <- cu_table %>%
  select(gene_id, all_of(codon_cols)) %>%
  pivot_longer(
    cols = all_of(codon_cols),
    names_to = "codon",
    values_to = "freq"
  )

# Add DESeq2 expression group to each gene/codon row
codon_long <- codon_long %>%
  left_join(
    merged %>% select(gene_id, translation_group),
    by = "gene_id"
  ) %>%
  mutate(
    status = case_when(
      translation_group == "Upregulated" ~ "up",
      translation_group == "Downregulated" ~ "down",
      TRUE ~ "other"
    ),
    status = factor(status, levels = c("down", "other", "up"))
  ) %>%
  filter(!is.na(freq))

# Save full codon-frequency table
write.table(
  codon_long,
  file = file.path(outdir, "codon_frequency_long_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================
# Statistical comparisons per codon
# =========================

my_comparisons <- list(
  c("down", "other"),
  c("up", "other"),
  c("down", "up")
)

codon_stats <- compare_means(
  freq ~ status,
  data = codon_long,
  comparisons = my_comparisons,
  group.by = "codon",
  method = "wilcox.test"
)

write.csv(
  codon_stats,
  file = file.path(outdir, "codon_frequency_group_stats.csv"),
  row.names = FALSE
)

# =========================
# Plot: codon frequency by group
# =========================

gg <- ggplot(codon_long, aes(x = status, y = freq, fill = status)) +
  geom_violin(trim = FALSE, na.rm = TRUE) +
  geom_boxplot(
    width = 0.07,
    show.legend = FALSE,
    outlier.size = 0.4,
    na.rm = TRUE
  ) +
  facet_wrap_paginate(
    ~ codon,
    ncol = 3,
    nrow = 3,
    page = 1,
    scales = "free_y"
  ) +
  geom_signif(
    test = "wilcox.test",
    comparisons = my_comparisons,
    map_signif_level = TRUE,
    step_increase = 0.06,
    size = 0.5,
    textsize = 2
  ) +
  ylab("Codon frequency") +
  xlab("Expression group") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

n_pages_total <- n_pages(gg)

pdf(
  file = file.path(outdir, "codon_frequency_by_DE_group_violinplot.pdf"),
  paper = "a4",
  width = 210 / 25.4,
  height = 297 / 25.4
)

for (i in seq_len(n_pages_total)) {
  print(
    ggplot(codon_long, aes(x = status, y = freq, fill = status)) +
      geom_violin(trim = FALSE, na.rm = TRUE) +
      geom_boxplot(
        width = 0.07,
        show.legend = FALSE,
        outlier.size = 0.4,
        na.rm = TRUE
      ) +
      facet_wrap_paginate(
        ~ codon,
        ncol = 3,
        nrow = 3,
        page = i,
        scales = "free_y"
      ) +
      geom_signif(
        test = "wilcox.test",
        comparisons = my_comparisons,
        map_signif_level = TRUE,
        step_increase = 0.06,
        size = 0.5,
        textsize = 2
      ) +
      ylab("Codon frequency") +
      xlab("Expression group") +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  )
}

dev.off()

# =========================
# Optional summary table
# =========================

codon_summary <- codon_long %>%
  group_by(codon, status) %>%
  summarise(
    Mean = mean(freq, na.rm = TRUE),
    Median = median(freq, na.rm = TRUE),
    SD = sd(freq, na.rm = TRUE),
    Min = min(freq, na.rm = TRUE),
    Max = max(freq, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

write.csv(
  codon_summary,
  file = file.path(outdir, "codon_frequency_by_group_summary.csv"),
  row.names = FALSE
)

message("Done.")
message("Results written to: ", outdir)
