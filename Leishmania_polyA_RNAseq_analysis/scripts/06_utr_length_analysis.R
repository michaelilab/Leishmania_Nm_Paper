#!/usr/bin/env Rscript

# ------------------------------------------------------------
# 06_utr_length_analysis.R
#
# Analyze 3' and 5' UTR lengths for different mRNA groups.
#
# Required input files:
#   1. 3' UTR table: gene_id length
#   2. 5' UTR table: gene_id length
#   3. Increased gene list
#   4. Decreased gene list
#
# Usage:
#   Rscript scripts/06_utr_length_analysis.R \
#     data/reference/LM_3pUTR.txt \
#     data/reference/LM_5pUTR.txt \
#     data/gene_lists/increased.txt \
#     data/gene_lists/decreased.txt \
#     results/utr
-----------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    "Usage:\n",
    "Rscript scripts/06_utr_length_analysis.R ",
    "<three_prime_utr_file> <five_prime_utr_file> ",
    "<increased_gene_list> <decreased_gene_list> <output_dir>\n\n",
    "Example:\n",
    "Rscript scripts/06_utr_length_analysis.R ",
    "data/reference/LM_3pUTR.txt ",
    "data/reference/LM_5pUTR.txt ",
    "data/gene_lists/increased.txt ",
    "data/gene_lists/decreased.txt ",
    "results/utr"
  )
}

three_prime_utr_file <- args[1]
five_prime_utr_file <- args[2]
increased_file <- args[3]
decreased_file <- args[4]
output_dir <- args[5]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Read gene lists
# ------------------------------------------------------------

# Original note:
# 11/04/2024 switched down/up
# increased.txt is treated as "down"
# decreased.txt is treated as "up"

down_genes <- scan(decreased_file, what = character(), quiet = TRUE)
up_genes <- scan(increased_file, what = character(), quiet = TRUE)

# ------------------------------------------------------------
# Function for UTR analysis
# ------------------------------------------------------------

analyze_utr_lengths <- function(
    utr_file,
    utr_label,
    y_limit,
    signif_y_positions,
    output_prefix
) {

  message("Analyzing ", utr_label, " UTR lengths")

  utr_table <- read.delim(
    utr_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(utr_table) < 2) {
    stop("UTR file must contain at least two columns: gene_id and length")
  }

  utr_table <- utr_table %>%
    select(1, 2) %>%
    rename(
      Gene = 1,
      Length = 2
    ) %>%
    mutate(
      Length = as.numeric(Length),
      status = case_when(
        Gene %in% up_genes ~ "up",
        Gene %in% down_genes ~ "down",
        TRUE ~ "other"
      ),
      status = factor(status, levels = c("down", "up", "other"))
    ) %>%
    filter(!is.na(Length))

  # Save annotated UTR table
  write.table(
    utr_table,
    file = file.path(output_dir, "tables", paste0(output_prefix, "_UTR_lengths_with_status.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Summary statistics
  summary_stats <- utr_table %>%
    group_by(status) %>%
    summarise(
      n = n(),
      mean_length = mean(Length, na.rm = TRUE),
      median_length = median(Length, na.rm = TRUE),
      sd_length = sd(Length, na.rm = TRUE),
      min_length = min(Length, na.rm = TRUE),
      max_length = max(Length, na.rm = TRUE),
      .groups = "drop"
    )

  write.table(
    summary_stats,
    file = file.path(output_dir, "tables", paste0(output_prefix, "_UTR_length_summary.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Statistical tests
  pairwise_wilcox <- pairwise.wilcox.test(
    utr_table$Length,
    utr_table$status,
    p.adjust.method = "bonferroni",
    exact = FALSE
  )

  kruskal_result <- kruskal.test(Length ~ status, data = utr_table)

  compare_means_result <- compare_means(
    Length ~ status,
    data = utr_table,
    method = "wilcox.test"
  )

  write.csv(
    compare_means_result,
    file = file.path(output_dir, "tables", paste0(output_prefix, "_UTR_compare_means.csv")),
    row.names = FALSE
  )

  kruskal_table <- data.frame(
    test = "Kruskal-Wallis",
    statistic = as.numeric(kruskal_result$statistic),
    df = as.numeric(kruskal_result$parameter),
    p_value = kruskal_result$p.value
  )

  write.table(
    kruskal_table,
    file = file.path(output_dir, "tables", paste0(output_prefix, "_UTR_kruskal_test.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  pairwise_table <- as.data.frame(as.table(pairwise_wilcox$p.value)) %>%
    filter(!is.na(Freq)) %>%
    rename(
      group1 = Var1,
      group2 = Var2,
      p_adj_bonferroni = Freq
    )

  write.table(
    pairwise_table,
    file = file.path(output_dir, "tables", paste0(output_prefix, "_UTR_pairwise_wilcox_bonferroni.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Boxplot
  boxplot_fig <- ggplot(
    utr_table,
    aes(x = status, y = Length, fill = status)
  ) +
    stat_boxplot(
      geom = "errorbar",
      position = "dodge2",
      width = 0.3
    ) +
    geom_boxplot(
      width = 0.5,
      outlier.shape = 1,
      position = "dodge2"
    ) +
    coord_cartesian(ylim = c(1, y_limit)) +
    ggtitle(paste0(utr_label, " UTR Lengths")) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    geom_signif(
      test = "wilcox.test",
      comparisons = list(
        c("up", "other"),
        c("down", "other"),
        c("down", "up")
      ),
      map_signif_level = TRUE,
      y_position = signif_y_positions
    ) +
    xlab("mRNA group") +
    ylab(paste0(utr_label, " UTR Length"))

  ggsave(
    filename = file.path(output_dir, "figures", paste0(output_prefix, "_UTR_length_boxplot.pdf")),
    plot = boxplot_fig,
    width = 4,
    height = 4
  )

  # Violin plot
  violin_fig <- ggplot(
    utr_table,
    aes(x = status, y = Length, fill = status)
  ) +
    geom_violin() +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    coord_cartesian(ylim = c(1, y_limit)) +
    theme_classic() +
    theme(
      legend.position = "none"
    ) +
    xlab("mRNA group") +
    ylab(paste0(utr_label, " UTR Length"))

  ggsave(
    filename = file.path(output_dir, "figures", paste0(output_prefix, "_UTR_length_violinplot.pdf")),
    plot = violin_fig,
    width = 4,
    height = 4
  )

  # Density plot
  density_fig <- ggplot(
    utr_table,
    aes(x = Length, color = status)
  ) +
    geom_density(linewidth = 1) +
    coord_cartesian(xlim = c(0, y_limit)) +
    theme_classic() +
    xlab(paste0(utr_label, " UTR Length")) +
    ylab("Density")

  ggsave(
    filename = file.path(output_dir, "figures", paste0(output_prefix, "_UTR_length_density.pdf")),
    plot = density_fig,
    width = 5,
    height = 4
  )

  message("Finished ", utr_label, " UTR analysis")

  return(
    list(
      utr_table = utr_table,
      summary_stats = summary_stats,
      compare_means = compare_means_result,
      kruskal = kruskal_table,
      pairwise_wilcox = pairwise_table
    )
  )
}

# ------------------------------------------------------------
# Run 3' UTR analysis
# ------------------------------------------------------------

three_prime_results <- analyze_utr_lengths(
  utr_file = three_prime_utr_file,
  utr_label = "3'",
  y_limit = 2500,
  signif_y_positions = c(1500, 1500, 2000),
  output_prefix = "3p"
)

# ------------------------------------------------------------
# Run 5' UTR analysis
# ------------------------------------------------------------

five_prime_results <- analyze_utr_lengths(
  utr_file = five_prime_utr_file,
  utr_label = "5'",
  y_limit = 9500,
  signif_y_positions = c(6000, 6000, 6500),
  output_prefix = "5p"
)

message("All UTR length analyses completed.")
