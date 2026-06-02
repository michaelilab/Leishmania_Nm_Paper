#!/usr/bin/env Rscript

library(tidyverse)

# ------------------------------------------------------------
# 04_make_count_matrix.R
#
# Combines individual HTSeq-count output files into one count
# matrix for DESeq2.
#
# Input:
#   results/counts/sample1.counts.txt
#   results/counts/sample2.counts.txt
#   ...
#
# Output:
#   data/counts/htseq_counts_matrix.tsv
# ------------------------------------------------------------

counts_dir <- "results/counts"
output_file <- "data/counts/htseq_counts_matrix.tsv"

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Find HTSeq-count files
count_files <- list.files(
  counts_dir,
  pattern = "\\.counts\\.txt$",
  full.names = TRUE
)

if (length(count_files) == 0) {
  stop("No count files found in: ", counts_dir)
}

# Function to read one HTSeq-count file
read_htseq_counts <- function(file) {
  sample_name <- basename(file) %>%
    str_replace("\\.counts\\.txt$", "")

  counts <- read.delim(
    file,
    header = FALSE,
    col.names = c("gene_id", sample_name),
    stringsAsFactors = FALSE
  )

  # Remove HTSeq summary rows, for example:
  # __no_feature, __ambiguous, __too_low_aQual, etc.
  counts <- counts %>%
    filter(!str_starts(gene_id, "__"))

  return(counts)
}

# Read and merge all count files
count_list <- lapply(count_files, read_htseq_counts)

count_matrix <- reduce(count_list, full_join, by = "gene_id")

# Replace missing counts with 0
count_matrix[is.na(count_matrix)] <- 0

# Make sure count columns are integers
count_matrix <- count_matrix %>%
  mutate(across(-gene_id, as.integer))

# Save count matrix
write.table(
  count_matrix,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Count matrix written to: ", output_file)
message("Number of genes: ", nrow(count_matrix))
message("Number of samples: ", ncol(count_matrix) - 1)
