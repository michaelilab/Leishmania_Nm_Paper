#!/usr/bin/env Rscript

##########################################################################
# 03_calculate_ribometh_scores.R
#
# This script calculates RiboMethSeq scores A, B, and C from preprocessed
# 5' and 3' count files.
#
# Expected input files:
#   <library>.sorted.init
#   <library>.sorted.3p
#
# Each file must contain the read counts in a selected column, default = 3.
#
# Output:
#   One CSV file per library containing:
#   5p, 3p, coverage, Sa, Sb, Sc, bp, and optional annotation columns
#
# Example:
#   Rscript calculate_ribometh_scores.R \
#     --input_dir Input \
#     --output_dir Output \
#     --fasta DB/LM_rRNA_whole_fasta.fa \
#     --annotation DB/LM_known_rRNA_Nms.txt \
#     --win_size 6 \
#     --count_col 3
##########################################################################

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  match_idx <- match(flag, args)
  if (!is.na(match_idx) && match_idx < length(args)) {
    return(args[match_idx + 1])
  }
  return(default)
}

input_dir  <- get_arg("--input_dir")
output_dir <- get_arg("--output_dir")
fasta_file <- get_arg("--fasta")
annot_file <- get_arg("--annotation", default = NULL)

win_size  <- as.numeric(get_arg("--win_size", default = 6))
count_col <- as.numeric(get_arg("--count_col", default = 3))

if (is.null(input_dir) || is.null(output_dir) || is.null(fasta_file)) {
  stop(
    paste(
      "Missing required arguments.\n\n",
      "Required:\n",
      "  --input_dir   Directory containing *.sorted.init and *.sorted.3p files\n",
      "  --output_dir  Directory for output CSV files\n",
      "  --fasta       FASTA file of the rRNA reference\n\n",
      "Optional:\n",
      "  --annotation  Annotation table matching rRNA positions\n",
      "  --win_size    Window size, default = 6\n",
      "  --count_col   Count column in input files, default = 3\n\n",
      "Example:\n",
      "  Rscript calculate_ribometh_scores.R --input_dir Input --output_dir Output ",
      "--fasta DB/LM_rRNA_whole_fasta.fa --annotation DB/LM_known_rRNA_Nms.txt\n"
    )
  )
}

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}

if (!file.exists(fasta_file)) {
  stop("FASTA file does not exist: ", fasta_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (win_size %% 2 != 0) {
  stop("win_size must be an even number. Current value: ", win_size)
}

if (!requireNamespace("seqinr", quietly = TRUE)) {
  install.packages("seqinr", repos = "https://cloud.r-project.org")
}

library(seqinr)

message("Reading FASTA file: ", fasta_file)

myfasta <- read.fasta(
  file = fasta_file,
  as.string = FALSE,
  set.attributes = FALSE
)

if (length(myfasta) < 1) {
  stop("No sequences found in FASTA file.")
}

# Use the first FASTA sequence as the rRNA reference
bp <- unlist(unname(myfasta[[1]]))
rRNA_length <- length(bp)

message("Reference length detected from FASTA: ", rRNA_length, " bp")

annot <- NULL

if (!is.null(annot_file)) {
  if (!file.exists(annot_file)) {
    stop("Annotation file does not exist: ", annot_file)
  }
  
  message("Reading annotation file: ", annot_file)
  
  annot <- read.table(
    annot_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (nrow(annot) != rRNA_length) {
    warning(
      "Annotation rows do not match FASTA length.\n",
      "FASTA length: ", rRNA_length, "\n",
      "Annotation rows: ", nrow(annot), "\n",
      "The script will continue, but the annotation may not align correctly."
    )
  }
}

init_files <- list.files(
  path = input_dir,
  pattern = "\\.sorted\\.init$",
  full.names = FALSE
)

if (length(init_files) == 0) {
  stop("No *.sorted.init files found in input directory: ", input_dir)
}

libraries <- sub("\\.sorted\\.init$", "", init_files)

message("Found ", length(libraries), " libraries.")

calculate_scores <- function(myinit, myends, win_size) {
  
  mycov <- myinit + myends
  mylength <- length(mycov)
  
  Sc <- rep(NA, mylength)
  Sa <- rep(NA, mylength)
  Sb <- rep(NA, mylength)
  
  W <- (1 + (1 - 0.1 * win_size)) * win_size / 2
  
  start_i <- win_size + 1
  end_i <- mylength - win_size - 1
  
  if (end_i <= start_i) {
    warning("Reference is too short for the selected window size.")
    return(data.frame(Sa = Sa, Sb = Sb, Sc = Sc))
  }
  
  for (i in start_i:end_i) {
    
    # Score A
    M_l <- mean(mycov[(i - win_size / 2):(i - 1)])
    S_l <- sd(mycov[(i - win_size / 2):(i - 1)])
    
    M_r <- mean(mycov[(i + 1):(i + win_size / 2)])
    S_r <- sd(mycov[(i + 1):(i + win_size / 2)])
    
    Sa[i] <- max(
      0,
      1 - (2 * mycov[i] + 1) /
        (0.5 * abs(M_l - S_l) + mycov[i] + 0.5 * abs(M_r - S_r) + 1)
    )
    
    # Scores B and C
    S1 <- 0
    for (j in 1:win_size) {
      S1 <- S1 + (1 - 0.1 * (j - 1)) * mycov[i - j]
    }
    S1 <- S1 / W
    
    S2 <- 0
    for (j in 1:win_size) {
      S2 <- S2 + (1 - 0.1 * (j - 1)) * mycov[i + j]
    }
    S2 <- S2 / W
    
    Sc[i] <- max(
      0,
      1 - 2 * mycov[i] / (S1 + S2)
    )
    
    Sb[i] <- abs(
      (mycov[i] - 0.5 * (S1 + S2)) /
        (mycov[i] + 1)
    )
  }
  
  return(data.frame(Sa = Sa, Sb = Sb, Sc = Sc))
}

pad_or_trim_counts <- function(x, target_length) {
  x <- as.numeric(x)
  
  if (length(x) < target_length) {
    x <- c(x, rep(0, target_length - length(x)))
  }
  
  if (length(x) > target_length) {
    warning("Count vector is longer than FASTA reference. Trimming to reference length.")
    x <- x[1:target_length]
  }
  
  return(x)
}

for (thisLib in libraries) {
  
  message("Processing library: ", thisLib)
  
  init_file <- file.path(input_dir, paste0(thisLib, ".sorted.init"))
  ends_file <- file.path(input_dir, paste0(thisLib, ".sorted.3p"))
  
  if (!file.exists(ends_file)) {
    warning("Missing 3p file for library ", thisLib, ": ", ends_file)
    next
  }
  
  init_table <- read.table(
    init_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  ends_table <- read.table(
    ends_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (ncol(init_table) < count_col) {
    warning("Skipping ", thisLib, ": init file has fewer columns than count_col.")
    next
  }
  
  if (ncol(ends_table) < count_col) {
    warning("Skipping ", thisLib, ": 3p file has fewer columns than count_col.")
    next
  }
  
  pre_myinit <- init_table[, count_col]
  pre_myends <- ends_table[, count_col]
  
  # Shift reads by 1 bp, following the original script logic
  myinit <- c(0, pre_myinit)
  myends <- c(pre_myends[-1], 0)
  
  myinit <- pad_or_trim_counts(myinit, rRNA_length)
  myends <- pad_or_trim_counts(myends, rRNA_length)
  
  score_df <- calculate_scores(
    myinit = myinit,
    myends = myends,
    win_size = win_size
  )
  
  mydf <- data.frame(
    `5p` = myinit,
    `3p` = myends,
    cov = myinit + myends,
    Sa = score_df$Sa,
    Sb = score_df$Sb,
    Sc = score_df$Sc,
    bp = bp,
    check.names = FALSE
  )
  
  if (!is.null(annot)) {
    mydf <- cbind(mydf, annot)
  }
  
  output_file <- file.path(output_dir, paste0(thisLib, ".csv"))
  
  write.csv(
    x = mydf,
    file = output_file,
    row.names = FALSE
  )
  
  message("Wrote: ", output_file)
}

message("Done.")
