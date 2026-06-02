#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# 02_align_reads_smalt.sh
#
# Align paired-end poly(A) RNA-seq reads to the L. major genome
# using SMALT, then create sorted and indexed BAM files.
#
# Usage:
#   bash scripts/02_align_reads_smalt.sh <fastq_dir> <smalt_genome_index> <output_dir>
#
# Example:
#   bash scripts/02_align_reads_smalt.sh \
#     data/raw \
#     data/reference/Lmajor_Friedlin_v9_smalt \
#     results/alignments
#
# Notes:
#   - Input R1 FASTQ files should match: *R1.fastq.gz
#   - R2 files are expected to have the same name with R1 replaced by R2.
#   - SMALT index should already exist.
# ------------------------------------------------------------

fastq_dir="$1"
genome_index="$2"
out_dir="$3"

mkdir -p "$out_dir"

for fastq_r1 in "$fastq_dir"/*R1.fastq.gz
do
    if [[ ! -e "$fastq_r1" ]]; then
        echo "No R1 FASTQ files found in: $fastq_dir"
        exit 1
    fi

    sample=$(basename "$fastq_r1" | cut -f1 -d.)

    fastq_r2="${fastq_r1/_R1/_R2}"

    if [[ ! -f "$fastq_r2" ]]; then
        echo "ERROR: Matching R2 file not found for:"
        echo "$fastq_r1"
        echo "Expected:"
        echo "$fastq_r2"
        exit 1
    fi

    echo "Processing sample: $sample"
    echo "R1: $fastq_r1"
    echo "R2: $fastq_r2"

    sam_file="$out_dir/${sample}_vs_genome.sam"
    sorted_bam="$out_dir/${sample}_vs_genome_sorted.bam"
   
    smalt map \
        "$genome_index" \
        "$fastq_r1" \
        "$fastq_r2" \
        > "$sam_file"

    
    samtools view \
        -bS \
        "$sam_file" \
        | samtools sort -o "$sorted_bam"

    samtools index "$sorted_bam"

    rm "$sam_file"

    echo "Finished sample: $sample"
done

echo "All SMALT alignments completed."
