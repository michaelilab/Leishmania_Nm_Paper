#!/usr/bin/env bash
set -euo pipefail

# Usage:
# bash scripts/03_htseq_count.sh <bam_dir> <annotation_gtf> <output_dir>
#
# Example:
# bash scripts/03_htseq_count.sh \
#   results/alignments \
#   data/reference/LM_mRNA_v9.gtf \
#   results/counts

bam_dir="$1"
annot="$2"
out_dir="$3"

mkdir -p "$out_dir"

for bam_file in "$bam_dir"/*sorted.bam
do
    if [[ ! -e "$bam_file" ]]; then
        echo "No sorted BAM files found in: $bam_dir"
        exit 1
    fi

    base=$(basename "$bam_file" | cut -f1,2 -d_)

    echo "Processing: $base"

    samtools index "$bam_file"

    htseq-count \
        -f bam \
        -s no \
        -r pos \
        "$bam_file" \
        "$annot" \
        > "$out_dir/${base}.htseqcount" &
done

wait

echo "HTSeq-count finished for all BAM files."
