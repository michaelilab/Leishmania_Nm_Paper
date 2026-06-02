# Poly(A) RNA-seq analysis in *Leishmania major*

This repository contains scripts used for poly(A) RNA-seq analysis in *Leishmania major* Friedlin v9.

The workflow includes:

1. Alignment of poly(A) RNA-seq reads to the *L. major* Friedlin v9 genome using SMALT.
2. Gene-level read counting using HTSeq-count.
3. Differential expression analysis using DESeq2.
4. UTR length analysis based on dominant 5′ spliced-leader and 3′ poly(A) sites.
5. Codon usage analysis using the R package `coRdon`.
6. Visualization of results using `ggplot2`.

## Software

The analysis was performed using:

- SMALT v0.7.5
- HTSeq-count
- R v4.2.2
- DESeq2 v1.38.3
- coRdon v1.6.0
- ggplot2

## Repository structure

```text
scripts/      Analysis scripts
data/         Input metadata, reference files, and count matrices
results/      Output tables and figures
environment/  Software versions and reproducibility files
docs/         Additional method descriptions
