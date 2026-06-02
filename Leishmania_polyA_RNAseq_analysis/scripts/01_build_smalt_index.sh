#!/bin/bash
set -euo pipefail

GENOME="data/reference/Lmajor_Friedlin_v9_genome.fa"
INDEX_PREFIX="data/reference/Lmajor_Friedlin_v9_smalt"

smalt index "$INDEX_PREFIX" "$GENOME"
