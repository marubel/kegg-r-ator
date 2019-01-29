#!/usr/bin/env Rscript

# Instead of this, we can use R() directly inside Snakemake.
source("kegg-r-ator.R")
KEGG_to_ko(snakemake@input[["sample"]], snakemake@input[["ko_mapping_file_fp"]], snakemake@output[["tsv"]])
