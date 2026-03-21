# Basic R Codes for Bioinformatics Analysis

This repository contains R scripts for:

- Differential expression analysis (DEG)
- Volcano plot visualization
- GSEA analysis
- Protein expression visualization

## Scripts

- `Volcano plot.R` – generate volcano plots
- `rank_plot_manual.R` – ranking visualization
- `protein_bubble_*.R` – protein expression plots

## Usage

Run scripts in R:

```r
source("Volcano plot.R")

## Reproducible figure check

This repository uses GitHub Actions to:
- install R and plotting packages
- check R scripts
- render an example figure automatically
