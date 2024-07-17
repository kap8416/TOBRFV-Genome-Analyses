# TOBRFV-Genome-Analyses

## Introduction

Tomato Brown Rugose Fruit Virus (ToBRFV) is a highly virulent virus that affects tomato and pepper plants, causing significant agricultural losses. The ToBRFV Mexican isolate is a specific strain of the virus identified in Mexico, known for its impact on local tomato crops. Understanding the phylogenetic relationships of this isolate with other species is crucial for developing effective management and control strategies.

## Requirements

Make sure you have the following packages installed:

- `ggtree`
- `treeio`
- `ape`
- `viridis`
- `ggplot2`
- `msa`
- `Biostrings`
- `stringr`
- `dplyr`
- `tidyr`

You can install them using the following commands:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ggtree")
BiocManager::install("treeio")
install.packages("ape")
install.packages("viridis")
install.packages("ggplot2")
install.packages("msa")
BiocManager::install("Biostrings")
install.packages("stringr")
install.packages("dplyr")
install.packages("tidyr"

Usage

Phylogenetic Tree Generation
Place your Newick formatted tree file (TOBRFV.nwk) in the same directory as the script.
Run the script generate_phylogenetic_tree.R.
The script will generate two files:
Phylogenetic_Tree_ToBRFV.pdf
Phylogenetic_Tree_ToBRFV.png
Genomic Differences Visualization
Place your FASTA file (lcl|Query_2456347.aln) in the same directory as the script.
Run the script visualize_genomic_differences.R.
The script will generate the file:
ToBRFV_Genomic_Differences_By_Host.png

Scripts

generate_phylogenetic_tree.R
visualize_genomic_differences.R

Summary

generate_phylogenetic_tree.R: This script generates a phylogenetic tree with highlighted species and saves the visualization as PDF and PNG files.
visualize_genomic_differences.R: This script visualizes genomic differences of ToBRFV grouped by host and saves the visualization as a PNG file.
README.md: This file explains how to install the necessary packages, how to use the scripts, and provides the script code for reference.
