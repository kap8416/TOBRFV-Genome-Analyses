# TOBRFV-Genome-Analyses
## Introduction

This repository contains a comprehensive bioinformatics workflow for analyzing the genome of the Mexican isolate of Tomato Brown Rugose Fruit Virus (ToBRFV), a highly virulent virus affecting tomato and pepper plants. The workflow includes RNA sequencing, quality control, genome assembly, and a series of analytical steps aimed at understanding the genetic diversity, phylogeny, and geographical distribution of ToBRFV isolates. The complete pipeline and code are provided for replicability and further exploration.

We then implemented a bioinformatics workflow for the Mexican ToBRFV isolate. The workflow includes genome retrieval, metadata extraction, genetic diversity calculations, phylogenetic modeling, and visualization. 

**Figure 1 provides an overview of the analysis steps. All scripts used are accessible in this GitHub repository**


![output-9](https://github.com/user-attachments/assets/0dfde311-73c6-4986-a50e-4aac691a96a6)


##The analysis pipeline includes:

**Sequence Retrieval**: Retrieval of full genome sequences of 100 ToBRFV isolates from GenBank.
**Metadata Extraction**: Extraction of geographical and host information for each isolate.
**Genetic Diversity Calculations**: Pairwise identity, haplotype diversity, and polymorphic site calculations to quantify genetic diversity.
**Nucleotide Diversity (π) and Average Differences (k)**: Measurement of nucleotide diversity (π) and average nucleotide differences (k) across isolates.
**Alignment and Phylogenetic Modeling**: Alignment with ClustalW, followed by IQ-TREE 2 phylogenetic modeling with 1000 bootstrap replicates.
**Phylogenetic Network and Tree Visualization**: Visualization of relationships among isolates using Neighbor-Net and ggtree.
**Geographical Distribution Mapping**: Mapping isolate regions using ggplot2 and sf to show distribution patterns.
**Genomic Variation Analysis**: Identification of nucleotide differences and mutation hotspots.
**Alignment Scores Analysis**: Calculation of alignment scores and scatter plot visualization by host and region.


##Phylogenetic Analysis

To assess genetic relationships among isolates:
101 ToBRFV genome sequences were aligned using ClustalW (in bash, with FASTA input). IQ-TREE2 was used for phylogenetic analysis, with Maximum Likelihood (ML) method and ModelFinder for optimal model selection. 1000 ultrafast bootstrap replicates were applied for clade support. Visualization was done in R using ggtree, with pairwise sequence identity matrices generated using ape and visualized as heatmaps with pheatmap. A Neighbor-Net phylogenetic network was created to examine relationships among isolates.
Key diversity metrics (haplotype diversity, haplotype count, nucleotide diversity, and average nucleotide differences) were calculated using pegas and ape, and compared visually with fmsb radar plots.

##Genomic Variation and Distribution Analysis of the Mexican ToBRFV Isolate
###Nucleotide Variation Analysis
Using the Biostrings package in RStudio, we analyzed nucleotide differences between the Mexican ToBRFV isolate and 100 reference genomes. Normalization was applied to ensure consistent measurement, and statistical summaries highlighted mutation rate variation. Visualization was done in ggplot2, revealing unique mutations and conserved regions in the Mexican isolate with potential functional or evolutionary relevance.
###Geographical Distribution of ToBRFV Isolates
Geographical metadata for each isolate was visualized using ggplot2 and sf in RStudio. A world map annotated with isolate counts highlighted the distribution and clustering of isolates genetically similar to the Mexican genome, illustrating global dissemination patterns.
###Alignment Scores by Host and Geographical Distribution
To examine alignment score variation across hosts and regions, alignment scores were calculated with msa and combined with host and geographical data. ggplot2 was used to generate scatter plots, with points colored by host and annotated by region.
Host-Specific Nucleotide Differences
###Nucleotide variations associated with different hosts were identified using Biostrings, linking genetic variation with host data. ggplot2 visualizations highlighted patterns along the genome, suggesting potential host-specific adaptations and evolutionary pressures.

**Requirements**

Ensure the following R packages are installed:

ggtree
treeio
ape
viridis
ggplot2
msa
Biostrings
stringr
dplyr
tidyr
pheatmap
phangorn
fmsb
sf

Installation Command:
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("ggtree", "treeio", "msa", "Biostrings"))
install.packages(c("ape", "viridis", "ggplot2", "stringr", "dplyr", "tidyr", "pheatmap", "phangorn", "fmsb", "sf"))
Usage
Phylogenetic Tree Generation
Place the Newick tree file (TOBRFV.nwk) in the working directory.
Run generate_phylogenetic_tree.R to generate the tree visualizations.
Genomic Differences Visualization
Place the FASTA file (lcl|Query_2456347.aln) in the working directory.
Run visualize_genomic_differences.R for genome-wide difference visualizations.
Genomic Differences by Host
Place the FASTA file (lcl|Query_2456347.aln) in the working directory.
Run GenomicDifferencesbyHost.R to analyze differences grouped by host.
Scripts
generate_phylogenetic_tree.R: Generates a phylogenetic tree with labeled clades, bootstrap values, and saves as PDF and PNG.
visualize_genomic_differences.R: Visualizes genomic differences across the ToBRFV genome.
GenomicDifferencesbyHost.R: Visualizes genomic differences categorized by host.
Summary
This repository offers a comprehensive workflow for analyzing ToBRFV’s genetic diversity, phylogenetic relationships, and geographical distribution patterns. It enables reproducible research and provides insights into the virus's transmission pathways, host adaptations, and global impact.
This README provides a structured and comprehensive overview of the ToBRFV Genome Analysis workflow, guiding users through the setup, usage, and functionalities of the scripts in the repository.

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
install.packages("tidyr")

Usage

Phylogenetic Tree Generation
Place your Newick formatted tree file (TOBRFV.nwk) in the same directory as the script.
Run the script generate_phylogenetic_tree.R.
The script will generate two files:
Phylogenetic_Tree_ToBRFV.pdf
Phylogenetic_Tree_ToBRFV.png

Genomic Differences Visualization
Genomic Differences Visualization
Place your FASTA file (lcl|Query_2456347.aln) in the same directory as the script.
Run the script visualize_genomic_differences.R.

Genomic Differences Visualization by host
Place your FASTA file (lcl|Query_2456347.aln) in the same directory as the script.
Run the script GenomicDifferencesbyHost.R.
The script will generate the file:
ToBRFV_Genomic_Differences_By_Host.png
Scripts

generate_phylogenetic_tree.R
visualize_genomic_differences.R
GenomicDifferencesbyHost.R

Summary

generate_phylogenetic_tree.R: This script generates a phylogenetic tree with highlighted species and saves the visualization as PDF and PNG files.
visualize_genomic_differences.R: This script visualizes genomic differences of ToBRFV grouped by host and saves the visualization as a PNG file.
GenomicDifferencesbyHost.R: This script visualizes genomic differences of ToBRFV grouped by host and saves the visualization as a PNG file.README.md: This file explains how to install the necessary packages, how to use the scripts, and provides the script code for reference.
