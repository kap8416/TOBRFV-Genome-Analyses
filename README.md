TOBRFV-Genome-Analyses

Introduction

This repository contains a comprehensive bioinformatics workflow for analyzing the genome of the Mexican isolate of Tomato Brown Rugose Fruit Virus (ToBRFV), a highly virulent virus that affects tomato and pepper plants. The workflow includes sequence retrieval, quality control, genome assembly, and a series of analyses aimed at understanding the genetic diversity, phylogeny, and geographical distribution of ToBRFV isolates. The complete pipeline and code are provided for replicability and further exploration.
This analysis pipeline includes genome retrieval, metadata extraction, genetic diversity calculations, phylogenetic modeling, and visualization to explore the relationships and variations within ToBRFV isolates.
Figure 1 provides an overview of the analysis steps, with all scripts accessible in this repository.


Analysis Pipeline Overview

Sequence Retrieval: Retrieve full genome sequences of 100 ToBRFV isolates from GenBank, including the Mexican isolate.
Metadata Extraction: Extract geographical and host information for each isolate.
Genetic Diversity Calculations: Compute pairwise identity, haplotype diversity, and polymorphic site statistics to quantify genetic diversity.
Nucleotide Diversity (π) and Average Differences (k): Measure nucleotide diversity (π) and average nucleotide differences (k) across isolates.
Alignment and Phylogenetic Modeling: Align sequences with ClustalW, followed by phylogenetic analysis with IQ-TREE 2, using 1000 bootstrap replicates.
Phylogenetic Network and Tree Visualization: Visualize relationships among isolates using Neighbor-Net and ggtree.
Geographical Distribution Mapping: Map isolate regions using ggplot2 and sf to visualize distribution patterns.
Genomic Variation Analysis: Identify nucleotide differences and mutation hotspots.
Alignment Scores Analysis: Calculate alignment scores and visualize them by host and region.
Phylogenetic Analysis

Steps
Alignment: Align 101 ToBRFV genome sequences using ClustalW (in bash, with FASTA input).
Phylogenetic Modeling: Use IQ-TREE2 for phylogenetic analysis, applying the Maximum Likelihood (ML) method with ModelFinder for optimal model selection, and 1000 ultrafast bootstrap replicates for clade support.
Visualization: Visualize the tree in R using ggtree. Generate pairwise sequence identity matrices using ape and visualize as heatmaps with pheatmap.
Network Analysis: Create a Neighbor-Net phylogenetic network to explore relationships among isolates.
Diversity Metrics: Calculate haplotype diversity, haplotype count, nucleotide diversity, and average nucleotide differences using pegas and ape. Visualize metrics with fmsb radar plots.
Genomic Variation and Distribution Analysis of the Mexican ToBRFV Isolate

Nucleotide Variation Analysis
Using the Biostrings package in RStudio, nucleotide differences were analyzed between the Mexican ToBRFV isolate and 100 reference genomes. Normalization was applied, and statistical summaries highlighted mutation rate variation. Visualization in ggplot2 revealed unique mutations and conserved regions in the Mexican isolate, potentially relevant for functional or evolutionary studies.
Geographical Distribution of ToBRFV Isolates
Geographical metadata for each isolate was visualized using ggplot2 and sf. A world map annotated with isolate counts shows the distribution and clustering of isolates genetically similar to the Mexican genome, indicating global dissemination patterns.
Alignment Scores by Host and Geographical Distribution
Alignment score variations across hosts and regions were analyzed using msa, with scores combined with host and geographical data for ggplot2 scatter plots. Points are colored by host and annotated by region.
Requirements

Make sure the following R packages are installed:
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
Install them using the following commands:
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("ggtree", "treeio", "msa", "Biostrings"))
install.packages(c("ape", "viridis", "ggplot2", "stringr", "dplyr", "tidyr", "pheatmap", "phangorn", "fmsb", "sf"))
Usage

Phylogenetic Tree Generation
Ensure the aligned sequences are saved as a Newick tree file (Resultado_Alignment.aln.contree) in the working directory.
Run generate_phylogenetic_tree.R to generate the phylogenetic tree visualization. This script will produce outputs in both PDF and PNG formats, displaying clade support values and branch length colors.
Genomic Differences Visualization
Place the aligned FASTA file (lcl|Query_2456347.aln) in the working directory.
Run visualize_genomic_differences.R to create visualizations of genome-wide nucleotide differences for the ToBRFV isolates. This script will generate a scatter plot indicating positions with significant variations.
Genomic Differences by Host
Place the aligned FASTA file (lcl|Query_2456347.aln) in the working directory.
Run GenomicDifferencesbyHost.R to analyze and visualize nucleotide differences grouped by host. The script will produce a plot highlighting host-specific genetic variations across the genome.


Scripts

generate_phylogenetic_tree.R: Generates a phylogenetic tree with labeled clades, bootstrap values, and saves the visualization as PDF and PNG.
visualize_genomic_differences.R: Visualizes genomic differences across the ToBRFV genome.
GenomicDifferencesbyHost.R: Visualizes genomic differences categorized by host.
