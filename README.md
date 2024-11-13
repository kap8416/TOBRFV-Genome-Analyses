TOBRFV-Genome-Analyses

Introduction

This repository provides a complete bioinformatics workflow to analyze the genome of the Mexican isolate of Tomato Brown Rugose Fruit Virus (ToBRFV), a highly virulent virus impacting tomato and pepper plants. The pipeline includes sequence retrieval, quality control, genome assembly, and various analyses to investigate the genetic diversity, phylogeny, and geographical distribution of ToBRFV isolates. This project is designed for reproducibility and further exploration, with all scripts available in this repository.
The workflow includes genome retrieval, metadata extraction, genetic diversity calculations, phylogenetic modeling, and visualization to uncover relationships and variations within ToBRFV isolates.

Figure 1 below illustrates the analysis steps, with all scripts accessible in this repository.




![output-9](https://github.com/user-attachments/assets/2c1938a2-6643-4d31-9780-73d232f3f748)





Analysis Pipeline Overview

Sequence Retrieval: Retrieve the full genome sequences of 100 ToBRFV isolates from GenBank, including the Mexican isolate.
Metadata Extraction: Extract geographical and host information for each isolate.
Genetic Diversity Calculations: Compute pairwise identity, haplotype diversity, and polymorphic site statistics to quantify genetic diversity.
Nucleotide Diversity (π) and Average Differences (k): Measure nucleotide diversity (π) and average nucleotide differences (k) across isolates.
Alignment and Phylogenetic Modeling: Align sequences with ClustalW and perform phylogenetic analysis with IQ-TREE 2 using 1000 bootstrap replicates.
Phylogenetic Network and Tree Visualization: Visualize isolate relationships using Neighbor-Net and ggtree.
Geographical Distribution Mapping: Use ggplot2 and sf to map isolate distribution.
Genomic Variation Analysis: Identify nucleotide differences and mutation hotspots.
Alignment Scores Analysis: Calculate alignment scores and visualize them by host and region.
Phylogenetic Analysis

Steps
Alignment: Align 101 ToBRFV genome sequences using ClustalW (in bash, with FASTA input).
Phylogenetic Modeling: Use IQ-TREE2 for phylogenetic analysis, applying the Maximum Likelihood (ML) method and ModelFinder for optimal model selection with 1000 ultrafast bootstrap replicates.
Visualization: Visualize the tree in R using ggtree. Generate pairwise sequence identity matrices using ape and display as heatmaps with pheatmap.
Network Analysis: Create a Neighbor-Net phylogenetic network to explore isolate relationships.
Diversity Metrics: Calculate key diversity metrics (haplotype diversity, haplotype count, nucleotide diversity, and average nucleotide differences) with pegas and ape. Visualize with fmsb radar plots.
Genomic Variation and Distribution Analysis of the Mexican ToBRFV Isolate

Nucleotide Variation Analysis
Using the Biostrings package in RStudio, we analyzed nucleotide differences between the Mexican ToBRFV isolate and 100 reference genomes. Normalization was applied, and statistical summaries highlighted mutation rate variation. Visualization with ggplot2 identified unique mutations and conserved regions in the Mexican isolate, with potential significance for functional or evolutionary studies.
Geographical Distribution of ToBRFV Isolates
Geographical metadata for each isolate was visualized using ggplot2 and sf. A world map with isolate counts shows the distribution and clustering of isolates genetically similar to the Mexican genome, revealing global dissemination patterns.
Alignment Scores by Host and Geographical Distribution
Using msa, alignment score variations were analyzed across hosts and regions, with scores combined with host and geographical data for ggplot2 scatter plots. Points are colored by host and annotated by region.
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
Ensure the Newick tree file (Resultado_Alignment.aln.contree) is in the working directory.
Run generate_phylogenetic_tree.R to generate the tree visualization, which will produce outputs in both PDF and PNG formats with clade support values.
Genomic Differences Visualization
Place the aligned FASTA file (lcl|Query_2456347.aln) in the working directory.
Run visualize_genomic_differences.R to create visualizations of genome-wide nucleotide differences.
Genomic Differences by Host
Place the aligned FASTA file (lcl|Query_2456347.aln) in the working directory.
Run GenomicDifferencesbyHost.R to analyze and visualize nucleotide differences grouped by host.
Scripts

generate_phylogenetic_tree.R: Generates a phylogenetic tree with labeled clades and bootstrap values, saving the visualization as PDF and PNG files.
visualize_genomic_differences.R: Visualizes genome-wide nucleotide differences across the ToBRFV genome.
GenomicDifferencesbyHost.R: Visualizes nucleotide differences grouped by host.
This README offers a structured and detailed guide for setting up, using, and interpreting results from the ToBRFV Genome Analysis workflow, enabling reproducible research and insights into the virus’s genetic diversity, transmission pathways, host adaptations, and global impact.
