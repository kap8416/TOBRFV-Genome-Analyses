# TOBRFV-Genome-Analyses

## Introduction

This repository provides a comprehensive bioinformatics workflow for analyzing the genome of the Mexican isolate of **Tomato Brown Rugose Fruit Virus (ToBRFV)**, a highly virulent virus affecting tomato and pepper plants. The workflow includes sequence retrieval, quality control, genome assembly, and analyses aimed at understanding the genetic diversity, phylogeny, and geographical distribution of ToBRFV isolates. This project is designed for reproducibility, with all scripts available in this repository.

The analysis pipeline includes genome retrieval, metadata extraction, genetic diversity calculations, phylogenetic modeling, and visualization to explore relationships and variations within ToBRFV isolates.

**Figure 1** below illustrates the analysis steps, with all scripts accessible in this repository.


![output-9](https://github.com/user-attachments/assets/6cf7e57a-57e3-4d01-b2c6-508247d209bc)

---

## Analysis Pipeline Overview

1. **Sequence Retrieval**: Retrieve full genome sequences of 100 ToBRFV isolates from GenBank, including the Mexican isolate.
2. **Metadata Extraction**: Extract geographical and host information for each isolate.
3. **Genetic Diversity Calculations**: Compute pairwise identity, haplotype diversity, and polymorphic site statistics to quantify genetic diversity.
4. **Nucleotide Diversity (π) and Average Differences (k)**: Measure nucleotide diversity (π) and average nucleotide differences (k) across isolates.
5. **Alignment and Phylogenetic Modeling**: Align sequences with ClustalW, followed by phylogenetic analysis with IQ-TREE 2, using 1000 bootstrap replicates.
6. **Phylogenetic Network and Tree Visualization**: Visualize isolate relationships using Neighbor-Net and `ggtree`.
7. **Geographical Distribution Mapping**: Map isolate regions using `ggplot2` and sf to visualize distribution patterns.
8. **Genomic Variation Analysis**: Identify nucleotide differences and mutation hotspots.
9. **Alignment Scores Analysis**: Calculate alignment scores and visualize them by host and region.

---

## Phylogenetic Analysis

### Steps

1. **Alignment**: Align 101 ToBRFV genome sequences using ClustalW (in bash, with FASTA input).
2. **Phylogenetic Modeling**: Use IQ-TREE2 for phylogenetic analysis, applying the Maximum Likelihood (ML) method with ModelFinder for optimal model selection, and 1000 ultrafast bootstrap replicates for clade support.
3. **Visualization**: Visualize the tree in R using `ggtree`. Generate pairwise sequence identity matrices using `ape` and display as heatmaps with `pheatmap`.
4. **Network Analysis**: Create a Neighbor-Net phylogenetic network to explore isolate relationships.
5. **Diversity Metrics**: Calculate haplotype diversity, haplotype count, nucleotide diversity, and average nucleotide differences using `pegas` and `ape`. Visualize metrics with `fmsb` radar plots.

---

## Genomic Variation and Distribution Analysis of the Mexican ToBRFV Isolate

### Nucleotide Variation Analysis

Using the `Biostrings` package in RStudio, nucleotide differences were analyzed between the Mexican ToBRFV isolate and 100 reference genomes. Normalization was applied, and statistical summaries highlighted mutation rate variation. Visualization in `ggplot2` revealed unique mutations and conserved regions in the Mexican isolate, with potential relevance for functional or evolutionary studies.

### Geographical Distribution of ToBRFV Isolates

Geographical metadata for each isolate was visualized using `ggplot2` and sf. A world map annotated with isolate counts shows the distribution and clustering of isolates genetically similar to the Mexican genome, indicating global dissemination patterns.

### Alignment Scores by Host and Geographical Distribution

Using `msa`, alignment score variations were analyzed across hosts and regions, with scores combined with host and geographical data for `ggplot2` scatter plots. Points are colored by host and annotated by region.

---

## Requirements

Ensure the following R packages are installed:

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
- `pheatmap`
- `phangorn`
- `fmsb`
- `sf`

Install them using the following commands:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("ggtree", "treeio", "msa", "Biostrings"))
install.packages(c("ape", "viridis", "ggplot2", "stringr", "dplyr", "tidyr", "pheatmap", "phangorn", "fmsb", "sf"))
```
---
## Usage

Phylogenetic Tree Generation

-Ensure the Newick tree file (Resultado_Alignment.aln.contree) is in the working directory.
-Run generate_phylogenetic_tree.R to generate the tree visualization, which will produce outputs in both PDF and PNG formats with clade support values.

Genomic Differences Visualization

-Place the aligned FASTA file (lcl|Query_2456347.aln) in the working directory.
-Run visualize_genomic_differences.R to create visualizations of genome-wide nucleotide differences.

Genomic Differences by Host
-Place the aligned FASTA file (lcl|Query_2456347.aln) in the working directory.
-Run GenomicDifferencesbyHost.R to analyze and visualize nucleotide differences grouped by host.

## Scripts

generate_phylogenetic_tree.R: Generates a phylogenetic tree with labeled clades and bootstrap values, saving the visualization as PDF and PNG files.
visualize_genomic_differences.R: Visualizes genome-wide nucleotide differences across the ToBRFV genome.
GenomicDifferencesbyHost.R: Visualizes nucleotide differences grouped by host.

