# Phylogenetics_and_Genetic_Diversity.R
# Comprehensive analysis of ToBRFV genetic diversity and phylogenetics.

# ==========================
# Load Required Libraries
# ==========================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio", "ape", "phangorn", "pheatmap", "fmsb"))
library(ggtree)
library(treeio)
library(ape)
library(phangorn)
library(pheatmap)
library(fmsb)

# ==========================
# Set Working Directory
# ==========================
# Update this path with the directory containing your data files
setwd("/Users/katiaavinapadilla/Desktop/SECUENCIASDEPSOMAGENERIKA/")

# ==========================
# Phylogenetic Tree Analysis
# ==========================

# Step 1: Perform alignment and phylogenetic analysis outside R (bash/terminal):
# Commands:
# clustalw -INFILE=101sequences.fasta -OUTFILE=Resultado_Alignment.aln -TYPE=DNA -OUTPUT=CLUSTAL
# iqtree2 -s Resultado_Alignment.aln -m GTR+G -bb 1000 -nt AUTO

# Step 2: Load the Newick tree generated by IQ-TREE2
phylo_tree <- read.tree("Resultado_Alignment.aln.contree")

# Step 3: Visualize the phylogenetic tree
# Create the tree visualization with bootstrap values using ggtree
p <- ggtree(phylo_tree, layout = "rectangular") +
  geom_tree(aes(color = branch.length), size = 0.5) + # Color branches by branch length
  scale_color_viridis_c() + # Professional color scale (viridis)
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "right") +
  geom_tiplab(aes(label = label), size = 2.5, color = "darkblue", hjust = -0.2) + # Tip labels
  geom_text2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 70, # Show bootstrap values >= 70
                 label = label), color = "red", hjust = -0.2, size = 3) + # Highlight high bootstrap values in red
  ggtitle("Phylogenetic Tree of ToBRFV Isolates") # Plot title

# Save and display the tree
ggsave("Phylogenetic_Tree_ToBRFV.pdf", plot = p, width = 12, height = 8)
print(p)

# ==========================
# Pairwise Sequence Identity Heatmap
# ==========================

# Load the alignment file
alignment <- read.dna("lcl|Query_2456347.aln", format = "fasta")

# Compute pairwise sequence identity
dist_matrix <- dist.dna(alignment, model = "raw")
identity_matrix <- (1 - dist_matrix) * 100

# Simplify labels
rownames(identity_matrix) <- sapply(rownames(identity_matrix), function(label) {
  match <- regmatches(label, regexpr("(?<=\\|)[A-Za-z0-9]+\\.[0-9]+", label, perl = TRUE))
  if (length(match) > 0) return(match) else return(label)
})
colnames(identity_matrix) <- rownames(identity_matrix)

# Create heatmap
heatmap <- pheatmap(identity_matrix,
                    color = colorRampPalette(c("blue", "green", "yellow", "red"))(50),
                    main = "Pairwise Sequence Identities (%) Between ToBRFV Isolates",
                    fontsize_row = 6,
                    fontsize_col = 6)
pdf("Pairwise_Sequence_Identities.pdf", width = 20, height = 20)
print(heatmap)
dev.off()

# ==========================
# Neighbor-Net Network
# ==========================
alignment <- read.dna("/Users/katiaavinapadilla/Desktop/SECUENCIASDEPSOMAGENERIKA/lcl|Query_2456347.aln", format = "fasta")
#Calculate the Distance Matrix
#We will use the Kimura 2-parameter model (K80) to compute the distance matrix between sequences.
dist_matrix <- dist.dna(alignment, model = "K80")
#Construct the Neighbor-Net Network
#Using the calculated distance matrix, you can construct a Neighbor-Net network using the phangorn package.
net <- neighborNet(dist_matrix)


#Review and Extract the Exact Identifier

#Extract the specific identifier after the first '|'
labels <- sapply(net$tip.label, function(label) {
  match <- regmatches(label, regexpr("(?<=\\|)[^|]+(?=\\|)", label, perl = TRUE))
  if (length(match) > 0) {
    return(match)  
  } else {
    return(label)  # If no match is found, return the original label
  }
})

# Assign the corrected labels to the network
net$tip.label <- labels
plot(net, tip.color = ifelse(net$tip.label == "Mexican_isolate", "red", "black"))
# Plot the network
plot(net)

# ==========================
# Genetic Diversity Metrics and Radar Plot
# ==========================
haplotypes <- haplotype(alignment)
H <- nrow(haplotypes)
Hd <- hap.div(haplotypes)
S <- length(seg.sites(alignment))
pi <- nuc.div(alignment)
k <- theta.k(alignment)

# Radar chart data
data_radar <- data.frame(
  "Number of Haplotypes (H)" = c(80, 0, H),
  "Haplotype Diversity (Hd)" = c(1, 0, Hd),
  "Polymorphic Sites (S)" = c(250, 0, S),
  "Nucleotide Diversity (π)" = c(0.01, 0, pi),
  "Average Nucleotide Differences (k)" = c(0.5, 0, k)
)

radarchart(data_radar, axistype = 1,
           pcol = rgb(0.2, 0.5, 0.5, 0.9), 
           pfcol = rgb(0.2, 0.5, 0.5, 0.5), 
           plwd = 2, 
           cglcol = "grey", 
           caxislabels = seq(0, 1, 0.2),
           vlcex = 0.8)
title(main = "Genetic Diversity Parameters of ToBRFV")


# ==========================
# End of Script
# ==========================

