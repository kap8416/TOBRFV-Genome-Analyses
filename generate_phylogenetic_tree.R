# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("treeio")
install.packages("ape")
install.packages("viridis")
install.packages("ggplot2")

library(ggtree)
library(treeio)
library(ape)
library(viridis)
library(ggplot2)

# Change the working directory to where the file is located, if necessary
setwd("/Users/katiaavinapadilla/Downloads")

# Verify the file exists in the directory
list.files()

# Read the tree from a Newick file
tree <- read.tree("TOBRFV.nwk")

# Create a data frame with labels and their colors
label_colors <- data.frame(label = tree$tip.label)
label_colors$color <- ifelse(label_colors$label == "TOBRFV Mexican isolate", "red", "darkblue")

# Customized visualization of the tree with adjusted labels
p <- ggtree(tree) + 
  geom_tree(aes(color=branch.length), size=0.5) +
  scale_color_viridis_c() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"),
        plot.title = element_text(size=20, face="bold", hjust = 0.5),
        legend.position="right") +
  geom_tiplab(aes(label=label, subset=label=="TOBRFV Mexican isolate"), size=2, color="red") +
  geom_tiplab(aes(label=label, subset=label!="TOBRFV Mexican isolate"), size=2, color="darkblue") +
  ggtitle("Phylogenetic Tree of ToBRFV")

# Display the plot
print(p)

# Save the plot as a PDF file
ggsave("Phylogenetic_Tree_ToBRFV.pdf", plot=p, width=10, height=8)

# Save the plot as a PNG file
ggsave("Phylogenetic_Tree_ToBRFV.png", plot=p, width=10, height=8, dpi=300)