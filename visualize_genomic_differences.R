# Install and load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("msa")
install.packages("ggplot2")
install.packages("viridis")

library(Biostrings)
library(msa)
library(ggplot2)
library(viridis)
library(viridisLite)
setwd("/Users/katiaavinapadilla/Desktop/SECUENCIASDEPSOMAGENERIKA/")

# Read the FASTA file
alignment <- readDNAStringSet("lcl|Query_2456347.aln", format = "fasta")

# Extract sequence names
seq_names <- names(alignment)

# Align sequences using 'msa'
aligned <- msa(alignment, method = "ClustalW")
aligned_matrix <- as.matrix(aligned)

# Identify the sequence of interest
query_index <- which(seq_names == "lcl|Query_2456347")
if (length(query_index) == 0) {
  stop("The sequence of interest 'lcl|Query_2456347' was not found in the alignment file.")
}

query_seq <- aligned_matrix[query_index, ]

# Calculate differences using vectorized operations
differences <- colSums(aligned_matrix != query_seq)

# Normalize differences
max_differences <- max(differences)
if (max_differences > 0) {
  differences <- differences / max_differences
} else {
  stop("No differences to normalize.")
}

# Create a data frame for ggplot
data <- data.frame(
  Position = 1:ncol(aligned_matrix),
  Differences = differences
)

# Create the plot
p <- ggplot(data, aes(x = Position, y = Differences, color = Differences)) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_color_viridis(option = "viridis") +
  labs(title = "Genomic Differences of ToBRFV Compared to 100 Species",
       x = "Nucleotide Position",
       y = "Normalized Differences") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# Save the figure
ggsave("ToBRFV_Genomic_Differences.png", plot = p, width = 15, height = 8, dpi = 300)

# Display the figure
print(p)

# Select the 5 positions with the highest differences
top_5_points <- data %>% arrange(desc(Differences)) %>% head(5)
print(top_5_points)

# Add labels to the plot
p <- p + 
  geom_text(data = top_5_points, aes(x = Position, y = Differences, label = Position), 
            vjust = -1, size = 4)

# Display the final plot with labels
print(p)
