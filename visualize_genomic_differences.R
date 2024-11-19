# Install and load the required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "msa"))
install.packages(c("ggplot2", "viridis"))

library(Biostrings)
library(msa)
library(ggplot2)
library(viridis)

# Read the FASTA file
alignment <- readDNAStringSet("lcl|Query_2456347.aln", format = "fasta")

# Extract the sequence names
seq_names <- names(alignment)
# Display the sequence names in the file
print(seq_names)

# Align the sequences using 'msa'
aligned <- msa(alignment, method = "ClustalW")

# Convert the alignment to a character matrix
aligned_matrix <- as.matrix(aligned)

# Identify the sequence of interest (|Mexican_isolate|)
query_index <- which(seq_names == "|Mexican_isolate| ")
if (length(query_index) == 0) {
  stop("The sequence of interest '|Mexican_isolate|' was not found in the alignment file.")
}

query_seq <- aligned_matrix[query_index, ]

# Calculate the differences
num_sequences <- nrow(aligned_matrix)
sequence_length <- ncol(aligned_matrix)

differences <- numeric(sequence_length)

for (i in 1:num_sequences) {
  if (i != query_index) {
    differences <- differences + as.numeric(aligned_matrix[i, ] != query_seq)
  }
}

# Normalize the differences
differences <- differences / max(differences)

# Create a data frame for ggplot
data <- data.frame(
  Position = 1:sequence_length,
  Differences = differences
)

# Create the plot
p <- ggplot(data, aes(x = Position, y = Differences, color = Differences)) +
  geom_point(size = 1.5) +
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

# Sort the data by the 'Differences' column in descending order
sorted_data <- data[order(-data$Differences), ]

# Select the 5 points with the highest values
top_5_points <- sorted_data[1:5, ]

# Display the positions and values of the top 5 points
cat("Top 5 positions with the highest genomic differences:\n")
print(top_5_points)
