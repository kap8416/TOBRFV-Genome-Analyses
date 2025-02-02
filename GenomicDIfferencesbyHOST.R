# Load required libraries
library(msa)
library(Biostrings)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggpubr)

# Convert the alignment to a character matrix
aligned_matrix <- as.matrix(aligned)

# Identify the query sequence (lcl|Query_2456347)
seq_names <- names(alignment)
query_index <- which(seq_names == "lcl|Query_2456347")
if (length(query_index) == 0) {
  stop("The query sequence 'lcl|Query_2456347' was not found in the alignment file.")
}

query_seq <- aligned_matrix[query_index, ]

# Calculate differences
num_sequences <- nrow(aligned_matrix)
sequence_length <- ncol(aligned_matrix)

differences <- numeric(sequence_length)

for (i in 1:num_sequences) {
  if (i != query_index) {
    differences <- differences + (aligned_matrix[i, ] != query_seq)
  }
}

# Normalize differences
max_differences <- max(differences)
differences <- differences / max_differences

# Create a data frame for ggplot
data <- data.frame(
  Position = 1:sequence_length,
  Differences = differences
)

# Identify highly variable positions
threshold <- 0.8  # Adjust this threshold to define "high variability"
high_variability_positions <- data[data$Differences >= threshold, ]

# Display highly variable positions
print(high_variability_positions)

# Extract host information
host_regex <- "\\[nat_host=([A-Za-z ]+)\\]"
hosts <- str_extract(seq_names, host_regex)
hosts <- str_replace_all(hosts, "\\[nat_host=|\\]", "")
hosts[hosts == "tomato"] <- "Solanum lycopersicum"

# Ensure both vectors have the same length by filling missing values with "Unknown"
if (length(seq_names) != length(hosts)) {
  missing_count <- length(seq_names) - length(hosts)
  hosts <- c(hosts, rep("Unknown", missing_count))
}

# Create a dataframe with sequence names and corresponding hosts
host_data <- data.frame(Sequence = seq_names, Host = hosts, stringsAsFactors = FALSE)

# Add host information to each sequence in the aligned matrix
aligned_df <- as.data.frame(aligned_matrix)
aligned_df$Sequence <- seq_names

# Find positions with differences
diff_positions <- which(differences > 0)

# Create a dataframe for differences and add host information
diff_data <- data.frame(
  Position = rep(diff_positions, each = num_sequences),
  Sequence = rep(seq_names, times = length(diff_positions)),
  Nucleotide = as.vector(aligned_matrix[, diff_positions])
)

# Merge difference data with host information
combined_data <- merge(diff_data, host_data, by = "Sequence")

# Filter differences for highly variable positions
high_variability_data <- combined_data %>%
  filter(Position %in% high_variability_positions$Position)

# Filter positions that have exactly two host levels
positions_with_two_levels <- high_variability_data %>%
  group_by(Position) %>%
  filter(n_distinct(Host) == 2) %>%
  pull(Position) %>%
  unique()

filtered_high_variability_data <- high_variability_data %>%
  filter(Position %in% positions_with_two_levels)

# Perform statistical significance analysis of differences by host
stat_tests <- filtered_high_variability_data %>%
  group_by(Position) %>%
  summarise(p_value = ifelse(n_distinct(Nucleotide) > 1, t.test(as.numeric(Nucleotide) ~ Host)$p.value, NA)) %>%
  filter(!is.na(p_value))

# Add statistical results to the dataframe
filtered_high_variability_data <- filtered_high_variability_data %>%
  left_join(stat_tests, by = "Position")

# Create a scatter plot highlighting regions of high variability and statistical significance
p <- ggplot(filtered_high_variability_data, aes(x = Position, y = as.numeric(Nucleotide), color = Host)) +
  geom_point(size = 3) +
  geom_vline(xintercept = high_variability_positions$Position, linetype = "dashed", color = "red", size = 1) +
  geom_text(data = stat_tests, aes(x = Position, y = Inf, label = sprintf("p=%.3f", p_value)),
            vjust = 1.5, color = "blue") +
  scale_color_viridis(discrete = TRUE) +
  labs(title = "Genomic Differences of ToBRFV Grouped by Host",
       x = "Nucleotide Position",
       y = "Differences",
       color = "Host") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))

# Save the figure
ggsave("ToBRFV_Genomic_Differences_By_Host.png", plot = p, width = 15, height = 8, dpi = 300)

# Display the figure
print(p)

# Extract and display highly variable positions with statistical significance
high_variability_significant <- filtered_high_variability_data %>%
  filter(p_value < 0.05)

print(high_variability_significant)

## Convert the alignment to a dataframe
alignment_df <- as.data.frame(as.matrix(alignment))

# Identify variable positions
variable_positions <- apply(alignment_df, 2, function(col) length(unique(col)) > 1)

# Filter only variable positions
variable_alignment <- alignment_df[, variable_positions]

# Extract host information
seq_names <- names(input_sequences)
host_regex <- "\\[nat_host=([A-Za-z ]+)\\]"
hosts <- str_extract(seq_names, host_regex)
hosts <- str_replace_all(hosts, "\\[nat_host=|\\]", "")
hosts <- ifelse(is.na(hosts), "Unknown", hosts)

# Add host information to the alignment
variable_alignment$Host <- hosts

# Convert the dataframe to long format for ggplot
variable_alignment_long <- melt(variable_alignment, id.vars = "Host")

# Create a visualization of differences by host
ggplot(variable_alignment_long, aes(x = variable, y = value, color = Host)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Variability in Alignments by Host",
       x = "Genomic Position",
       y = "Nucleotide",
       color = "Host") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
