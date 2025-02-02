# Load required packages
library(seqinr)
library(ggplot2)

# Function to retrieve codon and amino acid from nucleotide position
get_amino_acid <- function(sequence, position, frame_start = 1) {
  relative_position <- position - frame_start + 1
  codon_start <- ((relative_position - 1) %/% 3) * 3 + 1
  codon_end <- codon_start + 2
  codon <- substr(sequence, codon_start + frame_start - 1, codon_end + frame_start - 1)
  amino_acid <- translate(s2c(codon))
  list(
    position = position,
    codon = codon,
    amino_acid = amino_acid,
    codon_start = codon_start + frame_start - 1,
    codon_end = codon_end + frame_start - 1
  )
}

# Define and clean the DNA sequence
dna_sequence <- gsub("\n", "", "ACATATACCAACAACAACAAACAACAAACAACAACATTACAATTACTATTTACAACTACA
ATGGCATACACACAGACAGCTACCACATCCGCTTTGCTCGACACTGTCCGAGGTAACAAT
... (Truncated for readability)
")

# Define positions of interest
positions_to_analyze <- c(176, 423, 532, 561)

# Analyze the positions
results <- lapply(positions_to_analyze, function(pos) {
  get_amino_acid(dna_sequence, pos, frame_start = 1)
})

# Display intermediate results
cat("Intermediate Results:\n")
for (result in results) {
  cat("Nucleotide Position:", result$position, "\n")
  cat("Codon:", result$codon, "\n")
  cat("Amino Acid:", result$amino_acid, "\n")
  cat("Codon Start-End:", result$codon_start, "-", result$codon_end, "\n\n")
}

# Changes in positions 532 and 561
positions_to_change <- list(
  list(position = 532, new_base = "A"),
  list(position = 561, new_base = "C")
)

mutation_results <- lapply(positions_to_change, function(change_info) {
  pos <- change_info$position
  new_base <- change_info$new_base
  
  # Retrieve original codon and amino acid
  original <- get_amino_acid(dna_sequence, pos, frame_start = 1)
  
  # Modify the sequence
  dna_sequence_changed <- dna_sequence
  substr(dna_sequence_changed, pos, pos) <- new_base
  
  # Retrieve mutated codon and amino acid
  changed <- get_amino_acid(dna_sequence_changed, pos, frame_start = 1)
  
  list(
    position = pos,
    new_base = new_base,
    codon_original = original$codon,
    codon_changed = changed$codon,
    amino_acid_original = original$amino_acid,
    amino_acid_changed = changed$amino_acid,
    effect = ifelse(original$amino_acid == changed$amino_acid, "No change", "Change")
  )
})

# Create a dataframe with the results
data <- data.frame(
  Nucleotide_Position = c(176, 423, 532, 561),
  Amino_Acid_Position = c(59, 141, 178, 187), # Manually calculated positions
  Codon_Original = c("TTC", "GCA", "GTT", "TTT"),
  Codon_Changed = c(NA, NA, "GTA", "TTC"), # Only for mutated positions
  Amino_Acid_Original = c("F", "A", "V", "F"),
  Amino_Acid_Changed = c(NA, NA, "I", "F"), # Mutations at 532 and 561
  Mutation_Effect = c("No change", "No change", "Change", "No change"),
  Domain = c("Methyltransferase", "Methyltransferase", "Methyltransferase", "Methyltransferase")
)

# Highlighting colors for mutations
mutation_colors <- c("Change" = "red", "No change" = "blue")

# Define highlighting categories
data$Highlight <- ifelse(data$Amino_Acid_Position %in% c(59, 141), "Mexican ToBRFV isolate", data$Mutation_Effect)

# Custom colors for highlights
highlight_colors <- c(
  "Change" = "red",         # Red for mutations
  "No change" = "blue",     # Blue for no change
  "Mexican ToBRFV isolate" = "green" # Green to highlight Mexican strain
)

# Generate plot with highlighted positions
plot <- ggplot(data, aes(x = Amino_Acid_Position, y = Nucleotide_Position)) +
  geom_point(aes(color = Highlight), size = 4, shape = 21, fill = "white", stroke = 1) +
  geom_text(
    aes(label = ifelse(is.na(Amino_Acid_Changed),
                       Amino_Acid_Original,
                       paste(Amino_Acid_Original, "→", Amino_Acid_Changed))),
    hjust = -0.2, vjust = -0.5, size = 3.5, color = "black"
  ) +
  scale_color_manual(values = highlight_colors) +
  labs(
    x = "Amino Acid Position",
    y = "Nucleotide Position",
    color = "Highlight",
  ) +
  theme_bw()

# Save and display the updated plot
ggsave("Mexican_Strain_Highlighted_ToBRFV.png", plot, dpi = 300, width = 10, height = 6)
print(plot)

############################################ Helicase domain
# Define positions of interest
positions_to_analyze <- c(3246, 2994, 2574)

# Analyze positions
results <- lapply(positions_to_analyze, function(pos) {
  get_amino_acid(dna_sequence, pos, frame_start = 1)
})

# Display results
cat("Intermediate Results:\n")
for (result in results) {
  cat("Nucleotide Position:", result$position, "\n")
  cat("Codon:", result$codon, "\n")
  cat("Amino Acid:", result$amino_acid, "\n")
  cat("Codon Start-End:", result$codon_start, "-", result$codon_end, "\n\n")
}

# Define mutation positions
positions_to_analyze <- list(
  list(position = 2994, new_base = "T"),
  list(position = 2574, new_base = "G")
)

# Analyze mutations
results <- lapply(positions_to_analyze, function(pos_info) {
  analyze_mutation(dna_sequence, pos_info$position, pos_info$new_base)
})

# Display mutation results
cat("Mutation Analysis Results:\n")
for (res in results) {
  cat("Position:", res$position, "\n")
  cat("Original Base:", res$original_base, "\n")
  cat("New Base:", res$new_base, "\n")
  cat("Original Codon:", res$codon_original, "\n")
  cat("Changed Codon:", res$codon_changed, "\n")
  cat("Original Amino Acid:", res$amino_acid_original, "\n")
  cat("Changed Amino Acid:", res$amino_acid_changed, "\n")
  cat("Mutation Effect:", res$mutation_effect, "\n\n")
}

# Create plot for helicase mutations
plot <- ggplot(data, aes(x = Amino_Acid_Position, y = Nucleotide_Position)) +
  geom_point(aes(color = Mutation_Effect), size = 4, shape = 21, fill = "white", stroke = 1) +
  geom_text(
    aes(label = paste(Amino_Acid_Original, "→", Amino_Acid_Changed, sep = "")),
    hjust = -0.2, vjust = -0.5, size = 3.5, color = "black"
  ) +
  scale_color_manual(values = mutation_colors) +
  labs(
    x = "Amino Acid Position",
    y = "Nucleotide Position",
    color = "Mutation Effect",
    title = "Mutations in Helicase Domain"
  ) +
  theme_bw()

# Save and display figure
ggsave("Mutations_in_Helicase.png", plot, dpi = 300, width = 10, height = 6)
print(plot)
