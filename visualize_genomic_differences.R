# Install and load necessary packages
if (!requireNamespace("msa", quietly = TRUE)) install.packages("msa")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")

library(msa)
library(Biostrings)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# Change the working directory to where the file is located, if necessary
setwd("/Users/katiaavinapadilla/Downloads")

# Read the FASTA file
alignment <- readDNAStringSet("lcl|Query_2456347.aln", format="fasta")

# Perform the alignment
aligned <- msa(alignment, method = "ClustalW")
#Asegurar que los nombres de las columnas sean válidos
colnames(aligned_df) <- paste0("V", seq_len(ncol(aligned_df)))

# Extraer nombres de secuencias
seq_names <- names(alignment)

# Verificar la extracción de los nombres de secuencias
print(seq_names)

# Extraer información de hospedantes
host_regex <- "\\[nat_host=([A-Za-z ]+)\\]"
hosts <- str_extract(seq_names, host_regex)
hosts <- str_replace_all(hosts, "\\[nat_host=|\\]", "")
hosts[hosts == "tomato"] <- "Solanum lycopersicum"

# Verificar la extracción de hospedantes
print(hosts)

# Asegurar que ambos vectores tengan la misma longitud
if (length(seq_names) != length(hosts)) {
  missing_count <- length(seq_names) - length(hosts)
  hosts <- c(hosts, rep("Unknown", missing_count))
}

# Crear un data frame con los nombres de secuencias y sus respectivos hospedantes
host_data <- data.frame(Sequence = seq_names, Host = hosts, stringsAsFactors = FALSE)

# Combinar datos de alineación con información de hospedantes
combined_data <- cbind(aligned_df, host_data)

# Identificar posiciones con diferencias
diff_positions <- sapply(1:ncol(aligned_df), function(i) {
  length(unique(aligned_df[, i])) > 1
})

# Verificar las posiciones con diferencias
print(diff_positions)

# Extraer posiciones con diferencias
diff_positions <- which(diff_positions)

# Verificar las posiciones con diferencias extraídas
print(diff_positions)

# Filtrar datos combinados para incluir solo posiciones con diferencias
filtered_data <- combined_data[, c(diff_positions, ncol(combined_data))]

# Verificar los datos filtrados
print(head(filtered_data))

# Transformar datos para visualización
plot_data <- filtered_data %>%
  pivot_longer(cols = -Host, names_to = "Position", values_to = "Nucleotide") %>%
  mutate(Position = as.numeric(gsub("V", "", Position))) %>%
  filter(!is.na(Position))

# Verificar los datos transformados para visualización
print(head(plot_data))

# Calcular diferencias normalizadas
plot_data <- plot_data %>%
  group_by(Position, Host) %>%
  summarize(Differences = n_distinct(Nucleotide)) %>%
  ungroup() %>%
  mutate(Differences = scale(Differences))

# Verificar los datos calculados
print(head(plot_data))

# Guardar los datos de diferencias en un archivo CSV
write.csv(plot_data, file = "genomic_differences_by_host.csv", row.names = FALSE)

# Crear gráfico de dispersión
p <- ggplot(plot_data, aes(x = Position, y = Differences, color = Host)) +
  geom_point(size = 3) +
  scale_color_viridis(discrete = TRUE) +
  labs(title = "Genomic Differences of ToBRFV Grouped by Host",
       x = "Nucleotide Position",
       y = "Normalized Differences",
       color = "Host") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "right")

print(p)
Verificacio