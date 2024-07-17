# Crear el gráfico
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

# Guardar la figura
ggsave("ToBRFV_Genomic_Differences.png", plot = p, width = 15, height = 8, dpi = 300)

# Mostrar la figura
print(p)
#####
# Calcular las diferencias
num_sequences <- nrow(aligned_matrix)
sequence_length <- ncol(aligned_matrix)

differences <- numeric(sequence_length)

for (i in 1:num_sequences) {
  if (i != query_index) {
    differences <- differences + (aligned_matrix[i, ] != query_seq)
  }
}

# Normalizar las diferencias
max_differences <- max(differences)
differences <- differences / max_differences

# Crear un data frame para ggplot
data <- data.frame(
  Position = 1:sequence_length,
  Differences = differences
)

# Calcular la variabilidad (desviación estándar) en cada posición
variability <- apply(aligned_matrix, 2, function(column) sd(as.numeric(factor(column))))

# Identificar la región con mayor variabilidad
high_variability_region <- which.max(variability)

# Extraer información del hospedante
host_regex <- "\\[nat_host=([A-Za-z ]+)\\]"
hosts <- str_extract(seq_names, host_regex)
hosts <- str_replace_all(hosts, "\\[nat_host=|\\]", "")
hosts[hosts == "tomato"] <- "Solanum lycopersicum"

# Asegurar que ambos vectores tienen la misma longitud llenando la información faltante con "Unknown"
if (length(seq_names) != length(hosts)) {
  missing_count <- length(seq_names) - length(hosts)
  hosts <- c(hosts, rep("Unknown", missing_count))
}

# Crear un dataframe con los nombres de secuencia y sus respectivos hospedantes
host_data <- data.frame(Sequence = seq_names, Host = hosts, stringsAsFactors = FALSE)

# Añadir la información del hospedante a cada secuencia en la matriz alineada
aligned_df <- as.data.frame(aligned_matrix)
aligned_df$Sequence <- seq_names

# Encontrar las posiciones con diferencias
diff_positions <- which(differences > 0)

# Crear un dataframe para las diferencias y agregar la información del hospedante
diff_data <- data.frame(
  Position = rep(diff_positions, each = num_sequences),
  Sequence = rep(seq_names, times = length(diff_positions)),
  Nucleotide = as.vector(aligned_matrix[, diff_positions])
)

# Combinar los datos de diferencias con los datos de hospedantes
combined_data <- merge(diff_data, host_data, by = "Sequence")

# Identificar las posiciones más variables
threshold <- 0.8  # Este valor puede ajustarse para definir lo que consideras "alto"
high_variability_positions <- data[data$Differences >= threshold, ]

# Ver las posiciones con alta variabilidad
print(high_variability_positions)

# Filtrar las diferencias por hospedante
combined_data_filtered <- combined_data %>%
  filter(Position %in% high_variability_positions$Position)

# Crear el gráfico de puntos con resaltado de la región de alta variabilidad
p <- ggplot(combined_data_filtered, aes(x = Position, y = Nucleotide, color = Host)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = high_variability_positions$Position, linetype = "dashed", color = "red", size = 1) +
  scale_color_viridis(discrete = TRUE) +
  labs(title = "Genomic Differences of ToBRFV Grouped by Host",
       x = "Nucleotide Position",
       y = "Differences",
       color = "Host") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))

# Guardar la figura
ggsave("ToBRFV_Genomic_Differences_By_Host.png", plot = p, width = 15, height = 8, dpi = 300)

# Mostrar la figura
print(p)
###
Convertir el alineamiento a una matriz de caracteres
aligned_matrix <- as.matrix(aligned)

# Identificar la secuencia de interés (lcl|Query_2456347)
seq_names <- names(alignment)
query_index <- which(seq_names == "lcl|Query_2456347")
if (length(query_index) == 0) {
  stop("La secuencia de interés 'lcl|Query_2456347' no se encontró en el archivo de alineamiento.")
}

query_seq <- aligned_matrix[query_index, ]

# Calcular las diferencias
num_sequences <- nrow(aligned_matrix)
sequence_length <- ncol(aligned_matrix)

differences <- numeric(sequence_length)

for (i in 1:num_sequences) {
  if (i != query_index) {
    differences <- differences + (aligned_matrix[i, ] != query_seq)
  }
}

# Normalizar las diferencias
max_differences <- max(differences)
differences <- differences / max_differences

# Crear un data frame para ggplot
data <- data.frame(
  Position = 1:sequence_length,
  Differences = differences
)

# Identificar las posiciones más variables
threshold <- 0.8  # Este valor puede ajustarse para definir lo que consideras "alto"
high_variability_positions <- data[data$Differences >= threshold, ]

# Ver las posiciones con alta variabilidad
print(high_variability_positions)

# Extraer información del hospedante
host_regex <- "\\[nat_host=([A-Za-z ]+)\\]"
hosts <- str_extract(seq_names, host_regex)
hosts <- str_replace_all(hosts, "\\[nat_host=|\\]", "")
hosts[hosts == "tomato"] <- "Solanum lycopersicum"

# Asegurar que ambos vectores tienen la misma longitud llenando la información faltante con "Unknown"
if (length(seq_names) != length(hosts)) {
  missing_count <- length(seq_names) - length(hosts)
  hosts <- c(hosts, rep("Unknown", missing_count))
}

# Crear un dataframe con los nombres de secuencia y sus respectivos hospedantes
host_data <- data.frame(Sequence = seq_names, Host = hosts, stringsAsFactors = FALSE)

# Añadir la información del hospedante a cada secuencia en la matriz alineada
aligned_df <- as.data.frame(aligned_matrix)
aligned_df$Sequence <- seq_names

# Encontrar las posiciones con diferencias
diff_positions <- which(differences > 0)

# Crear un dataframe para las diferencias y agregar la información del hospedante
diff_data <- data.frame(
  Position = rep(diff_positions, each = num_sequences),
  Sequence = rep(seq_names, times = length(diff_positions)),
  Nucleotide = as.vector(aligned_matrix[, diff_positions])
)

# Combinar los datos de diferencias con los datos de hospedantes
combined_data <- merge(diff_data, host_data, by = "Sequence")

# Filtrar las diferencias para las posiciones de alta variabilidad
high_variability_data <- combined_data %>%
  filter(Position %in% high_variability_positions$Position)

# Analizar la significancia estadística de las diferencias por hospedante
stat_tests <- high_variability_data %>%
  group_by(Position) %>%
  summarise(p_value = t.test(Nucleotide ~ Host)$p.value)

# Agregar los resultados estadísticos al dataframe
high_variability_data <- high_variability_data %>%
  left_join(stat_tests, by = "Position")

# Crear el gráfico de puntos con resaltado de las regiones de alta variabilidad y significancia estadística
p <- ggplot(high_variability_data, aes(x = Position, y = Nucleotide, color = Host)) +
  geom_point(size = 1.5) +
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

# Guardar la figura
ggsave("ToBRFV_Genomic_Differences_By_Host.png", plot = p, width = 15, height = 8, dpi = 300)

# Mostrar la figura
print(p)

####
library(msa)
library(Biostrings)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(ggpubr)
# Convertir el alineamiento a una matriz de caracteres
aligned_matrix <- as.matrix(aligned)

# Identificar la secuencia de interés (lcl|Query_2456347)
seq_names <- names(alignment)
query_index <- which(seq_names == "lcl|Query_2456347")
if (length(query_index) == 0) {
  stop("La secuencia de interés 'lcl|Query_2456347' no se encontró en el archivo de alineamiento.")
}

query_seq <- aligned_matrix[query_index, ]

# Calcular las diferencias
num_sequences <- nrow(aligned_matrix)
sequence_length <- ncol(aligned_matrix)

differences <- numeric(sequence_length)

for (i in 1:num_sequences) {
  if (i != query_index) {
    differences <- differences + (aligned_matrix[i, ] != query_seq)
  }
}

# Normalizar las diferencias
max_differences <- max(differences)
differences <- differences / max_differences

# Crear un data frame para ggplot
data <- data.frame(
  Position = 1:sequence_length,
  Differences = differences
)

# Identificar las posiciones más variables
threshold <- 0.8  # Este valor puede ajustarse para definir lo que consideras "alto"
high_variability_positions <- data[data$Differences >= threshold, ]

# Ver las posiciones con alta variabilidad
print(high_variability_positions)

# Extraer información del hospedante
host_regex <- "\\[nat_host=([A-Za-z ]+)\\]"
hosts <- str_extract(seq_names, host_regex)
hosts <- str_replace_all(hosts, "\\[nat_host=|\\]", "")
hosts[hosts == "tomato"] <- "Solanum lycopersicum"

# Asegurar que ambos vectores tienen la misma longitud llenando la información faltante con "Unknown"
if (length(seq_names) != length(hosts)) {
  missing_count <- length(seq_names) - length(hosts)
  hosts <- c(hosts, rep("Unknown", missing_count))
}

# Crear un dataframe con los nombres de secuencia y sus respectivos hospedantes
host_data <- data.frame(Sequence = seq_names, Host = hosts, stringsAsFactors = FALSE)

# Añadir la información del hospedante a cada secuencia en la matriz alineada
aligned_df <- as.data.frame(aligned_matrix)
aligned_df$Sequence <- seq_names

# Encontrar las posiciones con diferencias
diff_positions <- which(differences > 0)

# Crear un dataframe para las diferencias y agregar la información del hospedante
diff_data <- data.frame(
  Position = rep(diff_positions, each = num_sequences),
  Sequence = rep(seq_names, times = length(diff_positions)),
  Nucleotide = as.vector(aligned_matrix[, diff_positions])
)

# Combinar los datos de diferencias con los datos de hospedantes
combined_data <- merge(diff_data, host_data, by = "Sequence")

# Filtrar las diferencias para las posiciones de alta variabilidad
high_variability_data <- combined_data
# Filtrar las diferencias para las posiciones de alta variabilidad
high_variability_data <- combined_data %>%
  filter(Position %in% high_variability_positions$Position)

# Filtrar posiciones que tengan exactamente dos niveles de hospedante
positions_with_two_levels <- high_variability_data %>%
  group_by(Position) %>%
  filter(n_distinct(Host) == 2) %>%
  pull(Position) %>%
  unique()

filtered_high_variability_data <- high_variability_data %>%
  filter(Position %in% positions_with_two_levels)

# Analizar la significancia estadística de las diferencias por hospedante
stat_tests <- filtered_high_variability_data %>%
  group_by(Position) %>%
  summarise(p_value = ifelse(n_distinct(Nucleotide) > 1, t.test(as.numeric(Nucleotide) ~ Host)$p.value, NA)) %>%
  filter(!is.na(p_value))

# Agregar los resultados estadísticos al dataframe
filtered_high_variability_data <- filtered_high_variability_data %>%
  left_join(stat_tests, by = "Position")

# Crear el gráfico de puntos con resaltado de las regiones de alta variabilidad y significancia estadística
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

# Guardar la figura
ggsave("ToBRFV_Genomic_Differences_By_Host.png", plot = p, width = 15, height = 8, dpi = 300)

# Mostrar la figura
print(p)

# Extraer y mostrar las posiciones de alta variabilidad con significancia estadística
high_variability_significant <- filtered_high_variability_data %>%
  filter(p_value < 0.05)

print(high_variability_significant)
##
# Convierte el alineamiento a un objeto de tipo data.frame
alignment_df <- as.data.frame(as.matrix(alignment))

# Extrae las posiciones variables
variable_positions <- apply(alignment_df, 2, function(col) length(unique(col)) > 1)

# Filtra solo las posiciones variables
variable_alignment <- alignment_df[, variable_positions]

# Extrae la información de hospedante
seq_names <- names(input_sequences)
host_regex <- "\\[nat_host=([A-Za-z ]+)\\]"
hosts <- str_extract(seq_names, host_regex)
hosts <- str_replace_all(hosts, "\\[nat_host=|\\]", "")
hosts <- ifelse(is.na(hosts), "Unknown", hosts)

# Agrega la información de hospedante al alineamiento
variable_alignment$Host <- hosts

# Convierte el data.frame a formato largo para ggplot
variable_alignment_long <- melt(variable_alignment, id.vars = "Host")

# Crea la visualización de las diferencias por hospedante
ggplot(variable_alignment_long, aes(x = variable, y = value, color = Host)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Variabilidad en Alineamientos por Hospedante",
       x = "Posición en el Genoma",
       y = "Nucleótido",
       color = "Hospedante") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
