# Install packages if they are not already installed
required_packages <- c("stringr", "dplyr", "ggplot2", "sf", "rnaturalearth", "rnaturalearthdata")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load the necessary libraries
library(stringr)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Extract sequence names from the alignment
seq_names <- names(alignment)

# Check if sequence names are not NULL
if (is.null(seq_names)) stop("Error: Sequence names are NULL. Please check the 'alignment' object.")

# Regular expression to extract the country from the sequence names
country_regex <- "\\[geo_loc_name=([A-Za-z ]+)\\]"

# Extract country names
countries <- str_extract(seq_names, country_regex)

# Clean country names and handle NA values
countries <- str_replace_all(countries, "\\[geo_loc_name=|\\]", "")
countries[is.na(countries)] <- "Unknown"  # Assign "Unknown" to NA values

# Create a dataframe with country names and their frequencies
country_data <- as.data.frame(table(countries))
colnames(country_data) <- c("Country", "Count")  # Rename columns

# Load the world map using rnaturalearth
world <- ne_countries(scale = "medium", returnclass = "sf")

# Merge country frequency data with the world map
world <- left_join(world, country_data, by = c("name" = "Country"))

# Create a map plot showing the geographical distribution
ggplot(data = world) +
  geom_sf(aes(fill = Count), color = "white", size = 0.1) +
  scale_fill_continuous(na.value = "lightgrey", low = "yellow", high = "red", name = "Frequency") +
  theme_minimal() +
  labs(title = "Geographical Distribution of ToBRFV Aligned Isolates",
       x = "Longitude",
       y = "Latitude") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
