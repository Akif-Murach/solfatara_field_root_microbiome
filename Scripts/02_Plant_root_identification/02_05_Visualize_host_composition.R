# 02_05_Visualize_host_composition.R

library(tidyverse)
library(here)

input<-here("Output", "02_Plant_root_identification", 
            "Metadata", "Plant")
output<-here("Output", "02_Plant_root_identification","Figures")
dir.create(output2,showWarnings = FALSE,recursive = TRUE)


host_data <- read.csv(here(input, "processed_host_metadata.csv"))

colors <- c(
  "#00005B", "#1B9E80", "#7F3C8D", "#1F78B4", "#E31A1C", "#FF7F00",
  "#A6761D", "#1B9E00", "#E45A8F", "#FFDC00", "#7570B3", "#7CAE50",
  "#F564E3", "#FFFF93", "#8DA0CB", "gray50", "gray70")

host_order <- c(
  "Enkianthus campanulatus", "Eubotryoides grayana", "Rhododendron spp.",
  "Vaccinium smallii", "Betula ermanii", "Pinus parviflora")

habitat_order <- c(
  "Arayu-Jigoku_Solfatara field",
  "Arayu-Jigoku_Forest edge",
  "Ofukasawa_Solfatara field",
  "Ofukasawa_Forest edge")

species_counts <- host_data |>
  count(site, habitat, host, name = "Frequency") |>
  mutate(
    site_habitat = factor(paste0(site, "_", habitat), levels = habitat_order),
    host = factor(host, levels = host_order))

p_bar <- ggplot(species_counts, aes(x = site_habitat, y = Frequency, fill = host)) +
  geom_col() +
  scale_fill_manual(name = "Host plant identity", values = colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15, face = "italic")
  )

ggsave(
  filename = here(output,"Plant_composition.pdf"),
  plot = p_bar,
  device = "pdf",  
  dpi = 300,
  width = 8,
  height = 6,
  units = "in")
