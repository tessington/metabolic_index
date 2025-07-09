source("lookup_taxa/lookup_taxa_helper_v2.R")
ylim <- c(0,0.42)
taxa.name <- "Teleostei"
p_teleost <- lookup_taxa(taxa.name, ylim)

taxa.name <- "Cottidae"
p_cottidae <- lookup_taxa(taxa.name, ylim)


taxa.name <- "Malacostraca"
p_malacostraca <- lookup_taxa(taxa.name, ylim)

taxa.name <- "Palinuridae"
p_palinuridae <- lookup_taxa(taxa.name, ylim)

p <- grid.arrange(p_teleost, p_malacostraca, p_cottidae, p_palinuridae, nrow = 2, ncol = 2)

ggsave(filename = "Figures/taxa_predictions.png",
      plot = p,
      width = 7,
      height = 6,
      units = "in")

taxa.name <- "Gadidae"
lookup_taxa(taxa.name, ylim)

taxa.name <- "Perciformes"
lookup_taxa(taxa.name, ylim)

taxa.name <- "Teleostei"
lookup_taxa(taxa.name, ylim)

taxa.name = "Elasmobranchii"
lookup_taxa(taxa.name, ylim)
