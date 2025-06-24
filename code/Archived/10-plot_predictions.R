source("code/fit_model_funs.R")
ylim <- c(0,0.42)
taxa.name <- "Teleostei"
p_teleost <- lookup_taxa(taxa.name, ylim)

taxa.name <- "Gadidae"
p_cottidae <- lookup_taxa(taxa.name, ylim)


taxa.name <- "Salmonidae"
p_malacostraca <- lookup_taxa(taxa.name, ylim)

taxa.name <- "Palinuridae"
p_palinuridae <- lookup_taxa(taxa.name, ylim)

p <- grid.arrange(p_palinuridae , p_cottidae, nrow = 1, ncol = 2)

ggsave(filename = "Figures/taxa_predictions.png",
       plot = p,
       width = 7,
       height = 6,
       units = "in")

