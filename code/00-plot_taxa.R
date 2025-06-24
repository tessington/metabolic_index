library(ape)
library(dplyr)
library(tidyverse)
library(forcats)
source("code/helper/fit_model_funs.R")
source("code/helper/plot.phylo2.funs.R")

all.dat <- load_data()
all.dat <- filter_data(all.dat)

all.taxa <- dplyr::select(all.dat, Phylum, Class, Order, Family, Genera, Species)
all.taxa <- dplyr::distinct(all.taxa)

# for each taxa, figure out if used oxyconform, SMR, or both
get_pcrit_method <- function(species, all.data) {
  species_data <- dplyr::filter(all.data, Species == species)
  unique_methods <- unique(species_data$EstMethod_Metric)
  if (length(unique(unique_methods)) == 1) return_method = unique_methods
  if (length(unique(unique_methods)) > 1) return_method = "both"
  return(return_method)
}



col_names <- names(all.taxa)
all.taxa <- all.taxa %>%
  mutate(across(all_of(col_names), as.factor))


frm <- ~Phylum/Class/Order/Family/Species

tr <- as.phylo(frm, data = all.taxa, use.labels = TRUE, collapse = FALSE)
ladderized_tree <- ladderize(tr)

tr$edge.length <- rep(1, nrow(tr$edge))


# Get methods for each tip label
tiplabels = tr$tip.label
method_by_species <- sapply(X = tiplabels, FUN = get_pcrit_method, all.data = all.dat)

# Set colors by species
color_by_species <- rep(NA, length(method_by_species))
color_by_species[method_by_species == "oxyconform"] <- "#fc8d62"
color_by_species[method_by_species == "smr"] <- "#8da0cb"
color_by_species[method_by_species == "both"] <- "#1b9e77"



par(mai = c(0.5, 0.5, 0.1, 0.1))
pdf(file = "figures/taxa_plot3.pdf",
    height = 7,
    width = 7)

plot.phylo2(tr, "fan", cex = 0.8, show.node.label = T, use.edge.length = T, show.tip.label =F, edge.width = 2, adj = 1)
#nodelabels()
tiplabels(pch = 21, col = "black", bg = color_by_species, cex = 1.5)
dev.off()

