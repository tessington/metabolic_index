library(ape)
library(dplyr)
library(tidyverse)
library(forcats)
source("code/fit_model_funs.R")
source("code/plot.phylo2.funs.R")

all.dat <- load_data()
rem_unknown <- T
all.dat$EstMethod_Metric <- tolower(all.dat$EstMethod_Metric)
if (rem_unknown) all.dat <- dplyr::filter(all.dat, !Method == "unknown", !EstMethod_Metric == "unknown")

rem_lethal <- T
if(rem_lethal) all.dat <- dplyr::filter(all.dat, !Method == "Lethal")
rem_escape <- T
if(rem_escape) all.dat <- dplyr::filter(all.dat, !Method == "Escape")


all.taxa <- dplyr::select(all.dat, Phylum, Class, Order, Family, Genera, Species)
all.taxa <- dplyr::distinct(all.taxa)



col_names <- names(all.taxa)
all.taxa <- all.taxa %>%
  mutate(across(all_of(col_names), as.factor))


frm <- ~Phylum/Class/Order/Family/Species

tr <- as.phylo(frm, data = all.taxa, collapse = F)
tr$edge.length <- rep(1, nrow(tr$edge))
par(mai = c(0.5, 0.5, 0.1, 0.1))
pdf(file = "figures/taxa_plot.pdf",
    height = 12,
    width = 11.5)
#plot(tr, "c", cex = 0.8)
plot.phylo2(tr, "c", cex = 0.875, show.node.label = T)
dev.off()
#plot(tr, "c", cex = 0.8, show.node.label = T)
