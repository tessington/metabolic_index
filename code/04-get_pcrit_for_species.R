#options(echo = FALSE)
rm(list = ls())
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)
library(Matrix)
library(worrms)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# Setup Code ####
## load functions ####
source("code/helper/fit_model_funs.R")

model.fit <- fit_model_augmented_taxa(fitnew = T)

# Do this for several types of scuplin
species.names <- c("Clinocottus globiceps", "Clinocottus analis", "Leptocottus armatus", "Scorpaenichthys marmoratus")

# using joint precision matrix of all parameters
result_df <- map_dfr(species.names, function(taxa) {
  result <- estimate_taxa_full(taxa.name = taxa,
                               w = 10,
                               temperature = 10,
                               method = "routine",
                               rep = model.fit$rep,
                               ParentChild_gz = model.fit$ParentChild_gz)
  tibble(species = taxa,
    log_pcrit = as.numeric(result$log_pcrit["logpcrit"]),
    se_pcrit = as.numeric(result$log_pcrit["se"] )
  )})
print(result_df)

# look at individual results
estimate_taxa_full(taxa.name = species.names[4],
                   w = 10,
                   temperature = 10,
                   method = "routine",
                   rep = model.fit$rep,
                   ParentChild_gz = model.fit$ParentChild_gz )


# using approximation based on conditional SE of random effects parameters (preserving var - covar structure)
result_df <- map_dfr(species.names, function(taxa) {
  result <- estimate_taxa(taxa.name = taxa,
                               w = 10,
                               temperature = 10,
                               method = "routine",
                               rep = model.fit$rep,
                               ParentChild_gz = model.fit$ParentChild_gz)
  tibble(species = taxa,
         log_pcrit = as.numeric(result$log_pcrit["logpcrit"]),
         se_pcrit = as.numeric(result$log_pcrit["se"] )
  )})
print(result_df)

# look at individual results
estimate_taxa(taxa.name = species.names[4],
                   w = 10,
                   temperature = 10,
                   method = "routine",
                   rep = model.fit$rep,
                   ParentChild_gz = model.fit$ParentChild_gz )


