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

# Do this for several types of sculplin in family cottidae
species.names <- c("Clinocottus globiceps", "Clinocottus analis", "Astrocottus leprops")

# do this for several types of >>>>
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
result_smr_df <- map_dfr(species.names, function(taxa) {
  result <- estimate_taxa_full(taxa.name = taxa,
                               w = 10,
                               temperature = 10,
                               method = "smr",
                               rep = model.fit$rep,
                               ParentChild_gz = model.fit$ParentChild_gz)
  tibble(species = taxa,
         log_pcrit = as.numeric(result$log_pcrit["logpcrit"]),
         se_pcrit = as.numeric(result$log_pcrit["se"] )
  )})
print(result_smr_df)

