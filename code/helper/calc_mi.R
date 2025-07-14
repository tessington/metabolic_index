# mock code for calculating 95 percent CI of metabolic index


# Code to generate MI estimate (confidence interval) given pO2, temp, and other info
calc_mi <- function(pO2, logw, inv_temperature, betas, var_covar, method = "smr", confidence_level = 0.95) {
  # code to calculate 95 % confidence interval of metabolic index given pO2 and temperature
  # logw is log(w / wref)
  # inv_temperature is (1 / kb) ( )1/T - 1 / Tref), where t is in kelvin and kb is boltzmann's constant
  if(method == "smr") x_predict <- as.vector(c(-1, logw, inv_temperature, -1))
  if(!method == "smr") x_predict <- c(-1, logw, inv_temperature, 0)
  
  log_mi_predict <- x_predict %*% betas + log(pO2)
  log_mi_se <-sqrt( t(x_predict) %*% var_covar %*% x_predict )
  log_lower_bound <- log_mi_predict - log_mi_se * qnorm(0.5 + confidence_level / 2) 
  log_upper_bound <- log_mi_predict + log_mi_se * qnorm(0.5 + confidence_level / 2)
  return(list ( mi = exp(log_mi_predict), lower_bound = exp(log_lower_bound), upper_bound = exp(log_upper_bound) ) )
}


# to use this
# calculate parameter values and variance - covariance
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)
library(Matrix)
library(worrms)
library(purrr)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# Setup Code ####
## load functions ####
source("code/helper/fit_model_funs.R")

taxa.name <- "Oncorhynchus tshawytscha"
model.fit <- fit_model_augmented_taxa(fitnew = F)
model_output <- estimate_taxa_full(taxa.name,
                                   w = 4000,
                                   temperature = 10,
                                   method = "routine",
                                   rep = model.fit$rep,
                                   ParentChild_gz = model.fit$ParentChild_gz)

betas <- model_output$parameters[,1]
print(betas)
print(model_output$var_covar)
calc_mi(5, log(500 / 10), inv_temperature = 0.75,  betas = betas, var_covar = taxa_estimates$var_covar, method = "smr", confidence_level = 0.95)
