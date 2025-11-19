rm(list = ls())
library(dplyr)
library(tidyr)
library(purrr)
library(INLA)
library(sdmTMB)
library(ggeffects)
library(purrr)
library(TMB)
library(worrms)
library(blockCV)
library(future)

#Load functions
source("code/helper/fit_model_funs.R")

#Load survey data data
files <- list.files(path = "data/survey_data", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

# keep only US WC and British Columbia
dat <- dplyr::filter(dat, region %in% c("cc", "bc") )

# create variable called catch_rate and assign values depending on survey
dat$catch_rate <- dat$catch_weight
dat$catch_rate[dat$survey == "iphc"] <- dat$cpue_count[dat$survey == "iphc"]

#remove iphc
remove_iphc <- F
if(remove_iphc){
  dat <- filter(dat, survey!="iphc")
}

#Remove other missing rows
dat <- dat  %>%
  drop_na(depth,year, po2, temperature_C, catch_rate, X, Y)


#Remove weird depths
dat <- filter(dat, depth>0)

# convert oxgyen to kPa
#Calculate o2 in kPa (data is in umol kg)
gas_const = 8.31
partial_molar_vol = 0.000032
boltz = 0.000086173324

SA = gsw_SA_from_SP(dat$salinity_psu,dat$depth,dat$longitude,dat$latitude) #absolute salinity for pot T calc
pt = gsw_pt_from_t(SA,dat$temp,dat$depth) #potential temp at a particular depth
O2_Sat0 = gsw_O2sol_SP_pt(dat$salinity_psu,pt)

press = exp(dat$depth*10000*partial_molar_vol/gas_const/kelvin(dat$temp))
O2_satdepth = O2_Sat0*press

#solubility at p=0
sol0 = O2_Sat0/0.209
sol_Dep = sol0*press
dat$po2 = dat$O2_umolkg/sol_Dep
dat$po2 <- dat$po2 * 101.325 # convert to kPa

# make present - absence
dat$present_absent <- 0
dat$present_absent[dat$catch_rate>0] <- 1

# make dataframes for each species
dover_sole <- dplyr::filter(dat, common_name == "dover sole")
pacific_hake <- dplyr::filter(dat, common_name == "pacific hake")
pacific_cod <- dplyr::filter(dat, common_name == "pacific cod")
pacific_halibut <- dplyr::filter(dat, common_name == "pacific halibut")
longspine_thornyhead <- dplyr::filter(dat, common_name == "longspine thornyhead")
canary_rockfish <- dplyr::filter(dat, common_name == "canary rockfish")
# function to streamline code.  Calculate P(po2 > pcrit) for each row of dataframe
update_df <- function(species_df, taxa.name, w, pcrit_type = "smr", rep, ParentChild_gz) {
# get parameter estimates
  taxa_estimates <- estimate_taxa(taxa.name, 
                                     w = w ,
                                     temperature = 10,
                                     method = pcrit_type)
  # use pmap to run calc_p_po2
  betas <- taxa_estimates$parameters[,1]
  var_covar <- taxa_estimates$var_covar
  species_df$p_pcrit <- pmap_dbl(list(po2 = species_df$po2, temp = species_df$temperature_C),
                                ~calc_p_po2(po2 = ..1, temp = ..2, w = w,  pcrit_type = pcrit_type, betas = betas, var_covar = var_covar)
                                )

return(species_df)
}

# get fitted model
model.fit <- fit_model_augmented_taxa(fitnew = F)

# calculate probability that each location po2 exceeds pcrit for each species
pcrit_type = "smr"
W = 1000

# Pacific cod
taxa.name <- "Gadus macrocephalus"
pacific_cod <- update_df(pacific_cod, taxa.name, pcrit_type = pcrit_type, w = W, model.fit$rep, model.fit$ParentChild_gz)
pacific_cod$cpue_weight <- with(pacific_cod, catch_weight / effort)

# Pacific Hake
taxa.name <- "Merluccius productus"
pacific_hake <- update_df(pacific_hake, taxa.name,pcrit_type = pcrit_type,w = W, model.fit$rep, model.fit$ParentChild_gz)

# Pacific halibut
taxa.name <- "Hippoglossus stenolepis"
pacific_halibut <- update_df(pacific_halibut, taxa.name, pcrit_type = pcrit_type,w = W, model.fit$rep, model.fit$ParentChild_gz)
pacific_halibut$cpue_weight <- with(pacific_halibut, catch_weight / effort)

# Dover sole
taxa.name <- "Microstomus pacificus"
dover_sole <- update_df(dover_sole, taxa.name, pcrit_type = pcrit_type, w = W,model.fit$rep, model.fit$ParentChild_gz)
dover_sole$cpue_weight <- with(dover_sole, catch_weight / effort)

# Longspine thornyhead
taxa.name <- "Sebastolobus altivelis"
longspine_thornyhead <- update_df(longspine_thornyhead, taxa.name, pcrit_type = pcrit_type, w = W,model.fit$rep, model.fit$ParentChild_gz)

taxa.name <- "Sebastes pinniger"
canary_rockfish <- update_df(canary_rockfish, taxa.name, pcrit_type = pcrit_type, w = W,model.fit$rep, model.fit$ParentChild_gz)

# filter data by depth of threshold percent of occurrences
# function
trim_depths <- function(data.2.use, threshold) {
  number_occurences <- sum(data.2.use$present_absent)
  sorted_data <- sort_by(data.2.use, ~ depth)
  sorted_data <- sorted_data %>% mutate(cumsum_occurence = cumsum(present_absent)) 
  return.data <- sorted_data[sorted_data$cumsum_occurence/number_occurences <= threshold,]
  return(return.data)
}
# trim data by depth ####
depth_threshold <- 0.975
pacific_cod <- trim_depths(data.2.use = pacific_cod,
                           threshold = depth_threshold)
pacific_halibut <- trim_depths(data.2.use = pacific_halibut,
                               threshold = depth_threshold)
pacific_hake <- trim_depths(data.2.use = pacific_hake,
                               threshold = depth_threshold)
dover_sole <- trim_depths(data.2.use = dover_sole,
                               threshold = depth_threshold)
longspine_thornyhead <- trim_depths(data.2.use = longspine_thornyhead,
                          threshold = depth_threshold)
canary_rockfish <- trim_depths(data.2.use = canary_rockfish,
                          threshold = depth_threshold)

# Functions to fit species distribution model via cross validation ####
build_mesh <- function(data.2.use, cutoff = 10) {
  bnd <- INLA::inla.nonconvex.hull(cbind(data.2.use$X, data.2.use$Y), convex = -0.05)
  inla_mesh <- INLA::inla.mesh.2d(boundary = bnd, max.edge = c(500, 1000),
                                  offset = -0.1, cutoff = NULL, min.angle = 5)
  sdmTMB::make_mesh(data.2.use, c("X","Y"), mesh = inla_mesh)
}

fit_model <- function(data.2.use, formula_string, spde, modeltype = "base", min_threshold = NULL) {
  
  # make sure formula doesn't carry huge environment
  formula_2_use <-  stats::as.formula(formula_string, env = globalenv())
  summary(formula_2_use)
  
  # setup sdmTMB control specifications
  if (modeltype == "base") {
    control <- sdmTMBcontrol(multiphase = TRUE) 
  }
  
  if (modeltype == "linear") {
    # Get the number of fixed effects parameters and set lower bound for the linear effects of oxygen (always second)
    nfixed_effects <- 4 + length(unique(data.2.use$survey)) -1 
    lower_bound_array <- rep(-Inf, nfixed_effects)
    upper_bound_array <- rep(Inf, nfixed_effects)
    lower_bound_array[2] <- 0 # Set lower bound for the ppcrit fixed effect to 0
    control = sdmTMBcontrol(lower = list(b_j = lower_bound_array), 
                            upper = list(b_j = upper_bound_array), 
                            newton_loops = 0,
                            multiphase = TRUE)
  } 
 
  if (modeltype == "threshold") {
    lower_bound_array = rep(-Inf, 2)
    upper_bound_array <- rep(Inf, 2)
    lower_bound_array[1] <- 0
    if(!is.null(min_threshold)) {
      lower_bound_array[2] <- min_threshold
    }
    
    control = sdmTMBcontrol(lower = list(b_threshold = lower_bound_array), 
                            upper = list(b_threshold = upper_bound_array), 
                            newton_loops = 0,
                            multiphase = TRUE)
  }
  
  future:: plan(multisession, workers = 6)
  cv_out <- sdmTMB_cv(
    formula = formula_2_use,
    data = data.2.use,
    family = binomial(),
    spatial = "on",
    anisotropy = FALSE,
    time = "year",
    spatiotemporal = "iid",
    extra_time = c(2016, 2020),
    control = control,
    fold_ids = data.2.use$fold,
    mesh = spde
  )

  return(cv_out)
}

make_spatial_blocks <- function(data.2.use, 
                                block_size = 100,
                                number_of_folds = 4,
                                type = "cells") {
  if (type == "cells") {
  data.2.use.adjusted <- sf::st_as_sf(data.2.use, coords = c("X", "Y"), crs = NA)
  spatial_blocking <- cv_spatial(
    x = data.2.use.adjusted,
    column = "present_absent",
    size = block_size ,
    k = number_of_folds,
    selection = "random",
    iteration = 50)
  }
  
  if (type == "rows") {
    data.vertical.range <- max(data.2.use$Y) - min(data.2.use$Y)
    target_width <- block_size
    
    number_rows <- ceiling(data.vertical.range / target_width)
    
    data.2.use.adjusted <- sf::st_as_sf(data.2.use, coords = c("X", "Y"), crs = NA)
    spatial_blocking <- cv_spatial(
      x = data.2.use.adjusted,
      column = "present_absent",
      rows_cols = c(number_of_folds,1),
      k = number_of_folds,
      hexagon = FALSE,
      selection = "systematic",
      extend = 0.5)
    
  }
  return_data <- data.2.use
  return_data$folds <- spatial_blocking$folds_ids
  return(return_data)
}

# Function to loop through models and save results

fit_all_models <- function(data.2.use, seed) {
  nfolds <- 5
  block_size <- 65
  # transform predictors
  data.2.use$log_depth_scaled <- scale(log(data.2.use$depth))
  data.2.use$temp_scaled <- scale(data.2.use$temperature_C)
  data.2.use$temp_scaled_squared <- with(data.2.use, scale(temperature_C^2) )
  data.2.use$po2_scaled <- scale(data.2.use$po2)
  data.2.use$p_pcrit_scaled <- scale(data.2.use$p_pcrit)
  
  # get spatial blocking
  set.seed(seed)
  data.2.use <- make_spatial_blocks(data.2.use, block_size = block_size, number_of_folds = nfolds, type = "cells")
  
  # make mesh
  #spde <- build_mesh(data.2.use)
  spde <- make_mesh(data.2.use, c("X", "Y"), n_knots = 250, type = "kmeans")
  
  # Run models
  base_model <- fit_model(data.2.use = data.2.use, 
                          formula_string = "present_absent ~ s(log_depth_scaled) + temp_scaled + temp_scaled_squared + survey",
                          modeltype = "base",
                          spde = spde)
  ppcrit_brkpoint <- fit_model(data.2.use = data.2.use,
                            formula_string = "present_absent ~ breakpt(p_pcrit_scaled)  + s(log_depth_scaled)  +  temp_scaled + temp_scaled_squared  + survey",
                            modeltype = "threshold",
                            spde = spde,
                            min_threshold = NULL# min(data.2.use$p_pcrit_scaled)
                            )
  ppcrit_linear <- fit_model(data.2.use = data.2.use,
                                   formula_string = "present_absent ~ p_pcrit_scaled  + s(log_depth_scaled)  +  temp_scaled + temp_scaled_squared  + survey",
                                   modeltype = "linear",
                                   spde = spde,
                                   min_threshold = NULL# min(data.2.use$p_pcrit_scaled)
  )
  po2_brkpoint <- fit_model(data.2.use = data.2.use,
                          formula_string = "present_absent ~ breakpt(po2_scaled)  + s(log_depth_scaled) + temp_scaled + temp_scaled_squared  + survey",
                         modeltype ="threshold",
                         spde = spde,
                         min_threshold = NULL# min(data.2.use$po2_scaled)
                         )
  po2_linear <- fit_model(data.2.use = data.2.use,
                         formula_string = "present_absent ~ (po2_scaled)  + s(log_depth_scaled) + temp_scaled + temp_scaled_squared  + survey",
                         modeltype ="linear",
                         spde = spde,
                         min_threshold = NULL# min(data.2.use$po2_scaled)
  )
  
  # calculate standard error of the cross validation score for model selection - following Yates et al. 2022
  
  #fold_loglik <- matrix(NA, nrow = nfolds, ncol = 3)
  ##fold_loglik[,1] <- base_model$fold_loglik
  #fold_loglik[,2] <- ppcrit_model$fold_loglik
  #fold_loglik[,3] <- po2_model$fold_loglik
  
  all_cv <- c(base_model$sum_loglik,
              ppcrit_brkpoint$sum_loglik,
              ppcrit_linear$sum_loglik,
              po2_brkpoint$sum_loglik,
              po2_linear$sum_loglik
              )
  
  best_model <- which(all_cv == max(all_cv))
  
  #sigma_best <- sd(fold_loglik[,best_model]) / sqrt(nfolds)
  
  #rho_best <- rep(0, 3)
  #for (i in 1:3) rho_best[i] <- cor(fold_loglik[,i], fold_loglik[, best_model])
  
  #sigma_m <- apply(fold_loglik, FUN = sd, MARGIN = 2) / sqrt(nfolds)
  
  #sigmas <- sigma_best * sqrt(1 - rho_best)
  
  Sm <- max(all_cv) - all_cv # pairwise differences
  #sigma_difference <- rep(0, 3)
  #for (i in 1:3) sigma_difference[i] <- sqrt(sigma_m[i]^2 + sigma_best^2 - 2 * rho_best[i] * sigma_m[i] * sigma_best)
  #sigma_difference[best_model] <- 0 # by definition
  
  return(list(base = base_model, ppcrit_brkpoint = ppcrit_brkpoint, ppcrit_linear = ppcrit_linear, po2_brkpoint= po2_brkpoint, po2_linear = po2_linear, cv_scores = all_cv) )
  
}
calc_cv_scores <- function(model_output, number_reps) {
  n_models <- 5
  # get table of CV scores across reps
  cv_scores <- matrix(NA, nrow = number_reps, ncol = n_models)
  for (i in 1:number_reps) cv_scores[i,1:n_models] <- model_output[[i]]$cv_scores
  best_model <- which(colMeans(cv_scores) == max(colMeans(cv_scores)) )
  sigma_best <- sd(cv_scores[,best_model]) / sqrt(number_reps)
  rho_best <- rep(0, n_models)
  for (i in 1:n_models) rho_best[i] <- cor(cv_scores[,i], cv_scores[, best_model])
  sigma_m <- apply(cv_scores, FUN = sd, MARGIN = 2) / sqrt(number_reps)
  sigmas <- sigma_best * sqrt(1 - rho_best)
  Sm <- max(colMeans(cv_scores)) - colMeans(cv_scores) # pairwise differences
  sigma_difference <- rep(0, n_models)
  for (i in 1:n_models) sigma_difference[i] <- sqrt(sigma_m[i]^2 + sigma_best^2 - 2 * rho_best[i] * sigma_m[i] * sigma_best)
  sigma_difference[best_model] <- 0 # by definition
  return(list(Sm = Sm, sigma = sigma_difference) ) 
}
run_replicate_cvs <- function(data.2.use, number_reps, seed) {
  model_fits <- lapply(seq_len(number_reps), function(i) {
    fit_all_models(data.2.use = data.2.use, seed = seed + i)
  })
  names(model_fits) <- paste0("rep_", seq_len(number_reps) )
  model_fits$CV_summary <- calc_cv_scores(model_fits, number_reps)
  return(model_fits)
}

# Fit models to each species ####
seed <- 707
number_reps <- 5

## Pacific cod ####
cat("Pacific Cod")
pacific_cod_fits <- run_replicate_cvs(pacific_cod, number_reps, seed)
## Pacific halibut ####
cat("Pacific Halibut")
pacific_halibut_fits <- run_replicate_cvs(data.2.use = pacific_halibut, number_reps = number_reps, seed = seed)
## Pacific hake ####
cat("Pacific Hake")
pacific_hake_fits <- run_replicate_cvs(data.2.use = pacific_hake, number_reps= number_reps, seed = seed)
## Dover sole ####
cat("Dover sole")
dover_sole_fits <- run_replicate_cvs(data.2.use = dover_sole, number_reps= number_reps, seed = seed)
## Longspine thornyhead ####
cat("Longspine thornyhead")
longspine_thornyhead_fits <- run_replicate_cvs(data.2.use = longspine_thornyhead, number_reps= number_reps, seed = seed)
## Canary rockfish ####
cat("Canary rockfish")
canary_rockfish_fits <- run_replicate_cvs(data.2.use = canary_rockfish, number_reps= number_reps, seed = seed)

# make master CV table ####
cv_table <- matrix(0,nrow = 6, ncol = 5)
rownames( cv_table ) <- c("Pacific cod", "Pacific halibut", "Canary rockfish", "Pacific hake", "Dover sole", "Longspine thornyhead")
colnames( cv_table ) <- c("base", "ppcrit-breakpoint", "ppcrit-linear",  "po2-breakpoint", "po2-linear")
cv_table["Pacific cod",] <- pacific_cod_fits$CV_summary$Sm
cv_table["Pacific halibut",] <-pacific_halibut_fits$CV_summary$Sm
cv_table["Canary rockfish",] <- canary_rockfish_fits$CV_summary$Sm
cv_table["Pacific hake",] <- pacific_hake_fits$CV_summary$Sm
cv_table["Dover sole",] <- dover_sole_fits$CV_summary$Sm
cv_table["Longspine thornyhead",] <- longspine_thornyhead_fits$CV_summary$Sm
print( cv_table )

sigma_table <- matrix(0,nrow = 6, ncol = 5)
rownames( sigma_table ) <- c("Pacific cod", "Pacific halibut", "Canary rockfish", "Pacific hake", "Dover sole", "Longspine thornyhead")
colnames( sigma_table ) <- c("base", "ppcrit-breakpoint", "ppcrit-linear",  "po2-breakpoint", "po2-linear")
sigma_table["Pacific cod",] <- pacific_cod_fits$CV_summary$sigma
sigma_table["Pacific halibut",] <-pacific_halibut_fits$CV_summary$sigma
sigma_table["Canary rockfish",] <- canary_rockfish_fits$CV_summary$sigma
sigma_table["Pacific hake",] <- pacific_hake_fits$CV_summary$sigma
sigma_table["Dover sole",] <- dover_sole_fits$CV_summary$sigma
sigma_table["Longspine thornyhead",] <- longspine_thornyhead_fits$CV_summary$sigma
print( sigma_table )

get_date <- Sys.Date()
save_filename <- paste0("analysis/cross_validation_results_", get_date)
save_list <- list(Sm = cv_table, se = sigma_table)
saveRDS(save_list, file = save_filename)
