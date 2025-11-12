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
  inla_mesh <- INLA::inla.mesh.2d(boundary = bnd, max.edge = c(150, 1000),
                                  offset = -0.1, cutoff = NULL, min.angle = 5)
  sdmTMB::make_mesh(data.2.use, c("X","Y"), mesh = inla_mesh)
}

fit_model <- function(data.2.use, formula_2_use, spde, modeltype = "base", max_threshold = NULL) {
  print(modeltype)
  # setup sdmTMB control specifications
  if (modeltype == "base") {
    control <- sdmTMBcontrol(multiphase = TRUE) 
  }
  
  if (modeltype == "po2") {
    # Get the number of fixed effects parameters and set lower bound for the po2 effect (always second)
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
    if(!is.null(max_threshold)) {
      upper_bound_array[2] <- max_threshold
    }
    
    control = sdmTMBcontrol(lower = list(b_threshold = lower_bound_array), 
                            upper = list(b_threshold = upper_bound_array), 
                            newton_loops = 0,
                            multiphase = TRUE)
  }
  
#  future::plan(multisession)
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
                                number_of_folds = 4) {
  data.2.use.adjusted <- sf::st_as_sf(data.2.use, coords = c("X", "Y"), crs = NA)
  spatial_blocking <- cv_spatial(
    x = data.2.use.adjusted,
    column = "present_absent",
    size = block_size ,
    k = number_of_folds,
    selection = "random",
    iteration = 50)
  return_data <- data.2.use
  return_data$folds <- spatial_blocking$folds_ids
  return(return_data)
}

# Function to loop through models and save results

fit_all_models <- function(data.2.use, seed) {
  # transform predictors
  data.2.use$log_depth_scaled <- scale(log(data.2.use$depth))
  data.2.use$temp_scaled <- scale(data.2.use$temperature_C)
  data.2.use$temp_scaled_squared <- with(data.2.use, scale(temperature_C^2) )
  data.2.use$po2_scaled <- scale(data.2.use$po2)
  data.2.use$p_pcrit_scaled <- scale(data.2.use$p_pcrit)
  
  # get spatial blocking
  set.seed(seed)
  data.2.use <- make_spatial_blocks(data.2.use, block_size = 65, number_of_folds = 5)
  
  # make mesh
  spde <- build_mesh(data.2.use)
  
  base_model <- fit_model(data.2.use = data.2.use, 
                          formula_2_use = formula(present_absent ~ s(log_depth_scaled) + temp_scaled + temp_scaled_squared + survey),
                          modeltype = "base",
                          spde = spde)
  ppcrit_model <- fit_model(data.2.use = data.2.use,
                            formula_2_use = formula(present_absent ~ breakpt(p_pcrit_scaled)  + s(log_depth_scaled)  +  temp_scaled + temp_scaled_squared  + survey),
                            modeltype = "threshold",
                            spde = spde,
                            max_threshold = max(data.2.use$p_pcrit_scaled)
                            )
  
  po2_model <- fit_model(data.2.use = data.2.use,
                          formula_2_use = formula(present_absent ~ breakpt(po2_scaled)  + s(log_depth_scaled) + temp_scaled + temp_scaled_squared  + survey),
                         modeltype ="threshold",
                         spde = spde,
                         max_threshold = max(data.2.use$po2_scaled))
  
  return(list(base = base_model, ppcrit = ppcrit_model, po2= po2_model, cv_scores = c(base_model$sum_loglik,
                                                                 ppcrit_model$sum_loglik,
                                                                 po2_model$sum_loglik) ) )
  
}
  
  

# Fit models to each species ####
## Pacific cod ####
pacific_cod_fits  <- fit_all_models(data.2.use = pacific_cod, seed = 1311)
## Pacific halibut ####
pacific_halibut_fits <- fit_all_models(data.2.use = pacific_halibut, seed = 1311)
## Pacific hake ####
pacific_hake_fits <- fit_all_models(data.2.use = pacific_hake, seed = 1311)
## Dover sole ####
dover_sole_fits <- fit_all_models(data.2.use = dover_sole, seed = 1311)
## Longspine thornyhead ####
longspine_thornyhead_fits <- fit_all_models(data.2.use = longspine_thornyhead, seed = 1311)
## Canary rockfish ####
canary_rockfish_fits <- fit_all_models(data.2.use = canary_rockfish, seed = 1311)

# make master CV table ####
cv_table <- matrix(0,nrow = 6, ncol = 3)
rownames( cv_table ) <- c("Pacific cod", "Pacific halibut", "Pacific hake", "Dover sole", "Longspine thornyhead", "Canary rockfish")
colnames( cv_table ) <- c("base", "ppcrit", "po2")
cv_table["Pacific cod",] <- pacific_cod_fits$cv_scores 
cv_table["Pacific halibut",] <-pacific_halibut_fits$cv_scores
cv_table["Pacific hake",] <- pacific_hake_fits$cv_scores
cv_table["Dover sole",] <- dover_sole_fits$cv_scores
cv_table["Longspine thornyhead",] <- longspine_thornyhead_fits$cv_scores
cv_table["Canary rockfish",] <- canary_rockfish_fits$cv_scores
print( cv_table )
