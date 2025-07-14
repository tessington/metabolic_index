rm(list = ls())
library(dplyr)
library(tidyr)
library(purrr)
library(INLA)
library(sdmTMB)
library(ggeffects)
library(purrr)
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

# function to streamline code.  Calculate P(po2 > pcrit) for each row of dataframe
update_df <- function(species_df, taxa.name, w, pcrit_type = "smr", rep, ParentChild_gz) {
# get parameter estimates
  taxa_estimates <- estimate_taxa(taxa.name, 
                                     w = w ,
                                     temperature = 10,
                                     method = pcrit_type,
                                     rep  = model.fit$rep, 
                                     ParentChild_gz = model.fit$ParentChild_gz)
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

# filter data by depth of 97.5 percent of occurences
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

# Function to fit species distribution model ####
fit_model <- function(data.2.use, formula_2_use, includepo2 = T) {
  data.2.use$log_depth_scaled <- scale(log(data.2.use$depth))
  
  #Make mesh
  library(INLA)
  library(sdmTMB)
  bnd <- INLA::inla.nonconvex.hull(cbind(data.2.use$X, data.2.use$Y), 
                                   convex = -0.05)
  inla_mesh <- INLA::inla.mesh.2d(
    boundary = bnd,
    max.edge = c(150, 1000),
    offset = -0.1, # default -0.1
    cutoff = NULL,
    min.angle = 5 # default 21
  )
  spde <- make_mesh(data.2.use, c("X", "Y"), mesh = inla_mesh)
  if (includepo2) {
    # Get the number of fixed effects parameters and set lower bound for the po2 effect (always second)
    nfixed_effects <- 4 + length(unique(data.2.use$survey)) -1 
    lower_bound_array <- rep(-Inf, nfixed_effects)
    upper_bound_array <- rep(Inf, nfixed_effects)
    lower_bound_array[2] <- 0 # Set lower bound for the ppcrit fixed effect to 0
    
    model <- sdmTMB(formula_2_use, 
                    mesh = spde,
                    family = binomial(),
                    data = data.2.use,
                    spatial = "on",
                    anisotropy = TRUE,
                    extra_time = c(2016, 2020),
                    time = "year",
                    spatiotemporal = "iid",
                    control = sdmTMBcontrol(lower = list(b_j = lower_bound_array), 
                                            upper = list(b_j = upper_bound_array), 
                                            newton_loops = 0,
                                            multiphase = TRUE)
                    )
  } else {
    model <- sdmTMB(formula_2_use, 
                    mesh = spde,
                    family = binomial(),
                    data = data.2.use,
                    spatial = "on",
                    anisotropy = TRUE,
                    extra_time = c(2016, 2020),
                    time = "year",
                    spatiotemporal = "iid",
                    control = sdmTMBcontrol(multiphase = TRUE)
                    )
    
  }

  return(model)
}

# Fit models to each species ####
## Pacific Cod ####

formula_2_use <- formula(present_absent ~ s(log_depth_scaled) + scale(temperature_C) + survey)

base_model_pacific_cod <- fit_model(data.2.use = pacific_cod, 
                               formula_2_use = formula_2_use,
                               includepo2 = FALSE)
formula_2_use <- formula(present_absent ~ scale(p_pcrit)  + s(log_depth_scaled)  + scale(temperature_C)  + survey)

ppcrit_model_pacific_cod <- fit_model(data.2.use = pacific_cod,
                                      formula_2_use = formula_2_use)
po2_model_pacific_cod <- fit_model(data.2.use = pacific_cod,
                                      formula_2_use = formula(present_absent ~ scale(po2)  + s(log_depth_scaled) + scale(temperature_C)  + survey))

aic_list_pacific_cod <- c(AIC(base_model_pacific_cod), AIC(ppcrit_model_pacific_cod), AIC(po2_model_pacific_cod))
## Pacific halibut ####

base_model_pacific_halibut <- fit_model(data.2.use = pacific_halibut, 
                                        formula_2_use = formula(present_absent ~ s(log_depth_scaled) + scale(temperature_C) + survey),
                                        includepo2 = F)
ppcrit_model_pacific_halibut <-fit_model(data.2.use = pacific_halibut, 
                                         formula_2_use = formula(present_absent ~ scale(p_pcrit) + s(log_depth_scaled) + scale(temperature_C) + survey),
                                         includepo2 = T)
po2_model_pacific_halibut <- fit_model(data.2.use = pacific_halibut, 
                                    formula_2_use = formula(present_absent ~ scale(po2)  + s(log_depth_scaled) + scale(temperature_C) + survey),
                                    includepo2 = T)
aic_list_pacific_halibut <- c(AIC(base_model_pacific_halibut), AIC(ppcrit_model_pacific_halibut), AIC(po2_model_pacific_halibut))
## Pacific hake ####
base_model_pacific_hake <- fit_model(data.2.use = pacific_hake, 
                                        formula_2_use = formula(present_absent ~ s(log_depth_scaled) + scale(temperature_C) + survey),
                                        includepo2 = F)
ppcrit_model_pacific_hake <-fit_model(data.2.use = pacific_hake, 
                                         formula_2_use = formula(present_absent ~ scale(p_pcrit) + s(log_depth_scaled) + scale(temperature_C) + survey),
                                         includepo2 = T)
po2_model_pacific_hake <- fit_model(data.2.use = pacific_hake, 
                                       formula_2_use = formula(present_absent ~ scale(po2)  + s(log_depth_scaled) + scale(temperature_C) + survey),
                                       includepo2 = T)
aic_list_pacific_hake <- c(AIC(base_model_pacific_hake), AIC(ppcrit_model_pacific_hake), AIC(po2_model_pacific_hake))

## Dover sole ####
base_model_dover_sole <- fit_model(data.2.use = dover_sole, 
                                   formula_2_use = formula(present_absent ~ s(log_depth_scaled) + scale(temperature_C) + survey),
                                   includepo2 = F)

ppcrit_model_dover_sole <- fit_model(data.2.use = dover_sole, 
                                      formula_2_use = formula(present_absent ~ scale(p_pcrit) + s(log_depth_scaled) + scale(temperature_C) + survey),
                                      includepo2 = T)
po2_model_dover_sole <- fit_model(data.2.use = dover_sole, 
                                  formula_2_use = formula(present_absent ~ scale(po2)  + s(log_depth_scaled) + scale(temperature_C) + survey),
                                  includepo2 = T)
aic_list_dover_sole <- c(AIC(base_model_dover_sole), AIC(ppcrit_model_dover_sole), AIC(po2_model_dover_sole))



# make master deltaAIC table ####
aic_table <- matrix(0,nrow = 3, ncol = 4)
colnames( aic_table ) <- c("Pacific cod", "Pacific halibut", "Pacific hake", "Dover sole")
rownames( aic_table ) <- c("base", "ppcrit", "po2")
aic_table[,"Pacific cod"] <- aic_list_pacific_cod - min(aic_list_pacific_cod)
aic_table[,"Pacific halibut"] <- aic_list_pacific_halibut - min(aic_list_pacific_halibut)
aic_table[,"Pacific hake"] <- aic_list_pacific_hake - min(aic_list_pacific_hake)
aic_table[,"Dover sole"] <- aic_list_dover_sole - min( aic_list_dover_sole )
print( aic_table )
