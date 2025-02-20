rm(list = ls())
# Code fit hierarchical model to Arrenhius Equation
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)
library(Matrix)
library(mvtnorm)
library(tmbstan)
library(shinystan)
library(egg)
library(cowplot)

# Setup Code ####
## ggplot setup ####
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             axis.text = element_text(color = "black") 
)

## load functions ####
source("code/fit_model_funs.R")
source("code/calculate_EDF_fn.R")

## load data ####
all.dat <- load_data()

# Switch to remove outlier  Styela Plicata
rem_styela <- T
if (rem_styela) all.dat <- all.dat %>%
  filter(Species != "Styela plicata")

# Remove all data where the estimation method could not be determined
rem_unknown <- T
all.dat$EstMethod_Metric <- tolower(all.dat$EstMethod_Metric)
if (rem_unknown) all.dat <- dplyr::filter(all.dat, !Method == "unknown", !EstMethod_Metric == "unknown")

# Keep only oxygen consumption methods
all.dat <- dplyr::filter(all.dat, Method == "OxygenConsumption")

# Create an array M, 1 or 0 depnding on method used
method_mat <- model.matrix(~ EstMethod_Metric , all.dat)
# number of method effect sizes that must be estimated
n_methods <- ncol(method_mat) -1

## Summarise data by taxa ####
class_summary <- all.dat %>%
  group_by(Class) %>%
  summarise(NoOrder = length(unique(Order)), NoFamily = length(unique(Family)),  NoSpecies = length(unique(Species)))

order_summary <- all.dat %>%
  group_by(Order, Class) %>%
  summarise(NoFamily = length(unique(Family)), NoSpecies = length(unique(Species)))

family_summary <- all.dat %>%
  group_by(Family, Order) %>%
  summarise(NoGenus = length(unique(Genera)), NoSpecies= length(unique(Species)))

genera_summary <- all.dat %>%
  group_by(Family, Order, Genera) %>%
  summarise(NoSpecies= length(unique(Species)))


## Sumarise data by source
source_summary <- all.dat %>%
  group_by(Source) %>%
  summarise(NoOrder = length(unique(Order)), NoFamily = length(unique(Family)), NoSpecies = length(unique(Species)))

sh_source_summary <- all.dat %>%
  group_by(SharedAuthor) %>%
  summarise(NoOrder = length(unique(Order)), NoFamily = length(unique(Family)), NoSpecies = length(unique(Species)))

## Setup data for TMB ####

kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W/wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$Pcrit_atm<- all.dat$Pcrit / 101.325 # convert from kPa to atm
all.dat$minuslogpo2 <- - log(all.dat$Pcrit)
taxa.list <- c("Class", "Order", "Family", "Genera", "Species")
all.dat$SourceNo <- as.numeric(as.factor(all.dat$Source))
all.dat$TeamNo <- as.numeric(as.factor(all.dat$SharedAuthor))
### Create new ParentChild matrix for reduced taxonomic structure ####
taxa.info <- make_taxa_tree(all.dat, taxa.list)
ParentChild_gz <- taxa.info$ParentChild_gz
PC_gz <- taxa.info$PC_gz
g_i <- taxa.info$g_i
g_i_i <- taxa.info$g_i_i
n_k <- taxa.info$n_k
n_j <- taxa.info$n_j
n_g <- taxa.info$n_g
n_i <- taxa.info$n_i
spc_in_PC_gz <- taxa.info$spc_in_PC_gz

# Fit base model ####
## Setup TMB ####
data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = all.dat$minuslogpo2,
             spc_in_PCgz = spc_in_PC_gz -1
)

parameters = list(alpha_j = rep(0,n_j),
                  L_z = rep(1, 6),
                  log_lambda = rep(-1, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  logsigma = 0
)

Random <- c("beta_gj")
model <- "hierarchical_mi_base"
compile(paste0("code/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("code/TMB/",model)))

## Run TUMB ####
obj_base <-
  MakeADFun(
    data = data,
    parameters = parameters,
    DLL = model,
    random = Random,
    silent = TRUE
  )
opt_base <- nlminb(obj_base$par, obj_base$fn, obj_base$gr)
rep = sdreport( obj_base,
                getReportCovariance = TRUE, 
                getJointPrecision=TRUE)


## Calculate cAIC
EDF <- calculate_EDF( obj_base,
                      opt_base,
                      nonvariance_fixed_effects = "",
                      prediction_name = "mu",
                      data_name = "minuslogpo2",
                      delta = 0.01,
                      show_progress = F,
                      refit = "random"
                      )

nll_data <- obj_base$report()[["nll_data"]]
cAIC_model1 <- 2 * nll_data + 2 * EDF

## Extract Fitted Parameters ####
re <- summary(rep, "random")
fixef <- summary(rep, "fixed")
beta_mle <- matrix(re[grep(rownames(re), pattern = "beta_gj"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta_gj"),2], nrow = n_g, ncol = 3, byrow = F)


## Summarize Estimates ####
sum_est <- summarize_estimates(beta_mle, beta_se, ParentChild_gz, taxa.list)

## Plot Estimates ####
plot_est <- F
if (plot_est) {
  plot_by_group(sum_est)
}

SpeciesEst <- sum_est$SpeciesEst
SpeciesEst$V = exp(SpeciesEst$logV)
## Plot model diagnostics ####
model_diagnostics <- T
if (model_diagnostics) {
  p_diagnostic <- plot_diagnostics(model= "base", 
                                   Pcrit = all.dat$Pcrit, 
                                   inv.temp = all.dat$inv.temp, 
                                   W = all.dat$W, 
                                   SpeciesEst = SpeciesEst)
  

}

# Fit model with method effects ####
## Setup TMB ####
data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = all.dat$minuslogpo2,
             spc_in_PCgz = spc_in_PC_gz -1,
             method_mat = method_mat[,-1]
             )

parameters = list(alpha_j = rep(0,n_j),
                  L_z = rep(1, 6),
                  log_lambda = rep(-1, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  beta_method = 0,
                  logsigma = 0
)

Random <- c("beta_gj")
model <- "hierarchical_mi"
compile(paste0("code/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("code/TMB/",model)))

## Run TUMB ####
obj_method <-
  MakeADFun(
    data = data,
    parameters = parameters,
    DLL = model,
    random = Random,
    silent = TRUE
  )
opt_method <- nlminb(obj_method$par, obj_method$fn, obj_method$gr)
rep = sdreport( obj_method,
                getReportCovariance = TRUE, 
                getJointPrecision=TRUE)


## Calculate cAIC
EDF <- calculate_EDF( obj= obj_method,
                      opt = opt_method,
                      nonvariance_fixed_effects = c("beta_method", "alpha_j", "L_z", "log_lambda"),
                      prediction_name = "mu",
                      data_name = "minuslogpo2",
                      delta = 0.01,
                      show_progress = F,
                      refit = "random"
)

nll_data <- obj_method$report()[["nll_data"]]
cAIC_model2 <- 2 * nll_data + 2 * EDF

## Extract Fitted Parameters ####
re <- summary(rep, "random")
fixef <- summary(rep, "fixed")
beta_mle <- matrix(re[grep(rownames(re), pattern = "beta_gj"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta_gj"),2], nrow = n_g, ncol = 3, byrow = F)
beta_method <- matrix(fixef[grep(rownames(fixef), pattern = "beta_method"),1], nrow = n_methods)
beta_method_se <- matrix(fixef[grep(rownames(fixef), pattern = "beta_method"),2], nrow = n_methods)

## Summarize Estimates ####
sum_est <- summarize_estimates(beta_mle, beta_se, ParentChild_gz, taxa.list)

## Plot Estimates ####
plot_est <- F
if (plot_est) {
  plot_by_group(sum_est)
}

SpeciesEst <- sum_est$SpeciesEst
SpeciesEst$V = exp(SpeciesEst$logV)
## Plot model diagnostics ####
model_diagnostics <- T
if (model_diagnostics) {
  p_diagnostic <- plot_diagnostics(model= "method", 
                                   Pcrit = all.dat$Pcrit, 
                                   inv.temp = all.dat$inv.temp, 
                                   W = all.dat$W, 
                                   SpeciesEst = SpeciesEst,
                                   beta_method = beta_method,
                                   method_mat = method_mat)
  
  
}


# Fit latent paper effects ####
## Setup TMB ####
n_sources <- length(unique(all.dat$Source))
data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = all.dat$minuslogpo2,
             spc_in_PCgz = spc_in_PC_gz -1,
             method_mat = method_mat[,-1],
             sourceno = all.dat$SourceNo -1
)

parameters = list(alpha_j = rep(0,n_j),
                  L_z = rep(1, 6),
                  log_lambda = rep(-1, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  beta_method = 0,
                  beta_source = rep(0, times = n_sources),
                  logsigma = 0,
                  logsigmap = 0
)
Random <- c("beta_gj", "beta_source")
model <- "hierarchical_mi_latent"
compile(paste0("code/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("code/TMB/",model)))

## Run TUMB ####
obj_source <-
  MakeADFun(
    data = data,
    parameters = parameters,
    DLL = model,
    random = Random,
    silent = TRUE
  )
opt_source <- nlminb(obj_source$par, obj_source$fn, obj_source$gr)
rep = sdreport( obj_source,
                getReportCovariance = TRUE, 
                getJointPrecision=TRUE)

# get EDF
EDF <- calculate_EDF(obj_source, opt_source)
nll_data <- obj_source$report()[["nll_data"]]
cAIC_model3 <- 2 * nll_data + 2 * EDF


re <- summary(rep, "random")
fixef <- summary(rep, "fixed")
beta_mle <- matrix(re[grep(rownames(re), pattern = "beta_gj"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta_gj"),2], nrow = n_g, ncol = 3, byrow = F)
beta_method <- matrix(fixef[grep(rownames(fixef), pattern = "beta_method"),1], nrow = n_methods)
beta_source <- matrix(re[grep(rownames(re), pattern = "beta_source"),1:2], nrow = n_sources, ncol = 2)

## Summarize Estimatess ####
sum_est <- summarize_estimates(beta_mle, beta_se, ParentChild_gz, taxa.list)

## Plot Estimates ####
plot_est <- F
if (plot_est) {
  plot_by_group(sum_est)
}

## Plot diagnostics ####
model_diagnostics <- T
if (model_diagnostics) {
  p_diagnostic <- plot_diagnostics(model= "paper", 
                                   Pcrit = all.dat$Pcrit, 
                                   inv.temp = all.dat$inv.temp, 
                                   W = all.dat$W, 
                                   SpeciesEst = sum_est$SpeciesEst,
                                   beta_method = beta_method,
                                   beta_source = beta_source,
                                   source_id = all.dat$SourceNo,
                                   method_mat = method_mat)
  
}

# Fit Latent effects by research group ####
n_sources <- length(unique(all.dat$SharedAuthor))
data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = all.dat$minuslogpo2,
             spc_in_PCgz = spc_in_PC_gz -1,
             method_mat = method_mat[,-1],
             sourceno = all.dat$TeamNo -1
)

parameters = list(alpha_j = rep(0,n_j),
                  L_z = rep(1, 6),
                  log_lambda = rep(-1, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  beta_method = 0,
                  beta_source = rep(0, times = n_sources),
                  logsigma = 0,
                  logsigmap = 0
)
Random <- c("beta_gj", "beta_source")
model <- "hierarchical_mi_latent"
compile(paste0("code/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("code/TMB/",model)))

## Run TUMB ####
obj_group <-
  MakeADFun(
    data = data,
    parameters = parameters,
    DLL = model,
    random = Random,
    silent = TRUE
  )
opt_group <- nlminb(obj_group$par, obj_group$fn, obj_group$gr)
rep = sdreport( obj_group,
                getReportCovariance = TRUE, 
                getJointPrecision=TRUE)
EDF <- calculate_EDF(obj_group, opt_group)
nll_data <- obj_group$report()[["nll_data"]]
cAIC_model4 <- 2 * nll_data + 2 * EDF


re <- summary(rep, "random")
fixef <- summary(rep, "fixed")
beta_mle <- matrix(re[grep(rownames(re), pattern = "beta_gj"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta_gj"),2], nrow = n_g, ncol = 3, byrow = F)
beta_method <- matrix(fixef[grep(rownames(fixef), pattern = "beta_method"),1], nrow = n_methods)
beta_source <- matrix(re[grep(rownames(re), pattern = "beta_source"),1:2], nrow = n_sources, ncol = 2)

## Summarize Estimatess ####
sum_est <- summarize_estimates(beta_mle, beta_se, ParentChild_gz, taxa.list)

## Plot Estimates ####
plot_est <- F
if (plot_est) {
  plot_by_group(sum_est)
}

SpeciesEst = sum_est$SpeciesEst
## Plot diagnostics ####
model_diagnostics <- T
if (model_diagnostics) {
  p_diagnostic <- plot_diagnostics(model= "group", 
                                   Pcrit = all.dat$Pcrit, 
                                   inv.temp = all.dat$inv.temp, 
                                   W = all.dat$W, 
                                   SpeciesEst = sum_est$SpeciesEst,
                                   beta_method = beta_method,
                                   beta_source = beta_source[,1],
                                   source_id = all.dat$TeamNo,
                                   method_mat = method_mat)
  
}

cAIC <- c(cAIC_model1,
          cAIC_model2,
          cAIC_model3,
          cAIC_model4)
print(cAIC - min(cAIC))
