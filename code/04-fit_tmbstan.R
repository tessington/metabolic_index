# Run tmbstan on a fitted tmb object #####
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

source("code/helper/fit_model_funs.R")

# load and fit TMB model
all.dat <- load_data()
all.dat <- filter_data(all.dat)


method_mat <- model.matrix(~ EstMethod_Metric , all.dat)
n_methods <- ncol(method_mat) -1

#### Create new ParentChild matrix for reduced taxonomic structure ####
kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W/wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$minuslogpo2 <- - log(all.dat$Pcrit) # fit using pO2 in atm
taxa.list <- c("Class","Order","Family", "Genera", "Species")
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

# Run TUMB ####
  # Setup TMB ####
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
                  log_lambda = rep(0, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  beta_method = 0,
                  logsigma = 0
)
  
  model <- "hierarchical_mi_stan"
  compile(paste0("code/TMB/", model, ".cpp"))
  dyn.load(dynlib(paste0("code/TMB/",model)))
  Random <- c("beta_gj")
  obj <-
    MakeADFun(
      data = data,
      parameters = parameters,
      DLL = model,
      random = Random,
      silent = TRUE,
      hessian = TRUE
    )
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  rep = sdreport( obj,
                  getReportCovariance = TRUE, 
                  getJointPrecision=TRUE)
  saveRDS(file = "analysis/modelfit_stan.RDS",list(obj = obj, opt = opt, rep = rep))

  

nchain <- 8
niter <- 6000
nwarm <- 1000
nthin <- 40
maxtree <- 12
rstan_options(auto_write = TRUE)  # this option stops Stan from re-compiling if not necessary
options(mc.cores = parallel::detectCores())
sims <- tmbstan(
  obj,
  laplace = F,
  chains = nchain,
  iter = niter,
  init = obj$env$last.par.best,
  cores = nchain,
  thin = nthin,
  warmup = nwarm,
  control = list(
    adapt_engaged = TRUE,
    max_treedepth = maxtree
  )
)

saveRDS(sims, file = "analysis/mcmcoutput_class.RDS")
sims <- readRDS(file = "analysis/mcmcoutput_class.RDS")
launch_shinystan(sims)
