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

# ggplot setup ####
theme_set(theme_bw(base_size = 16))
theme_update(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text = element_text(color = "black")
)

# load functions ####
source("code/fit_model_funs.R")

# Get data with  taxonomy tree ####
all.dat <- load_data()

# Remove cancer irroratus?
rem_canc_irroratus <- T
if (rem_canc_irroratus)  all.dat <- dplyr::filter(all.dat, !Species == "Cancer irroratus")

rem_gardiner <- T
if (rem_gardiner)  all.dat <- dplyr::filter(all.dat, !Source == "Gardiner et al. (2010)")

all.dat$Source <- factor(all.dat$Source)
all.dat$SourceNo <- as.numeric(all.dat$Source)
n_p <- length(unique(all.dat$SourceNo))

# load fitted model
fitted <- readRDS(file = "analysis//modelfit.RDS")
opt <- fitted$opt
rep <- fitted$rep
obj <- fitted$obj
re <- summary(rep, "random")
fixef <- summary(rep, "fixed")



# Setup TMB data and parameters ####

kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W / wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1 / (tref + 273.15))
all.dat$minuslogpo2 <- -log(all.dat$Pcrit)
taxa.list <- c("Order", "Family", "Genera", "Species")

## Create new Parent Child matrix for  taxonomic structure ####
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

## get fitted beta_gj ####
real_beta_gj <- matrix(re[grep(rownames(re), pattern = "beta_gj"), 1],
                       nrow = n_g,
                       ncol = 3,
                       byrow = F)
real_beta_p <- re[grep(rownames(re), pattern = "beta_p"), 1]

## get species beta_gj ####
real_species_beta_gj <- real_beta_gj[which(ParentChild_gz[, 4] == 4), ]

# function to extract parameters from fitted model
extract_pars <- function(opt) {
  # Extract Parameters ####
  L_z <- opt$par[grep(names(opt$par), pattern = "L_z")]
  lambda <- exp(opt$par[grep(names(opt$par), pattern = "log_lambda")])
  alpha_j <- opt$par[grep(names(opt$par), pattern = "alpha_j")]
  sigma <- exp(opt$par[grep(names(opt$par), pattern = "\\blogsigma\\b")])
  
  # Create var covar function
  L <- matrix(0, nrow = n_j, ncol = n_j)
  #### Fill in Cholesky Matrix ####
  Count = 1
  D <- 0.00001
  Count = 1
  
  for (r in 1:3) {
    for (c in 1:3) {
      if (r == c) {
        L[r, c] = L_z[Count]
        Count <- Count + 1
      }
      if (r > c) {
        L[r, c] = L_z[Count]
        Count <- Count + 1
      }
    }
  }
  # variance in logV across all taxa (not including estimation uncertainty)
  Cov_jj <- L %*% t(L) + D
  
  logV_sigma_sum <- sqrt(Cov_jj[1, 1] * (1 +  sum(lambda)))
  
  
  
  return(
    list(
      lambda = lambda,
      sigma = sigma,
      alpha_j = alpha_j,
      Cov_jj = Cov_jj,
      logV_sigma_sum = logV_sigma_sum
    )
  )
}
## extract and assign fixed effect parameters ####
pars <- extract_pars(opt)
lambda <- pars$lambda
alpha_j <- pars$alpha_j
Cov_jj <- pars$Cov_jj
logV_sigma_sum = pars$logV_sigma_sum
sigma <- pars$sigma

logsigma_p <- fixef[grep(pattern = "logsigma_p", x = rownames(fixef)), 1]
sigma_p <- exp(logsigma_p)

# Begin Simulation ####
nsims <- 100
vmult <- 5 # multiplier to increase variance of V 
V_rsquare <- V_rmse <- V_bias <- rep(NA, nsims)
beta_p_rsquare <- beta_p_rmse <- beta_p_bias <- rep(NA, nsims)
for (sim in 1:nsims) {
  ## Generate Taxa-level parameters
  sim_beta_gj = matrix(0, nrow = n_g, ncol = n_j)
  for (g in 1:n_g) {
    Child_num = PC_gz[g, 2]
    Parent_row = PC_gz[g, 1] + 1
    lambda_num = Child_num
    if (Child_num == 0)
      covmult = 1
    if (Child_num >= 1)
      covmult = lambda[lambda_num]
    tmpCov_jj = Cov_jj * covmult
    tmpCov_jj[1,1] <- vmult * tmpCov_jj[1,1]
    
    # get means
    if (PC_gz[g, 2] == 0)
      Parent_j = alpha_j
    if (PC_gz[g, 2] >= 1)
      Parent_j = sim_beta_gj[PC_gz[g, 1] + 1, ]
    # Generate parameters
    sim_beta_gj[g, ] <- rmvnorm(n = 1, mean = Parent_j, tmpCov_jj)
    
  }
  # for testing, override with actual estimated beta_gj ####
  #sim_beta_gj <- real_beta_gj
  
  # get species beta_gj
  sim_species_beta_gj <- sim_beta_gj[which(ParentChild_gz[, 4] == 4), ]
  
  ## simulate the data
  V <- exp(sim_beta_gj[which(ParentChild_gz[, 4] == 4), 1])
  n_pow <- sim_beta_gj[which(ParentChild_gz[, 4] == 4), 2]
  Eo <- sim_beta_gj[which(ParentChild_gz[, 4] == 4), 3]
  
  n_d <- nrow(all.dat)
  invtemp = all.dat$inv.temp
  logW = log(all.dat$W)
  taxa_id = g_i_i
  paper = all.dat$SourceNo
  mu = rep(NA, n_d)
  ## simulate paper effects
  sim_beta_p <- rnorm(n = n_p, mean = 0, sd = sigma_p)
  # for testing, override with actual betap
  #sim_beta_p <- real_beta_p
  for (id in 1:n_d) {
    mu[id] =  Eo[taxa_id[id]] * invtemp[id] + n_pow[taxa_id[id]] * logW[id] - log(V[taxa_id[id]]) + sim_beta_p[paper[id]]
  }
  
  minuslogpo2 <- rnorm(n = n_d, mean = mu, sd = sigma)
  
  # fit hierarchical model
  data <- list(
    PC_gz = PC_gz,
    g_i = g_i - 1,
    invtemp = all.dat$inv.temp,
    logW = log(all.dat$W),
    taxa_id = g_i_i - 1,
    minuslogpo2 = minuslogpo2,
    spc_in_PCgz = spc_in_PC_gz - 1,
    paper = all.dat$SourceNo - 1
    
  )
  
  parameters = list(
    alpha_j = rep(0, n_j),
    L_z = rep(1, 6),
    log_lambda = rep(0, length(unique(PC_gz[, 2])) - 1),
    beta_gj = matrix(0, nrow = n_g, ncol = n_j),
    beta_p = rep(0, times = n_p),
    logsigma = 0,
    logsigma_p = 0
  )
  Random <- c("beta_gj", "beta_p")
  model <- "hierarchical_mi"
  compile(paste0("code/TMB/", model, ".cpp"))
  dyn.load(dynlib(paste0("code/TMB/", model)))
  
  ## Run TUMB ####
  obj_sim <-
    MakeADFun(
      data = data,
      parameters = parameters,
      DLL = model,
      random = Random,
      silent = TRUE
    )
  opt_sim <- nlminb(obj_sim$par, obj_sim$fn, obj_sim$gr)
  rep_sim <- sdreport(obj_sim)
  re_sim <- summary(rep_sim, "random")
  
  # get fitted beta_gj
  beta_gj_fitted <- matrix(re_sim[grep(rownames(re_sim), pattern = "beta_gj"), 1],
                           nrow = n_g,
                           ncol = 3,
                           byrow = F)
  # get species beta_gj
  fitted_species_beta_gj <- beta_gj_fitted[which(ParentChild_gz[, 4] == 4), ]
  fitted_V <- exp(fitted_species_beta_gj[, 1])
  # get betap
  beta_p_fitted <- re_sim[grep(rownames(re_sim), pattern = "beta_p"), 1]
  
  # Save results
  V_rsquare[sim] <- cor(fitted_V, (V)) ^ 2
  V_rmse[sim] <- sqrt(sum((V - fitted_V) ^ 2) / length(V))
  V_bias[sim] <- mean(fitted_V - V)
  
  beta_p_rsquare[sim] <- cor(beta_p_fitted, sim_beta_p) ^ 2
  beta_p_rmse[sim] <- sqrt(sum((beta_p_fitted - sim_beta_p) ^ 2) / length(sim_beta_p))
  beta_p_bias[sim] <- mean(sim_beta_p - beta_p_fitted)
  
}
