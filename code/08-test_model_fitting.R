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
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             axis.text = element_text(color = "black")
             )

# get MLE estimate for betamethod
mle <- readRDS("analysis/modelfit.RDS")
fixef <- summary(mle$rep, "fixed")
beta_method_mle <- fixef[grep(x = rownames(fixef), pattern = "beta_method"),]

# load functions ####
source("code/fit_model_funs.R")

# Get data with  taxonomy tree ####
all.dat <- load_data()

rem_styela <- T
if (rem_styela) all.dat <- all.dat %>%
  filter(Species != "Styela plicata")

rem_unknown <- T
all.dat$EstMethod_Metric <- tolower(all.dat$EstMethod_Metric)
if (rem_unknown) all.dat <- dplyr::filter(all.dat, !Method == "unknown", !EstMethod_Metric == "unknown")

# Keep only oxygen consumption methods
all.dat <- dplyr::filter(all.dat, Method == "OxygenConsumption")


# Setup TMB data and parameters ####

kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W/wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$minuslogpo2 <- - log(all.dat$Pcrit)
taxa.list <- c("Class", "Order", "Family", "Genera", "Species")

## Create new ParentChild matrix for reduced taxonomic structure ####
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
method_mat <- model.matrix(~ EstMethod_Metric , all.dat)

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
                  log_lambda = rep(0, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  logsigma = 0
)
Random <- c("beta_gj")
model <- "hierarchical_mi_no_method"
compile(paste0("code/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("code/TMB/",model)))

## Run TUMB ####
obj <-
  MakeADFun(
    data = data,
    parameters = parameters,
    DLL = model,
    random = Random,
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = TRUE, 
                getJointPrecision=TRUE)
# get fitted beta_gj
beta_gj_real <- matrix(rep$par.random, nrow = n_g, ncol = n_j)

# get species beta_gj
species_beta_gj <- beta_gj_real[ which(ParentChild_gz[,4] == 5),]
sd_beta_gj <- apply(X = species_beta_gj, FUN = sd, MAR = 2)


extract_pars <- function(opt){
  # Extract Parameters ####
  L_z <- opt$par[grep(names(opt$par), pattern = "L_z")]
  lambda <- exp(opt$par[grep(names(opt$par), pattern = "log_lambda")])
  alpha_j <- opt$par[grep(names(opt$par), pattern = "alpha_j")]
  sigma <- exp(opt$par[grep(names(opt$par), pattern = "logsigma")])
  
  # Create var covar function
  L <- matrix(0, nrow = n_j, ncol = n_j)
  #### Fill in Cholesky Matrix ####
  Count = 1
  D <- 0.00001
  Count = 1;
  for(r in 1:3){
    for(c in 1:3){
      if(r == c) {
        L[r,c] = L_z[Count]
        Count<- Count + 1
      }
      if(r>c){
        L[r,c] = L_z[Count]
        Count <- Count + 1
      }
    }
  }
  # variance in logV across all taxa (not including estimation uncertainty)
  Cov_jj <- L %*% t(L) + D
  
  logV_sigma_sum <- sqrt(Cov_jj[1,1] *(  1 +  sum(lambda)))
  
  
  
  return(list(lambda = lambda,
              sigma = sigma,
              alpha_j = alpha_j,
              Cov_jj = Cov_jj,
              logV_sigma_sum = logV_sigma_sum))
}

pars <- extract_pars(opt)
lambda <- pars$lambda
alpha_j <- pars$alpha_j
Cov_jj <- pars$Cov_jj
logV_sigma_sum = pars$logV_sigma_sum
sigma <- pars$sigma

nsims <- 100
sim_beta_method <- est_beta_method <- est_V_sd <- rep(NA, nsims)

for (sim in 1:nsims) {
## Generate Taxa-level parameters
sim_beta_gj = matrix(0, nrow = n_g, ncol = n_j)
for (g in 1:n_g) {
  Child_num = PC_gz[g,2] 
  Parent_row = PC_gz[g,1] + 1
  lambda_num = Child_num
  if( Child_num==0 ) covmult = 1
  if( Child_num>=1 ) covmult = lambda[ lambda_num ]
  tmpCov_jj = Cov_jj * covmult
  
  # get means
  if( PC_gz[g,2]==0 ) Parent_j = alpha_j
  if( PC_gz[g,2]>=1 ) Parent_j = sim_beta_gj[PC_gz[g,1] + 1,]
    # Generate parameters
  sim_beta_gj[g,] <- rmvnorm(n = 1, mean = Parent_j, tmpCov_jj)
  
  
}

# get species beta_gj
sim_species_beta_gj <- sim_beta_gj[ which(ParentChild_gz[,4] == 5),]
sim_sd_beta_gj <- apply(X = sim_species_beta_gj, FUN = sd, MAR = 2)

## simulate the data
V <- exp(sim_beta_gj[ which(ParentChild_gz[,4]==5),1])
n_pow <- sim_beta_gj[ which(ParentChild_gz[,4]==5),2]
Eo <- sim_beta_gj[ which(ParentChild_gz[,4]==5),3]
n_d <- nrow(all.dat)
invtemp = all.dat$inv.temp
logW = log(all.dat$W)
taxa_id = g_i_i
mu = rep(NA, n_d)
## simulate method effects
beta_method <- runif(1, min = -1, max = 0)


for(id in 1:n_d){
  mu[ id ] =  Eo[ taxa_id[ id ] ] * invtemp[ id ] + n_pow[ taxa_id[ id ] ]* logW[ id  ] - log(V[ taxa_id[ id ] ]) - beta_method * method_mat[ id, 2]
}

minuslogpo2 <- rnorm(n = n_d, 
                     mean = mu,
                     sd = sigma)

# fit hierarchical model
data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = minuslogpo2,
             spc_in_PCgz = spc_in_PC_gz -1,
             method_mat =  method_mat[,-1]
             
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
obj <-
  MakeADFun(
    data = data,
    parameters = parameters,
    DLL = model,
    random = Random,
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)

# get estimated sd of species
# get fitted beta_gj
beta_gj_fitted <- matrix(rep$par.random, nrow = n_g, ncol = n_j)

# get species beta_gj
fitted_species_beta_gj <- beta_gj_fitted[ which(ParentChild_gz[,4] == 5),]
fitted_sd_beta_gj <- apply(X = fitted_species_beta_gj, FUN = sd, MAR = 2)

est_pars <- extract_pars(opt)
est_V_sd[sim] <- est_pars$logV_sigma_sum
est_beta_method[sim] <- opt$par[grep(names(opt$par), pattern = "beta_method")]
sim_beta_method[sim] <- beta_method
}



# compare to fit to real data

sim_df <- tibble(simulated = sim_beta_method,
                 estimated = est_beta_method)

ggplot(data = sim_df, aes(x = simulated, y = estimated)) +
  geom_point(size = 2) +
  xlab("Simulated Parameter Value") +
  ylab("Estimated Parameter Value") +
  geom_abline(slope = 1, linewidth = 1)

