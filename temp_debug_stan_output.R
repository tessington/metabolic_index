# Debug tmbstan on a fitted tmb object #####
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

source("code/fit_model_funs.R")
make_Sigma <- function(L_z) {
  L <- matrix(0, nrow = n_j, ncol = n_j)
  #### Fill iin Cholesky Matrix ####
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
  sigma <- L %*% t(L) + D
  return(sigma)
}
# load and fit TMB model
all.dat <- load_data()


rem_styela <- T
if (rem_styela) all.dat <- all.dat %>%
  filter(Species != "Styela plicata")
rem_unknown <- T
all.dat$EstMethod_Metric <- tolower(all.dat$EstMethod_Metric)
if (rem_unknown) all.dat <- dplyr::filter(all.dat, !Method == "unknown", !EstMethod_Metric == "unknown")

# Keep only oxygen consumption methods
all.dat <- dplyr::filter(all.dat, Method == "OxygenConsumption")
method_mat <- model.matrix(~ EstMethod_Metric , all.dat)
n_methods <- ncol(method_mat) -1

#### Create new ParentChild matrix for reduced taxonomic structure ####
kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W/wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$minuslogpo2 <- - log(all.dat$Pcrit) # fit using pO2 in atm
taxa.list <- c("Class", "Order","Family", "Genera", "Species")
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
nll_mle_fit <- opt$objective


# get ML-estimate of class-level variance
fixef <- summary(rep, "fixed")
L_z <- fixef[grep(x = rownames(fixef),
                    pattern = "L_z"), ]

sigma_base <- make_Sigma(L_z[ , 1])
               
  
trans <- summary(rep, "report")
log_lambda <- fixef[grep(rownames(fixef), pattern = "log_lambda"), ]
lambda <- exp(log_lambda)
sigma_class <- sigma_base * sum(lambda[2:4,1])
sigma_order <- sigma_class * sum(lambda[3:4,1])
sqrt(diag(sigma_class))
sqrt(diag(sigma_order))
# load in stan fit


sims <- readRDS(file = "analysis/mcmcoutput.RDS")

alpha_sim <- rstan::extract(sims, "alpha_j")$alpha_j
beta_gj_sim <- rstan::extract(sims, "beta_gj")$beta_gj
lambda_sim<- exp(rstan::extract(sims, "log_lambda")$log_lambda)

log_lambda_sim <- rstan::extract(sims, "log_lambda")$log_lambda
log_lambda_median <-  apply(X = log_lambda_sim, MARGIN = 2, FUN = median)
L_z_sim <- rstan::extract(sims, "L_z")$L_z
# posterior medians for L_z
L_z_median <- apply(L_z_sim, FUN = median, MARGIN = 2)

## setup simulation ####
parnames <- names(obj$env$last.par.best)
n.sims <- nrow(alpha_sim)

# create empty list to store output
sigma.list <- sigma.class.list <- sigma.order.list <- list()



for (sim in 1:n.sims) {
  L_z_random <- L_z_sim[sim,]
  lambda_random <- lambda_sim[sim,]
  Sigma <- make_Sigma(L_z_random) # base variance - variance
  Sigma_class <- make_Sigma( L_z_random ) *  ( sum(lambda_random[2:4]) )
  Sigma_order <- make_Sigma( L_z_random ) *  ( sum(lambda_random[3:4]) )
  sigma.list[[sim]] <- Sigma
  sigma.class.list[[sim]] <- Sigma_class
  sigma.order.list[[sim]] <- Sigma_order
}

# get posterior median  
sigma_array <-  array(unlist(sigma.list), c(dim(sigma.list[[1]]), length(sigma.list)))
sigma_median <- apply(sigma_array, 1:2, median)
print(sigma_median)
print(sigma_base)


sigma_class_array <-  array(unlist(sigma.class.list), c(dim(sigma.class.list[[1]]), length(sigma.class.list)))
sigma_class_median <- apply(sigma_class_array, 1:2, median)

print(sigma_class_median)
print(sigma_class)

sigma_order_array <- array(unlist(sigma.order.list), c(dim(sigma.order.list[[1]]), length(sigma.order.list)))
sigma_order_median <- apply(sigma_order_array, 1:2, median)

print(sigma_order_median)
print(sigma_order)


print(sqrt(diag(sigma_class)))
print(sqrt(diag(sigma_class_median)))

print(sqrt(diag(sigma_order)))
print(sqrt(diag(sigma_order_median)))

# get the correlation among parameters in bayesan fit
corr <- matrix(0, nrow = 3, ncol = 3)
corr[2,1] <- sigma_median[2,1] / (sqrt(sigma_median[1,1]) * sqrt(sigma_median[2,2]) )
corr[3,1] <- sigma_median[3,1] / (sqrt(sigma_median[1,1]) * sqrt(sigma_median[3,3]) )
corr[3,2] <- sigma_median[2,3] / (sqrt(sigma_median[2,2]) * sqrt(sigma_median[3,3]) )
corr^2

# refit ML model using posterior medians
data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = all.dat$minuslogpo2,
             spc_in_PCgz = spc_in_PC_gz -1,
             method_mat = method_mat[,-1],
             L_z = L_z_median
             
)

parameters = list(alpha_j = rep(0,n_j),
                  log_lambda = rep(0, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  beta_method = 0,
                  logsigma = 0
)

model <- "hierarchical_mi_stan_fixLz"
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
nll_bayes_med <- opt$objective

print(c(nll_mle_fit, nll_bayes_med))
