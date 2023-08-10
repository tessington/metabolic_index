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

## Set the number of random simulations #####
n.sims <- 10000
### ggplot setup ####
theme_set(theme_bw(base_size = 14))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

# load functions ####
source("code/fit_model_funs.R")

# Generate Evolutionary Trait Structure ####

## Get data with  taxonomy tree ####
all.dat <- load_data()

# describe the data a bit - how many unique  order per class, families per order, etc. ####
class_summary <- all.dat %>%
  group_by(Class) %>%
  summarise(NoOrder = length(unique(Order)), NoFamily = length(unique(Family)),  NoSpecies = length(unique(Species)))
print(class_summary, n = 30)

order_summary <- all.dat %>%
  group_by(Order, Class) %>%
  summarise(NoFamily = length(unique(Family)), NoSpecies = length(unique(Species)))
print(order_summary, n = 50)

family_summary <- all.dat %>%
  group_by(Family, Order) %>%
  summarise(NoGenus = length(unique(Genera)), NoSpecies= length(unique(Species)))
print(family_summary, n = 70)

# Setup TMB data and parameters ####
#### Create new ParentChild matrix for reduced taxonomic structure ####
kb <-  8.617333262145E-5
tref <- 15
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$Pcrit_atm<- all.dat$Pcrit / 101.325 # convert from kPa to atm
all.dat$minuslogpo2 <- - log(all.dat$Pcrit_atm) # fit using pO2 in atm
taxa.list <- c("Class", "Order", "Family", "Species")


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

# Setup TMB ####
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

model <- "hierarchical_mi"
compile(paste0("code/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("code/TMB/",model)))

# Run TUMB ####
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

saveRDS(list(obj = obj, opt = opt, rep = rep), "analysis/modelfit.RDS")

re <- summary(rep, "random")
fixef <- summary(rep, "fixed")
beta_mle <- matrix(re[grep(rownames(re), pattern = "beta"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta"),2], nrow = n_g, ncol = 3, byrow = F)


### Plot Estimates ####
#### Plot Class Estimates ####

ClassEst <- make_df_plot(level = 1, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

Est.2.plot <- merge(ClassEst, class_summary)
plot_est <- T
if (plot_est) {
Aoplot <- plotest(Est.2.plot, logAo, Class, logAomin, logAomax)
Eoplot <- plotest(Est.2.plot, Eo, Class, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Class, nmin, nmax)
grid.arrange(Aoplot, Eoplot, ncol = 2)
}
#### Plot Orders #####

OrderEst <- make_df_plot(level = 2, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

OrderEst <- merge(OrderEst, order_summary)
Est.2.plot <- dplyr::filter(OrderEst, NoSpecies >=3)
if(plot_est) {
Aoplot <- plotest(Est.2.plot, logAo, Order, logAomin, logAomax)
Eoplot <- plotest(Est.2.plot, Eo, Order, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Order, nmin, nmax)
grid.arrange(Aoplot, Eoplot, ncol = 2)
}

#### Plot Families ####
FamilyEst <- make_df_plot(level = 3, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)
FamilyEst <- merge(FamilyEst, family_summary)
Est.2.plot <- dplyr::filter(FamilyEst, NoSpecies >=2)
if (plot_est) {
Aoplot <- plotest(Est.2.plot, logAo, Family, logAomin, logAomax)
Eoplot <- plotest(Est.2.plot, Eo, Family, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Family, nmin, nmax)
grid.arrange(Aoplot, Eoplot, ncol = 2)
}
### Compare fits to individual species ####
#### Add to species the Ao and se
SpeciesEst <- make_df_plot(level = 4, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

SpeciesEst$Ao_atm = exp(SpeciesEst$logAo)
SpeciesEst$Aoind <- NA
SpeciesEst$Eoind <- NA


#### Load previous fits ####
spc.fits <- readRDS(file = "Data/species_estimates.RDS")
# reload the all.dat to be able to retrieve APHIAID for matching

for (i in 1:nrow(spc.fits)) {
  spc.2.use <- tolower(spc.fits$Species[i])
  spc.index <- which(tolower(SpeciesEst$Species) == spc.2.use)
  SpeciesEst$Aoind[spc.index] <- exp(spc.fits$logAo[i])
  SpeciesEst$Eoind[spc.index] <- spc.fits$Eo[i]
  }

ggplot(SpeciesEst, aes(x = Aoind, y = Ao_atm)) +
  geom_point() + 
  geom_abline(a = 0, b = 1) +
  xlab("Ao estimated independently") +
  ylab("Ao estimated hierarchically")
