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

# load functions ####
source("code/fit_model_funs.R")

# Generate Evolutionary Trait Structure ####

# Get data with  taxonomy tree ####
all.dat <- load_data()


all.dat$Source <- factor(all.dat$Source)
all.dat$SourceNo <- as.numeric(all.dat$Source)
n_p <- length(unique(all.dat$SourceNo))



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


genera_summary <- all.dat %>%
  group_by(Family, Order, Genera) %>%
  summarise(NoSpecies= length(unique(Species)))
print(genera_summary, n = 70)

# Setup TMB data and parameters ####

kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W/wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$Pcrit_atm<- all.dat$Pcrit / 101.325 # convert from kPa to atm
all.dat$minuslogpo2 <- - log(all.dat$Pcrit)
taxa.list <- c("Order", "Family", "Genera", "Species")

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

## Setup TMB ####
data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = all.dat$minuslogpo2,
             spc_in_PCgz = spc_in_PC_gz -1,
             paper = all.dat$SourceNo - 1
             
)

parameters = list(alpha_j = rep(0,n_j),
                  L_z = rep(1, 6),
                  log_lambda = rep(0, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  beta_p = rep(0, times = n_p),
                  logsigma = 0,
                  logsigma_p = 0
)
Random <- c("beta_gj", "beta_p")
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
rep = sdreport( obj,
                getReportCovariance = TRUE, 
                getJointPrecision=TRUE)

saveRDS(list(obj = obj, opt = opt, rep = rep), "analysis/modelfit.RDS")

re <- summary(rep, "random")
fixef <- summary(rep, "fixed")
beta_mle <- matrix(re[grep(rownames(re), pattern = "beta_gj"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta_gj"),2], nrow = n_g, ncol = 3, byrow = F)



# Summarize Estimatess ####

## By Order #####
OrderEst <- make_df_plot(level = 1, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

OrderEst <- merge(OrderEst, order_summary)

## By Family ####
FamilyEst <- make_df_plot(level = 2, 
                          beta_mle,
                          beta_se,
                          ParentChild_gz,
                          groups = taxa.list)
FamilyEst <- merge(FamilyEst, family_summary)

GeneraEst <- make_df_plot(level = 3, 
                          beta_mle,
                          beta_se,
                          ParentChild_gz,
                          groups = taxa.list)
GeneraEst <- merge(GeneraEst, genera_summary)

# Plot Estimates ####
plot_est <- T
if (plot_est) {
  
## Plot Orders #####
  Est.2.plot <- dplyr::filter(OrderEst, NoSpecies >=2)
  Vploto <- plotest(Est.2.plot, logV, Order, logVmin, logVmax)
  Vploto <- Vploto + xlab("log(V)") + xlim(c(0.8, 1.75))
  Eoploto <- plotest(Est.2.plot, Eo, Order, Eomin, Eomax)
  Eoploto <- Eoploto + theme(axis.text.y = element_blank(), 
                             axis.title.y = element_blank()
                             )  +
    xlab(expression("E"[o])) +
    xlim(c(0, 0.65))
  
  nploto <- plotest(Est.2.plot, n, Order, nmin, nmax)
  nploto <- nploto + xlim(c(-0.175, 0.05)) + theme(axis.text.y = element_blank(), 
                                                           axis.title.y = element_blank()
  )
  
  order_plot <- ggarrange(Vploto, 
                          nploto,
                          Eoploto,
                          nrow = 1)
  
  ##Plot Families #####
  Est.2.plot <- dplyr::filter(FamilyEst, NoSpecies >=2)
  Vplotf <- plotest(Est.2.plot, logV, Family, logVmin, logVmax)
  Vplotf <- Vplotf + xlab("log(V)") + xlim(c(0.8, 1.75))
  Eoplotf <- plotest(Est.2.plot, Eo, Family, Eomin, Eomax)
  Eoplotf <- Eoplotf + theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank()
                             )  +
    xlab(expression("E"[o]))+
    xlim(c(0, 0.66))
  nplotf <- plotest(Est.2.plot, n, Family, nmin, nmax)
  nplotf <- nplotf + xlim(c(-0.175, 0.05)) + theme(axis.text.y = element_blank(),
                                                 axis.title.y = element_blank()
  )  
family_plot <-  ggarrange(Vplotf, 
                          nplotf,
                          Eoplotf,
                          nrow = 1)


## Multi plot of order 
ggsave(filename = "figures/order_plot_mle.png",
       plot = order_plot,
    units = "px",
    scale = 3,
    width = 1029,
    height = 629)

ggsave(filename = "figures/family_plot_mle.png",
       plot = family_plot,
       units = "px",
       scale = 3,
       width = 1029,
       height = 629)
}
# Compare fits to individual species ####
## Add to species the Ao and se ####
SpeciesEst <- make_species_df(level = 4, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

SpeciesEst$V = exp(SpeciesEst$logV)
saveRDS(SpeciesEst, file = "analysis/hierarchical_species_estimates.RDS")
