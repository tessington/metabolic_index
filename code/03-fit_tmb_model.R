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

## Setup data for TMB ####

kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W/wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$Pcrit_atm<- all.dat$Pcrit / 101.325 # convert from kPa to atm
all.dat$minuslogpo2 <- - log(all.dat$Pcrit)
taxa.list <- c("Class","Order", "Family", "Genera", "Species")
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
beta_method <- matrix(fixef[grep(rownames(fixef), pattern = "beta_method"),1], nrow = n_methods)
trans <- summary(rep, "report")
lambda <- trans[grep(rownames(trans), pattern = "lambda"), ]

# Summarize Estimatess ####
sum_est <- summarize_estimates(beta_mle, beta_se, ParentChild_gz, taxa.list)

# Plot Estimates ####
plot_est <- T
if (plot_est) {
  add_number <- function(txt, num) paste0(txt, " (", num,")")
  
if (taxa.list[1] == "Class") {
  Est.2.plot <- dplyr::filter(sum_est$ClassEst, NoSpecies >=2)
  Est.2.plot <- Est.2.plot %>%
    mutate(Class_num = add_number(Class, NoSpecies))
  
  Vploto <- plotest(Est.2.plot, logV, Class_num, logVmin, logVmax)
  Vploto <- Vploto + xlab("log(V)") + xlim(c(0.5, 2.5)) + ylab("Class")
  Eoploto <- plotest(Est.2.plot, Eo, Class, Eomin, Eomax)
  Eoploto <- Eoploto + theme(axis.text.y = element_blank(), 
                             axis.title.y = element_blank()
  )  +
    xlab(expression("E"[o])) +
    xlim(c(0, 0.75)) 
  
  
  nploto <- plotest(Est.2.plot, n, Class, nmin, nmax)
  nploto <- nploto + xlim(c(-0.175, 0.05)) + theme(axis.text.y = element_blank(), 
                                                   axis.title.y = element_blank()
  )
  
  class_plot <- ggarrange(Vploto, 
                          nploto,
                          Eoploto,
                          nrow = 1)
  print(class_plot)
  }
  
  ## Plot Orders #####
  
  
  Est.2.plot <- dplyr::filter(sum_est$OrderEst, NoSpecies >=2)
  Est.2.plot <- Est.2.plot %>%
    mutate(Order_num = add_number(Order, NoSpecies))
  
  Vploto <- plotest(Est.2.plot, logV, Order_num, logVmin, logVmax)
  Vploto <- Vploto + xlab("log(V)") + xlim(c(0.5, 2.5)) + ylab("Order")
  Eoploto <- plotest(Est.2.plot, Eo, Order, Eomin, Eomax)
  Eoploto <- Eoploto + theme(axis.text.y = element_blank(), 
                             axis.title.y = element_blank()
                             )  +
    xlab(expression("E"[o])) +
    xlim(c(0, 0.75)) +
    ggtitle("Figure 2")

  
  nploto <- plotest(Est.2.plot, n, Order, nmin, nmax)
  nploto <- nploto + xlim(c(-0.175, 0.05)) + theme(axis.text.y = element_blank(), 
                                                           axis.title.y = element_blank()
  )
  
  order_plot <- ggarrange(Vploto, 
                          nploto,
                          Eoploto,
                          nrow = 1)
  

  
  ##Plot Families #####
  Est.2.plot <- dplyr::filter(sum_est$FamilyEst, NoSpecies >=2)
  Est.2.plot <- Est.2.plot %>%
    mutate(Family_num = add_number(Family, NoSpecies))
  Vplotf <- plotest(Est.2.plot, logV, Family_num, logVmin, logVmax)
  Vplotf <- Vplotf + xlab("log(V)") + xlim(c(0.5, 2.5)) + ylab("Family")
  Eoplotf <- plotest(Est.2.plot, Eo, Family, Eomin, Eomax)
  Eoplotf <- Eoplotf + theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank()
                             )  +
    xlab(expression("E"[o]))+
    xlim(c(0, 0.7)) +
    ggtitle("Figure 3")
  
  nplotf <- plotest(Est.2.plot, n, Family, nmin, nmax)
  nplotf <- nplotf + xlim(c(-0.2, 0.05)) + theme(axis.text.y = element_blank(),
                                                 axis.title.y = element_blank()
  )  
family_plot <-  ggarrange(Vplotf, 
                          nplotf,
                          Eoplotf,
                          nrow = 1)


## Multi plot of order 
ggsave(filename = "figures/order_plot_mle.pdf",
       plot = order_plot,
       device = "pdf",
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
SpeciesEst <- make_species_df(level = which(taxa.list == "Species"), 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

SpeciesEst$V = exp(SpeciesEst$logV)
saveRDS(SpeciesEst, file = "analysis/hierarchical_species_estimates.RDS")

p_diagnostic <- plot_diagnostics(model= "method", 
                                 Pcrit = all.dat$Pcrit, 
                                 inv.temp = all.dat$inv.temp, 
                                 W = all.dat$W, 
                                 SpeciesEst = SpeciesEst,
                                 beta_method = beta_method,
                                 method_mat = method_mat)
ggsave(filename= "figures/diagnostic.png",
       plot = p_diagnostic,
       units = "px",
       scale = 3,
       width = 1029,
       height = 1029)

