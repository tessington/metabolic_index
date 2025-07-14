# Code fit hierarchical model to Arrenhius Equation
rm(list = ls())
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
library(rlang)

# Setup Code ####


## load functions ####
source("code/helper/fit_model_funs.R")
source("code/helper/calculate_EDF_fn.R")

# ggplot setup ####
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             axis.text = element_text(color = "black"),
             axis.text.y = element_text(size = 16)
)
# ploting and summary functions used below ####
make_group_plot <- function(group_to_plot,sum_est, saveplot = T) {
  group_sym <- enquo(group_to_plot)
  dataset_name <- paste0(as_name(group_sym), "Est")
  
  Est.2.plot <- sum_est[[dataset_name]] %>%
    filter(NoSpecies  >= 2) %>%
    mutate(group_num = paste0(!!group_sym, " (",  NoSpecies , ")"))
  
  Vplot <- plotest(Est.2.plot, logV,  logVse, group_num)
  Vplot <- Vplot + xlab("log(V)") + xlim(x_v_lims) + ylab("")
  Eoplot <- plotest(Est.2.plot, Eo, Eose, group_num)
  Eoplot <- Eoplot + theme(axis.text.y = element_blank(), axis.title.y = element_blank())  + 
    xlab(expression("E"[o])) + xlim(x_eo_lims) 
  
  nplot <- plotest(Est.2.plot, n, nse, group_num)
  nplot <- nplot + xlim(x_n_lims) + theme(axis.text.y = element_blank(), 
                                          axis.title.y = element_blank()
  )
  
  group_plot <- ggarrange(Vplot, 
                          nplot,
                          Eoplot,
                          nrow = 1)
  print(group_plot)
  if (saveplot) {
    saveplotname <- paste0("figures/", as_name(group_sym), "_plot_mle.png")
    ggsave(filename = saveplotname,
           plot = group_plot,
           device = "png",
           units = "px",
           scale = 3,
           width = 1029,
           height = 629)
  }
  return(group_plot)
}
# Make a dataframe for plotting group-level estimates #
make_df_plot <- function(level, beta_mle, beta_se, ParentChild_gz, groups) {
  
  group_index <- which(ParentChild_gz$ChildTaxon==level)
  groupname <- groups[level]
  GroupNames <- gsub(".*_","",ParentChild_gz$ChildName[group_index])
  
  Est <- tibble("{groupname}" := GroupNames,
                logV = beta_mle[group_index,1],
                logVse =  beta_se[group_index,1],
                Eo = beta_mle[group_index,3],
                Eose = beta_se[group_index,3],
                n = beta_mle[group_index,2],
                nse =  beta_se[group_index,2]
  )
  return(Est)
}
# Function for plotting estimates by specified group ####
plotest <- function(dataest, trait, xse, groupname) {
  trait <- enquo(trait)
  groupname <- enquo(groupname)
  xse <- enquo(xse)
  
  groupplot <- ggplot(data = dataest, aes(x = !!trait, y = !!groupname)) +
    geom_point(size = 2) + 
    scale_y_discrete(limits = rev) +
    geom_errorbar(aes(y = !!groupname,
                      xmin = !!trait - !!xse,
                      xmax = !!trait + !!xse)
    )
  
  return(groupplot)
}
# Function to make a species-level dataframe of estimates ####
make_species_df <- function(level, beta_mle, beta_se, ParentChild_gz, groups) {
  group_index <- which(ParentChild_gz$ChildTaxon==level)
  groupname <- groups[level]
  GroupNames <- gsub(".*_","",ParentChild_gz$ChildName[group_index])
  
  
  Est <- tibble("{groupname}" := GroupNames,
                logV = beta_mle[group_index,1],
                logVSE = beta_se[group_index,1],
                Eo = beta_mle[group_index,3],
                EoSE = beta_se[group_index,3],
                n = beta_mle[group_index,2],
                nSE = beta_se[group_index,2],
                
  )
  return(Est)
}
summarize_estimates <- function (beta_mle, beta_se, ParentChild_gz, taxa.list){
  
  if (taxa.list[1] == "Family") baselevel = 0
  if (taxa.list[1] == "Order") baselevel =1
  if (taxa.list[1] == "Class") baselevel =2
  if (taxa.list[1] == "Phylum") baselevel =3
  
  ## Phylum ####
  if (taxa.list[1] == "Phylum") {
    PhylumEst <- make_df_plot(level = 1,
                              beta_mle,
                              beta_se,
                              ParentChild_gz,
                              groups = taxa.list)
    PhylumEst <- merge(PhylumEst, phylum_summary)
  }
  ## Class ####
  if (taxa.list[1] %in% c("Phylum", "Class")) {
    ClassEst <- make_df_plot(level = 2,
                             beta_mle,
                             beta_se,
                             ParentChild_gz,
                             groups = taxa.list)
    ClassEst <- merge(ClassEst, class_summary)
  }
  if (taxa.list[1] %in% c("Order", "Class", "Phylum") ) {
    ## By Order #####
    OrderEst <- make_df_plot(level = baselevel, 
                             beta_mle,
                             beta_se,
                             ParentChild_gz,
                             groups = taxa.list)
    
    OrderEst <- merge(OrderEst, order_summary)
  }
  if (taxa.list[1] %in% c("Family","Order", "Class", "Phylum")) {
    ## By Family ####
    FamilyEst <- make_df_plot(level = baselevel + 1, 
                              beta_mle,
                              beta_se,
                              ParentChild_gz,
                              groups = taxa.list)
    FamilyEst <- merge(FamilyEst, family_summary)
  }
  
  
  SpeciesEst <- make_species_df(level = baselevel + 3, 
                                beta_mle,
                                beta_se,
                                ParentChild_gz,
                                groups = taxa.list)
  
  
  output <- list(FamilyEst = FamilyEst,
                 SpeciesEst = SpeciesEst)
  if (taxa.list[1] == "Phylum") output$PhylumEst = PhylumEst
  if (taxa.list[1] %in% c("Phylum", "Class") ) output$ClassEst = ClassEst
  if (taxa.list[1] %in% c("Phylum", "Class", "Order") ) output$OrderEst = OrderEst
  return(output)
}

# load data ####
all.dat <- load_data()
all.dat <- filter_data(all.dat)
method_mat <- model.matrix(~ EstMethod_Metric , all.dat)
n_methods <- ncol(method_mat) -1

 # Summarise data by taxa ####
phylum_summary <- all.dat %>%
  group_by(Phylum) %>%
  summarise(NoClass = n_distinct(Class), NoOrder = n_distinct(Order), NoFamily = n_distinct(Family),  NoGenus = n_distinct(Genus), NoSpecies = n_distinct(Species))

class_summary <- all.dat %>%
  group_by(Class) %>%
  summarise(NoOrder = length(unique(Order)), NoFamily = length(unique(Family)),  NoSpecies = length(unique(Species)))

order_summary <- all.dat %>%
  group_by(Order, Class) %>%
  summarise(NoFamily = length(unique(Family)), NoSpecies = length(unique(Species)))

family_summary <- all.dat %>%
  group_by(Family, Order) %>%
  summarise(NoGenus = length(unique(Genus)), NoSpecies= length(unique(Species)))

genus_summary <- all.dat %>%
  group_by(Family, Order, Genus) %>%
  summarise(NoSpecies= length(unique(Species)))

# Setup data for TMB ####

kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W/wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$Pcrit_atm<- all.dat$Pcrit / 101.325 # convert from kPa to atm
all.dat$minuslogpo2 <- - log(all.dat$Pcrit)
taxa.list <- c("Phylum", "Class","Order", "Family", "Genus", "Species")
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
             minuslogpo2 = -log(all.dat$Pcrit),
             spc_in_PCgz = spc_in_PC_gz -1,
             method_mat = method_mat[,-1]
)

parameters = list(alpha_j = c(0, 0, 0),
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
# set axis limits
x_v_lims <- c(0.75, 2.2)
x_eo_lims <- c(-0.125, 0.6)
x_n_lims <- c(-0.2, 0.075)

plot_est <- T
if (plot_est) {
  if (taxa.list[1] == "Phylum") make_group_plot(Phylum, sum_est, saveplot = T)

  if (taxa.list[1] %in% c("Phylum", "Class") ) make_group_plot(Class, sum_est, saveplot = T)
  make_group_plot(Order, sum_est, saveplot = T)
  make_group_plot(Family, sum_est,  saveplot = T)
}
# Compare fits to individual species ####
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


# Get species-level trait variance and make table
sampled_species_standard_deviation <- apply(X = SpeciesEst[,c("logV","n", "Eo")], MAR = 2, FUN = sd)

# get variance-covariance from cholesky matrix
L_z <- fixef[grep(rownames(fixef), pattern = "L_z"),1]
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

# Get all species sigma

## a plot of all species pcrit, vs. temperature and body mass
alldata_plot_temperature <- ggplot(data = all.dat, aes(x = inv.temp, y = log(Pcrit), col = as.factor(Phylum) )) + 
  scale_colour_viridis_d(option = "turbo") +
  geom_point(size = 2) +
  xlab("Inverse Temperature") + 
  ylab(bquote( log( p["crit"] ) ) ) +
  labs(col = "Phylum")

alldata_plot_w<- ggplot(data = all.dat, aes(x =log(W), y = log(Pcrit), col = as.factor(Phylum) )) + 
  scale_colour_viridis_d(option = "turbo") +
  geom_point(size = 2) +
  xlab ("log(W)" ) + 
  ylab(bquote( log( p["crit"] ) ) ) +
  labs(col = "Phylum")

alldata_plot <- 
  cowplot::plot_grid(alldata_plot_temperature + theme(legend.position="none"),
             alldata_plot_w + theme(legend.position="none"),
             nrow = 1,
             align = "v")
  
legend <- get_legend(
    # create some space to the left of the legend
    alldata_plot_temperature + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
alldata_plot <- plot_grid(alldata_plot, legend, rel_widths = c(3, 0.85) )

ggsave(file = "figures/alldata_plot.png",
       plot = alldata_plot,
       units = "px",
       scale = 2,
       height = 500,
       width = 1200,
       bg = "white")

# print out parameter values reported in manuscript

## levels of variance - the log_lambdas
print(exp(fixef[rownames(fixef) == "log_lambda", "Estimate"]))

## print out alpha_j - the mean trait values
fixef[rownames(fixef) == "alpha_j", ]
### First element is log(V).  Convert to V and get SE using delta method
c(exp(fixef[rownames(fixef) == "alpha_j", "Estimate"][1])  , se =  exp(fixef[rownames(fixef) == "alpha_j", "Estimate"][1]) * fixef[rownames(fixef) == "alpha_j", "Std. Error"][1] )

# for each species, document how many unique body sizes and unique temperatures

number_bodysize_temperature <- all.dat %>%
  group_by(Species) %>%
  summarise(n_sizes = length(unique(W)), n_temps = length(unique(Temp)))
n_species_2_temperature <- length(which(number_bodysize_temperature$n_temps > 2))
n_species_2_sizes <- length(which(number_bodysize_temperature$n_sizes > 2))
print(c(n_species_2_sizes, n_species_2_temperature))

