### simulate_taxa.R #####
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

# ggplot setup ####
colpal <- colorRampPalette(c("white", "#deebf7", "#9ecae1","#3182bd"))
theme_set(theme_bw(base_size = 14))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

# load functions ####
source("code/fit_model_funs.R")
# load data ###
all.dat <- load_data()
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

# load TMB fitted object
fit <- readRDS("analysis/modelfit_stan.RDS")
obj <- fit$obj
opt <- fit$opt
rep <- fit$rep

# get information by taxonomic level ####
class_summary <- all.dat %>%
  group_by(Class) %>%
  summarise(NoOrder = length(unique(Order)), NoFamily = length(unique(Family)),  NoSpecies = length(unique(Species)))

order_summary <- all.dat %>%
  group_by(Order, Class) %>%
  summarise(NoFamily = length(unique(Family)), NoSpecies = length(unique(Species)))

family_summary <- all.dat %>%
  group_by(Family, Order) %>%
  summarise(NoGenus = length(unique(Genera)), NoSpecies= length(unique(Species)))

# simulate trait values for out of sample species calling sim_taxa() ####
sim_betas <- sim_taxa(obj = obj,
                      ParentChild_gz = ParentChild_gz,
                      Groups = c(unique(all.dat$Class),
                                 unique(all.dat$Order),
                                 unique(all.dat$Family),
                                 unique(all.dat$Species)
                                 )
                      )

# Plot by class ####
# include in this plot "Other Class" for out of sample

 dclass<-  ggplot(data =dplyr::filter(sim_betas, level == 1),
                  aes( x = logAo, y = Eo)) + 
   geom_density_2d_filled( stat = "density_2d_filled", h = c(1, 0.2),
                           show.legend = F) +
   xlim(1.5, 5.5) + 
   ylim(-0.10, 0.8) + 
   stat_ellipse(type = "norm",
                level = 0.9,
                linewidth = 1.5,
                col = "black") +
   scale_fill_manual(palette = colpal) +
   labs(x = expression(log(A[o])), y = expression(E[o])) +
   facet_wrap(vars(Group), nrow = 4, ncol = 4)
 dclass
 
 # plot by order ####
 
 groups.2.use <- dplyr::filter(order_summary, NoFamily >=2)$Order
  dorder<-  ggplot(data =dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use),
                  aes( x = logAo, y = Eo)) + 
   geom_density_2d_filled( stat = "density_2d_filled", h = c(1, 0.2),
                           show.legend = F) +
   xlim(1.5, 5.5) + 
   ylim(-0.10, 0.8) + 
   stat_ellipse(type = "norm",
                level = 0.9,
                linewidth = 1.5,
                col = "black") +
   scale_fill_manual(palette = colpal) +
   labs(x = expression(log(A[o])), y = expression(E[o])) +
   facet_wrap(vars(Group), nrow = 3, ncol = 4)
 dorder
 
 # Plot by family ####
 groups.2.use <- dplyr::filter(family_summary, NoSpecies >=2)$Family
 dfamily<-  ggplot(data =dplyr::filter(sim_betas, level == 3, Group %in% groups.2.use),
                  aes( x = logAo, y = Eo)) + 
   geom_density_2d_filled( stat = "density_2d_filled", h = c(1, 0.2),
                           show.legend = F) +
   xlim(1.5, 5.5) + 
   ylim(-0.10, 0.8) + 
   stat_ellipse(type = "norm",
                level = 0.9,
                linewidth = 1.5,
                col = "black") +
   scale_fill_manual(palette = colpal) +
   labs(x = expression(log(A[o])), y = expression(E[o])) +
   facet_wrap(vars(Group), nrow = 4, ncol = 4)
 dfamily
 
 saveRDS(sim_betas, "analysis/taxa_sims.RDS")
 