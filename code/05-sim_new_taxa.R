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
             strip.background = element_blank()
)


# load functions ####
source("code/fit_model_funs.R")
# load data ###
all.dat <- load_data()
taxa.list <- c("Order", "Family", "Genera", "Species")
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
order_summary <- all.dat %>%
  group_by(Order, Class) %>%
  summarise(NoFamily = length(unique(Family)), NoSpecies = length(unique(Species)))

family_summary <- all.dat %>%
  group_by(Family, Order) %>%
  summarise(NoGenus = length(unique(Genera)), NoSpecies= length(unique(Species)))

genera_summary <- all.dat %>%
  group_by(Family, Order, Genera) %>%
  summarise(NoSpecies= length(unique(Species)))

species_summary <- all.dat %>%
  group_by(Species) %>%
  summarise(NoSpecies= length(unique(Species)))

# simulate trait values for out of sample species calling sim_taxa() ####
sim_betas <- sim_taxa(obj = obj,
                      ParentChild_gz = ParentChild_gz
                      )

 # Plot by order ####

xlims <- c(-2, -0.5)
ylims <- c(-0.05, 0.85)
vticks <- c(1, 2, 3, 4, 5, 6, 7, 8)

 groups.2.use <- dplyr::filter(order_summary, NoFamily >=2)$Order
  dorder<-  ggplot(data =dplyr::filter(sim_betas, level == 1, Group %in% groups.2.use),
                  aes( x = logAo, y = Eo)) + 
   geom_density_2d_filled( stat = "density_2d_filled", h = NULL,
                           show.legend = F) +
    scale_x_continuous(limits = c(xlims), sec.axis = sec_axis(~exp(-.), name="V (kPa)", 
                                                              breaks = vticks)) +
    ylim(ylims) +  
    stat_ellipse(type = "norm",
                level = 0.8,
                linewidth = 1.5,
                col = "black") +
   scale_fill_manual(palette = colpal) +
   labs(x = expression(log(A[o])), y = expression(E[o])) +
   facet_wrap(vars(Group), nrow = 2, ncol = 3) + 
    theme(axis.line.x.top = element_line(color = "black"),
          axis.ticks.x.top = element_line(color = "black"),
          axis.text.x.top = element_text(color = "black"))
 dorder
 
 # Plot by family ####
 
 groups.2.use <- dplyr::filter(family_summary, NoSpecies >=2)$Family
 dfamily<-  ggplot(data =dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use),
                  aes( x = logAo, y = Eo)) + 
   geom_density_2d_filled( stat = "density_2d_filled", h = NULL,
                           show.legend = F) +
   ylim(ylims) +  
   stat_ellipse(type = "norm",
                level = 0.8,
                linewidth = 1.5,
                col = "black") +
   scale_x_continuous(limits = c(xlims), sec.axis = sec_axis(~exp(-.), name="V (atm)", 
                                                                 breaks = vticks)) +
   scale_fill_manual(palette = colpal) +
   labs(x = expression(log(A[o])), y = expression(E[o])) +
   facet_wrap(vars(Group), nrow = 4, ncol = 4) +
   theme(axis.line.x.top = element_line(color = "black"),
         axis.ticks.x.top = element_line(color = "black"),
         axis.text.x.top = element_text(color = "black"))
 dfamily
 
 
 nfamily<-  ggplot(data =dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use),
                   aes( x = n)) + 
   geom_density(fill = "#9ecae1") +
   ylab("Density") + 
   scale_y_continuous(n.breaks = 3) +
   scale_x_continuous(limits = c(-0.5, 0.5), n.breaks = 4) +
   facet_wrap(vars(Group), nrow = 4, ncol = 4) 
 
 
 dfamily
 nfamily
 ## Plot by Genera ####
 groups.2.use <- dplyr::filter(genera_summary, NoSpecies >=2)$Genera
 dgenera<-  ggplot(data =dplyr::filter(sim_betas, level == 3, Group %in% groups.2.use),
                   aes( x = logAo, y = Eo)) + 
   geom_density_2d_filled( stat = "density_2d_filled", h = NULL,
                           show.legend = F) +
   xlim(xlims) + 
   ylim(ylims) + 
   stat_ellipse(type = "norm",
                level = 0.9,
                linewidth = 1.5,
                col = "black") +
   scale_fill_manual(palette = colpal) +
   labs(x = expression(log(A[o])), y = expression(E[o])) +
   facet_wrap(vars(Group), nrow = 4, ncol = 4)
 dgenera
 saveRDS(sim_betas, "analysis/taxa_sims.RDS")
 
 savefiles <- T
 if (savefiles) {
 
   ggsave(plot = dorder,
          filename = "figures/order_plot.png",
          width = 1184*2,
          height = 745*2,
          units = "px")
   
   ggsave(plot = dfamily,
          filename = "figures/family_plot.png",
          width = 1184*2,
          height = 1184*2,
          units = "px")
   ggsave(plot = dgenera,
          filename = "figures/genera_plot.png",
          width = 1184*2,
          height = 745*2,
          units = "px")
 }
   
   