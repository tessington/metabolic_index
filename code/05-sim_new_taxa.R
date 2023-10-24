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
library(egg)
library(cowplot)
library(KernSmooth)
library(ggridges)

# ggplot setup ####

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
sim_betas <- readRDS("analysis/taxa_sims.RDS")

 # Plot by family ####
 xlims <- c(0.75, 2)
 ylims <- c(-0.25, 0.2)
 vticks <- c(1, 2, 3, 4, 5, 6, 7, 8)
### Plot usingdensity ridgeline plot #### 
 groups.2.use <- dplyr::filter(family_summary, NoSpecies >=2)$Family
 logVfamily <- ggplot(data = dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use),
                      aes(x = logV, y = Group,fill = Group)) +
   geom_density_ridges() +
   scale_x_continuous(limits = c(xlims), sec.axis = sec_axis(~exp(.), name = "kPa",
                                                             breaks = vticks)) +
   scale_y_discrete(limits = rev) + 
   ylab("") +
   scale_fill_viridis_d() +
   theme(legend.position = "none")
 
 ylabs<-get_plot_component(logVfamily,pattern = "axis-l", return_all = T)
 xlims <- c(0.75, 2)
 groups.2.use <- dplyr::filter(family_summary, NoSpecies >=2)$Family
 logVfamily_nolab <- ggplot(data = dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use),
                      aes(x = logV, y = Group,fill = Group)) +
   geom_density_ridges() +
   scale_x_continuous(limits = c(xlims), sec.axis = sec_axis(~exp(.), name = "kPa",
                                                             breaks = vticks)) +
   scale_y_discrete(limits = rev) + 
   ylab("") +
   scale_fill_viridis_d() +
   theme(legend.position = "none",
         axis.text.y = element_blank())
 
 
 xlims <-  c(-0.25, 0.2)
 nfamily <- ggplot(data = dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use),
                      aes(x = n, y = Group,fill = Group)) +
   geom_density_ridges() +
   scale_x_continuous(limits = (xlims)) +
   scale_y_discrete(limits = rev) +
   ylab("") +
   scale_fill_viridis_d() +
   theme(legend.position = "none",
         axis.text.y = element_blank())
 
 
 Eofamily <- ggplot(data = dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use),
                      aes(x = Eo, y = Group, fill = Group)) +
   geom_density_ridges() +
   scale_x_continuous(limits = c(-0.1, 1.0), n.breaks = 4) +
   scale_y_discrete(limits = rev) +
   scale_fill_viridis_d() +
   ylab("") +
   theme(legend.position = "none",
         axis.text.y = element_blank())
 
 # make a goofy plot that only has y axis labels
 yaxis <- ggplot( data = dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use),
 aes(x = 1, y = Group)) +
   scale_y_discrete(limits = rev) + 
   scale_x_continuous(limits = c(0,1), sec.axis = sec_axis(~.+1, name = "kPa", breaks = c(1.95,2.05))) +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 100)
 
 p_row<- plot_grid(logVfamily_nolab, nfamily, Eofamily,
           align = "hv", axis = "tblr", nrow = 1)
 plot_grid(yaxis, p_row, axis = "tblr", rel_widths = c(1,3))
 
 ggsave(file = "figures/unaligned_family.pdf",
        width = 10,
        height = 8,
        units= "in")

 
 ### plot using inner 90 quantiles ####
 groups.2.use <- dplyr::filter(family_summary, NoSpecies >=2)$Family
 make_bar <- function(i, tmp.data) {
   b <- quantile(tmp.data, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
   rect(ybottom = i - .2, 
        ytop =i + 0.2,
        xleft = b[2],
        xright = b[4],
        lwd = 1.5,
        col = "lightblue"
   )
   arrows(y0 = i,
          y1 = i,
          x0 = b[2],
          x1 = b[1],
          angle = 90,
          length = 0.075,
          lwd = 1.5)
   
   arrows(y0 = i,
          y1 = i,
          x0 = b[4],
          x1 = b[5],
          angle = 90, 
          length = 0.075,
          lwd = 1.5)
   lines(x = rep(b[3],2),
         y = c(i-0.2, i + 0.2 ),
         lwd = 1.5)
   
 }
 # set layout
 # relative size of log(v) plot to others
 pdf(file= "figures/trait_by_family.pdf",
     height = 6,
     width = 8)
 ngroups <- length(groups.2.use)
 groups.2.use <- rev(groups.2.use)
 
 rel.size <- 1.6
 # set size of axis tick labels
 cex.mult <- 1.4
 # make layout matrix so that first plot has relsize
 #
 layout(mat = matrix(c(1,2,3), nrow = 1), width = c(rel.size, 1,1))
 
 # Make Empty plot for logV
 xlims <- c(0.69, 2)
 ylims <- c(-0.25, 0.2)
 vticks <- c(1, 2, 3, 4, 5, 6, 7, 8)
 par(mar = c(5,12,5,1))
 # blank plot for logV
 plot(0, type = "n",
      axes = F,
      xlab = "log(V)",
      ylab= "",
      xlim = xlims,
      ylim = c(0.5, ngroups + 0.5),
      xaxs = "i",
      yaxs = "i",
      cex.lab = 1.75)
 axis(side = 1, at = seq(0.8, 2.0, by = 0.4), cex.axis = cex.mult)
 axis(side = 2, at = 1:ngroups, labels = groups.2.use, las = 1, cex.axis= cex.mult)
 axis(side = 3, at = log(vticks), labels = vticks, line, cex.axis = cex.mult)
 box()
 mtext(side = 3, "kPa", line = 2, cex = 1.25)
 
 for (i in 1:length(groups.2.use)) {
   tmp.data <-  dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use[i])
   make_bar(i, tmp.data$logV)
 }

 # blank plot for n
 par(mar = c(5,1,5,1))
 plot(0, type = "n",
      axes = F,
      xlab = "n",
      ylab= "",
      xlim = c(-0.3, 0.21),
      ylim = c(0.5, ngroups + 0.5),
      xaxs = "i",
      yaxs = "i",
      cex.lab = 1.75)
 axis(side = 1, at = seq(-0.2, 0.2, by = 0.2), cex.axis = cex.mult)
 axis(side = 2, at = 1:ngroups, labels = F)
 box()
 for (i in 1:length(groups.2.use)) {
   tmp.data <-  dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use[i])
   make_bar(i, tmp.data$n)
 }
 
 # blank plot for Eo
 par(mar = c(5,1,5,1))
 plot(0, type = "n",
      axes = F,
      xlab = expression(E[o]),
      ylab= "",
      xlim = c(-0.1, 0.9),
      ylim = c(0.5, ngroups + 0.5),
      xaxs = "i",
      yaxs = "i",
      cex.lab = 1.75)
 axis(side = 1, at = seq(-0.0, 0.8, by = 0.4),cex.axis = cex.mult)
 axis(side = 2, at = 1:ngroups, labels = F)
 box()
 for (i in 1:length(groups.2.use)) {
   tmp.data <-  dplyr::filter(sim_betas, level == 2, Group %in% groups.2.use[i])
   make_bar(i, tmp.data$Eo)
 }
 
 dev.off()
 
 saveRDS(sim_betas, "analysis/taxa_sims.RDS")
 
 
   