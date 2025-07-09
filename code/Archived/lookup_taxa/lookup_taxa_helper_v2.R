library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(tmbstan)
library(knitr)
library(egg)
library(cowplot)
library(KernSmooth)
library(grid)

source("code/fit_model_funs.R")
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

# Functions for calculating credibility interval
# function to return the probability in the tails below a theshold density
getprob <- function(dens, x) {
  smooth <- bkde (x)
  deltax <- diff(smooth$x)[1]
  # turn into probability
  probs <- smooth$y * deltax
  index <- which(smooth$y <dens)
  return(sum(probs[index]))
}
# function to get the differrence between tail probability and target probability
fit_dens <- function(dens, x, p) {
  pstar <- getprob(dens, x)
  target_p <- 1 - p
  return(target_p - pstar)
}
# master function to solve for the density that produces a tail probability equal
# to target, and to give the limits of that range (the inner p percentile interval)
get_x_range <- function(x, p) {
  dens <- uniroot(f = fit_dens, interval = c(0.01, 8), x = x, p = 0.9)
  smooth <- bkde (x)
  index <- which(smooth$y >=dens$root)
  xbounds = c(min(smooth$x[index]), max(smooth$x[index]))
  return(xbounds)
}



lookup_taxa <- function(taxa.name, ylim) {
  all.dat <- load_data()
  # do the usual filtering of data
  all.dat <- filter_dat(all.dat)
  wref <- 5
  tref <- 15
  kb <-  8.617333262145E-5
  sim_beta <- readRDS("analysis/taxa_sims.RDS")
  ltaxa.name <- tolower(taxa.name)
  lookup.taxa <- ltaxa.name %in% tolower(sim_beta$Group)
  
  if(lookup.taxa) {
    options(warn = -1)
    sims <- dplyr::filter(sim_beta, Group == taxa.name)
    if (sims$level[1] == 1) taxa_dat <- dplyr::filter(all.dat, tolower(Class) == ltaxa.name)
    if (sims$level[1] == 2) taxa_dat <- dplyr::filter(all.dat, tolower(Order) == ltaxa.name)
    if (sims$level[1] == 3) taxa_dat <- dplyr::filter(all.dat, tolower(Family) == ltaxa.name)
    
    Wmed <- median(all.dat$W/ wref)
    
    # get prediction interval for V, Eo, n
    Eo_range <- get_x_range(x = sims$Eo, p = 0.9)
    V_range <- get_x_range(x = exp(sims$logV), p = 0.9)
    n_range <- get_x_range(x = sims$n, p = 0.9)
    
    cat(" Eo: ", Eo_range[1], " - ", median(sims$Eo), " - ", Eo_range[2], "\n")
    cat(" V: ", V_range[1], " - ", median(sims$logV)," - ", V_range[2], "\n")
    cat(" n: ", n_range[1], " - ",median(sims$Eo), " - ", n_range[2], "\n")
    
    # calculate Pcrit for 10degrees
    t_est1 <- 10
    inv.temp <- (1 / kb) * (1 / kelvin(t_est1) - 1 / kelvin(tref) )
    pcrit_1 <- sims$logV - sims$n * log(Wmed) - sims$Eo * inv.temp
    
    # calculate Pcrit for 20 degrees
    t_est2 <- 20
    inv.temp <- (1 / kb) * (1 / kelvin(t_est2) - 1 / kelvin(tref) )
    pcrit_2 <- sims$logV - sims$n * log(Wmed) - sims$Eo * inv.temp
    
    # put in a tibble
    pcrit_df <- tibble(pcrit = c(pcrit_1, pcrit_2),
                       temp = c(rep(t_est1, times = length(pcrit_1)),
                                rep(t_est2, times = length(pcrit_2))
                                )
                       )
    
    # calculate posterior prediction interval defined so that each tail has the same density
    # and print out results
    x_1_range <- get_x_range(x = exp(pcrit_1), p = 0.9)
    x_2_range <- get_x_range(x = exp(pcrit_2), p = 0.9)
    cat(" 10 degrees: ", x_1_range[1], " - ", x_1_range[2], "\n")
    cat(" 20 degrees: ", x_2_range[1], " - ", x_2_range[2] )
    
    # Add taxonomic name to upper right hand column
    grob <- grobTree(textGrob(taxa.name, x=0.9,  y=0.95, hjust=1,
                              gp=gpar(fontsize=18)))
    
    p <- ggplot(data = pcrit_df, aes(x = exp(pcrit), fill = as.factor(temp))) +
      geom_density(alpha = 0.75) + 
      ylab("Density") + 
      xlab(bquote(pO[2~crit])) + 
      scale_fill_manual(values = c("#67a9cf", "#ef8a62")) +
      theme(legend.position = "none") +
      scale_x_continuous(expand = c(0,0), limits = c(0, 17) ) + 
      scale_y_continuous(expand = c(0,0), limits = ylim) +
      annotation_custom(grob)
    
    return(p)
  }
  
  
} 

