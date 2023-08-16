library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(tmbstan)
library(knitr)

source("code/fit_model_funs.R")
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

make_plot<- function(sims) {
  xlim <- c(1, 6)
  ylim <- c(-0.25, 0.9)
  colpal <- colorRampPalette(c("white", "#deebf7", "#9ecae1","#3182bd"))
  
  p2d <-  ggplot(data =sims,
                 aes( x = logAo, y = Eo)) + 
    geom_density_2d_filled( stat = "density_2d_filled", h = c(1, 0.2),
                            show.legend = F) +
    scale_x_continuous(limits = xlim, sec.axis = sec_axis(~exp(-.), name="V (atm)", 
                                                          breaks = c(0.3, 0.1, 0.03, 0.01, 0.003))) +
    ylim(ylim) + 
    stat_ellipse(type = "norm",
                 level = 0.8,
                 linewidth = 1.5,
                 col = "black") +
    geom_point() + 
    scale_fill_manual(palette = colpal) +
    labs(x = expression(log(A[o])), y = expression(E[o])) + 
    theme(axis.line.x.top = element_blank(),
          axis.ticks.x.top = element_line(color = "black"),
          axis.text.x.top = element_blank(),
          axis.title.x.top = element_blank())
    
  
  pao <- ggplot(data =sims,
                 aes( x =logAo)) + 
    geom_density(fill = "#9ecae1") +
    ylab("Density") + 
    scale_y_continuous(n.breaks = 3) +
    scale_x_continuous(limits = xlim, sec.axis = sec_axis(~exp(-.), name="V (atm)", 
                                                               breaks = c(0.3, 0.1, 0.03, 0.01, 0.003))) +
    
    theme(axis.line.x.top = element_blank(),
          axis.ticks.x.top = element_line(color = "black"),
          axis.text.x.top = element_text(color = "black"),
          axis.text.x.bottom = element_blank(),
          axis.title.x.bottom = element_blank())
  
  
  peo <- ggplot(data =sims,
                 aes( x =Eo)) + 
    geom_density(fill = "#9ecae1") +
    xlim(ylim) +
    ylab("Density") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    coord_flip()
  
  pn <- ggplot(data =sims,
                 aes( x =n)) + 
    geom_density(fill = "#9ecae1") +
    ylab("Density") + 
    scale_y_continuous(n.breaks = 3)
  
  pblank <- ggplot() + geom_blank() + theme_void()
  aoeo <- ggarrange(pao, pblank,p2d, peo,
            widths = c(4,1),
            heights = c(1,4)
            )
  allplot <- plot_grid(aoeo, pn, labels = "auto", nrow = 2, rel_heights = c(5,2))
  return(allplot)
}


lookup_taxa <- function(taxa.name) {
  all.dat <- load_data()
  sim_beta <- readRDS("analysis/taxa_sims.RDS")
  ltaxa.name <- tolower(taxa.name)
  lookup.taxa <- ltaxa.name %in% tolower(sim_beta$Group)
  
  if(lookup.taxa) {
    options(warn = -1)
    sims <- dplyr::filter(sim_beta, Group == taxa.name)
    
    # retrieve posterior medians
    logao.taxa <- sims$logAo
    n.taxa <- sims$n
    eo.taxa <- sims$Eo
    
    ao.med <- median(logao.taxa)
    n.med <- median(n.taxa)
    eo.med <- median(eo.taxa)
    
    cov_jj <- cov(cbind(logao.taxa, n.taxa, eo.taxa))
    cat(paste("Simulation medians for", taxa.name, "\n"))
    
    medlist <- data.frame(logAo = c(ao.med, sd(logao.taxa)),
                          n = c(n.med, sd(n.taxa)),
                          Eo = c(eo.med, sd(eo.taxa))
    )
    
    rownames(medlist) <- c("Median", "SD")
    
    print(knitr::kable(medlist))
    cat("\n")
    cov_jj <- as.data.frame(cov_jj)
    colnames(cov_jj) <- rownames(cov_jj) <- c("log Ao", "n", "Eo")
    cat("Covariance Matrix \n")
    print(knitr::kable(cov_jj))
    make_plot(sims)
    options(warn = 0)
  }
  if(!lookup.taxa) {
    options(warn = -1)
    cat("Taxonomic group not in sample, showing distribution for unknown Class \n")
    cat("Run print.order(),  print.family() or print.genera() to see list \n of taxonomic groups")
    options(warn = -1)
    sims <- dplyr::filter(sim_beta, is.na(Group))
    # retrieve posterior medians
    logao.taxa <- sims$logAo
    n.taxa <- sims$n
    eo.taxa <- sims$Eo
    
    ao.med <- median(logao.taxa)
    n.med <- median(n.taxa)
    eo.med <- median(eo.taxa)
    
    cov_jj <- cov(cbind(logao.taxa, n.taxa, eo.taxa))
    cat(paste("Simulation medians for all Classes \n"))
    
    medlist <- data.frame(logAo = c(ao.med, sd(logao.taxa)),
                          n = c(n.med, sd(n.taxa)),
                          Eo = c(eo.med, sd(eo.taxa))
    )
    
    
    print(knitr::kable(medlist))
    cat("\n")
    cov_jj <- as.data.frame(cov_jj)
    colnames(cov_jj) <- rownames(cov_jj) <- c("log Ao", "n", "Eo")
    cat("Covariance Matrix \n")
    print(knitr::kable(cov_jj))
    allplot <- make_plot(sims)
    options(warn = 0)
    return(allplot)
  }
}

print.genera <- function() {
  all.dat <- load_data()
  sprintf(unique(all.dat$Genera))
}

print.order <- function() {
  all.dat <- load_data()
  print(unique(all.dat$Order))
}

print.family <- function() {
  all.dat <- load_data()
  print(unique(all.dat$Order))
}
