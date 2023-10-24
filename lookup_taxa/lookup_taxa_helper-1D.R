library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(tmbstan)
library(knitr)
library(egg)
library(cowplot)

source("code/fit_model_funs.R")
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

make_plot<- function(sims) {
  xlim <- c(0.5, 2.5)
  ylim <- c(-0.35, 0.2)
  vticks <- c(2, 3, 4, 5, 6, 7, 9, 11, 13)
  colpal <- colorRampPalette(c("white", "#deebf7", "#9ecae1","#3182bd"))
  
  #p2d <-  ggplot(data =sims,
  #               aes( x = logV, y = n)) + 
  #  geom_density_2d_filled( stat = "density_2d_filled", h = NULL,
  #                          show.legend = F) +
  #  scale_x_continuous(limits = xlim, sec.axis = sec_axis(~exp(.), name="V (kPa)", 
  #                                                        breaks = vticks)) +
  #  ylim(ylim) + 
  #  stat_ellipse(type = "norm",
  #               level = 0.8,
  #               linewidth = 1.5,
  #               col = "black") +
  #  geom_point(size = 0.5, alpha = 0.25) + 
  #  scale_fill_manual(palette = colpal) +
  #  labs(x = "log(V)", y = "n") + 
  #  theme(axis.line.x.top = element_blank(),
  #        axis.ticks.x.top = element_line(color = "black"),
  #        axis.text.x.top = element_blank(),
  #        axis.title.x.top = element_blank())
    
  
  pV <- ggplot(data =sims,
                 aes( x =logV)) + 
    geom_density(fill = "#9ecae1") +
    ylab("Density") + 
    xlab("log(V)") +
    scale_y_continuous(n.breaks = 3) +
    scale_x_continuous(limits = xlim, sec.axis = sec_axis(~exp(.), name="V (kPa)", 
                                                               breaks = vticks)) +
    theme(axis.line.x.top = element_blank(),
          axis.ticks.x.top = element_line(color = "black"),
          axis.text.x.top = element_text(color = "black"))#,
          #axis.text.x.bottom = element_blank(),
          #axis.title.x.bottom = element_blank())
  
  
  pn <- ggplot(data =sims,
                 aes( x =n)) + 
    geom_density(fill = "#9ecae1") +
    xlim(ylim) +
    ylab("Density") +
    xlab("n") #+
    #theme(axis.text.y = element_blank(),
    #      axis.title.y = element_blank()) #+
    #coord_flip()
  
  pEo <- ggplot(data =sims,
                 aes( x =Eo)) + 
    geom_density(fill = "#9ecae1", adjust = 1.2) +
    ylab("Density") + 
    xlab(expression(E[o])) +
    scale_y_continuous(n.breaks = 3) +
    scale_x_continuous(limits = c(-0.2, 0.9))
  
  #pblank <- ggplot() + geom_blank() + theme_void()
  #vn <- ggarrange(pV, pblank,p2d, pn,
  #          widths = c(4,1),
  #          heights = c(1,4)
  #          )
  allplot<-plot_grid(pV,pn, pEo, nrow = 3)
  #allplot <- plot_grid(vn, pEo, labels = "auto", nrow = 2, rel_heights = c(5,2))
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
    logV.taxa <- sims$logV
    n.taxa <- sims$n
    eo.taxa <- sims$Eo
    
    V.med <- median(logV.taxa)
    n.med <- median(n.taxa)
    eo.med <- median(eo.taxa)
    
    cov_jj <- cov(cbind(logV.taxa, n.taxa, eo.taxa))
    cat(paste("Simulation medians for", taxa.name, "\n"))
    
    medlist <- data.frame(logV = c(V.med, sd(logV.taxa)),
                          n = c(n.med, sd(n.taxa)),
                          Eo = c(eo.med, sd(eo.taxa))
    )
    
    rownames(medlist) <- c("Median", "SD")
    
    print(knitr::kable(medlist))
    cat("\n")
    cov_jj <- as.data.frame(cov_jj)
    colnames(cov_jj) <- rownames(cov_jj) <- c("log V", "n", "Eo")
    cat("Covariance Matrix \n")
    print(knitr::kable(cov_jj))
    allplot <- make_plot(sims)
    options(warn = 0)
  }
  if(!lookup.taxa) {
    options(warn = -1)
    cat("Taxonomic group not in sample, showing distribution for unknown Class \n")
    cat("Run print.order(),  print.family() or print.genera() to see list \n of taxonomic groups")
    options(warn = -1)
    sims <- dplyr::filter(sim_beta, is.na(Group))
    # retrieve posterior medians
    logV.taxa <- sims$logV
    n.taxa <- sims$n
    eo.taxa <- sims$Eo
    
    V.med <- median(logV.taxa)
    n.med <- median(n.taxa)
    eo.med <- median(eo.taxa)
    
    cov_jj <- cov(cbind(logV.taxa, n.taxa, eo.taxa))
    cat(paste("Simulation medians for all Classes \n"))
    
    medlist <- data.frame(logV = c(V.med, sd(logV.taxa)),
                          n = c(n.med, sd(n.taxa)),
                          Eo = c(eo.med, sd(eo.taxa))
    )
    
    
    print(knitr::kable(medlist))
    cat("\n")
    cov_jj <- as.data.frame(cov_jj)
    colnames(cov_jj) <- rownames(cov_jj) <- c("log V", "n", "Eo")
    cat("Covariance Matrix \n")
    print(knitr::kable(cov_jj))
    allplot <- make_plot(sims)
  
    options(warn = 0)
  
  }
  # print inner 90% intervals
  logVint <- quantile(logV.taxa, probs = c(0.05, 0.95))
  Vint <- exp(logVint)
  Eoint <- quantile(eo.taxa, probs = c(0.05, 0.95))
  nint <- quantile(n.taxa, probs = c(0.05, 0.95))
  cat(paste0("Inner 90% interval for V is (", 
             round(Vint, 3)[1], 
             " ",
             round(Vint, 3)[2],
             ")"
  ))
  cat("\n")
  cat(paste0("Inner 90% interval for Eo is ", 
             round(Eoint, 3)[1],
             " ",
             round(Eoint, 3)[2],
             ")"
  ))
  
  cat("\n")
  cat(paste0("Inner 90% interval for n is ", 
             round(nint, 3)[1],
             " ",
             round(nint, 3)[2],
             ")"
  ))
  return(allplot)
  
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
