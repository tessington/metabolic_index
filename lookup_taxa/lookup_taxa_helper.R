library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(tmbstan)
library(knitr)
library(egg)
library(cowplot)
library(KernSmooth)

source("code/fit_model_funs.R")
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank())
make_poly <- function(simdata) {
  smoothed <- bkde(simdata)
  polygon(x = c(smoothed$x, rev(smoothed$x)),
          y= c(rep(0, length(smoothed$x)), rev(smoothed$y)),
          col = "gray50")
  b <- quantile(simdata, probs = c(0.05, 0.95))
  lower.index <- which(smoothed$x <= b[1])
  lowersmoothed <- list(x = smoothed$x[lower.index],
                        y = smoothed$y[lower.index])
  polygon(x = c(lowersmoothed$x, rev(lowersmoothed$x)),
          y= c(rep(0, length(lowersmoothed$x)), rev(lowersmoothed$y)),
          col = "gray90",
          border = NA)
  upper.index <- which(smoothed$x >= b[2])
  uppersmoothed <- list(x = smoothed$x[upper.index],
                        y = smoothed$y[upper.index])
  polygon(x = c(uppersmoothed$x, rev(uppersmoothed$x)),
          y= c(rep(0, length(uppersmoothed$x)), rev(uppersmoothed$y)),
          col = "gray90",
          border = NA)
}

make_plot<- function(sims) {
  # set size of axis tick labels
  rel.size = 1.2
  cex.mult <- 1.75
  # make layout matrix so that first plot has relsize
  #
  allplot <- layout(mat = matrix(c(1,2,3), nrow = 3), height = c(rel.size, 1,1))
  
  # Make Empty plot for logV
  xlims <- c(0.69, 2)
  ylims <- c(-0.25, 0.2)
  vticks <- c(1, 2, 3, 4, 5, 6, 7, 8)
  par(mar = c(5,5,4,1))
  # blank plot for logV
  plot(0, type = "n",
       axes = F,
       xlab = "log(V)",
       ylab= "",
       ylim = c(0,2),
       xlim = xlims,
       xaxs = "i",
       yaxs = "i",
       cex.lab = 1.75)
  axis(side = 1, at = seq(0.8, 2.0, by = 0.4), cex.axis = cex.mult)
  axis(side = 2, cex.axis= cex.mult, las = 1)
  axis(side = 3, at = log(vticks), labels = vticks, line, cex.axis = cex.mult)
  box()
  mtext(side = 3, "kPa", line = 2, cex = 1.25)
  make_poly(simdata = sims$logV)
  par(mar = c(5,5,1,1))
 # make blank plot for n
  plot(0, type = "n",
       axes = F,
       xlab = "n",
       ylab= "",
       ylim = c(0,6),
       xlim = c(-0.3, 0.2),
       xaxs = "i",
       yaxs = "i",
       cex.lab = 1.75)
  axis(side = 1, at = seq(-0.2, 0.2, by = 0.2), cex.axis = cex.mult)
  axis(side = 2, cex.axis= cex.mult, las = 1)
  box()
  make_poly(simdata = sims$n)
  mtext(side = 2, text = "Density", line = 3, cex = 2) 
  
  
  # blank plot for Eo
  plot(0, type = "n",
       axes = F,
       xlab = expression(E[o]),
       ylab= "",
       ylim = c(0,2),
       xlim = c(-0.2, 1),
       xaxs = "i",
       yaxs = "i",
       cex.lab = 1.75)
  axis(side = 1, at = seq(0, 1, by = 0.25), cex.axis = cex.mult)
  axis(side = 2, cex.axis= cex.mult, las = 1)
  box()
  make_poly(simdata = sims$Eo)
  
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
    make_plot(sims)
  
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
