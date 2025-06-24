library(KernSmooth)
library(mgcv)
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
library(gsw)

# Functions ####
kelvin <- function(x) x+273.15
find_index <- function(x,y) y <- which(y == x)
filter_data <- function(data.df) {
  
  data.df <- data.df %>%
      filter(Species != "Styela plicata")
  
  data.df$EstMethod_Metric <- tolower(data.df$EstMethod_Metric)
  data.df <- dplyr::filter(data.df, !Method == "unknown", !EstMethod_Metric == "unknown")
  ## Keep only oxygen consumption methods
  data.df <- dplyr::filter(data.df, Method == "OxygenConsumption")
  return(data.df)
}
  
  
make_df_plot <- function(level, beta_mle, beta_se, ParentChild_gz, groups) {
  
  group_index <- which(ParentChild_gz$ChildTaxon==level)
  groupname <- groups[level]
  GroupNames <- gsub(".*_","",ParentChild_gz$ChildName[group_index])
  
  
  Est <- tibble("{groupname}" := GroupNames,
                logV = beta_mle[group_index,1],
                logVmin = beta_mle[group_index,1] - beta_se[group_index,1],
                logVmax = beta_mle[group_index,1]  + beta_se[group_index,1],
                Eo = beta_mle[group_index,3],
                Eomin = beta_mle[group_index,3] - beta_se[group_index,3],
                Eomax = beta_mle[group_index,3] + beta_se[group_index,3],
                n = beta_mle[group_index,2],
                nmin = beta_mle[group_index,2] - beta_se[group_index,2],
                nmax = beta_mle[group_index,2] + beta_se[group_index,2],
  )
  return(Est)
}

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

plotest <- function(dataest, trait, groupname, xmin, xmax) {
  trait <- enquo(trait)
  groupname <- enquo(groupname)
  xmin <- enquo(xmin)
  xmax <- enquo(xmax)
  
  groupplot <- ggplot(data = dataest, aes(x = !!trait, y = !!groupname)) +
    geom_point(size = 2) + 
    scale_y_discrete(limits = rev) +
    geom_errorbar(aes(y = !!groupname,
                      xmin = !!xmin,
                      xmax = !!xmax)
    )
  
  return(groupplot)
}

# Make taxonomic matrix ####
make_taxa_tree <- function(all.dat, taxa.list) {
  Z_ik_main <- dplyr::select(all.dat, all_of(taxa.list))
  Z_ik <- unique(Z_ik_main, MARGIN = 1)
  ParentChild_gz = NULL
  # 1st column: child taxon name
  # 2nd column: parent taxon name
  # 3rd column: parent row-number in ParentChild_gz
  # 4th column: Taxon level
  # Loop through
  for (colI in 1:ncol(Z_ik)) {
    Taxa_Names = apply(Z_ik[, 1:colI, drop = FALSE],
                       MARGIN = 1,
                       FUN = paste,
                       collapse = "_")
    Unique_Taxa = unique(Taxa_Names)
    for (uniqueI in 1:length(Unique_Taxa)) {
      Which = which(Taxa_Names == Unique_Taxa[uniqueI])
      if (colI == 1) {
        ParentChild_gz = rbind(ParentChild_gz, c(Unique_Taxa[uniqueI], NA, NA, colI))
      } else{
        if (length(unique(Z_ik[Which, colI - 1])) > 1)
          stop("Taxa has multiple parents")
        ChildName = Unique_Taxa[uniqueI]
        ParentName = paste(rev(rev(strsplit(
          ChildName, "_"
        )[[1]])[-1]), collapse = "_")
        ParentChild_gz = rbind(ParentChild_gz, c(
          ChildName,
          ParentName,
          match(ParentName, ParentChild_gz[, 1]),
          colI
        ))
      }
    }
  }
  
  
  
  # Relabel
  ParentChild_gz = data.frame(ParentChild_gz)
  colnames(ParentChild_gz) = c("ChildName", "ParentName", "ParentRowNumber", "ChildTaxon")
  ParentChild_gz[, 'ParentRowNumber'] = as.numeric(as.character(ParentChild_gz[, 'ParentRowNumber']))
  ParentChild_gz[, 'ChildTaxon'] = as.numeric(as.character(ParentChild_gz[, 'ChildTaxon']))
  PC_gz <- as.matrix(ParentChild_gz[, c('ParentRowNumber', 'ChildTaxon')]) - 1
  # Identify location for every observation
  Taxa_Names = apply(Z_ik,
                     MARGIN = 1,
                     FUN = paste,
                     collapse = "_")
  g_i = match(Taxa_Names, ParentChild_gz[, 'ChildName'])
  n_k = ncol(Z_ik)
  n_j = 3 # three traits
  n_g = nrow(ParentChild_gz)
  n_i <- length(g_i)
  
  ## Create index of data to Parent - Child ####
  #Z_ik_dat <- dplyr::select(all.dat, Class, Order, Family, Species)
  Taxa_Names_dat <-  apply(Z_ik_main,
                           MARGIN = 1,
                           FUN = paste,
                           collapse = "_")
  g_i_dat = match(Taxa_Names_dat, ParentChild_gz[, 'ChildName'])
  g_i_i <- sapply(FUN = find_index, X = g_i_dat, y = g_i)
  
  # Create index of species to Parent  - Child
  spc_in_PC_gz <- which(PC_gz[, 2] == max(PC_gz[, 2]))
  return(
    list(
      ParentChild_gz = ParentChild_gz,
      PC_gz = PC_gz,
      g_i = g_i,
      g_i_i = g_i_i,
      n_k = n_k,
      n_j = n_j,
      n_g = n_g,
      n_i = n_i,
      spc_in_PC_gz = spc_in_PC_gz
    )
  )
}

load_data <- function() {
  all.dat <- readRDS(file = "data/alldata_taxonomy.RDS")
  
  # remove data with missing body size
  all.dat <- all.dat %>%
    filter(!is.na(W))

  # change name of incertae sedis orders to one name
  tmp.index <- which(all.dat$Order == "Eupercaria incertae sedis")
  all.dat$Order[tmp.index] <- "Eupercaria"
  tmp.index <- which(all.dat$Order == "Ovalentaria incertae sedis")
  all.dat$Order[tmp.index] <- "Ovalentaria"
  # remove midwater crustacean 'Pleuroncodes planipes', as authors believed it lived anaerobically
  all.dat <- dplyr::filter(all.dat, Species !="Pleuroncodes planipes")
  
  return(all.dat)
}

# Function to simulate taxa within each taxonomic group ####

sim_taxa <- function(obj, ParentChild_gz) {
  ## Load mcmc output ####
  sims <- readRDS(file = "analysis/mcmcoutput.RDS")
  
  ## extract mcmc output ####
  alpha_sim <- rstan::extract(sims, "alpha_j")$alpha_j
  beta_gj_sim <- rstan::extract(sims, "beta_gj")$beta_gj
  lambda_sim<- exp(rstan::extract(sims, "log_lambda")$log_lambda)
  L_z_sim <- rstan::extract(sims, "L_z")$L_z
  beta_method <- rstan::extract(sims, "beta_method")$beta_method
  ## setup simulation ####
  parnames <- names(obj$env$last.par.best)
  n.sims <- nrow(alpha_sim)
  ### get number of class, order, family and species ####
  class_index <-  which(ParentChild_gz$ChildTaxon == 1)
  order_index <-  which(ParentChild_gz$ChildTaxon == 2)
  family_index <-  which(ParentChild_gz$ChildTaxon == 3)
  genera_index <- which(ParentChild_gz$ChildTaxon == 4)
  species_index <- which(ParentChild_gz$ChildTaxon == 5)
  n_gcj <- length(class_index)
  n_goj <- length(order_index)
  n_gfj <- length(family_index)
  n_ggj <- length(genera_index)
  n_gj <- length(species_index)
  nbeta <- ncol(beta_gj_sim)
  n_j <- 3 # number of traits
  allgroups <- c("class", "order", "family", "genera", "species")
  ngroups_array <- c(n_gcj, n_goj, n_gfj, n_ggj)
  nreps <-  sum(ngroups_array) + n_gj
  nlevels <- max(ParentChild_gz$ChildTaxon) 
  beta_sim_list <- list()
  Groups <- gsub(".*_","",ParentChild_gz$ChildName)
  # Iterate through each mcmc iteration
  
  for (sim in 1:n.sims) {
    #### Extract the "sim"th MCMC simulation
   
    alpha_j_random <- alpha_sim[sim, ]
    beta_gj_random <- beta_gj_sim[sim, ]
    beta_gj_random_matrix <- matrix(beta_gj_random, ncol = 3, byrow = F)
    L_z_random <- L_z_sim[sim,]
    L_z_random[1] <- exp(L_z_random[1]) # exponentiate the first to deal with log transformation
    lambda_random <- lambda_sim[sim,]
    beta_method_random <- beta_method[sim]
    
    # Place in single matrix with mean trait values as rows, including groupname as a fourth column ####
    # This code sets logV assuming method = SMR
    group_sims <- tibble(
      logV = beta_gj_random_matrix[,1] + beta_method_random,
      n = beta_gj_random_matrix[,2],
      Eo = beta_gj_random_matrix[,3],
      Group = Groups,
      level = c(rep(1, n_gcj), rep(2, n_goj), rep(3, n_gfj), rep(4, n_ggj), rep(5, n_gj))
    )
  
    # Make covar matrix
    L <- matrix(0, nrow = n_j, ncol = n_j)
    #### Fill iin Cholesky Matrix ####
    Count = 1
    D <- 0.00001
    Count = 1;
    for(r in 1:3){
      for(c in 1:3){
        if(r == c) {
          L[r,c] = L_z_random[Count]
          Count<- Count + 1
        }
        if(r>c){
          L[r,c] = L_z_random[Count]
          Count <- Count + 1
        }
      }
    }
    
    sigma <- L %*% t(L) + D
    
    ### using random draw for taxonomic level mean simulate value for a random species for each known taxa at that level
    beta_sim_i <- list()
    # loop through all levels but omiting species
    for (l in 1:(nlevels - 1) )  {
      level_sims <- dplyr::filter(group_sims, level == l)
      lambdasum <- sum(lambda_random[l:(nlevels - 1)])
      tmpsims <- level_sims[,1:3] + rmvnorm(n =  ngroups_array[l],
                                            mean = rep(0, 3),
                                            sigma = lambdasum * sigma) 
      tmp_df <- tibble(logV =tmpsims[,1],
                       n = tmpsims[,2],
                       Eo = tmpsims[,3],
                       Group = dplyr::filter(group_sims, level==l)$Group,
                       level = l)
      beta_sim_i[[l]] <- tmp_df
    }
    
    # add species - level posterior 
    beta_sim_i[[nlevels]] <- dplyr::filter(group_sims, level == nlevels)
    
    
    # add out of sample classes
    lambdasum <- sum(lambda_random) + 1
    tmpsims <- alpha_j_random + rmvnorm(n =  1,
                                        mean = rep(0, 3),
                                        sigma = lambdasum * sigma) 
    tmp_df <- tibble(logV =tmpsims[,1],
                     n = tmpsims[,2],
                     Eo = tmpsims[,3],
                     Group = NA,
                     level = 0)
    beta_sim_i[[nlevels + 1]] <- tmp_df
    # combine all simulated taxa into a single tibble
    beta_sim_list[[sim]] <- do.call("rbind",beta_sim_i)
  }
  # combine all "sim" results into a single tibble
  sim_betas <- do.call("rbind", beta_sim_list)
  return(sim_betas)
}

# Function to extract predictions
  
  lookup_taxa <- function(taxa.name) {
    all.dat <- load_data()
    sim_beta <- readRDS("analysis/taxa_sims.RDS")
    ltaxa.name <- tolower(taxa.name)
    lookup.taxa <- ltaxa.name %in% tolower(sim_beta$Group)
    sim_beta$Group <- tolower(sim_beta$Group)
    if(lookup.taxa) {
      options(warn = -1)
      sims <- dplyr::filter(sim_beta, Group == ltaxa.name)
    }
    if(!lookup.taxa) {
      options(warn = -1)
      cat("Taxonomic group not in sample, showing distribution for unknown Class \n")
      cat("Run print.order(),  print.family() or print.genera() to see list \n of taxonomic groups")
      options(warn = -1)
      sims <- dplyr::filter(sim_beta, is.na(Group))
    }
    
    # retrieve posterior medians
    logV.taxa <- sims$logV
    n.taxa <- sims$n
    eo.taxa <- sims$Eo
    return_obj <- list(logV = logV.taxa, n = n.taxa, Eo = eo.taxa)
    return(return_obj)
  }
  
  # Function to make residual plots ####
  plot_diagnostics <- function(model, Pcrit, inv.temp, W, SpeciesEst, method_mat=NULL, beta_method = NULL,
                               beta_source= NULL, source_id = NULL) {
    if (!model %in% c("base", "method", "paper", "group")) stop("input model 
                                                              must be either base, method, paper or group"
    )
    
    # Calculate expected values ####
    ndata <- length(Pcrit)
    est <- rep(NA, times = length(Pcrit))
    for (i in 1:ndata) {
      species.index <- which(SpeciesEst$Species == all.dat$Species[i])
      mipars <-cbind(SpeciesEst$logV, SpeciesEst$n, SpeciesEst$Eo)[species.index,]
      #  paper.index <- all.dat$SourceNo[i]
      est[i] <- mipars[1] - log(W[i]) * mipars[2] - inv.temp[i] * mipars[3]
    }
    if (model == "method")
      est <- est + method_mat[,-1] %*% matrix(beta_method, ncol = 1)
    if (model %in% c("paper, group"))
      est <- est + method_mat[,-1] %*% matrix(beta_method, ncol = 1) + 
      beta_source[source_id]
    
    res <- log(Pcrit) - est
    
    plot.dat <-tibble(Pcrit = Pcrit,
                      est = est,
                      res = res,
                      )
    
    pred_plot <- ggplot(plot.dat, aes(x = est, y = log(Pcrit) ) ) + 
      geom_point(size = 2) + 
      xlab("Predicted log(Pcrit)") +
      geom_abline(slope = 1, linewidth = 1)
    
    
    pred_resid <- ggplot(plot.dat, aes(x = est, y = res)) + 
      geom_point(size = 2) + 
      xlab("Predicted") +
      ylab("Residual") 
    
    
    qq_plot <- ggplot(plot.dat, aes(sample = res)) +
      stat_qq() + stat_qq_line(linewidth = 1.0) +
      xlab("Theoretical Quantiles") +
      ylab("Sample Quantiles")
    
    density_plot <- ggplot(plot.dat, aes(x = res)) +
      geom_density(adjust = 1.5, linewidth = 1.0) +
      xlab("Residual") +
      ylab("Density")
    
    
    return(grid.arrange(pred_plot, pred_resid, qq_plot, density_plot, nrow = 2, ncol = 2))
  }

  
  summarize_estimates <- function (beta_mle, beta_se, ParentChild_gz, taxa.list){
    
    if (taxa.list[1] == "Order") baselevel =1
    if (taxa.list[1] == "Class") baselevel =2
    if (taxa.list[1] == "Family") baselevel = 0
    
    ## Class ####
    if (taxa.list[1] == "Class") {
      ClassEst <- make_df_plot(level = 1,
                  beta_mle,
                  beta_se,
                  ParentChild_gz,
                  groups = taxa.list)
      ClassEst <- merge(ClassEst, class_summary)
    }
    if (taxa.list[1] %in% c("Order", "Class") ) {
      ## By Order #####
      OrderEst <- make_df_plot(level = baselevel, 
                               beta_mle,
                               beta_se,
                               ParentChild_gz,
                               groups = taxa.list)
      
      OrderEst <- merge(OrderEst, order_summary)
    }
    if (taxa.list[1] %in% c("Family","Order", "Class")) {
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
    if (taxa.list[1] == "Class") output$ClassEst = ClassEst
      
    
    if (taxa.list[1] %in% c("Class", "Order")) output$OrderEst = OrderEst
    return(output)
  }
  
  
  
  plot_by_group <- function(sum_est) {
    
    ### Plot Orders #####
    add_number <- function(txt, num) paste0(txt, " (", num,")")
    
    Est.2.plot <- dplyr::filter(sum_est$OrderEst, NoSpecies >=2)
    Est.2.plot <- Est.2.plot %>%
      mutate(Order_num = add_number(Order, NoSpecies))
    
    Vploto <- plotest(Est.2.plot, logV, Order_num, logVmin, logVmax)
    Vploto <- Vploto + xlab("log(V)") + xlim(c(0.5, 2)) + ylab("Order")
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
    
    
    
    ## Plot Families #####
    Est.2.plot <- dplyr::filter(sum_est$FamilyEst, NoSpecies >=2)
    Est.2.plot <- Est.2.plot %>%
      mutate(Family_num = add_number(Family, NoSpecies))
    Vplotf <- plotest(Est.2.plot, logV, Family_num, logVmin, logVmax)
    Vplotf <- Vplotf + xlab("log(V)") + xlim(c(0.5, 2.25)) + ylab("Family")
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
    
    print(order_plot)
    print(family_plot)
    
  }
  
  filter_dat <- function(dat) {
    rem_styela <- T
    if (rem_styela) dat <- dat %>%
        filter(Species != "Styela plicata")
    rem_unknown <- T
    dat$EstMethod_Metric <- tolower(dat$EstMethod_Metric)
    if (rem_unknown) dat <- dplyr::filter(dat, !Method == "unknown", !EstMethod_Metric == "unknown")
    
    # Keep only oxygen consumption methods
    dat <- dplyr::filter(dat, Method == "OxygenConsumption")
    return(dat)
  }
  
  
    
  
  
  
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
  # function to get the difference between tail probability and target probability
  fit_dens <- function(dens, x, p) {
    pstar <- getprob(dens, x)
    target_p <- 1 - p
    return(target_p - pstar)
  }
  # master function to solve for the density that produces a tail probability equal
  # to target, and to give the limits of that range (the inner p percentile interval)
  get_x_range <- function(x, p) {
    dens <- uniroot(f = fit_dens, interval = c(0.01, 1), x = x, p = p)
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
      
      # calculate Pcrit for 10degrees
      t_est1 <- 8
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
      # Print SD of simulated species
      print(c( sd ( exp ( pcrit_1 ) ), sd( exp ( pcrit_2 ) ) ) )
      
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
        xlab(bquote(p[crit])) + 
        scale_fill_manual(values = c("#67a9cf", "#ef8a62")) +
        theme(legend.position = "none") +
        scale_x_continuous(expand = c(0,0), limits = c(0, 17) ) + 
        scale_y_continuous(expand = c(0,0), limits = ylim) +
        annotation_custom(grob)
      
      return(p)
    }
    
    
  } 
  
  
  lookup_taxa_t <- function(taxa.name, t.range, w.2.use) {
    
    all.dat <- load_data()
    # remove data as needed
    all.dat <- filter_data(all.dat)
    
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
      if (sims$level[1] == 4) taxa_dat <- dplyr::filter(all.dat, tolower(Genus) == ltaxa.name)
      if (sims$level[1] == 5) taxa_dat <- dplyr::filter(all.dat, tolower(Species) == ltaxa.name)
    }
      
      Wmed <- median(w.2.use/ wref)
      
      pcrit_df <- tibble(Temp = NA, 
                         Pcritlower90 = NA,
                         Pcritupper90= NA,
                         Pcritlower50 = NA,
                         Pcritupper50= NA)
      
      for (i in seq_along(t.range)) {
      # calculate Pcrit for given temperature
      inv.temp <- (1 / kb) * (1 / kelvin(t.range[i]) - 1 / kelvin(tref) )
      pcrit <- sims$logV - sims$n * log(Wmed) - sims$Eo * inv.temp
      pcrit_90 <- get_x_range(x = exp(pcrit), p = 0.9)
      pcrit_50 = get_x_range(x = exp(pcrit), p = 0.5)
      pcrit_df <- pcrit_df %>% add_row(Temp = t.range[i],
                            Pcritlower90 = pcrit_90[1], 
                            Pcritupper90 = pcrit_90[2],
                            Pcritlower50 = pcrit_50[1],
                            Pcritupper50 = pcrit_50[2]
                           )
      
      }
      # delete placeholder NA row
      pcrit_df <- na.omit(pcrit_df)
      # smooth lower and upper bounnds
      lower90mod = gam(Pcritlower90~ s(Temp), data = pcrit_df)
      pcrit_df$lower90s = predict.gam(lower90mod, data = pcrit_df)
      upper90mod = gam(Pcritupper90~ s(Temp), data = pcrit_df)
      pcrit_df$upper90s = predict.gam(upper90mod, data = pcrit_df)
      
      lower50mod = gam(Pcritlower50~ s(Temp), data = pcrit_df)
      pcrit_df$lower50s = predict.gam(lower50mod, data = pcrit_df)
      upper50mod = gam(Pcritupper50~ s(Temp), data = pcrit_df)
      pcrit_df$upper50s = predict.gam(upper50mod, data = pcrit_df)
      
      
      return(pcrit_df)
  } 
  
  # calc o2 solubility, relies on o2 in umol/kg
  gsw_O2sol_SP_pt <- function(sal,pt) {
    library(gsw)
    x = sal
    pt68 = pt*1.00024
    y = log((298.15 - pt68)/(273.15 + pt68))
    
    a0 =  5.80871
    a1 =  3.20291
    a2 =  4.17887
    a3 =  5.10006
    a4 = -9.86643e-2
    a5 =  3.80369
    b0 = -7.01577e-3
    b1 = -7.70028e-3
    b2 = -1.13864e-2
    b3 = -9.51519e-3
    c0 = -2.75915e-7
    
    O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
    return(O2sol)
  }
  
  
  calc_p_po2 <- function(temp, po2, taxa.name,  w) {
    # function to calculate the proportion of MCMC simulations with pcrit less than po2,
    # given a temperature and po2
    
    wref <- 5
    tref <- 15
    kb <-  8.617333262145E-5
    sim_beta <- readRDS("analysis/taxa_sims.RDS")
    ltaxa.name <- tolower(taxa.name)
    lookup.taxa <- ltaxa.name %in% tolower(sim_beta$Group)
    
    if(lookup.taxa) {
      options(warn = -1)
      sims <- dplyr::filter(sim_beta, Group == taxa.name)
    }
    
      # calculate Pcrit for given temperature
      inv.temp <- (1 / kb) * (1 / kelvin(temp) - 1 / kelvin(tref) )
      pcrit <- exp(sims$logV - sims$n * log(w/wref) - sims$Eo * inv.temp)
      # get the proportion of mcmc sims with pcrit less than the po2
      p_po2_exceeds_pcrit <- length(which(pcrit < po2)) / nrow(sims)
      return(p_po2_exceeds_pcrit)
  }

  calc_many_p_po2 <- function(temp, po2, taxa.name,  w) {
    # function to calculate the proportion of MCMC simulations with pcrit less than po2,
    # given an array of temperature and po2
    
    wref <- 5
    tref <- 15
    kb <-  8.617333262145E-5
    sim_beta <- readRDS("analysis/taxa_sims.RDS")
    ltaxa.name <- tolower(taxa.name)
    lookup.taxa <- ltaxa.name %in% tolower(sim_beta$Group)
    
    if(lookup.taxa) {
      options(warn = -1)
      sims <- dplyr::filter(sim_beta, Group == taxa.name)
  
    }
    ndata <- length(temp)
    p_po2_exceeds_pcrit <- vector(length = ndata)
    for (i in 1:ndata) {
    # calculate Pcrit for given temperature
    inv.temp <- (1 / kb) * (1 / kelvin(temp[i]) - 1 / kelvin(tref) )
    pcrit <- exp(sims$logV - sims$n * log(w/wref) - sims$Eo * inv.temp)
    # get the proportion of mcmc sims with pcrit less than the po2
    p_po2_exceeds_pcrit[i] <- length(which(pcrit < po2[i])) / nrow(sims)
    }
    return(p_po2_exceeds_pcrit)
  }
  
      