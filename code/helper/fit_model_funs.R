# functions used throughout project
library(mgcv)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(egg)
library(cowplot)
library(KernSmooth)
library(grid)
library(gsw)

# Simple Functions ####
kelvin <- function(x) x+273.15
find_index <- function(x,y) y <- which(y == x)
filter_data <- function(data.df) {
  data.df$EstMethod_Metric <- tolower(data.df$EstMethod_Metric)
  data.df <- dplyr::filter(data.df, !Method == "unknown", !EstMethod_Metric == "unknown")
  ## Keep only oxygen consumption methods
  data.df <- dplyr::filter(data.df, Method == "OxygenConsumption")
  return(data.df)
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
  
  ## Create index of data to Parent - Child 
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
    dplyr::filter(!is.na(W))

  # change name of incertae sedis orders to one name
  tmp.index <- which(all.dat$Order == "Eupercaria incertae sedis")
  all.dat$Order[tmp.index] <- "Eupercaria"
  tmp.index <- which(all.dat$Order == "Ovalentaria incertae sedis")
  all.dat$Order[tmp.index] <- "Ovalentaria"
  # remove midwater crustacean 'Pleuroncodes planipes', as authors believed it lived anaerobically
  all.dat <- dplyr::filter(all.dat, Species !="Pleuroncodes planipes")
  
  return(all.dat)
}

  # Function to make diagnostic plots ####
  plot_diagnostics <- function(model, Pcrit, inv.temp, W, SpeciesEst, method_mat=NULL, beta_method = NULL,
                               beta_source= NULL, source_id = NULL) {
    if (!model %in% c("base", "method")) stop("input model must be either base or method")
    
    # Calculate expected values ##
    ndata <- length(Pcrit)
    est <- rep(NA, times = length(Pcrit))
    for (i in 1:ndata) {
      species.index <- which(SpeciesEst$Species == all.dat$Species[i])
      mipars <-cbind(SpeciesEst$logV, SpeciesEst$n, SpeciesEst$Eo)[species.index,]
      est[i] <- mipars[1] - log(W[i]) * mipars[2] - inv.temp[i] * mipars[3]
    }
    if (model == "method")
      est <- est + method_mat[,-1] %*% matrix(beta_method, ncol = 1)
    
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
  
  # Calculate the probability that an input value po2 is greater than pcrit ####
  calc_p_po2 <- function(po2, temp, w, pcrit_type = "smr", betas, var_covar) {
    
    # Setup fixed values
    wref <- 5
    tref <- 15
    kb <-  8.617333262145E-5
    # Transform w and temperature
    inv_temperature <- (1/ kb) * (1 / kelvin(temp) - 1 / kelvin(tref) )
    logw <- log(w / wref )
    # calculate prediction x's
    if(pcrit_type == "smr") x_predict <- as.vector(c(1, -logw, -inv_temperature, 1))
    if(!pcrit_type == "smr") x_predict <- c(1, -logw, -inv_temperature, 0)
    # generate prediction and standard error
    log_pcrit_predict <- x_predict %*% betas
    log_pcrit_se <-as.numeric(sqrt( t(x_predict) %*% var_covar %*% x_predict ))
    # calculate probability that po2 exceeds pcrit
    p_po2_exceeds_pcrit <- pnorm(q = log(po2), mean = log_pcrit_predict, sd = log_pcrit_se)
    
    return(p_po2_exceeds_pcrit)
  }
  
  ## Functions for model advancement using blank taxa names for
  # out of sample prediction
  
  # Augment dataframe with blank taxa ####
  augment_taxa <- function(all.dat) {
    ## Create "empty" species for each genera, each family, each order, each class, and one
    ## out of sample
    
    # 1. Blank species per Genus
    genus_placeholders <- all.dat %>%
      distinct(Phylum, Class, Order, Family, Genus) %>%
      mutate(Species = 'sp' ,
             source = "genus_placeholder")
    
    
    # 2. Blank species + genus per Family
    family_placeholders <- all.dat %>%
      distinct(Phylum, Class, Order, Family) %>%
      mutate(Genus = "blankgenusinfamily",
             Species = "blankgenusinfamily sp",
             source = "family_placeholder")
    
    # 3. Blank species + genus + family per Order
    order_placeholders <- all.dat %>%
      distinct(Phylum, Class, Order) %>%
      mutate(Family = "blankfamilyinorder",
             Genus = "blankgenusinorder",
             Species = "blankgenusinorder sp",
             source = "order_placeholder")
    
    # 4. Blank species + genus + family + order per Class
    class_placeholders <- all.dat %>%
      distinct(Phylum, Class) %>%
      mutate(Order = "blankorderinclass",
             Family = "blankfamilyinclass",
             Genus = "blankgenusinclass",
             Species = "blankgenusinclass sp",
             source = "class_placeholder")
    
    phylum_placeholders <- all.dat %>%
      distinct(Phylum) %>%
      mutate(Class = "blankclassinphylum",
             Order = "blankorderinphylum",
             Family = "blankfamilyinphylum",
             Genus = "blankgenusinphylum",
             Species = "blankgenusinphylum sp",
             source = "class_placeholder")
    
    # Combine all placeholders
    all_placeholders <- bind_rows(genus_placeholders,
                                  family_placeholders,
                                  order_placeholders,
                                  class_placeholders,
                                  phylum_placeholders)
    
    # Add NA to other columns  Pcrit, Temp, W, and EstMethod_Metrix
    all_placeholders <- all_placeholders %>%
      mutate(Temp = NA_real_, inv.temp = NA_real_, W = NA_real_, EstMethod_Metric = NA_character_, Pcrit = NA_real_)
    
    # Combine with  original data 
    augmented_df <- bind_rows(all.dat %>% mutate(source = "real_data"),
                              all_placeholders)
    return(augmented_df)
  }
  
  # Create new ParentChild matrix for reduced taxonomic structure ####
  make_taxa_tree_augmented <- function(augmented_df, taxa.list) {
    # pull only those columns in taxa.list
    Z_ik_main <- dplyr::select(augmented_df, all_of(taxa.list))
    # get only unique taxa
    Z_ik <- unique(Z_ik_main, MARGIN = 1)
    # set up ParentChild_gz matrix
    ParentChild_gz = NULL
    # 1st column: child taxon name
    # 2nd column: parent taxon name
    # 3rd column: parent row-number in ParentChild_gz
    # 4th column: Taxon level
    # Loop through all unique taxa
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
    
    # Create index of data to Parent - Child ###
    #Z_ik_dat <- dplyr::select(all.dat, Class, Order, Family, Species)
    Taxa_Names_dat <-  apply(Z_ik_main,
                             MARGIN = 1,
                             FUN = paste,
                             collapse = "_")
    g_i_dat = match(Taxa_Names_dat, ParentChild_gz[, 'ChildName'])
    g_i_i <- sapply(FUN = find_index, X = g_i_dat, y = g_i)
    
    # Create index of species to Parent  - Child, but only for those
    # that have data (chatGPT code assistance, confirmed through testing)
    # Step 1: Identify species-level nodes in the parent-child graph
    species_rows <- which(PC_gz[, 2] == max(PC_gz[, 2]))
    
    # Step 2: Identify rows from real data
    real_rows <- which(!is.na(augmented_df$Pcrit == "real_data"))
    
    # Step 3: Get taxa names from real data
    real_taxa_names <- apply(Z_ik_main[real_rows, , drop = FALSE], 1, paste, collapse = "_")
    
    # Step 4: Match these to species-level taxon names in the parent-child table
    species_names <- ParentChild_gz$ChildName[species_rows]
    
    # Step 5: Filter only those species-level entries that are in the real data
    spc_in_PC_gz <- species_rows[species_names %in% real_taxa_names]
    
    
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
  
  # lookup lowest taxonomic group represntation ####
  lookup_taxonomic_group <- function(taxa.name, all.dat, ParentChild_gz) {
    # lookup taxa in sealifebase
    # look to see if species is in database
    taxa_in_data <- taxa.name %in% all.dat$Species
    # do this if taxa is within the database 'alldat'
    if (taxa_in_data) {
      # lookup row
      row_in_PC <- grep(x = ParentChild_gz$ChildName, pattern = taxa.name)
      cat("This taxa is in the dataset")
      cat("\n")
      return(row_in_PC)
    }
    
    # do this if taxa is not in data (lookup lowest group represented)
    if (!taxa_in_data) {
      # lookup full taxonomy of taxa
      worms_lookup <- try(as.data.frame(wm_records_name(taxa.name)), silent = T)
      if (class(worms_lookup) == "try-error") {
        stop("taxa name not found in WoRMS database.  Check spelling")
        
      } else {
        # reduce df to include only Phylum, Class, etc. Remove unaccepted taxonomy
        taxa_lookup <- as.data.frame(worms_lookup) %>%
          filter(status == "accepted") %>%
          select(phylum, class, order, family, genus, scientificname)
        
        # Identify the lowest taxonomic order that is represented in the df all.dat
        taxa_levels <- c("phylum", "class", "order", "family", "genus")
        
        # make taxa_df from all.dat that only has taxonomic information
        taxa_df <- all.dat %>%
          dplyr::select(Phylum, Class, Order, Family, Genus)
        colnames(taxa_df) = tolower(colnames(taxa_df))
        
        
        # Start from most specific and work upward to find lowest shared rank
        lowest_shared_rank <- NULL
        lowest_shared_value <- NULL
        
        for (i in seq_along(rev(taxa_levels))) {
          level <- rev(taxa_levels)[i]
          value <- taxa_lookup[[level]]
          
          if (!is.null(value) && level %in% colnames(taxa_df)) {
            if (value %in% taxa_df[[level]]) {
              lowest_shared_rank <- level
              lowest_shared_value <- value
              break
            }
          }
        }
        
        # If match was found, construct the placeholder name accordingly
        if (!is.null(lowest_shared_rank)) {
          # Create the placeholder row string based on the level
          ranks_to_use <- taxa_levels[1:which(taxa_levels == lowest_shared_rank)]
          base_taxa <- sapply(ranks_to_use, function(r)
            taxa_lookup[[r]])
          
          # Add placeholder at appropriate level
          placeholder_suffix <- list(
            genus = "sp",
            family = "blankgenusinfamily_blankgenusinfamily sp",
            order = "blankfamilyinorder_blankgenusinorder_blankgenusinorder sp",
            class = "blankorderinclass_blankfamilyinclass_blankgenusinclass_blankgenusinclass sp",
            phylum = "blankclassinplyum_blankorderinphylum_blankfamilyinphylum_blankgenusinphylum_blankgenusinphylum sp"
          )
          
          placeholder <- placeholder_suffix[[lowest_shared_rank]]
          
          # Combine taxon string
          collapsed_taxon <- paste(c(base_taxa, placeholder), collapse = "_")
          
          # Lookup row in ParentChild_gz
          row_in_PC <- which(ParentChild_gz$ChildName == collapsed_taxon)
          
          if (length(row_in_PC) > 0) {
            cat( "using estimate from the following taxa: \n" )
            cat( collapsed_taxon )
            cat("\n")
            
            return(row_in_PC)
          } else {
            message("no matching taxonomic rank found in data")
            return(NULL)
          }
        }
      }
    }
  }
  
  # Estimate pcrit for taxa using full covariance ####
  estimate_taxa <- function(taxa.name, 
                            w ,
                            temperature,
                            method = "smr",
                            rep,
                            ParentChild_gz,
                            psigma = 0) {
    
    # first load taxa
    all.dat <- load_data()
    all.dat <- filter_data(all.dat)
    
    # lookup rownumber for listed taxa within ParentChild matrix
    lookup_row_number <- lookup_taxonomic_group(taxa.name = taxa.name,
                                                all.dat = all.dat,
                                                ParentChild_gz = ParentChild_gz)
    
    # Extract fitted parameters from fitted object
    re <- summary(rep, "random")
    fixef <- summary(rep, "fixed")
    names <- rownames(summary(rep, "all"))
    
    # Lookup estimated logsigma
    logsigma <- fixef["logsigma", "Estimate"]
    sigma <- exp(logsigma)
    
    n_j <- 3 # number of traits
    n_g <- nrow(ParentChild_gz) # number of taxagroups
    
    beta_mle <- matrix(re[grep(rownames(re), pattern = "beta_gj"), 1],
                       nrow = n_g,
                       ncol = 3,
                       byrow = F)
    
    beta_se <- matrix(re[grep(rownames(re), pattern = "beta_gj"), 2],
                      nrow = n_g,
                      ncol = 3,
                      byrow = F)
    
    # get all rownumbers for beta_gj indices
    beta_gj_positions <- which(names == "beta_gj")
    all_beta_gj_matrix_indices <- matrix(beta_gj_positions,
                                         nrow = n_g,
                                         ncol = n_j,
                                         byrow = FALSE)
    
    # Extract taxonomic-specific values
    beta_mle_taxa <- beta_mle[lookup_row_number, ]
    beta_se_taxa <- beta_se[lookup_row_number, ]
    
    
    # Lookup beta_method
    beta_method_index <- which(names == "beta_method")
    beta_method <- fixef["beta_method", 1]
    beta_method_se <- fixef["beta_method", 2]
    
    # Do calculations to approximate the variance- covariance matrix
    # I assume no covariance between traits and beta_method
    
    # get sparse precision matrix
    J <- rep$jointPrecision
    # convert to var-covar
    
    V_joint <- solve(J)
    
    beta_gj_indices <- all_beta_gj_matrix_indices[lookup_row_number,]
    V <- V_joint[c(beta_gj_indices, beta_method_index), c(beta_gj_indices, beta_method_index)]
    
    
    # some placeholder settings
    wref <- 5
    tref <- 15
    kb <-  8.617333262145E-5
    
    betas <- as.vector(c(beta_mle_taxa, beta_method))
    logw <- log(w / wref)
    inv.temp <- 1 / kb * (1 / kelvin(temperature) - 1 / kelvin(tref))
    
    # Set up x_predict
    if (method == "smr") x_predict <- as.vector(c(1, -logw, -inv.temp, 1))
    if (!method == "smr") x_predict <- as.vector(c(1, -logw, -inv.temp, 0))
    
    # calculate prediction and SE of prediction
    log_pcrit_predict <-  x_predict %*% betas
    log_pcrit_se <- as.numeric(sqrt(t(x_predict) %*% V %*% x_predict + psigma * sigma^2))
    
    # format results nicely
    parameter_estimates <- cbind(c(beta_mle_taxa, beta_method),
                                 c(beta_se_taxa, beta_method_se))
    rownames(parameter_estimates) <- c("log(V)", "n", "Eo", "beta_method")
    colnames(parameter_estimates) <- c("Estimate", "SE")
    
    log_pcrit_estimate <- c(logpcrit = log_pcrit_predict, se = log_pcrit_se)
    return(list(parameters = parameter_estimates, log_pcrit = log_pcrit_estimate, var_covar = V))
  }
  # Fit TMB model using augmented taxa dataframe ####
  fit_model_augmented_taxa <- function(fitnew = F) {
    if (fitnew) {
      library(dplyr)
      library(MASS)
      library(TMB)
      library(ggplot2)
      library(gridExtra)
      library(readxl)
      library(Matrix)
      
      conflicted::conflict_prefer("select", "dplyr")
      conflicted::conflict_prefer("filter", "dplyr")
      
      # Setup Code ####
      ## load functions ####
      source("code/helper/fit_model_funs.R")
      
      ## load data ####
      all.dat <- load_data()
      all.dat <- filter_data(all.dat)
      augmented_df <- augment_taxa(all.dat)
      
      # Do all of the setup stuff on data for TMB
      method_mat <- model.matrix( ~ EstMethod_Metric , all.dat)
      n_methods <- ncol(method_mat) - 1
      kb <-  8.617333262145E-5
      tref <- 15
      wref <- 5
      all.dat$W <- all.dat$W / wref
      all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1 / (tref + 273.15))
      all.dat$minuslogpo2 <- -log(all.dat$Pcrit)
      
      taxa.list <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
      
      taxa.info <- make_taxa_tree_augmented(augmented_df, taxa.list)
      
      ParentChild_gz <- taxa.info$ParentChild_gz
      PC_gz <- taxa.info$PC_gz
      g_i <- taxa.info$g_i
      g_i_i <- taxa.info$g_i_i
      n_k <- taxa.info$n_k
      n_j <- taxa.info$n_j
      n_g <- taxa.info$n_g
      n_i <- taxa.info$n_i
      # which rows of taxa ParentChild_gz are attached to "real" species in data
      spc_in_PC_gz <- taxa.info$spc_in_PC_gz
      ## Setup TMB ####
      data <- list(
        PC_gz = PC_gz,
        g_i = g_i - 1,
        invtemp = all.dat$inv.temp,
        logW = log(all.dat$W),
        taxa_id = g_i_i - 1,
        minuslogpo2 = -log(all.dat$Pcrit),
        spc_in_PCgz = spc_in_PC_gz - 1,
        method_mat = method_mat[, -1]
      )
      
      parameters = list(
        alpha_j = c(0, 0, 0),
        L_z = rep(1, 6),
        log_lambda = rep(-1, length(unique(PC_gz[, 2])) - 1),
        beta_gj = matrix(0, nrow = n_g, ncol = n_j),
        beta_method = 0,
        logsigma = 0
      )
      Random <- c("beta_gj")
      model <- "hierarchical_mi"
      compile(paste0("code/TMB/", model, ".cpp"))
      dyn.load(dynlib(paste0("code/TMB/", model)))
      
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
      rep = sdreport(obj,
                     getReportCovariance = TRUE,
                     getJointPrecision = TRUE)
      # save rep and ParentChild_gz in analysis/
      saveRDS(object = list(rep = rep, ParentChild_gz = ParentChild_gz),
              file = "analysis/model_fit_augmented.RDS")
      return(list(rep = rep, ParentChild_gz = ParentChild_gz) )
    }
    if (!fitnew) {
      # check to see if file exists in analysis subfolder
      check.file <- file.exists("analysis/model_fit_augmented.RDS")
      if (!check.file)
        stop("model_fit_augmented.RDS not found. Re-run setting fitnew = TRUE")
      if (check.file) {
        model.fit <- readRDS("analysis/model_fit_augmented.RDS")
        return(model.fit)
      }
    }
  }
  

  
  