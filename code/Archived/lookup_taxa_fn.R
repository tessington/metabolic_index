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
          return(row_in_PC)
        } else {
          message("no matching taxonomic rank found in data")
          return(NULL)
        }
      }
    }
  }
}
# Estimate pcrit for taxa ####
estimate_taxa <- function(taxa.name, 
                          w ,
                          temperature,
                          method = "smr",
                          rep,
                          ParentChild_gz) {
  
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
  names <- rownames(rep)
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
  
  # Extract taxonomic-specific values
  beta_mle_taxa <- beta_mle[lookup_row_number, ]
  beta_se_taxa <- beta_se[lookup_row_number, ]
  
  
  # Lookup beta_method
  beta_method_index <- which(names == "beta_method")
  beta_method <- fixef["beta_method", 1]
  beta_method_se <- fixef["beta_method", 2]
  
  # Do calculations to approximate the variance- covariance matrix
  # I assume no covariance between traits and beta_method
  Sigma_hat <- calc_sigma(model.fit$rep)
  scaling <- beta_se_taxa / sqrt(diag(Sigma_hat) )
  V_approx_beta_j <- diag(scaling) %*% Sigma_hat %*% diag(scaling)
  # now add 0s for beta_method
  V_approx <- matrix(0, nrow = 4, ncol = 4)
  V_approx[1:3, 1:3] <- V_approx_beta_j
  V_approx[4, 4] <- beta_method_se^2

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
  log_pcrit_se <- as.numeric(sqrt(t(x_predict) %*% V_approx %*% x_predict))
  
  # format results nicely
  parameter_estimates <- cbind(c(beta_mle_taxa, beta_method),
                               c(beta_se_taxa, beta_method_se))
  rownames(parameter_estimates) <- c("log(V)", "n", "Eo", "beta_method")
  colnames(parameter_estimates) <- c("Estimate", "SE")
  
  log_pcrit_estimate <- c(logpcrit = log_pcrit_predict, se = log_pcrit_se)
  
  return(list(parameters = parameter_estimates, log_pcrit = log_pcrit_estimate, var_covar = V_approx))
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

# Calculated variance- covariance matrix ####
calc_sigma <- function(rep) {
  n_j <- 3
  fixef <- rep$par.fixed
  L_z <- fixef[grep(names(fixef), pattern = "L_z")]
  
  # Make covar matrix
  L <- matrix(0, nrow = n_j, ncol = n_j)
  #### Fill iin Cholesky Matrix ####
  Count = 1
  D <- 0.00001
  Count = 1
  
  for (r in 1:3) {
    for (c in 1:3) {
      if (r == c) {
        L[r, c] = L_z[Count]
        Count <- Count + 1
      }
      if (r > c) {
        L[r, c] = L_z[Count]
        Count <- Count + 1
      }
    }
  }
  
  sigma <- L %*% t(L) + D
  return(sigma)
}

