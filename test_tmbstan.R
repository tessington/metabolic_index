sim_taxa <-
  function(obj,
           ParentChild_gz,
           runMCMC = F,
           Groups,
           level) {
    if (runMCMC) {
      rstan_options(auto_write = TRUE)  # this option stops Stan from re-compiling if not necessary
      options(mc.cores = parallel::detectCores())
      sims <- tmbstan(
        obj,
        chains = 5,
        iter = 500,
        init = obj$env$last.par.best,
        cores = 5,
        warmup = floor(500 / 2),
        control = list(
          adapt_engaged = TRUE,
          adapt_delta = 0.9,
          max_treedepth = 15
        )
      )
      
      saveRDS(sims, file = "analysis/mcmcoutput.RDS")
    }
    if (!runMCMC)
      sims <- readRDS(file = "analysis/mcmcoutput.RDS")
    
    
    alpha_sim <- extract(sims, "alpha_j")$alpha_j
    beta_gj_sim <- extract(sims, "beta_gj")$beta_gj
    loglambda_sim <- extract(sims, "log_lambda")$log_lambda
    L_sim <- extract(sims, "L_z")$L_z
    parnames <- names(obj$env$last.par.best)
    n.sims <- nrow(alpha_sim)
    ### assign beta_gj to logAo, no, and Eo matrixes ####
    class_index <- which(ParentChild_gz$ChildTaxon == 1)
    order_index <-  which(ParentChild_gz$ChildTaxon == 2)
    family_index <-  which(ParentChild_gz$ChildTaxon == 3)
    species_index <- which(ParentChild_gz$ChildTaxon == 4)
    n_gcj <- length(class_index)
    n_goj <- length(order_index)
    n_gfj <- length(family_index)
    n_gj <- length(species_index)
    nbeta <- ncol(beta_gj_sim)
    n_j <- 3 # number of traits
    allgroups <- c("class", "order", "family")
    ngroups_array <- c(n_gcj, n_goj, n_gfj)
    nreps <- ifelse(level == 0, 1,  ngroups_array[level])
    # matrix to store simulated species drawn within taxomic groups
    
    
    for (sim in 1:n.sims) {
      #### Extract the ith MCMC simulation
      alpha_j_random <- alpha_sim[sim, ]
      beta_gj_random <- beta_gj_sim[sim, ]
      log_lambda_random <- loglambda_sim[sim, ]
      L_z_random <- L_sim[sim, ]
      L <- matrix(0, nrow = n_j, ncol = n_j)
      #### Fill iin Cholesky Matrix ####
      Count = 1
      D <- 0.001
      for (i in 1:n_j) {
        for (j in 1:n_j) {
          if (i >= j) {
            L[i, j] <- L_z_random[Count]
            Count <- Count + 1
          }
        }
      }
      sigma <- L %*% t(L) + D
      # Make dataframe with Variable and Level so later can look up corresponding elements in beta_gj_random
      # make an array for looking things up
      beta_df <-
        tibble(
          Var = rep(c("logAo", "n", "Eo"), each = (n_gcj + n_goj + n_gfj + n_gj)),
          Level = rep(c(
            rep("class", n_gcj),
            rep("order", n_goj),
            rep("family", n_gfj),
            rep("species", n_gj)
          ), n_j),
          Est = beta_gj_random
        )
      # Save just the betas corresponding to the taxonomic level
      beta_df_group <-
        dplyr::filter(beta_df, Level == allgroups[level])
      
      
      # Extract out vectors for each parameter
      logAos <-
        matrix(
          dplyr::filter(beta_df_group, Var == "logAo",)$Est,
          ncol = 1,
          nrow = ngroups_array[level]
        )
      ns <-
        matrix(dplyr::filter(beta_df_group, Var == "n")$Est,
               ncol = 1,
               nrow = ngroups_array[level])
      Eos <-
        matrix(
          dplyr::filter(beta_df_group, Var == "Eo")$Est,
          ncol = 1,
          nrow = ngroups_array[level]
        )
      # Place in single matrix with mean trait values as rows, including groupname as a fourth column ####
      group_sims <- tibble(
        logAo = logAos,
        n = ns,
        Eo = Eos,
        Group = Groups
      )
      
      
      ### using random draw for taxonomic level mean simulate value for a random species for each known taxa at that level
      
      if (level > 0)
        lambda_sum <- sum(exp(log_lambda_random[level:3]))
      if (level == 0)
        lambda_sum <- sum(exp(log_lambda_random[1:3])) + 1
      
      sim_beta_species <-
        group_sims[, 1:3] + rmvnorm(n =  nreps,
                                    mean = rep(0, 3),
                                    sigma = lambda_sum * sigma)
      ### Save result of ith simulation in a tibble ####
      sim_beta_df_i <- tibble(
        logAo = sim_beta_species[, 1],
        Eo = sim_beta_species[, 3],
        n = sim_beta_species[, 2],
        Group = Groups
      )
      ### Combine this tibble into a larger tibble by adding rows ####
      if (sim == 1)
        sim_beta_df <- sim_beta_df_i
      if (sim > 1)
        sim_beta_df <- bind_rows(sim_beta_df, sim_beta_df_i)
      
    }
    return(sim_beta_df)
  }
