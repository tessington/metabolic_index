# Functions ####
find_index <- function(x,y) y <- which(y == x)
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
for( colI in 1:ncol(Z_ik)){
  Taxa_Names = apply( Z_ik[,1:colI,drop=FALSE], MARGIN=1, FUN=paste, collapse="_")
  Unique_Taxa = unique(Taxa_Names)
  for( uniqueI in 1:length(Unique_Taxa) ){
    Which = which( Taxa_Names == Unique_Taxa[uniqueI] )
    if( colI==1 ){
      ParentChild_gz = rbind( ParentChild_gz, c(Unique_Taxa[uniqueI], NA, NA, colI) )
    }else{
      if( length(unique(Z_ik[Which,colI-1]))>1 ) stop("Taxa has multiple parents")
      ChildName = Unique_Taxa[uniqueI]
      ParentName = paste(rev(rev(strsplit(ChildName,"_")[[1]])[-1]),collapse="_")
      ParentChild_gz = rbind( ParentChild_gz, c(ChildName, ParentName, match(ParentName,ParentChild_gz[,1]), colI) )
    }
  }
}



# Relabel
ParentChild_gz = data.frame( ParentChild_gz )
colnames(ParentChild_gz) = c("ChildName", "ParentName", "ParentRowNumber", "ChildTaxon")
ParentChild_gz[,'ParentRowNumber'] = as.numeric(as.character(ParentChild_gz[,'ParentRowNumber']))
ParentChild_gz[,'ChildTaxon'] = as.numeric(as.character(ParentChild_gz[,'ChildTaxon']))
PC_gz<- as.matrix(ParentChild_gz[, c('ParentRowNumber', 'ChildTaxon')]) - 1
# Identify location for every observation
Taxa_Names = apply( Z_ik, MARGIN=1, FUN=paste, collapse="_")
g_i = match( Taxa_Names, ParentChild_gz[,'ChildName'] )
n_k = ncol(Z_ik)
n_j = 3 # three traits
n_g = nrow(ParentChild_gz)
n_i <- length(g_i)

## Create index of data to Parent - Child ####
#Z_ik_dat <- dplyr::select(all.dat, Class, Order, Family, Species)
Taxa_Names_dat <-  apply( Z_ik_main, MARGIN=1, FUN=paste, collapse="_")
g_i_dat = match( Taxa_Names_dat, ParentChild_gz[,'ChildName'] )
g_i_i <- sapply(FUN = find_index, X = g_i_dat, y = g_i)

# Create index of species to Parent  - Child
spc_in_PC_gz <- which(PC_gz[,2] == max(PC_gz[,2]))
return(list(ParentChild_gz = ParentChild_gz,
            PC_gz = PC_gz,
            g_i = g_i,
            g_i_i = g_i_i,
            n_k = n_k,
            n_j = n_j,
            n_g = n_g,
            n_i = n_i,
            spc_in_PC_gz = spc_in_PC_gz)
       )
}

load_data <- function() {
  all.dat <- readRDS(file = "data/alldata_taxonomy.RDS")
  
  # remove data with missing body size
  all.dat <- all.dat %>%
    filter(!is.na(W))
# when a species isn't listed, make a genus_spp.
    naIndex <- which(is.na(all.dat$Species))
  for (i in 1:length(naIndex)) all.dat$Species[naIndex[i]] <- paste0(all.dat$Genera[naIndex[i]], " spc")
  
  # get median mass for each species
  #species.median.mass <- all.dat %>%
  #  group_by(Species) %>%
  #  summarize(Wmed = median(W))
  # get median T for each species
  #species.median.temp <- all.dat %>%
  #  group_by(Species) %>%
  #  summarize(Tmed = median(Temp))
  
  #Tref <- median(species.median.temp$Tmed)
  # divide Actual mass by median mass for that species
  #for (i in 1:nrow(species.median.mass)) {
  #  spc.index <- which(all.dat$Species == species.median.mass$Species[i])
  #  all.dat$W[spc.index] <- all.dat$W[spc.index] / species.median.mass$Wmed[i]
  #}
  # for species with only 1 mass, replace the NA with 1
  #all.dat$W[is.na(all.dat$W)] = 1
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
  ## setup simulation ####
  parnames <- names(obj$env$last.par.best)
  n.sims <- nrow(alpha_sim)
  ### get number of class, order, family and species ####
  order_index <-  which(ParentChild_gz$ChildTaxon == 1)
  family_index <-  which(ParentChild_gz$ChildTaxon == 2)
  genera_index <- which(ParentChild_gz$ChildTaxon == 3)
  species_index <- which(ParentChild_gz$ChildTaxon == 4)
  n_goj <- length(order_index)
  n_gfj <- length(family_index)
  n_ggj <- length(genera_index)
  n_gj <- length(species_index)
  nbeta <- ncol(beta_gj_sim)
  n_j <- 3 # number of traits
  allgroups <- c("order", "family", "genera")
  ngroups_array <- c(n_goj, n_gfj, n_ggj)
  nreps <-  sum(ngroups_array) + n_gj
  nlevels <- max(ParentChild_gz$ChildTaxon) - 1
  beta_sim_list <- list()
  Groups <- gsub(".*_","",ParentChild_gz$ChildName)
  # Iterate through each mcmc iteration
  for (sim in 1:n.sims) {
    #### Extract the "sim"th MCMC simulation
    
    alpha_j_random <- alpha_sim[sim, ]
    beta_gj_random <- beta_gj_sim[sim, ]
    L_z_random <- L_z_sim[sim,]
    lambda_random <- lambda_sim[sim,]
    
    # Place in single matrix with mean trait values as rows, including groupname as a fourth column ####
    group_sims <- tibble(
      logV = beta_gj_random[1:nreps],
      n = beta_gj_random[(nreps +1):(2* nreps)],
      Eo = beta_gj_random[(2 * nreps + 1):(3*nreps)],
      Group = Groups,
      level = c(rep(1, n_goj), rep(2, n_gfj), rep(3, n_ggj), rep(4, n_gj))
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
    
    for (l in 1:nlevels) {
      level_sims <- dplyr::filter(group_sims, level == l)
      lambdasum <- sum(lambda_random[l:nlevels])
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
    
    # add out of sample orders
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