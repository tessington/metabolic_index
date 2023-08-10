# Code fit hierarchical model to Arrenhius Equation
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)
library(Matrix)
library(mvtnorm)

## Set the number of random simulations #####
n.sims <- 10000
### ggplot setup ####
theme_set(theme_bw(base_size = 14))
theme_update(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

### Functions ####
find_index <- function(x,y) y <- which(y == x)
make_df_plot <- function(level, beta_mle, beta_se, ParentChild_gz, groups) {
  
  group_index <- which(ParentChild_gz$ChildTaxon==level)
  groupname <- groups[level]
  GroupNames <- gsub(".*_","",ParentChild_gz$ChildName[group_index])
  
  
  Est <- tibble("{groupname}" := GroupNames,
                logAo = beta_mle[group_index,1],
                logAomin = beta_mle[group_index,1] - beta_se[group_index,1],
                logAomax = beta_mle[group_index,1]  + beta_se[group_index,1],
                Eo = beta_mle[group_index,3],
                Eomin = beta_mle[group_index,3] - beta_se[group_index,3],
                Eomax = beta_mle[group_index,3] + beta_se[group_index,3],
                n = beta_mle[group_index,2],
                nmin = beta_mle[group_index,2] - beta_se[group_index,2],
                nmax = beta_mle[group_index,2] + beta_se[group_index,2],
  )
  return(Est)
}

plotest <- function(dataest, trait, groupname, xmin, xmax) {
  trait <- enquo(trait)
  groupname <- enquo(groupname)
  xmin <- enquo(xmin)
  xmax <- enquo(xmax)
  
  
  groupplot <- ggplot(data = dataest, aes(x = !!trait, y = !!groupname)) +
    geom_point() + 
    scale_y_discrete(limits = rev) +
    geom_errorbar(aes(y = !!groupname,
                      xmin = !!xmin,
                      xmax = !!xmax)
    ) 
  
  return(groupplot)
}
# Function to simulate draws from mvnormal given precision matrix ####
rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
  set.seed( random_seed )
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L = Matrix::Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  return(mu + z)
}


# Function to simulate species within a specified grouping ####
sim_spc_in_group <- function(obj, 
                             rep,
                             n.sims, 
                             ParentChild_gz, 
                             Groups,
                             level) {
  
  # get random draws for all variables
  newpar = rmvnorm_prec( mu = obj$env$last.par.best,
                         prec = rep$jointPrecision, 
                         n.sims = n.sims,
                         random_seed = sample(1:1000000, 1))
  
  parnames <- names(obj$env$last.par.best)
  # extract the simulated beta_gj
  beta_gj_random <-newpar[grep(parnames, pattern = "beta_gj"),]
  
  # calculate sigmae
  L_z <- obj$env$last.par.best[grep(parnames, pattern = "L_z")]
  L <- matrix(0, nrow = n_j, ncol = n_j)
  #### Fill iin Cholesky Matrix ####
  Count = 1
  for (i in 1 : n_j) {
    for (j in 1 : n_j) {
      if (i>=j) {
        L[i,j] <- L_z[Count]
        Count <- Count + 1
      }
    }
  }
  sigma <- L %*% t(L)
  ### extract log_lambda ####
  log_lambda <- obj$env$last.par.best[grep(parnames, pattern = "log_lambda")]
  
  ### assign beta_gj to Ao, no, and Eo matrixes ####
  class_index <- which(ParentChild_gz$ChildTaxon == 1)
  order_index <-  which(ParentChild_gz$ChildTaxon == 2)
  family_index <-  which(ParentChild_gz$ChildTaxon == 3)
  species_index <- which(ParentChild_gz$ChildTaxon == 4)
  n_gcj <- length(class_index)
  n_goj <- length(order_index)
  n_gfj <- length(family_index)
  n_gj <- length(species_index)
  nbeta <- nrow(beta_gj_random)
  n_j <- 3 # number of traits
  # make an array for looking things up
  
  # Make datafraem with Variable and Level so later can look up corresponding elements in beta_gj_random
  beta_df <- tibble(Var = rep(c("Ao", "n", "Eo"), each = (n_gcj + n_goj + n_gfj + n_gj)),
                      Level = rep(c(rep("class", n_gcj), rep("order", n_goj), rep("family", n_gfj), rep("species", n_gj)), n_j),
                      Est = beta_gj_random)
  
  allgroups <- c("class", "order", "family")
  ngroups_array <- c(n_gcj, n_goj, n_gfj)
  # Filter out means corresponding the the taxonomic level ####
  beta_df_group <- dplyr::filter(beta_df, Level == allgroups[level])
  # Extract out vectors for each parameter
  Aos <- matrix(dplyr::filter(beta_df_group, Var == "Ao")$Est, ncol = 1, nrow = n.sims * ngroups_array[level])
  ns <- matrix(dplyr::filter(beta_df_group, Var == "n")$Est, ncol = 1, nrow = n.sims * ngroups_array[level])
  Eos <- matrix(dplyr::filter(beta_df_group, Var == "Eo")$Est, ncol = 1, nrow = n.sims * ngroups_array[level])
  
  # Place in single matrix with traits as rows, including groupname as a fourth colun ####
  group_sims <- tibble(Ao = Aos,
                       n = ns,
                       Eo = Eos,
                       Group = rep(Groups, each = n.sims)
  )
  
  ### using random draw for taxonomic level mean simulate value for random species within that level
  lambda_sum <- sum(exp(log_lambda[level:3]))
  sim_beta_species <- group_sims[,1:3]+ rmvnorm(n = n.sims * ngroups_array[level], mean = rep(0,3), sigma = lambda_sum * sigma)
  
  ### Save result in a tibble ####
  sim_beta_df <- tibble(Ao = sim_beta_species[,1],
                        Eo = sim_beta_species[,3],
                        n = sim_beta_species[,2],
                        Group = rep(Groups,n.sims)
  )
  return(sim_beta_df)
}

### Generate Evolutionary Trait Structure ####

#### Get taxonomy tree ####
all.dat <- readRDS(file = "data/alldata_taxonomy.RDS")
# remove data where there is no body size data
all.dat <- dplyr::filter(all.dat, !is.na(W))
naIndex <- which(is.na(all.dat$Species))
for (i in 1:length(naIndex)) all.dat$Species[naIndex[i]] <- paste0(all.dat$Genera[naIndex[i]], " spc")

# get median mass for each species
species.median.mass <- all.dat %>%
  group_by(Species) %>%
  summarize(Wmed = median(W))
# divide Actual mass by median mass for that species
for (i in 1:nrow(species.median.mass)) {
  spc.index <- which(all.dat$Species == species.median.mass$Species[i])
  all.dat$W[spc.index] <- all.dat$W[spc.index] / species.median.mass$Wmed[i]
}
# for species with only 1 mass, replace the NA with 1
all.dat$W[is.na(all.dat$W)] = 1
# describe the data a bit - how many unique  order per class, families per order, etc.
class_summary <- all.dat %>%
  group_by(Class) %>%
  summarise(NoOrder = length(unique(Order)), NoFamily = length(unique(Family)),  NoSpecies = length(unique(Species)))
print(class_summary, n = 30)

order_summary <- all.dat %>%
  group_by(Order, Class) %>%
  summarise(NoFamily = length(unique(Family)), NoSpecies = length(unique(Species)))
print(order_summary, n = 50)

family_summary <- all.dat %>%
  group_by(Family, Order) %>%
  summarise(NoGenus = length(unique(Genera)), NoSpecies= length(unique(Species)))
print(family_summary, n = 70)

rem.classes <- FALSE
# Remove class with a single order
if (rem.classes) {
keep.class <- class_summary$Class[class_summary$NoOrder>1]
all.dat <- dplyr::filter(all.dat, Class %in% keep.class)
}

fish.only <- FALSE
if (fish.only) all.dat <- dplyr::filter(all.dat, Class == "Actinopteri")

### Setup TMB data and parameters ####
#### Create new ParentChild matrix for reduced taxonomic structure ####
kb <-  8.617333262145E-5
tref <- 15
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$Pcrit_atm<- all.dat$Pcrit / 101.325 # convert from kPa to atm
all.dat$minuslogpo2 <- - log(all.dat$Pcrit_atm) # fit using pO2 in atm
taxa.list <- c("Class", "Order", "Family", "Species")


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

#### Create index of data to Parent - Child ####
#Z_ik_dat <- dplyr::select(all.dat, Class, Order, Family, Species)
Taxa_Names_dat <-  apply( Z_ik_main, MARGIN=1, FUN=paste, collapse="_")
g_i_dat = match( Taxa_Names_dat, ParentChild_gz[,'ChildName'] )
g_i_i <- sapply(FUN = find_index, X = g_i_dat, y = g_i)

# Create index of species to Parent  - Child
spc_in_PC_gz <- which(PC_gz[,2] == max(PC_gz[,2]))

# Setup TMB ####
data <- list(PC_gz = PC_gz,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = all.dat$minuslogpo2,
             spc_in_PCgz = spc_in_PC_gz -1
             
)

parameters = list(alpha_j = rep(0,n_j),
                  L_z = rep(1, 6),
                  log_lambda = rep(0, length(unique(PC_gz[,2])) -1),
                  beta_gj = matrix(0, nrow = n_g, ncol = n_j),
                  logsigma = 0
)
Random <- c("beta_gj")

model <- "hierarchical_mi"
compile(paste0("code/TMB/", model, ".cpp"))
dyn.load(dynlib(paste0("code/TMB/",model)))


obj <-
  MakeADFun(
    data = data,
    parameters = parameters,
    DLL = model,
    random = Random,
    silent = TRUE
  )
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep = sdreport( obj,
                getReportCovariance = TRUE, 
                getJointPrecision=TRUE)
summary(rep, "fixed")
re <- summary(rep, "report")


spc_parameters <- re[grep(rownames(re), pattern  = "spc_ij"),1]
spc_parameters_se <-re[grep(rownames(re), pattern  = "spc_ij"),2]
spc_ij_mle <- matrix(spc_parameters, nrow = n_i, ncol = 3, byrow = F)
spc_ij_se <- matrix(spc_parameters_se, nrow = n_i, ncol = 3, byrow = F)


re <- summary(rep, "random")
fixef <- summary(rep, "fixed")
beta_mle <- matrix(re[grep(rownames(re), pattern = "beta"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta"),2], nrow = n_g, ncol = 3, byrow = F)


### Plot Estimates ####
#### Plot Class Estimates ####

ClassEst <- make_df_plot(level = 1, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

Est.2.plot <- merge(ClassEst, class_summary)
#Est.2.plot <- dplyr::filter(Est.2.plot, NoSpecies >=3)
Aoplot <- plotest(Est.2.plot, logAo, Class, logAomin, logAomax)
Eoplot <- plotest(Est.2.plot, Eo, Class, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Class, nmin, nmax)
grid.arrange(Aoplot, Eoplot, ncol = 2)

#### Plot Orders #####

OrderEst <- make_df_plot(level = 2, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

OrderEst <- merge(OrderEst, order_summary)
Est.2.plot <- dplyr::filter(OrderEst, NoSpecies >=3)
Aoplot <- plotest(Est.2.plot, logAo, Order, logAomin, logAomax)
Eoplot <- plotest(Est.2.plot, Eo, Order, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Order, nmin, nmax)
grid.arrange(Aoplot, Eoplot, ncol = 2)


#### Plot Families ####
FamilyEst <- make_df_plot(level = 3, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)
FamilyEst <- merge(FamilyEst, family_summary)
Est.2.plot <- dplyr::filter(FamilyEst, NoSpecies >=2)
Aoplot <- plotest(Est.2.plot, logAo, Family, logAomin, logAomax)
Eoplot <- plotest(Est.2.plot, Eo, Family, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Family, nmin, nmax)
grid.arrange(Aoplot, Eoplot, ncol = 2)


#### Plot Species ####
SpeciesEst <- make_df_plot(level = 4, 
                          beta_mle,
                          beta_se,
                          ParentChild_gz,
                          groups = taxa.list)
Est.2.plot <- SpeciesEst
Aoplot <- plotest(Est.2.plot, logAo, Species, logAomin, logAomax)
Eoplot <- plotest(Est.2.plot, Eo, Species, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Species, nmin, nmax)
grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)



### Compare fits to individual species ####
#### Add to species the Ao and se
SpeciesEst$Ao_atm = spc_ij_mle[,1] 
SpeciesEst$Ao_atm_SE = spc_ij_se[,1] 
SpeciesEst$Aoind <- NA
SpeciesEst$Eoind <- NA


#### Load previous fits ####
spc.fits <- readRDS(file = "Data/species_estimates.RDS")
# reload the all.dat to be able to retrieve APHIAID for matching

for (i in 1:nrow(spc.fits)) {
  spc.2.use <- tolower(spc.fits$Species[i])
  spc.index <- which(tolower(SpeciesEst$Species) == spc.2.use)
  SpeciesEst$Aoind[spc.index] <- exp(spc.fits$logAo[i])
  SpeciesEst$Eoind[spc.index] <- spc.fits$Eo[i]
  }

ggplot(SpeciesEst, aes(x = Aoind, y = Ao_atm)) +
  geom_point() + 
  geom_abline(a = 0, b = 1) +
  xlab("Ao estimated independently") +
  ylab("Ao estimated hierarchically")


ggplot(SpeciesEst, aes(x = Eoind, y = Eo)) +
  geom_point() + 
  geom_abline(a = 0, b = 1) + 
  xlab("Eo estimated independently") +
  ylab("Eo estimated hierarchically")


#### Plot All Classess ####

Group.2.use <- ClassEst$Class
  sim_beta_df_group <- sim_spc_in_group(obj  = obj,
                                        rep = rep,
                                        n.sims = n.sims,
                                        ParentChild_gz = ParentChild_gz,
                                        Group = Group.2.use,
                                        level = 1
  )
  

d <- ggplot(data = sim_beta_df_group, aes( x = Ao, y = Eo)) + 
  geom_density_2d_filled( stat = "density_2d_filled", h = c(1, 0.2),
                          show.legend = F) +
  xlim(1.5, 5.5) + 
  ylim(0, 0.8) + 
  facet_wrap(vars(Group), nrow = 4, ncol = 3)
print(d)  


#### Plot Orders with at least 3 species ####

redOrderEst <- dplyr::filter(OrderEst, NoSpecies >=3)$Order
Groups.2.use <- OrderEst$Order

sim_beta_df_group <- sim_spc_in_group(obj  = obj,
                                      rep = rep,
                                      n.sims = n.sims,
                                      ParentChild_gz = ParentChild_gz,
                                      Groups = Groups.2.use,
                                      level = 2
)

  

d <- ggplot(data = dplyr::filter(sim_beta_df_group, Group %in% redOrderEst), 
            aes( x = Ao, y = Eo)) + 
  geom_density_2d_filled( stat = "density_2d_filled", h = c(1, 0.2),
                          show.legend = F) +
  xlim(1.5, 5.5) + 
  ylim(0, 0.8) + 
  facet_wrap(vars(Group), nrow = 2, ncol = 2)
print(d)  

#### Plot families with at least 2 species ####

redFamilyNames <- dplyr::filter(FamilyEst, NoSpecies >=2)$Family
allFamilyNames <- FamilyEst$Family

  Groups.2.use <- allFamilyNames
  sim_beta_df_group <- sim_spc_in_group(obj  = obj,
                                        rep = rep,
                                        n.sims = n.sims,
                                        ParentChild_gz = ParentChild_gz,
                                        Groups = Groups.2.use,
                                        level = 3
  )
  
sim_beta_df_group$V <- 1 / sim_beta_df_group$Ao

d <- ggplot(data = dplyr::filter(sim_beta_df_group, Group %in% redFamilyNames), 
            aes( x = Ao, y = Eo)) + 
  geom_density_2d_filled( stat = "density_2d_filled", h = c(1, 0.2),
                          show.legend = F) +
  xlim(1.25, 5.5) + 
  ylim(-0.1, 0.8) + 
  facet_wrap(vars(Group), nrow = 3, ncol = 5)
print(d)  

# fit MVN using method of moments
getquantile <- function(Ao, Eo) {
  mu <- c(mean(Ao), mean(Eo))
  sigma <- matrix(c(var(Ao), cov(Ao, Eo), cov(Ao, Eo), var(Eo)), nrow = 2,ncol = 2, byrow = T)
  quant <- qmvnorm(p = 0.05, tail = "upper.tail",
                   mean = mu,
                   sigma = sigma)
  return(quant$quantile)
}

ggplot(data = dplyr::filter(sim_beta_df_group, Group %in% redFamilyNames),
       aes(x = V, y = Eo, col = Group)) +
  stat_ellipse(type = "norm",
               level = 0.9,
               linewidth = 1.5) +
  scale_fill_viridis_d(option = "turbo")
  


