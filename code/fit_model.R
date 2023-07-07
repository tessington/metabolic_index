# Code to test TMB fitting
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)
library(Matrix)
library(mvtnorm)

### Functions ####
find_index <- function(x,y) y <- which(y == x)
make_df_plot <- function(level, beta_mle, beta_se, ParentChild_gz, groups) {
  
  group_index <- which(ParentChild_gz$ChildTaxon==level)
  groupname <- groups[level]
  GroupNames <- gsub(".*_","",ParentChild_gz$ChildName[group_index])
  
  
  Est <- tibble("{groupname}" := GroupNames,
                Ao = beta_mle[group_index,1],
                Aomin = beta_mle[group_index,1] - beta_se[group_index,1],
                Aomax = beta_mle[group_index,1]  + beta_se[group_index,1],
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
    geom_errorbar(aes(y = !!groupname,
                      xmin = !!xmin,
                      xmax = !!xmax)
    )
  return(groupplot)
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
  summarise(NoOrder = length(unique(Order)))
print(class_summary, n = 30)

order_summary <- all.dat %>%
  group_by(Order) %>%
  summarise(NoFamily = length(unique(Family)))
print(order_summary, n = 50)

family_summary <- all.dat %>%
  group_by(Family) %>%
  summarise(NoGenus = length(unique(Genera)))
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
all.dat$Pcrit_atm <- all.dat$Pcrit / 101.325 # convert from KPa to atm
all.dat$minuslogpo2 <- - log(all.dat$Pcrit_atm) 
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


ClassEst <- make_df_plot(level = 1, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

Est.2.plot <- ClassEst
Aoplot <- plotest(Est.2.plot, Ao, Class, Aomin, Aomax)
Eoplot <- plotest(Est.2.plot, Eo, Class, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Class, nmin, nmax)
grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)

## Plot Orders #####
OrderEst <- make_df_plot(level = 2, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)

Est.2.plot <- OrderEst
Aoplot <- plotest(Est.2.plot, Ao, Order, Aomin, Aomax)
Eoplot <- plotest(Est.2.plot, Eo, Order, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Order, nmin, nmax)
grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)


### Plot Families ####
FamilyEst <- make_df_plot(level = 3, 
                         beta_mle,
                         beta_se,
                         ParentChild_gz,
                         groups = taxa.list)
Est.2.plot <- FamilyEst
Aoplot <- plotest(Est.2.plot, Ao, Family, Aomin, Aomax)
Eoplot <- plotest(Est.2.plot, Eo, Family, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Family, nmin, nmax)
grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)


### Plot Species ####
SpeciesEst <- make_df_plot(level = 4, 
                          beta_mle,
                          beta_se,
                          ParentChild_gz,
                          groups = taxa.list)
Est.2.plot <- SpeciesEst
Aoplot <- plotest(Est.2.plot, Ao, Species, Aomin, Aomax)
Eoplot <- plotest(Est.2.plot, Eo, Species, Eomin, Eomax)
nplot <- plotest(Est.2.plot, n, Species, nmin, nmax)
grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)



### Compare fits to Penn et al. ####
#### Convert Ao from kPa to atm
SpeciesEst$Ao_atm = spc_ij_mle[,1] 
SpeciesEst$Ao_atm_SE = spc_ij_se[,1] 


#### Load previous fits ####
penn.fits <- read_xlsx(path = "data/MI_traits_Est_by_Curtis.xlsx") 

# reload the all.dat to be able to retrieve APHIAID for matching
all.dat <- readRDS(file = "data/alldata_taxonomy.RDS")
SpeciesEst$PennAo <- NA
SpeciesEst$PennEo <- NA

unique.species <- unique(SpeciesEst$Species)
for (i in 1:length(unique.species)) {
  spc.2.use <- tolower(unique.species[i])
  all.dat.index <- which(all.dat$scientific.name == spc.2.use) [ 1 ]
  aphiaID <- all.dat$AphiaID[all.dat.index]
  penn.index <- which(penn.fits$AphiaID == aphiaID)
  if (length(penn.index) ==1 ) {
    SpeciesEst$PennAo[i] <- 1 / penn.fits$`Vh (atm)`[penn.index]
    SpeciesEst$PennEo[i] <- penn.fits$`Eo (eV)`[penn.index]
  }
}
ggplot(SpeciesEst, aes(x = PennAo, y = Ao)) +
  geom_point()

ggplot(SpeciesEst, aes(x = PennEo, y = Eo)) +
  geom_point()

### Make Sigma ####
#### Make empty Cholesky Matrix ####
L <- matrix(0, nrow = n_j, ncol = n_j)
# extract L_z
fixed_names <- rownames(fixef)

L_z <- fixef[grep("L_z", fixed_names),1]
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

##### Calculate Sigma and extract log_lambdas
sigma <- L %*% t(L)
log_lambda <- fixef[grep("log_lambda",fixed_names),]

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

n.sims <- 10000
sim_beta_actinop <- matrix(NA, nrow = n.sims, ncol = 3)
sim_beta_species <- matrix(NA, nrow = n.sims, ncol = 3)
actinop_index <- which(ClassEst$Class == "Actinopteri")
newpar = rmvnorm_prec( mu = obj$env$last.par.best,
                       prec = rep$jointPrecision, 
                       n.sims = n.sims,
                       random_seed = sample(1:1000000, 1))

for (i in 1:n.sims) {
  # simulate new parameters 
  newpar_i <- newpar[,i] 
  parnames <- names(obj$env$last.par.best)
  # extract the simulated beta_gj
  beta_gj_random <-newpar_i[grep(parnames, pattern = "beta_gj")]
  # assign beta_gj to Ao, no, and Eo matrixes
  Ao_sim <- beta_gj_random[1:n_g]
  no_sim <- beta_gj_random[(n_g +1): (2 * n_g) ]
  Eo_sim <- beta_gj_random[(2 * n_g + 1) : (3 * n_g)]
  
  # save results to sim_beta_actinop
  sim_beta_actinop[i,1] <- Ao_sim[actinop_index]
  sim_beta_actinop[i,2] <- no_sim[actinop_index]
  sim_beta_actinop[i,3] <- Eo_sim[actinop_index]
  
  # using random draw for Actinopteri, simulate mean for a random order, random family, and then the value for a given species
  sim_beta_actinop_ord <- sim_beta_actinop[i,] + rmvnorm(n = 1, mean = rep(0,3), sigma = exp(log_lambda[1]) * sigma)
  sim_beta_actinop_fam <- sim_beta_actinop_ord + rmvnorm(n = 1, mean = rep(0,3), sigma = exp(log_lambda[2]) * sigma)
  sim_beta_actinop_spc <- sim_beta_actinop_fam + rmvnorm(n = 1, mean = rep(0,3), sigma = exp(log_lambda[3]) * sigma)
  # Save the species result into matrix
  sim_beta_species[i,] <- sim_beta_actinop_spc
  #sim_beta_species[i,] <- c(Ao_sim[actinop_index], no_sim[actinop_index], Eo_sim[actinop_index])
}


hist(sim_beta_species[,3], breaks = 20)

sim_beta_df <- tibble(Ao = sim_beta_species[,1],
                           Eo = sim_beta_species[,3],
                           n = sim_beta_species[,2])


d <- ggplot(data = sim_beta_df, aes( x = Ao, y = Eo)) + 
  geom_density_2d_filled( stat = "density_2d_filled", h = c(0.75, 0.15),
                          show.legend = F) +
  xlim(1.5, 5.5) + 
  ylim(0, 0.8)
print(d)  
