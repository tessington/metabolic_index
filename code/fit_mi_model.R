### Code to fit simplest bayesian model using data compiled by Deutsch et al 2020

# Install Libraries
library(tidyr)
library(dplyr)
library(rstan)
library(shinystan)

# Load data
thedata <- readRDS(file = "data/alldata_taxonomy.RDS")

Z_ik <- dplyr::select(thedata, Class,Order)
Z_ik <- unique(Z_ik, MARGIN = 1)

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


# Add top predictive
ParentChild_gz = rbind( ParentChild_gz, c("predictive", NA, NA, 1) )
for( colI in 2:ncol(Z_ik)) ParentChild_gz = rbind( ParentChild_gz, c(paste(rep("predictive",colI),collapse="_"), paste(rep("predictive",colI-1),collapse="_"), match(paste(rep("predictive",colI-1),collapse="_"),ParentChild_gz[,1]), colI) )
# Relabel
ParentChild_gz = data.frame( ParentChild_gz )
colnames(ParentChild_gz) = c("ChildName", "ParentName", "ParentRowNumber", "ChildTaxon")
ParentChild_gz[,'ParentRowNumber'] = as.numeric(as.character(ParentChild_gz[,'ParentRowNumber']))
ParentChild_gz[,'ChildTaxon'] = as.numeric(as.character(ParentChild_gz[,'ChildTaxon']))

# Identify location for every observation
Taxa_Names = apply( Z_ik, MARGIN=1, FUN=paste, collapse="_")
g_i = match( Taxa_Names, ParentChild_gz[,'ChildName'] )
n_k = ncol(Z_ik)
n_g = nrow(ParentChild_gz)



summary.table <- Z_ik %>%
   group_by(Class, Order, Family) %>%
   summarise(n = n())
View(summary.table)
### Remove data with missing bits and convert temperature

thedata <- thedata %>%
  filter(!is.na(W), !is.na(Temp))

 kb <-  8.617333262145E-5
 tref <- 15
 thedata$inv.temp <- (1 / kb)  * ( 1 / (thedata$Temp + 273.15) - 1 / (tref + 273.15))

 ### Prepare data for Stan ####
 N <- nrow(thedata)
 unique.species <- unique(thedata$AphiaID)
 n <- length(unique.species)
 spc <- rep(NA, times = N)
 for (i in 1:n)  spc[which(thedata$AphiaID==unique.species[i])] <- i
 logb <- log(thedata$W)
 y <- -log(thedata$Pcrit)
 stan.data <- list(N = N,
                   nspc = n,
                   y = y,
                   logb = logb,
                   inv_temp = thedata$inv.temp,
                   spc = spc
                   )
 
 # Function to initialize chains 
 init.chain <- function(chain_id = 1, n) {
   list(
     logmu_Ao = runif(1, -1, 1),
     logsigma_Ao = runif(1, 0, 2),
     logAo_raw = rnorm(n, 0,1),
     mu_Eo = runif(1, -1, 1),
     sigma_Eo = runif(1, 0, 2),
     Eo_raw = rnorm(n, 0,1),
     mu_n = runif(1, -1, 1),
     sigma_n = runif(1, 0,1),
     n_raw = rnorm(n, 0,1),
     logsigma = runif(1, -1,1)
   )
 }
 n_chains <- 3
 init_ll <- lapply(1:n_chains, function(id)
   init.chain(chain_id = id, n))
   
### Which parameters to monitor
 stan.pars <- c("logmu_Ao","mu_Eo", "mu_n", "logsigma_Ao", "sigma_Eo", "sigma_n","logAo", "Eo", "n")
 
 stan.model <- "code/fit_mi_varEo.stan"
 
 ### Run the Stan Model
 rstan_options(auto_write = TRUE)  # this option stops Stan from re-compiling if not necessary
 options(mc.cores = parallel::detectCores()) # this is nice because it allows each chain to be run in parallel on a separate core of your processor
 niters <- 1000 # how long should each chain be.  1000 is probably fine
 
 fit <- stan(
   file = stan.model,
   data = stan.data,
   iter = niters,
   pars = stan.pars,
   warmup = floor(niters / 2),
   chains = n_chains,
   thin = 1,
   algorithm = 'NUTS',
   init = init_ll,
   verbose = FALSE,
   control = list(adapt_engaged = TRUE, adapt_delta = 0.8, max_treedepth = 10)
 )
 
                  
   