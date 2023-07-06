# Code to test TMB fitting
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)

### Generate Evolutionary Trait Structure ####
set.seed(123)
#### Get taxonomy tree ####
all.dat <- readRDS(file = "data/alldata_taxonomy.RDS")
# remove data where there is no body size data
all.dat <- dplyr::filter(all.dat, !is.na(W))

all.dat <- dplyr::select(all.dat, Temp, W, Pcrit, Class, Order, Family, Genera, Species)
naIndex <- which(is.na(all.dat$Species))
for (i in 1:length(naIndex)) all.dat$Species[naIndex[i]] <- paste0(all.dat$Genera[naIndex[i]], " spc")
Z_ik <- dplyr::select(all.dat, Class, Order, Family, Genera, Species)
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
#ParentChild_gz = rbind( ParentChild_gz, c("predictive", NA, NA, 1) )
#for( colI in 2:ncol(Z_ik)) ParentChild_gz = rbind( ParentChild_gz, c(paste(rep("predictive",colI),collapse="_"), paste(rep("predictive",colI-1),collapse="_"), match(paste(rep("predictive",colI-1),collapse="_"),ParentChild_gz[,1]), colI) )
# Relabel
ParentChild_gz = data.frame( ParentChild_gz )
colnames(ParentChild_gz) = c("ChildName", "ParentName", "ParentRowNumber", "ChildTaxon")
ParentChild_gz[,'ParentRowNumber'] = as.numeric(as.character(ParentChild_gz[,'ParentRowNumber']))
ParentChild_gz[,'ChildTaxon'] = as.numeric(as.character(ParentChild_gz[,'ChildTaxon']))
PC_gz = as.matrix(ParentChild_gz[, c('ParentRowNumber', 'ChildTaxon')]) 

# Identify location for every observation
Taxa_Names = apply( Z_ik, MARGIN=1, FUN=paste, collapse="_")
g_i = match( Taxa_Names, ParentChild_gz[,'ChildName'] )
n_k = ncol(Z_ik)
n_j = 3 # three traits
n_g = nrow(ParentChild_gz)
n_i <- length(g_i)


#### Specify all parameters ####
alpha_j <- c(0, 2, 0.8) # means above class
# lambda - relative covariance w/in groups
lambda_g <- c(1, .8, 0.6, 0.6, 0.4)
# Specify  variance covariance and generate L_z
sigma1 <- 0.1
sigma2 <- .5
sigma3 <- 0.2
rho12 <- 0.5
rho13 <- 0.25
rho23 <- 0.5
L_z <- c(sigma1, 
         sigma1 * sigma2 * rho12, 
         sigma2,
         sigma1 * sigma3 * rho13, 
         sigma2 * sigma3 * rho23, 
         sigma3)


### Create evolutionary covariance matrix
cov_matrix <- function(L_val, logmult_col, min_var, n_rows, n_cols, invertTF ){
  L_rc<- matrix(0, n_rows, abs(n_cols))
  Cov_rr <- matrix(0, n_rows, n_rows);
  Return_rrN <- matrix(0, n_rows, n_rows);
  
  # Loadings matrix with zero upper-diagonal
  Count = 1
    for(r in 0:(n_rows -1)){
      for(c in 0:(n_cols -1)){
        if(r>=c){
          L_rc[r+1,c+1] = L_val[Count]
          Count <- Count +1 
        }else{
          L_rc[r+1,c+1] = 0.0
        }
      }
      
    }

  
  ## Additive constant on diagonal
  for(r in 0:(n_rows -1)){
    Cov_rr[r+1,r+1] = Cov_rr[r+1,r+1] + min_var # Necesary to prevent crashes during innner optimizer when using SR data
  }
  ## Combine and return
  Cov_rr = Cov_rr + L_rc %*% t(L_rc)
  if(invertTF==FALSE) Return_rr = Cov_rr
  if(invertTF==TRUE) Return_rr = atomic::matinv( Cov_rr )
  return (Return_rr)
}


beta_gj <- matrix(0, n_g, n_j)
Cov_z <- cov_matrix(L_z, logmult_col = c(0,0,0), min_var = 0.01,  n_rows = 3, n_cols = 3, invertTF = FALSE)
for(g in 1:n_g){
  if( PC_gz[g,2]==1 ) {
    beta_gj[g,] <- alpha_j + mvrnorm(n = 1, mu = rep(0, 3), Sigma = Cov_z)
  }
  if( PC_gz[g,2]>=2 ) {
    tmpCov_z <- Cov_z * lambda_g[PC_gz[g,2]]
    beta_gj[g,] = beta_gj[PC_gz[g,1],] + mvrnorm(n = 1, mu = rep(0, 3), Sigma = tmpCov_z)
  }
}

### Extract trait values for each taxa ####
Y_ij <- Yhat_ij <- matrix(0,  n_i, n_j )
for(i in 1:n_i)  Y_ij[i,] <- beta_gj[g_i[i],]




### Now make experimental data from traits
spc.na <- which(is.na(all.dat$species))
for (i in 1:length(spc.na)) all.dat$Species[spc.na[i]] <- paste(all.dat$Genera[spc.na[i]], "spc")

Z_ik_dat <- dplyr::select(all.dat, Class, Order, Family, Genera, Species)
Taxa_Names_dat <-  apply( Z_ik_dat, MARGIN=1, FUN=paste, collapse="_")
g_i_dat = match( Taxa_Names_dat, ParentChild_gz[,'ChildName'] )

find_index <- function(x,y) y <- which(y == x)

# this is the index of which element of Y_ij to use for each row of the data i
g_i_i <- sapply(FUN = find_index, X = g_i_dat, y = g_i)

obs_sigma <- 0.05 # dunno, why not
kb <-  8.617333262145E-5
tref <- 15

all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
n_data <- nrow(all.dat)

mu <- rep(NA, times = n_data)

for (i in 1:n_data) mu[i] <- (Y_ij[g_i_i[i],1] * log(all.dat$W[i]) +
  Y_ij[g_i_i[i],3] * all.dat$inv.temp[i] + 
    Y_ij[g_i_i[i],2] )

all.dat$minusLogPo2 <- rnorm(n_data, mu, obs_sigma)

# for debugging, how many unique sizes per taxa

nsizes <- all.dat %>%
  group_by(Species) %>%
  summarise(nsize = length(unique(W)), sizerange = log(max(W) / min(W)))
          
  



#### Setup TMB data and parameters ####
### Create new ParentChild matrix for reduced taxonomic structure

Z_ik_est <- dplyr::select(all.dat, Class, Order, Family, Species)
Z_ik_est <- unique(Z_ik_est, MARGIN = 1)
ParentChild_gz_est = NULL
# 1st column: child taxon name
# 2nd column: parent taxon name
# 3rd column: parent row-number in ParentChild_gz
# 4th column: Taxon level
# Loop through
for( colI in 1:ncol(Z_ik_est)){
  Taxa_Names = apply( Z_ik_est[,1:colI,drop=FALSE], MARGIN=1, FUN=paste, collapse="_")
  Unique_Taxa = unique(Taxa_Names)
  for( uniqueI in 1:length(Unique_Taxa) ){
    Which = which( Taxa_Names == Unique_Taxa[uniqueI] )
    if( colI==1 ){
      ParentChild_gz_est = rbind( ParentChild_gz_est, c(Unique_Taxa[uniqueI], NA, NA, colI) )
    }else{
      if( length(unique(Z_ik_est[Which,colI-1]))>1 ) stop("Taxa has multiple parents")
      ChildName = Unique_Taxa[uniqueI]
      ParentName = paste(rev(rev(strsplit(ChildName,"_")[[1]])[-1]),collapse="_")
      ParentChild_gz_est = rbind( ParentChild_gz_est, c(ChildName, ParentName, match(ParentName,ParentChild_gz_est[,1]), colI) )
    }
  }
}

# Add top predictive

#ParentChild_gz_est = rbind( ParentChild_gz_est, c("predictive", NA, NA, 1) )
#for( colI in 2:ncol(Z_ik_est)) ParentChild_gz_est = rbind( ParentChild_gz_est, c(paste(rep("predictive",colI),collapse="_"), paste(rep("predictive",colI-1),collapse="_"), match(paste(rep("predictive",colI-1),collapse="_"),ParentChild_gz_est[,1]), colI) )
# Relabel
ParentChild_gz_est = data.frame( ParentChild_gz_est )
colnames(ParentChild_gz_est) = c("ChildName", "ParentName", "ParentRowNumber", "ChildTaxon")
ParentChild_gz_est[,'ParentRowNumber'] = as.numeric(as.character(ParentChild_gz_est[,'ParentRowNumber']))
ParentChild_gz_est[,'ChildTaxon'] = as.numeric(as.character(ParentChild_gz_est[,'ChildTaxon']))
PC_gz_tmb <- as.matrix(ParentChild_gz_est[, c('ParentRowNumber', 'ChildTaxon')]) - 1
# Identify location for every observation
Taxa_Names = apply( Z_ik_est, MARGIN=1, FUN=paste, collapse="_")
g_i = match( Taxa_Names, ParentChild_gz_est[,'ChildName'] )
n_k = ncol(Z_ik_est)
n_j = 3 # three traits
n_g = nrow(ParentChild_gz_est)
n_i <- length(g_i)

#### Redo index of data to Parent - Child ####
Z_ik_dat <- dplyr::select(all.dat, Class, Order, Family, Species)
Taxa_Names_dat <-  apply( Z_ik_dat, MARGIN=1, FUN=paste, collapse="_")
g_i_dat = match( Taxa_Names_dat, ParentChild_gz_est[,'ChildName'] )


g_i_i <- sapply(FUN = find_index, X = g_i_dat, y = g_i)

# make index that identifies which of the rows of PC_gz correspond to species
spc_in_PC_gz <- which(PC_gz_tmb[,2] == max(PC_gz_tmb[,2]))

data <- list(PC_gz = PC_gz_tmb,
             g_i = g_i - 1,
             invtemp = all.dat$inv.temp,
             logW = log(all.dat$W),
             taxa_id = g_i_i -1,
             minuslogpo2 = all.dat$minusLogPo2,
             spc_in_PCgz = spc_in_PC_gz -1
             
)

parameters = list(alpha_j = rep(0,n_j),
                  L_z = rep(1, length(L_z)),
                  log_lambda = rep(0, length(unique(PC_gz_tmb[,2])) -1),
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
                getReportCovariance = FALSE)
re <- summary(rep, "report")
spc_parameters <- re[grep(rownames(re), pattern  = "spc_ij"),1]
spc_parameters_se <-re[grep(rownames(re), pattern  = "spc_ij"),2]
re <- summary(rep, "random")
spc_ij_mle <- matrix(spc_parameters, nrow = n_i, ncol = 3, byrow = F)
spc_ij_se <- matrix(spc_parameters_se, nrow = n_i, ncol = 3, byrow = F)

beta_mle <- matrix(re[grep(rownames(re), pattern = "beta"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta"),2], nrow = n_g, ncol = 3, byrow = F)


# make dataframe to hold everything
datatest <- tibble(trueAo = Y_ij[,3],
                  truen = Y_ij[,1],
                  trueEo = Y_ij[,2],
                  estAo = spc_ij_mle[,3],
                  estn = spc_ij_mle[,1],
                  estEo = spc_ij_mle[,2],
                  Aose = spc_ij_se[,3],
                  nse = spc_ij_se[,1],
                  Eose = spc_ij_se[,2],
                  nsize = nsizes$nsize,
                  sizerange = nsizes$sizerange
                  )


aoplot <- ggplot(datatest, aes(x = trueAo, y = estAo)) +
  geom_point() +
  geom_errorbar(aes(x = trueAo,
                ymin = estAo - Aose,
                ymax = estAo + Aose)) + 
  geom_abline(slope = 1,
              intercept = 0)

nplot <- ggplot(datatest, aes(x = truen, y = estn, col = sizerange)) +
  geom_point() + 
  scale_color_continuous() +
  geom_errorbar(aes(x = truen,
                    ymin = estn - nse,
                    ymax = estn + nse)) + 
  geom_abline(slope = 1,
              intercept = 0)
eoplot <- ggplot(datatest, aes(x = trueEo, y = estEo)) +
  geom_point() +
  geom_errorbar(aes(x = trueEo,
                    ymin = estEo - Eose,
                    ymax = estEo + Eose)) + 
  geom_abline(slope = 1,
              intercept = 0)
grid.arrange(aoplot, nplot, eoplot, nrow = 1)


ClassEst <- tibble(Class = ParentChild_gz$ChildName[1:12],
                   Aomle = beta_mle[1:12,2],
                   Aose = beta_se[1:12,2],
                   Eomle = beta_mle[1:12,3],
                   Eose = beta_se[1:12,3],
                   nmle = beta_mle[1:12,1],
                   nse = beta_se[1:12, 1])

plotest <- function(dataest, trait, groupname) {
  # make min and max
  eval(parse(text = paste0("dataest$min <- dataest$", trait,"mle - 
  dataest$",trait,"se")))
  eval(parse(text = paste0("dataest$max <- dataest$", trait,"mle + 
  dataest$",trait,"se")))
  
                   
groupplot <- ggplot(data = dataest, aes_string(x = paste0(trait, "mle"), y = groupname)) +
  geom_point() +
  geom_errorbar(aes_string(y = groupname,
                    xmin = "min",
                    xmax = "max")
  )
return(groupplot)
}

Aoplot <- plotest(ClassEst, "Ao", "Class")
Eoplot <- plotest(ClassEst, "Eo", "Class")
nplot <- plotest(ClassEst, "n", "Class")

grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)


### Plot actinopterygii ####
actinIndex <- grep(x = ParentChild_gz$ParentName, pattern = "\\bActinopteri\\b")
longOrderNames <- ParentChild_gz$ChildName[actinIndex]

# Remove all before and up to ":":
OrderNames <- gsub(".*_","",longOrderNames)

ActinEst <- tibble(Order = OrderNames,
                   Aomle = beta_mle[actinIndex,2],
                   Aose = beta_se[actinIndex,2],
                   Eomle = beta_mle[actinIndex,3],
                   Eose = beta_se[actinIndex,3],
                   nmle = beta_mle[actinIndex,1],
                   nse = beta_se[actinIndex,1]
)

Aoplot <- plotest(ActinEst, "Ao", "Order")
Eoplot <- plotest(ActinEst, "Eo", "Order")
nplot <- plotest(ActinEst, "n", "Order")

grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)

summary(rep, "fixed")
