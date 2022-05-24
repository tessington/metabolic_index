# Code to test TMB fitting
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)

### Generate Evolutionary Trait Structure ####

#### Get taxonomy tree ####
all.dat <- readRDS(file = "data/alldata_taxonomy.RDS")
# remove data where there is no body size data
all.dat <- dplyr::filter(all.dat, !is.na(W))
naIndex <- which(is.na(all.dat$Species))
for (i in 1:length(naIndex)) all.dat$Species[naIndex[i]] <- paste0(all.dat$Genera[naIndex[i]], " spc")



# for debugging, how many unique sizes per taxa

nsizes <- all.dat %>%
  group_by(Species) %>%
  summarise(nsize = length(unique(W)), sizerange = log(max(W) / min(W)))
          
  

#### Setup TMB data and parameters ####
### Create new ParentChild matrix for reduced taxonomic structure
kb <-  8.617333262145E-5
tref <- 15
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$minuslogpo2 <- - log(all.dat$Pcrit)

Z_ik_est <- dplyr::select(all.dat, Class, Order, Species)
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
Z_ik_dat <- dplyr::select(all.dat, Class, Order, Species)
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
             minuslogpo2 = all.dat$minuslogpo2,
             spc_in_PCgz = spc_in_PC_gz -1
             
)

parameters = list(alpha_j = rep(0,n_j),
                  L_z = rep(1, length(L_z)),
                  cov_logmult_z = rep(0, length(unique(PC_gz_tmb[,2])) -1),
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
spc_ij_mle <- matrix(spc_parameters, nrow = n_i, ncol = 3, byrow = F)
spc_ij_se <- matrix(spc_parameters_se, nrow = n_i, ncol = 3, byrow = F)

re <- summary(rep, "random")
beta_mle <- matrix(re[grep(rownames(re), pattern = "beta"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta"),2], nrow = n_g, ncol = 3, byrow = F)


# make dataframe to hold everything
datatest <- tibble(Species = NA,
                   estAo = spc_ij_mle[,3],
                  estn = spc_ij_mle[,1],
                  estEo = spc_ij_mle[,2],
                  Aose = spc_ij_se[,3],
                  nse = spc_ij_se[,1],
                  Eose = spc_ij_se[,2],
                  nsize = nsizes$nsize,
                  sizerange = nsizes$sizerange
                  )




ClassEst <- tibble(Class = ParentChild_gz_est$ChildName[1:11],
                   Aomle = -beta_mle[1:11,2],
                   Aose = beta_se[1:11,2],
                   Eomle = -beta_mle[1:11,3],
                   Eose = beta_se[1:11,3],
                   nmle = -beta_mle[1:11,1],
                   nse = beta_se[1:11, 1])

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
actinIndex <- grep(x = ParentChild_gz_est$ParentName, pattern = "\\bActinopteri\\b")
longOrderNames <- ParentChild_gz_est$ChildName[actinIndex]

# Remove all before and up to ":":
OrderNames <- gsub(".*_","",longOrderNames)

ActinEst <- tibble(Order = OrderNames,
                   Aomle = -beta_mle[actinIndex,2],
                   Aose = beta_se[actinIndex,2],
                   Eomle = -beta_mle[actinIndex,3],
                   Eose = beta_se[actinIndex,3],
                   nmle = -beta_mle[actinIndex,1],
                   nse = beta_se[actinIndex,1]
)

Aoplot <- plotest(ActinEst, "Ao", "Order")
Eoplot <- plotest(ActinEst, "Eo", "Order")
nplot <- plotest(ActinEst, "n", "Order")

grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)

# plot all species
longSpeciesNames <- ParentChild_gz_est$ChildName[ParentChild_gz_est[,'ChildTaxon']==3]
SpeciesNames <- gsub(".*_","",longSpeciesNames)
SpeciesEst <- tibble(Species = SpeciesNames,
                 Aomle =-spc_ij_mle[,2],
                 Aose = spc_ij_se[,2],
                 Eomle = -spc_ij_mle[,3],
                 Eose = spc_ij_se[,3],
                 nmle = -spc_ij_mle[,1],
                 nse = spc_ij_se[,1]
)
Aoplot <- plotest(SpeciesEst, "Ao", "Species")
Eoplot <- plotest(SpeciesEst, "Eo", "Species")
nplot <- plotest(SpeciesEst, "n", "Species")

grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)

