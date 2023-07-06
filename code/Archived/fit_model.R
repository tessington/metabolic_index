# Code to test TMB fitting
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)


### Functions ####
find_index <- function(x,y) y <- which(y == x)
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
                getReportCovariance = FALSE)
summary(rep, "fixed")
re <- summary(rep, "report")


spc_parameters <- re[grep(rownames(re), pattern  = "spc_ij"),1]
spc_parameters_se <-re[grep(rownames(re), pattern  = "spc_ij"),2]
spc_ij_mle <- matrix(spc_parameters, nrow = n_i, ncol = 3, byrow = F)
spc_ij_se <- matrix(spc_parameters_se, nrow = n_i, ncol = 3, byrow = F)


re <- summary(rep, "random")
beta_mle <- matrix(re[grep(rownames(re), pattern = "beta"),1], nrow = n_g, ncol = 3, byrow = F)
beta_se <- matrix(re[grep(rownames(re), pattern = "beta"),2], nrow = n_g, ncol = 3, byrow = F)


### Plot Estimates ####

order_index <- which(ParentChild_gz$ChildTaxon==1)
OrderEst <- tibble(Order = ParentChild_gz$ChildName[order_index],
                   Ao = beta_mle[order_index,1],
                   Aomin = beta_mle[order_index,1] - beta_se[order_index,1],
                   Aomax = beta_mle[order_index,1]  + beta_se[order_index,1],
                   Eo = beta_mle[order_index,3],
                   Eomin = beta_mle[order_index,3] - beta_se[order_index,3],
                   Eomax = beta_mle[order_index,3] + beta_se[order_index,3],
                   n = beta_mle[order_index,2],
                   nmin = beta_mle[order_index,2] - beta_se[order_index,2],
                   nmax = beta_mle[order_index,2] + beta_se[order_index,2],
                   )



Aoplot <- plotest(OrderEst, Ao, Order, Aomin, Aomax)
Eoplot <- plotest(OrderEst, Eo, Order, Eomin, Eomax)
nplot <- plotest(OrderEst, n, Order, nmin, nmax)

grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)


### Plot Families ####
#actinIndex <- grep(x = ParentChild_gz$ParentName, pattern = "\\bActinopteri\\b")
longFamilyNames <- ParentChild_gz$ChildName[ParentChild_gz[,'ChildTaxon']==2]

# Remove all before and up to ":":
FamilyNames <- gsub(".*_","",longFamilyNames)
f_index <- which(ParentChild_gz$ChildTaxon ==2)
Est <- tibble(Family = FamilyNames,
              Ao = beta_mle[f_index,1],
              Aomin = beta_mle[f_index,1] - beta_se[f_index,1],
              Aomax = beta_mle[f_index,1]  + beta_se[f_index,1],
              Eo = beta_mle[f_index,3],
              Eomin = beta_mle[f_index,3] - beta_se[f_index,3],
              Eomax = beta_mle[f_index,3] + beta_se[f_index,3],
              n = beta_mle[f_index,2],
              nmin = beta_mle[f_index,2] - beta_se[f_index,2],
              nmax = beta_mle[f_index,2] + beta_se[f_index,2],
)

Aoplot <- plotest(Est, Ao, Family, Aomin, Aomax)
Eoplot <- plotest(Est, Eo, Family, Eomin, Eomax)
nplot <- plotest(Est, n, Family, nmin, nmax)

grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)

# plot all species
longSpeciesNames <- ParentChild_gz$ChildName[ParentChild_gz[,'ChildTaxon']==3]
SpeciesNames <- gsub(".*_","",longSpeciesNames)
Est <- tibble(Species = SpeciesNames,
                 Ao =spc_ij_mle[,1],
                 Aomin = spc_ij_mle[,1] - spc_ij_se[,1] ,
                 Aomax = spc_ij_mle[,1] + spc_ij_se[,1] ,
                 Eo = spc_ij_mle[,3],
                 Eomin = spc_ij_mle[,3] - spc_ij_se[,3],
                 Eomax = spc_ij_mle[,3] + spc_ij_se[,3],
                 n = spc_ij_mle[,2],
                 nmin = spc_ij_mle[,2] - spc_ij_se[,2],
                 nmax = spc_ij_mle[,2] + spc_ij_se[,2],
)
Aoplot <- plotest(Est, Ao, Species, Aomin, Aomax)
Eoplot <- plotest(Est, Eo, Species, Eomin, Eomax)
nplot <- plotest(Est, n, Species, nmin, nmax)


grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)


### Compare fits to Penn et al. ####
#### Convert Ao from kPa to atm
SpeciesEst <- Est
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
# Make Sigma
L <- matrix(0, nrow = n_k, ncol = n_k)
# extract L_z
fixef <- summary(rep, "fixed")
fixed_names <- rownames(fixef)

L_z <- fixef[L_z <- grep("L_z", fixed_names),]
Count = 1
for (i in 1 : n_k) {
  for (j in 1 : n_k) {
    if (i>=j) {
      L[i,j] <- L_z[Count,1]
      Count <- Count + 1
    }
  }
}

sigma <- L %*% t(L)
