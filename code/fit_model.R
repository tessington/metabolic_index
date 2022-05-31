# Code to test TMB fitting
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)


### Functions ####
find_index <- function(x,y) y <- which(y == x)
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

### Generate Evolutionary Trait Structure ####

#### Get taxonomy tree ####
all.dat <- readRDS(file = "data/alldata_taxonomy.RDS")
# remove data where there is no body size data
all.dat <- dplyr::filter(all.dat, !is.na(W))
naIndex <- which(is.na(all.dat$Species))
for (i in 1:length(naIndex)) all.dat$Species[naIndex[i]] <- paste0(all.dat$Genera[naIndex[i]], " spc")


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
### Setup TMB data and parameters ####
#### Create new ParentChild matrix for reduced taxonomic structure ####
kb <-  8.617333262145E-5
tref <- 15
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$minuslogpo2 <- - log(all.dat$Pcrit)
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

ClassEst <- tibble(Class = ParentChild_gz$ChildName[1:11],
                   Aomle = beta_mle[1:11,2],
                   Aose = beta_se[1:11,2],
                   Eomle = beta_mle[1:11,3],
                   Eose = beta_se[1:11,3],
                   nmle = beta_mle[1:11,1],
                   nse = beta_se[1:11, 1])



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

# plot all species
longSpeciesNames <- ParentChild_gz$ChildName[ParentChild_gz[,'ChildTaxon']==4]
SpeciesNames <- gsub(".*_","",longSpeciesNames)
SpeciesEst <- tibble(Species = SpeciesNames,
                 Aomle =spc_ij_mle[,2],
                 Aose = spc_ij_se[,2],
                 Eomle = spc_ij_mle[,3],
                 Eose = spc_ij_se[,3],
                 nmle = spc_ij_mle[,1],
                 nse = spc_ij_se[,1]
)
Aoplot <- plotest(SpeciesEst, "Ao", "Species")
Eoplot <- plotest(SpeciesEst, "Eo", "Species")
nplot <- plotest(SpeciesEst, "n", "Species")

grid.arrange(Aoplot, nplot, Eoplot, ncol = 3)

