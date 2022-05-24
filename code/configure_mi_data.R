library(tidyr)
library(R.matlab)
library(rMR)


rm(list = ls())



milist  <- readMat('data/Pdat.mat')$Pdat
speciesNames <- tolower(unlist(milist[[4]]))
allTemp <- milist[[1]]
allW <- milist[[2]]
allpcrit <- milist[[3]]
stripNa <- function(x) x[!is.nan(x)&x>0]
for (i in 1:length(speciesNames)) {
  species.index <- which(speciesNames == speciesNames[i])
  tempT <- stripNa(allTemp[species.index,])
  tempW <- stripNa(allW[species.index,])
  temppcrit <- stripNa(allpcrit[species.index,])
  if (length(tempT)>0) {
    # deal with issues of missing data
    maxN <- max(c(length(tempT), length(tempW), length(temppcrit)))
    if(length(tempT)==0) tempT <- rep(NA, maxN)
    if(length(tempW)==0) tempW <- rep(NA, maxN)
    if(length(temppcrit)==0) temppcrit <- rep(NA, maxN)
  if (!exists("all.dat")) {
    all.dat <- tibble(species = speciesNames[i],
                          Temp = tempT,
                          W = tempW,
                          Pcrit = temppcrit)
  } else {
    tmp.dat <- tibble(species = speciesNames[i],
                          Temp = tempT,
                          W = tempW,
                          Pcrit = temppcrit)
    all.dat <- rbind(all.dat, tmp.dat)
  }
  }
}

### Add in data that I compiled for Essington et al. 2022

newDat <- read.csv("data/metabolic_index_data_from_Essington.csv", header = T)

newDat$po2 <- NA
#### cycle by species and measurement type, and get po2 for each ####
spc.list <- unique(newDat$spc)
c.to.K <- function(x) x + 273.15
for (i in 1:length(spc.list)) {
  spc.index <- which(newDat$spc == spc.list[i])
  tmp.dat <- dplyr::filter(newDat, spc == spc.list[i])
  po2.tmp <- DO.unit.convert(tmp.dat$lc50,
                             DO.units.in = tmp.dat$units[1],
                             DO.units.out = "PP",
                             bar.press = 1,
                             bar.units.in= "atm",
                             temp.C = tmp.dat$temp,
                             bar.units.out = "kpa",
                             salinity = tmp.dat$salinity[1],
                             salinity.units = "uS")
  newDat$po2[spc.index] = po2.tmp
}

#### Merge into larger tibble ####
updated.dat <- tibble(species = tolower(newDat$spc),
                      Temp = newDat$temp,
                      W = newDat$b,
                      Pcrit = newDat$po2)
all.dat <- rbind(all.dat, updated.dat)

# save to CSV file
write.csv(x = all.dat, file = "data/allmidata.csv")

### Read in adjusted file that includes Alphia IDs for each species / genera
all.dat <- read.csv(file = "data/allmidata_alphia.csv", header = T)


install.packages("jsonlite", repos="http://cran.r-project.org")
install.packages("httr")

#Use the libraries
library(jsonlite) #https://cran.r-project.org/web/packages/jsonlite/
library(httr)

#Fill in the AphiaID you need

unique.taxa <- unique(all.dat$scientific.name)
all.dat$Class <- NA
all.dat$Subclass <- NA
all.dat$Superorder <- NA
all.dat$Order <- NA
all.dat$Suborder <- NA
all.dat$Infraorder <- NA
all.dat$Superfamily <- NA
all.dat$Family <- NA
all.dat$Genera <- NA
all.dat$Species <- NA

lookup.taxa <- function(x, level) {
  xunlist <- unlist(x)
  level.index <- which(xunlist == level)
  if(length(level.index)>0) {
    return(unname(xunlist[level.index + 1]))
  } else {return(NA)
  }
}

for (i in 1:length(unique.taxa)) {
  taxa <- unique.taxa[i]
  taxa.index <- which(all.dat$scientific.name == taxa)
  aphiaID <- all.dat$AphiaID[taxa.index[1]]
  #Build the URL to get the data from
  url <- sprintf("https://www.marinespecies.org/rest/AphiaClassificationByAphiaID/%d", aphiaID);
  #Get the actual data from the URL
  classificationTree <- fromJSON(url)
  all.dat$Class[taxa.index] <- lookup.taxa(classificationTree, "Class")
  all.dat$Subclass[taxa.index] <- lookup.taxa(classificationTree, "Subclass")
  all.dat$Superorder[taxa.index] <- lookup.taxa(classificationTree, "Superorder")
  all.dat$Order[taxa.index] <- lookup.taxa(classificationTree, "Order")
  all.dat$Suborder[taxa.index] <- lookup.taxa(classificationTree, "Suborder")
  all.dat$Infraorder[taxa.index] <- lookup.taxa(classificationTree, "Infraorder")
  all.dat$Superfamily[taxa.index] <- lookup.taxa(classificationTree, "Superfamily")
  all.dat$Family[taxa.index] <- lookup.taxa(classificationTree, "Family")
  all.dat$Genera[taxa.index] <- lookup.taxa(classificationTree, "Genus")
  if(all.dat$lowest.taxon[taxa.index[1]] =="species")  all.dat$Species[taxa.index] <- lookup.taxa(classificationTree, "Species")

}

saveRDS(all.dat, file = "data/alldata_taxonomy.RDS")
all.dat <- readRDS(file = "data/alldata_taxonomy.RDS")

