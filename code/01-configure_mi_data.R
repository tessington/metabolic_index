library(tidyr)
library(dplyr)

rm(list = ls())
# save to CSV file  
all.dat <- read.csv(file = "data/allmidata_alphia.csv", header = T)



#Use the libraries
library(jsonlite) #https://cran.r-project.org/web/packages/jsonlite/
library(httr)

#Fill in the AphiaID you need

unique.taxa <- unique(all.dat$scientific.name)
all.dat$Phylum <- NA
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
  all.dat$Phylum[taxa.index] <- lookup.taxa(classificationTree, "Phylum")
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
nrow(all.dat)
# Add uncertain order for Capitellidae
all.dat$Order[which(all.dat$Family == "Capitellidae")] = "Capitellidae"

saveRDS(all.dat, file = "data/alldata_taxonomy.RDS")


# Get summary statistics ####
sum_by_species <- all.dat %>%
  filter(lowest.taxon == "species") %>%
  group_by(Species) %>%
  summarise(nstudy = length(unique(Source))) %>%
  group_by(nstudy) %>%
  summarise(nspecies = n())

sum_by_study <- all.dat %>%
  filter(lowest.taxon == "species") %>%
  group_by(Source) %>%
  summarise(nspc = length(unique(Species))) %>%
  group_by(nspc) %>%
  summarise(nstudies = n())

sum_by_temp_w <- all.dat %>%
  filter(lowest.taxon == 'species') %>%
  group_by(Species) %>%
  summarise(ntemp = length(unique(Temp)), nw = length(unique(W)))

hist_temp <- ggplot(sum_by_temp_w, aes(x = ntemp)) + 
  geom_histogram(bins = 30)
print(hist_temp)

hist_w <- ggplot(sum_by_temp_w, aes(x = nw)) + 
  geom_histogram()
print(hist_w)
