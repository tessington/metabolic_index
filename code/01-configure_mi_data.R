library(tidyverse)

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
all.dat$Order[which(all.dat$Family == "Capitellidae")] = "Capitellida"
# create Genus spp. names for taxa identified to genus
no_spc <- which(is.na(all.dat$Species))
for (i in no_spc) {
  all.dat$Species[i] <- paste(all.dat$Genera[i], "spp.")
}


saveRDS(all.dat, file = "data/alldata_taxonomy.RDS")

all.dat <- readRDS(file ="data/alldata_taxonomy.RDS")

# Remove unknown methods
rem_unknown <- T
if (rem_unknown) all.dat <- dplyr::filter(all.dat, !Method == "unknown")


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
  summarise(nspc = length(unique(Species)), 
            ngenus = length(unique(Genera)),
            nfamily = length(unique(Family)),
            norder = length(unique(Order)),
            nclass = length(unique(Class))
            )

sum_by_author <- all.dat %>%
  filter(lowest.taxon == "species") %>%
  group_by(SharedAuthor) %>%
  summarise(nspc = length(unique(Species)), 
            ngenus = length(unique(Genera)),
            nfamily = length(unique(Family)),
            norder = length(unique(Order)),
            nclass = length(unique(Class))
  )

spc_by_study <- sum_by_study %>%
  group_by(nspc) %>%
  summarise(nstudies = n())

genus_by_study <- sum_by_study %>%
  group_by(ngenus) %>%
  summarise(nstudies = n())

family_by_study <- sum_by_study %>%
  group_by(nfamily) %>%
  summarise(nstudies = n())

order_by_study <- sum_by_study %>%
  group_by(norder) %>%
  summarise(nstudies = n())

class_by_study <- sum_by_study %>%
  group_by(nclass) %>%
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


# output for supplemental table
taxa_summary <- dplyr::select(all.dat, Source, Species, Phylum, Order, Family)
taxa_summary <- dplyr::distinct(taxa_summary)
  

# write to csv.file,
write.csv(taxa_summary, file = "analysis/taxa_source.csv", row.names = F)
