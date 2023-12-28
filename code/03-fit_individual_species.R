### Script to run through all species in dataset, and fit linear model if there is sufficient data

source("code/fit_model_funs.R")
#### Get all data ####
all.dat <- load_data()

#### Load hierarchical parameters from bigger fit ####
SpeciesEst <- readRDS(file = "analysis/hierarchical_species_estimates.RDS")

# make change to species name to match dataframe
SpeciesEst$Species[SpeciesEst$Species == "Doryteuthis (Amerigo) pealeii"] <- "Doryteuthis pealeii"

kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W / wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$minuslogpo2 <- - log(all.dat$Pcrit) # fit using pO2 in atm



spc.list <- unique(all.dat$scientific.name)
spc.est <- tibble(Species = spc.list,
                        logAo = NA,
                        logAose = NA,
                        Eo = NA,
                        Eose = NA)


# generally, the experiments are not well designed to separate out the effects of temp and body size, so 
# add an offset so that the Ao are comparable and body-size adjusted, but otherwise don't estimate n

for (i in 1:length(spc.list)) {
  red.data <- dplyr::filter(all.dat, scientific.name == spc.list[i])
  n_w <- length(unique(red.data$W))
  n_t <- length(unique(red.data$inv.temp))
  lowest.taxon <- unique(red.data$lowest.taxon)
  
  # only estimate if there is more than 1 temperature and is a species
  if(n_t >2 & lowest.taxon == "species") {
          n.2.use <- SpeciesEst$n[tolower(SpeciesEst$Species) == tolower(spc.list[i])]
          nlogw <- n.2.use * log(red.data$W)
          m <- lm(minuslogpo2 ~ inv.temp + offset(nlogw), data = red.data)
          msum <- summary(m)
          coefs <- msum$coefficients
          spc.est$logAo[i] <- coefs[1,1]
          spc.est$logAose[i] <- coefs[1,2]
          spc.est$Eo[i] <- coefs[2,1]
          spc.est$Eose[i] <- coefs[2,2]
  }
}


# Save result into RDS file
saveRDS(object = spc.est, file = "analysis/ind_species_estimates.RDS")

  
  