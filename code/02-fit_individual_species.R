### Script to run through all species in dataset, and fit linear model if there is sufficient data


#### Get all data ####
all.dat <- load_data()
# set W to 1 if NA


kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W / wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$Pcrit_atm<- all.dat$Pcrit / 101.325 # convert from kPa to atm
all.dat$minuslogpo2 <- - log(all.dat$Pcrit) # fit using pO2 in atm



spc.list <- unique(all.dat$scientific.name)
spc.est <- tibble(Species = spc.list,
                        logAo = NA,
                        logAose = NA,
                        Eo = NA,
                        Eose = NA,
                  ntemp = NA,
                  nw = NA)


# generally, the experiments are not well designed to separate out the effects of temp and body size, so 
# add an offset so that the Ao are comparable and body-size adjusted, but otherwise don't estimate n

for (i in 1:length(spc.list)) {
  red.data <- dplyr::filter(all.dat, scientific.name == spc.list[i])
  n_w <- length(unique(red.data$W))
  n_t <- length(unique(red.data$inv.temp))
  
  if(n_t >2) {
      # if no W are present, fit without body size
      if (any(is.na(unique(red.data$W)))) {
        m <- lm(minuslogpo2 ~ inv.temp , data = red.data)
      } else {
      # otherwise, use an offset
      red.data$nlogw <- -7.987253e-02 * log(red.data$W)
      m <- lm(minuslogpo2 ~ inv.temp + offset(nlogw), data = red.data)
      }
    msum <- summary(m)
    coefs <- msum$coefficients
    spc.est$logAo[i] <- coefs[1,1]
    spc.est$logAose[i] <- coefs[1,2]
    spc.est$Eo[i] <- coefs[2,1]
    spc.est$Eose[i] <- coefs[2,2]
    spc.est$ntemp[i] <- n_t
    spc.est$nw[i] <- n_w
  }
}

# Save result into RDS file
saveRDS(object = spc.est, file = "analysis/ind_species_estimates.RDS")

  
  