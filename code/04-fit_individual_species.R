### Script to run through all species in dataset, and fit linear model if there is sufficient data
library(dplyr)
source("code/fit_model_funs.R")
#### Get all data ####
all.dat <- load_data()

rem_styela <- T
if (rem_styela) all.dat <- all.dat %>%
  filter(Species != "Styela plicata")

rem_unknown <- T
all.dat$EstMethod_Metric <- tolower(all.dat$EstMethod_Metric)
if (rem_unknown) all.dat <- dplyr::filter(all.dat, !Method == "unknown", !EstMethod_Metric == "unknown")

# Keep only oxygen consumption methods
all.dat <- dplyr::filter(all.dat, Method == "OxygenConsumption")

#### Load hierarchical parameters from bigger fit ####
SpeciesEst <- readRDS(file = "analysis/hierarchical_species_estimates.RDS")


kb <-  8.617333262145E-5
tref <- 15
wref <- 5
all.dat$W <- all.dat$W / wref
all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1/(tref + 273.15))
all.dat$minuslogpo2 <- - log(all.dat$Pcrit) # fit using pO2 in atm



spc.list <- unique(all.dat$Species)
spc.est <- tibble(Species = spc.list,
                        logAo = NA,
                        logAose = NA,
                        Eo = NA,
                        Eose = NA)



# generally, the experiments are not well designed to separate out the effects of temp and body size, so 
# add an offset so that the Ao are comparable and body-size adjusted, but otherwise don't estimate n

for (i in 1:length(spc.list)) {
  red.data <- dplyr::filter(all.dat, Species == spc.list[i])
  n_w <- length(unique(red.data$W))
  n_t <- length(unique(red.data$inv.temp))
  lowest.taxon <- unique(red.data$lowest.taxon)
  
  # only estimate if there is more than 2 temperatures and is a species
  if(n_t >2 & lowest.taxon == "species") {
          n.2.use <- SpeciesEst$n[tolower(SpeciesEst$Species) == tolower(spc.list[i])]
          red.data$nlogw <- n.2.use * log(red.data$W)
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

### Alternative method : fit as random effects

# assign n for each species, even if data are poor
spc.index <- match(all.dat$Species, SpeciesEst$Species)
ns <- SpeciesEst$n[spc.index]
all.dat$nlogw = ns * log(all.dat$W)

library(lme4)
mod <- lmer(log(Pcrit) ~ 1 + EstMethod_Metric + offset(nlogw) +  inv.temp +  (1 + inv.temp|Species) , data = all.dat)
fixed_effect <- fixef(mod)
species_pars <- ranef(mod)
species_V = exp(fixed_effect["(Intercept)"] + species_pars$Species["(Intercept)"] )
species_Eo = - fixed_effect["inv.temp"] - species_pars$Species["inv.temp"]

species.index <- match(SpeciesEst$Species, rownames(species_Eo))

# plot using ggplot
library(ggplot2)

SpeciesEst$Vind = species_V[species.index,1]
SpeciesEst$Eoind = species_Eo[species.index,1]

## Make Plot ####
Vplot <- ggplot(SpeciesEst, aes(x = Vind, y = V)) +
  geom_point(size = 2.5) + 
  scale_y_continuous(limits = c(0, 8.5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 17), expand = c(0, 0)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 1.25) +
  xlab(("V (kPa) from random effect model")) +
  ylab("V (kPa) from taxonomic hierarchical model")

eoplot <- ggplot(SpeciesEst, aes(x = Eoind, y = Eo)) +
  geom_point(size = 2.5) + 
  geom_abline(intercept = 0, slope = 1, linewidth = 1.25) +
  xlab(expression("E"[o]~ "from random effect model")) +
  ylab(expression("E"[o]~ "from taxonomic hierarchical model"))
grid.arrange(Vplot, eoplot, ncol = 2)

