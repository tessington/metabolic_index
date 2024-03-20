# get fits from hierarchical model
SpeciesEst <- readRDS(file = "analysis/hierarchical_species_estimates.RDS")

SpeciesEst$Vind <- NA
SpeciesEst$Eoind <- NA
SpeciesEst$nind <- NA

## Load individual species fits ####
spc.fits <- readRDS(file = "analysis/ind_species_estimates.RDS")

## Combine into one df ####
for (i in 1:nrow(spc.fits)) {
  spc.2.use <- tolower(spc.fits$Species[i])
  spc.index <- which(tolower(SpeciesEst$Species) == spc.2.use)
  SpeciesEst$Vind[spc.index] <- exp(-spc.fits$logAo[i])
  SpeciesEst$Eoind[spc.index] <- spc.fits$Eo[i]
}
## Make Plot ####
Vplot <- ggplot(SpeciesEst, aes(x = Vind, y = V)) +
  geom_point(size = 2.5) + 
  geom_abline(intercept = 0, slope = 1, linewidth = 1.25) +
  xlab(("V (kPa) estimated independently")) +
  ylab("V (kPa) estimated hierarchically")

eoplot <- ggplot(SpeciesEst, aes(x = Eoind, y = Eo)) +
  geom_point(size = 2.5) + 
  geom_abline(intercept = 0, slope = 1, linewidth = 1.25) +
  xlab(expression("E"[o]~ "estimated independently")) +
  ylab(expression("E"[o]~ "estimated hierarchically"))
grid.arrange(Vplot, eoplot, ncol = 2)

# save species - level estimates 
write.csv(x = SpeciesEst, "analysis/species_estimates.csv", row.names = F)

# plot study effects
beta_p <- data.frame(re[grep(rownames(re), pattern = "beta_p"),])

ggplot(beta_p, aes(x = Estimate, weight = Std..Error)) +
  geom_histogram(bins = 20)

# make a dataframe of study names and estimated effect sizes
studydf <- data.frame(Source = levels(as.factor(all.dat$Source)),
                      effect_size = beta_p$Estimate,
                      effect_SE = beta_p$Std..Error)

nspecies <- all.dat %>%
  group_by(Source) %>%
  summarize(n = length(unique(AphiaID)))

studydf <- left_join(studydf, nspecies)
View(studydf)

# Plot Species fits with source effects ####
# lookup fits in hierarchical model
species.2.use <- c("Menidia menidia", "Oplophorus gracilirostris")
study.2.use <- c("Hoff 1967", "Cowles et al 1991 MarBiol 110:75:83")

### Make blank plot ####
par(mar = c(5,5,1,1))
plot(0, 0,
     type = "n",
     ylab = expression(log(pO[2]~W^n)),
     xlab = expression(kb^-1~(T^-1 - T[ref]^-1)),
     las = 1,
     ylim = c(0.0, 2),
     xlim = c(-2, 2),
     xaxs = "i",
     yaxs = "i",
     cex.axis = 1.5,
     cex.lab = 1.5
)


cols <- c("darkred", "darkblue")
for (i in 1:length(species.2.use)) {
  
  ## Lookup species estimates ####
  spc.index <- which(SpeciesEst$Species==species.2.use[i])
  study.index <- which(studydf$Source==study.2.use[i])
  
  Vind <- SpeciesEst$Vind[spc.index]
  Eoind <- SpeciesEst$Eoind[spc.index]
  
  V <- SpeciesEst$V[spc.index]
  Eo <- SpeciesEst$Eo[spc.index]
  n <- SpeciesEst$n[spc.index]
  
  beta_p_source <- studydf$effect_size[study.index]
  
  ## Lookup data for that species ####
  species_data <- filter(all.dat, Species == species.2.use[i])
  min.inv.t <- min(min(species_data$inv.temp),0)
  max.inv.t <- max(max(species_data$inv.temp), 0)
  Wbar <- mean(species_data$W)
  
  inv.t.line <- c(min.inv.t, max.inv.t)
  logpcrit_hatind <- log(Vind) - Eoind * inv.t.line 
  
  points(species_data$inv.temp,  log(species_data$Pcrit * Wbar^n),
         type = "p",
         pch = 21,
         bg = cols[i],
  )
  lines(inv.t.line, logpcrit_hatind,
        lwd = 2,
        col = cols[i],
        lty = "dashed")
  
  #### Make fitted line with source effect removed ####
  logpcrit_hat <- log(V) - Eo * inv.t.line 
  lines(inv.t.line, logpcrit_hat,
        lwd = 2,
        col = cols[i])
  arrows(x0 = 0,
         x1 = 0,
         y0 = log(Vind),
         y1 = log(Vind) + beta_p_source,
         lwd = 2)
}
