library(tidyr)
library(dplyr)
library(mvtnorm)
library(respR)
library(sp)
library(gsw)
library(sdmTMB)
library(ggplot2)
library(gridExtra)

source("code/fit_model_funs.R")

survey.name = "SYN WCVI"
env_data <- readRDS("data/bc-synoptic-env.rds")
trawl_data <- readRDS("data/bc-synoptic-trawls.rds")
dat_all <- left_join(trawl_data, env_data, by = "fishing_event_id")
species.list <- unique(dat_all$species)

dat_all <- rename(dat_all, o2 = do, temp = temperature, cpue_kg_m2 = density_kgpm2, sal = salinity, 
                  depth = depth_m, longitude_dd = longitude, latitude_dd = latitude)

species.name <- "english sole"
dat = filter(dat_all, species == species.name, survey == survey.name,
             !is.na(temp), !is.na(o2), !is.na(sal), !is.na(depth),
             year != c(2016))

dat$o2_pressure <- convert_DO(dat$o2,
                              from = "ml/l",
                              to = "kPa", 
                              t = dat$temp,
                              S = dat$sal)
kelvin = 273.15
kb = 0.000086173324
Tref <- 15 + kelvin

# get traits

mi.traits <- lookup_taxa(taxa.name = "Pleuronectidae")


Ao <-exp( -mi.traits$logV)
n <- mi.traits$n
Eo <- mi.traits$Eo
nsims <- length(Ao)

dat$pphi1 <- NA
dat$pphi2 <- NA
dat$pphi3 <- NA
# loop through data
wref <- 5
W <- 250/wref  # assuming 500 g individual
x1 <- 1
x2 <- 2
x3 <- 3
for (i in 1:nrow(dat)) {
  Temp <- dat$temp[i] + kelvin
  po2 <- dat$o2_pressure[i]
  mi.tmp <- rep(NA, nsims)
  for (j in 1:nsims) {
    mi.tmp[j] = po2 * Ao[j] * W^n[j] * exp(Eo[j] / kb * (1 / Temp - 1 / Tref))
  }
  dat$pphi1[i] <- length(which(mi.tmp >= x1)) / nsims
  dat$pphi2[i] <- length(which(mi.tmp >= x2)) / nsims
  dat$pphi3[i] <- length(which(mi.tmp >= x3)) / nsims
}


dat_ll = dat
coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
# convert to utm with spTransform
dat_utm = spTransform(dat_ll, 
                      CRS("+proj=utm +zone=9 +datum=WGS84 +units=km")) # check on zone
# convert back from sp object to data frame
dat = as.data.frame(dat_utm)
dat = dplyr::rename(dat, longitude = coords.x1, 
                   latitude = coords.x2)


spde <- make_mesh(dat, xy_cols = c("longitude", "latitude"), n_knots = 250) # choose # knots
pa <- function(x) {
if (x>0) y = 1
if (x ==0) y = 0
return(y)
}
dat$cpue_kg_km2 <- 1E6 * dat$cpue_kg_m2
dat$present <- sapply(X = dat$cpue_kg_m2, FUN = pa)

## abundance SDMS ####
mod <- sdmTMB(cpue_kg_km2 ~ scale(pphi1) +  (scale((depth))) + as.factor(year),
              data = dat,
              mesh = spde,
              family = tweedie(link = "logit"),
              anisotropy = TRUE,
              spatial = "on",
              spatiotemporal = "off"
)

mod2 <- sdmTMB(cpue_kg_km2 ~ scale(o2_pressure) + s(scale((depth))) + as.factor(year),
              data = dat,
              mesh = spde,
              family = tweedie(link = "logit"),
              anisotropy = TRUE
)
mod3 <- sdmTMB(cpue_kg_km2 ~ scale(temp) + s(scale((depth))) + as.factor(year),
               data = dat,
               mesh = spde,
               family = tweedie(link = "logit"),
               anisotropy = TRUE
)
summary(mod)
summary(mod2)
summary(mod3)

## Presence - Absence SDMS ####
mod0 <- sdmTMB(present ~  s(scale((depth))) + as.factor(year),
              data = dat,
              mesh = spde,
              family = binomial(link = "logit"),
              anisotropy = TRUE,
              spatial = "on",
              spatiotemporal = "off"
)

mod <- sdmTMB(present ~ scale(pphi1) +  s(scale((depth))) + as.factor(year),
              data = dat,
              mesh = spde,
              family = binomial(link = "logit"),
              anisotropy = TRUE,
              spatial = "on",
              spatiotemporal = "off"
)

mod2 <- sdmTMB(present ~ scale(o2_pressure) + s(scale((depth))) + as.factor(year),
               data = dat,
               mesh = spde,
               family = binomial(link = "logit"),
               anisotropy = TRUE,
               spatial = "on",
               spatiotemporal = "off"
)
mod3 <- sdmTMB(present ~ scale(temp) + s(scale((depth))) + as.factor(year),
               data = dat,
               mesh = spde,
               family = binomial(link = "logit"),
               anisotropy = TRUE,
               spatial = "on",
               spatiotemporal = "off"
)
summary(mod0)
summary(mod)
summary(mod2)
summary(mod3)

phi1 <- ggplot(dat, aes(x= pphi1, y = depth, col = factor(present))) + 
  geom_point(show.legend = T) +
  scale_y_reverse(limits = c(400, 0)) + 
  scale_color_viridis_d(alpha = 1, begin = 0, end = 0.4 ) +
  xlab(expression(P(phi)>1))

phi2 <- ggplot(dat, aes(x= pphi2, y = depth, col = factor(present))) + 
  geom_point(show.legend = F) +
  scale_y_reverse(limits = c(400, 0)) + 
  scale_color_viridis_d(alpha = 0.6, begin = 0, end = .4 ) +
  xlab(expression(P(phi)>2))

phi3 <- ggplot(dat, aes(x= pphi3, y = depth, col = factor(present))) + 
  geom_point(show.legend = F) +
  scale_y_reverse(limits = c(400, 0)) + 
  scale_color_viridis_d(alpha = 0.6, begin = 0, end = 0.4 ) +
  xlab(expression(P(phi)>3)) 


grid.arrange(phi1, phi2, phi3, ncol = 3)
print(phi1)
