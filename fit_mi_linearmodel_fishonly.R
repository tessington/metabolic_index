# Code to test TMB fitting
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)
library(glmmTMB)


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

all.dat <- dplyr::filter(all.dat, Class %in% c("Actinopteri", "Elasmobranchii"))


kb <-  8.617333262145E-5


all.dat$spc <- as.factor(all.dat$Species)
nspc <- length(levels(all.dat$spc))


tref <- 12
all.dat$inv.temp <- (1 / kb)  * ( 1 / (all.dat$Temp + 273.15) - 1 / (tref + 273.15))

### Fit Models ####
fit <- glmmTMB(-log(Pcrit) ~  diag(1 + inv.temp|spc) + inv.temp + log(W), data = all.dat,
               family = gaussian(link = "identity"),
               se = TRUE
)
fit2 <- glmmTMB(-log(Pcrit) ~ (1|spc) + log(W) + inv.temp, data = all.dat,
                family = gaussian(link = "identity"),
                se = TRUE
)

summary(fit)
summary(fit2)
AICmat <- c(AIC(fit), AIC(fit2))
print(AICmat)
## model with random effects of temp by species fits much better
# Extract parameter estimates
spc <- levels(all.dat$spc)
n <- fixef(fit)$cond[2]
Eo <- ranef(fit)$cond$spc[,2]
Aospc.effects <- exp(fixef(fit2)$cond[1] + ranef(fit2)$cond$spc$`(Intercept)`)

all.dat$bnpo2 <- log(all.dat$W^n * all.dat$Pcrit)


