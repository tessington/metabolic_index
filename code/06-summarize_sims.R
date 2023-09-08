## Make table of all taxanomic group median and SD
library(dplyr)
library(tidyr)

source("code/fit_model_funs.R")

all.dat <- load_data()


order_summary <- all.dat %>%
  group_by(Order) %>%
  summarise(
    NSpecies = length(unique(Species))
  )

sim_beta <- readRDS("analysis/taxa_sims.RDS")
taxa_summary <- list()
Count <- 1
for (i in 1:nrow(order_summary)) {
  # get class median and SD
  taxa.name <- order_summary$Order[i]
  
  sims <- dplyr::filter(sim_beta, Group == taxa.name)
  order.medians <-
    apply(
      X = cbind(sims$logV, sims$n, sims$Eo),
      MAR = 2,
      FUN = median
    )
  order.sds <-
    apply(
      X = cbind(sims$logV, sims$n, sims$Eo),
      MAR = 2,
      FUN = sd
    )
  order.low <- apply(
    X = cbind(sims$logV, sims$n, sims$Eo),
    MAR = 2,
    FUN = quantile,
    probs = 0.05
  )
  order.high <- apply(
    X = cbind(sims$logV, sims$n, sims$Eo),
    MAR = 2,
    FUN = quantile,
    probs = 0.95
  )
  tmp_df <- tibble(
    Group = taxa.name,
    logV = order.medians[1],
    logVsd = order.sds[1],
    logVlow = order.low[1],
    logVhigh = order.high[1],
    n = order.medians[2],
    nsd = order.sds[2],
    nlow = order.low[2],
    nhigh = order.high[2],
    Eo = order.medians[3],
    Eosd = order.sds[3],
    Eolow = order.low[3],
    Eohigh = order.high[3],
    Nspecies = dplyr::filter(order_summary, Order == taxa.name)$NSpecies,
    level = "Order"
  )
  
  taxa_summary[[Count]] <- tmp_df
  Count = Count + 1
  
  # setup loop of families
  family_summary <- all.dat %>%
    filter(Order == order_summary$Order[i]) %>%
    group_by(Family) %>%
    summarise(NSpecies = length(unique(Species)))
  n.family <- nrow(family_summary)
  for (j in 1:n.family) {
    family.name <- family_summary$Family[j]
    sims <- dplyr::filter(sim_beta, Group == family.name)
    family.medians <-
      apply(
        X = cbind(sims$logV, sims$n, sims$Eo),
        MAR = 2,
        FUN = median
      )
    family.sds <-
      apply(
        X = cbind(sims$logV, sims$n, sims$Eo),
        MAR = 2,
        FUN = sd
      )
    family.low <- apply(
      X = cbind(sims$logV, sims$n, sims$Eo),
      MAR = 2,
      FUN = quantile,
      probs = 0.05
    )
    family.high <- apply(
      X = cbind(sims$logV, sims$n, sims$Eo),
      MAR = 2,
      FUN = quantile,
      probs = 0.95
    )
    
    tmp_df <- tibble(
      Group = family.name,
      logV = family.medians[1],
      logVsd = family.sds[1],
      logVlow = family.low[1],
      logVhigh = family.high[1],
      n = family.medians[2],
      nsd = family.sds[2],
      nlow = family.low[2],
      nhigh = family.high[2],
      Eo =  family.medians[3],
      Eosd = family.sds[3],
      Eolow = family.low[3],
      Eohigh = family.high[3],
      Nspecies = dplyr::filter(family_summary, Family == family.name)$NSpecies,
      level = "Family"
    )
    taxa_summary[[Count]] <- tmp_df
    Count = Count + 1
}
}

taxa_table <- do.call("rbind", taxa_summary)  

# save as csv file
write.csv(x = taxa_table, file = "analysis/taxa_table.csv")



