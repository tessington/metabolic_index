## Make table of all taxanomic group median and SD
library(dplyr)
library(tidyr)

source("code/fit_model_funs.R")

all.dat <- load_data()


class_summary <- all.dat %>%
  group_by(Class) %>%
  summarise(
    NSpecies = length(unique(Species))
  )

sim_beta <- readRDS("analysis/taxa_sims.RDS")
class_list <- list()
taxa_summary <- list()
Count <- 1
for (i in 1:nrow(class_summary)) {
  # get class median and SD
  taxa.name <- class_summary$Class[i]
  
  sims <- dplyr::filter(sim_beta, Group == taxa.name)
  class.medians <-
    apply(
      X = cbind(sims$logAo, sims$n, sims$Eo),
      MAR = 2,
      FUN = median
    )
  class.sds <-
    apply(
      X = cbind(sims$logAo, sims$n, sims$Eo),
      MAR = 2,
      FUN = sd
    )
  tmp_df <- tibble(
    Group = taxa.name,
    logAo = class.medians[1],
    logAosd = class.sds[1],
    n = class.medians[2],
    nsd = class.sds[2],
    Eo = class.medians[3],
    Eosd = class.sds[3],
    Nspecies = dplyr::filter(class_summary, Class == taxa.name)$NSpecies,
    level = "Class"
  )
  
  taxa_summary[[Count]] <- tmp_df
  Count = Count + 1
  order_summary <- all.dat %>%
    filter(Class == class_summary$Class[i]) %>%
    group_by(Order) %>%
    summarise(NSpecies = length(unique(Species)))
  n.order <- nrow(order_summary)
  for (j in 1:n.order) {
    order.name <- order_summary$Order[j]
    sims <- dplyr::filter(sim_beta, Group == order.name)
    order.medians <-
      apply(
        X = cbind(sims$logAo, sims$n, sims$Eo),
        MAR = 2,
        FUN = median
      )
    order.sds <-
      apply(
        X = cbind(sims$logAo, sims$n, sims$Eo),
        MAR = 2,
        FUN = sd
      )
    tmp_df <- tibble(
      Group = order.name,
      logAo = order.medians[1],
      logAosd = order.sds[1],
      n = order.medians[2],
      nsd = order.sds[2],
      Eo = order.medians[3],
      Eosd = order.sds[3],
      Nspecies = dplyr::filter(order_summary, Order == order.name)$NSpecies,
      level = "Order"
    )
    
    family_summary <- all.dat %>%
      filter(Order == order_summary$Order[j]) %>%
      group_by(Family) %>%
      summarise(NSpecies = length(unique(Species)))
    n.family <- nrow(family_summary)
    taxa_summary[[Count]] <- tmp_df
    Count = Count + 1
    
    for (k in 1:n.family) {
      family.name <- family_summary$Family[k]
      sims <- dplyr::filter(sim_beta, Group == family.name)
      family.medians <-
        apply(
          X = cbind(sims$logAo, sims$n, sims$Eo),
          MAR = 2,
          FUN = median
        )
      family.sds <-
        apply(
          X = cbind(sims$logAo, sims$n, sims$Eo),
          MAR = 2,
          FUN = sd
        )
      tmp_df <- tibble(
        Group = family.name,
        logAo = family.medians[1],
        logAosd = family.sds[1],
        n = family.medians[2],
        nsd = family.sds[2],
        Eo = family.medians[3],
        Eosd = family.sds[3],
        Nspecies = dplyr::filter(family_summary, Family == family.name)$NSpecies,
        level = "Family"
      )
      taxa_summary[[Count]] <- tmp_df
      Count = Count + 1
    }
  }
}

taxa_table <- do.call("rbind", taxa_summary)  

# save as csv file
write.csv(x = taxa_table, file = "analysis/taxa_table.csv")
