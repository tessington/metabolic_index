library(ggplot2)
library(gridExtra)
library(gtable)
rm(list = ls())
source("code/helper/fit_model_funs.R")

# turn off warnings
options( warn = -1)

# calcualte Pcrit for each taxa and temperature
pcrit_calculation <- function(sims, all.dat) {
  
  Wmed <- median(all.dat$W/ wref)
  tref <- 15
  
  # calculate Pcrit for 10degrees
  t_est1 <- 10
  inv.temp <- (1 / kb) * (1 / kelvin(t_est1) - 1 / kelvin(tref) )
  pcrit_1 <- sims$logV - sims$n * log(Wmed) - sims$Eo * inv.temp
  
  # calculate Pcrit for 20 degrees
  t_est2 <- 20
  inv.temp <- (1 / kb) * (1 / kelvin(t_est2) - 1 / kelvin(tref) )
  pcrit_2 <- sims$logV - sims$n * log(Wmed) - sims$Eo * inv.temp
  
  # put in a tibble
  pcrit_df <- tibble(pcrit = c(pcrit_1, pcrit_2),
                     temp = c(rep(t_est1, times = length(pcrit_1)),
                              rep(t_est2, times = length(pcrit_2))
                     )
  )
}
taxa.names <- c("Cottidae", "Enophrys bison")

ylims <- list( ylim_V = c(0, 2.5), ylim_n = c(0, 10), ylim_Eo = c(0, 4)) # y limits for metabolic index traits
ylim <- c(0, .75) # y limits for pcrit

# load data and filter
all.dat <- load_data()
all.dat <- filter_dat(all.dat)
# set up default parameters

wref <- 5
tref <- 15
kb <-  8.617333262145E-5
# load simulated  traits
sim_beta <- readRDS("analysis/taxa_sims.RDS")
# look up taxa names in the object "sim_beta"
ltaxa.names <- tolower(taxa.names)
lookup.taxa <- ltaxa.names %in% tolower(sim_beta$Group) # returns T if provided taxa names are in database
taxa_sims <- setNames(vector("list", length(taxa.names)), taxa.names) # create a list with elements named after "taxa.names"

# check to see if all taxa in taxa.names are in the MCMC output
if(!any(lookup.taxa)) cat("at least one of the taxa names are not in dataset")
if(all(lookup.taxa)) {

  # lookup simulated traits for each taxa and assign to list
  for (i in 1:length(taxa.names)) taxa_sims[[i]] <- dplyr::filter(sim_beta, Group == taxa.names[i])

# convert to df
taxa_sims_df <- bind_rows( taxa_sims)
# plot colors
colors_2_plot <- c("#67a9cf", "#ef8a62")

# begin plotting, first the metabolic index traits, and then pcrit
logV_plot <- ggplot(data = taxa_sims_df, aes(x = logV, fill = as.factor(Group))) +
  geom_density(alpha = 0.75, adjust = 1.5) + 
  ylab("Density") + 
  xlab("log(V)") + 
  scale_fill_manual(values = colors_2_plot) +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 2.5) ) + 
  scale_y_continuous(expand = c(0,0), limits = ylims[[1]])

n_plot <- ggplot(data = taxa_sims_df, aes(x = n, fill = as.factor(Group))) +
  geom_density(alpha = 0.75, adjust = 1.5) + 
  ylab("Density") + 
  xlab("n") + 
  scale_fill_manual(values = colors_2_plot) +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(-0.3, 0.1) ) + 
  scale_y_continuous(expand = c(0,0), limits = ylims[[2]])

Eo_plot <- ggplot(data = taxa_sims_df, aes(x = Eo, fill = as.factor(Group))) +
  geom_density(alpha = 0.75, adjust = 1.5) + 
  ylab("Density") + 
  xlab(bquote(E[o])) + 
  scale_fill_manual(values = colors_2_plot) +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(-0.2, 0.8) ) + 
  scale_y_continuous(expand = c(0,0), limits = ylims[[3]])

# make a dataframe of pcrit calculated for each simulated set of traits
pcrit_list <- lapply(X = taxa_sims, FUN = pcrit_calculation, all.dat = all.dat)
pcrit_df <- bind_rows(pcrit_list, .id = "Groups")
# plot results for each temperature
temperature_list <- c(10, 20)

pcrit_plot_cool <- ggplot(data = dplyr::filter(pcrit_df, temp == temperature_list[1]), aes(x = exp(pcrit), fill = as.factor(Groups))) +
  geom_density(alpha = 0.75, adjust = 1.5) + 
  ylab("Density") + 
  xlab(bquote(p[crit])) + 
  scale_fill_manual(values = c("#67a9cf", "#ef8a62")) +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 12) ) + 
  scale_y_continuous(expand = c(0,0), limits = ylim) +
  annotate("text", 
           x = 0.75 * 12 , y = 0.9 * ylim[2], 
           label = "10 \u00B0C", 
           hjust = 0, vjust = 1,
           size = 5)
pcrit_plot_warm <- ggplot(data = dplyr::filter(pcrit_df, temp == temperature_list[2]), aes(x = exp(pcrit), fill = as.factor(Groups))) +
  geom_density(alpha = 0.75, adjust = 1.5) + 
  ylab("Density") + 
  xlab(bquote(p[crit])) + 
  scale_fill_manual(values = c("#67a9cf", "#ef8a62")) +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 12) ) + 
  scale_y_continuous(expand = c(0,0), limits = ylim) + 
  annotate("text", 
           x = 0.75 * 12 , y = 0.9 * ylim[2], 
           label = "20 \u00B0C", 
           hjust = 0, vjust = 1,
           size = 5)


# get legend object from first plot and make a final grob object that is just a legend - code via chatGPT
tmp <- ggplotGrob(
  logV_plot +
    theme(
      legend.position = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = 18)
    )
)

legend_gtable <- gtable_filter(tmp, "guide-box")

grobs_list <- list(
  logV_plot,
  n_plot,
  Eo_plot,
  pcrit_plot_cool,
  pcrit_plot_warm,
  legend_gtable  # now this is compatible
)

# put all of the plots together
gridExtra::grid.arrange(grobs = grobs_list, ncol = 3, nrow = 2)

# get credibility intervals for pcrit and posterior medians
# turn pcrit_list into a nested list (temperature nested under Groups)
pcrit_nested_list <- lapply(pcrit_list, function(df) {
  split(df, df$temp)
})

intervals <- lapply(X = pcrit_nested_list,FUN = function(temp_list, p = 0.9) {
  lapply(temp_list, function(df, p = 0.9) {
    # get p-level credibility interval
    get_x_range(df$pcrit, p )  
  })
})

medians <- lapply(X = pcrit_nested_list,FUN = function(temp_list) {
  lapply(temp_list, function(df) {
    # get medians
    median(df$pcrit,)  
  })
})

print(intervals)
print(medians)

}