library(ggplot2)
library(dplyr)
library(tidyr)


source("code/fit_model_funs.R")
# Load O2 Data
dat <- as.data.frame(readRDS("data/all_o2_dat_filtered.rds"))

#Remove any rows with missing data, outliers
dat <- dat %>%
drop_na(o2, temp, sigma0, doy, X, Y, year)
dat <- filter(dat, o2<1500)

#Remove AI, EBS
dat<- dat %>%
  filter(region %in% c("cc", "bc", "goa"))

#Remove weird depths
dat <- filter(dat, depth>0)

# remove super deep depths
dat <- filter(dat, depth <1000)

#Remove oxygen outliers
dat <- filter(dat, o2<1500)

#Just trawl survey data and IPHC long line data
dat<- dat %>%
  filter(survey %in% c("nwfsc", "dfo", "goa", "iphc", "wcoa", "codap", "calCOFI", "NewportLine"))

#Calculate o2 in kPa (data is in umol kg)
gas_const = 8.31
partial_molar_vol = 0.000032
boltz = 0.000086173324

SA = gsw_SA_from_SP(dat$salinity_psu,dat$depth,dat$longitude,dat$latitude) #absolute salinity for pot T calc
pt = gsw_pt_from_t(SA,dat$temp,dat$depth) #potential temp at a particular depth
O2_Sat0 = gsw_O2sol_SP_pt(dat$salinity_psu,pt)

press = exp(dat$depth*10000*partial_molar_vol/gas_const/kelvin(dat$temp))
O2_satdepth = O2_Sat0*press
  
#solubility at p=0
sol0 = O2_Sat0/0.209
sol_Dep = sol0*press
dat$po2 = dat$o2/sol_Dep
dat$po2 <- dat$po2 * 101.325 # convert to kPa
dat$temp_K <- kelvin(dat$temp)

#Region labels
labs <- c("British Columbia", "California Current", "Gulf of Alaska")
names(labs) <- c("bc", "cc","goa")
dat$region <- factor(dat$region, levels=c("cc", "bc", "goa")) 

#Set themes
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

# Pick a taxonomic group to simulate pcrit based on each temperature

taxa.2.use <-  "Myctophidae"
W <- 500

t.range <- seq(min(dat$temp), max(dat$temp), length.out = 100)

pcrit_df <- lookup_taxa_t(taxa.2.use, t.range = t.range, w.2.use = W)

# make an x_y positions df based on inner 50% range
positions50 <- data.frame(
  Temp = c(pcrit_df$Temp, rev(pcrit_df$Temp)),
  ys = c(pcrit_df$lower50s, rev(pcrit_df$upper50s))
)
positions90 <- data.frame(
  temp = c(pcrit_df$Temp, rev(pcrit_df$Temp)),
  ys = c(pcrit_df$lower90s, rev(pcrit_df$upper90s))
)

# depth legend labels
depthlabel = c(75, 200, 500, 1000)

#Separate
ggplot(data = dat, aes(x = temp, y=po2, col = log(depth))) +
  facet_wrap("region", ncol=3, labeller = labeller(region=labs))+
  geom_point(alpha = 0.75)+
  labs(col = "Depth(m)") +
  xlab("Temperature (C)") +
  ylab(bquote(pO[2]~"(kPa)")) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 21) ) +
  scale_y_continuous(expand = c(0,0), limits = c(0,25.5) ) +
  theme(panel.spacing = unit(2, "lines")) +
  scale_colour_viridis_c(option = "mako", direction = -1, limits = log(c(50, 1500)), end = 0.9,
                         oob = scales::squish, breaks = log(depthlabel),
                         labels = depthlabel)  + 
  geom_polygon(data = positions50, aes(x = Temp, y = ys),fill = "grey1", alpha = 0.5, color = NA) +
  geom_polygon(data = positions90, aes(x = temp, y = ys), fill = "lightgrey", color = NA, alpha = 0.5) 

  
savefilename <- paste0("figures/",
                       gsub(x = taxa.2.use, pattern = " ", replacement = "_"),
                       "_temp_o2_multi_plot.png"
                       )
ggsave(filename = savefilename,
       plot = last_plot(),
       width = 8,
       height = 3.5,
       units = "in")



