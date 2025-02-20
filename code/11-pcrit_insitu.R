library(ggplot2)
library(dplyr)
library(tidyr)
library(gsw)
library(wesanderson)

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

#Remove oxygen outliers
dat <- filter(dat, o2<1500)

#Just trawl survey data and IPHC long line data
dat<- dat %>%
  filter(survey %in% c("nwfsc", "dfo", "goa", "iphc"))

#Calculate o2 in kPa (data is in umol kg)
gas_const = 8.31
partial_molar_vol = 0.000032
kelvin = 273.15
boltz = 0.000086173324

SA = gsw_SA_from_SP(dat$salinity_psu,dat$depth,dat$longitude,dat$latitude) #absolute salinity for pot T calc
pt = gsw_pt_from_t(SA,dat$temp,dat$depth) #potential temp at a particular depth
O2_Sat0 = gsw_O2sol_SP_pt(dat$salinity_psu,pt)

press = exp(dat$depth*10000*partial_molar_vol/gas_const/(dat$temp+kelvin))
O2_satdepth = O2_Sat0*press
  
#solubility at p=0
sol0 = O2_Sat0/0.209
sol_Dep = sol0*press
dat$po2 = dat$o2/sol_Dep
dat$po2 <- dat$po2 * 101.325 # convert to kPa
dat$temp_K <- dat$temp+kelvin

#Region labels
labs <- c("British Columbia", "California Current", "Gulf of Alaska")
names(labs) <- c("bc", "cc","goa")
dat$region <- factor(dat$region, levels=c("cc", "bc", "goa")) 

#Set themes
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

#Plot all regions together
ggplot(data = dat, aes(x = temp, y=po2, colour=region)) +
  geom_point(alpha=0.3)+
  xlab("Temperature (C)") +
  ylab(bquote(pO[2]~"(kPa)")) +
  theme(legend.position = "top") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 15) ) + #or truncate to 15
  scale_y_continuous(expand = c(0,0), limits = c(0,30) )+
  scale_color_manual(values = c("#9970AB", "#ACD39E", "#FFEE99"), labels=labs, name="Region")+
   guides(colour = guide_legend(override.aes = list(size=6, alpha=1)))

#ggsave(filename = "temp_o2_single_plot.png",
#       plot = last_plot(),
#       width = 7,
#       height = 6,
#      units = "in")


# Pick a taxonmic group to simulate pcrit based on each temperature

family = "Myctophidae"


#Separate
ggplot(data = dat, aes(x = temp, y=po2)) +
  facet_wrap("region", ncol=3, labeller = labeller(region=labs))+
  geom_point(alpha=1)+
  xlab("Temperature (C)") +
  ylab(bquote(pO[2]~"(kPa)")) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 15) ) +
  scale_y_continuous(expand = c(0,0), limits = c(0,30) )+
  theme(panel.spacing = unit(2, "lines"))

ggsave(filename = "temp_o2_multi_plot.png",
       plot = last_plot(),
       width = 7,
       height = 3,
       units = "in")



