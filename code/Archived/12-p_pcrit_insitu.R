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

taxa.2.use <-  "Gadidae"
W <- 500

# specify values of temperature and po2 to use for pcrit calculations
t.range <- seq(min(dat$temp), 15, length.out = 20)
po2.range <- seq(min(dat$po2), 15, length.out = 20)

temp_po2_grid <- expand.grid(t.range, po2.range)
names(temp_po2_grid) <- c("temp", "po2")

# calculate the probability that each pair of temp, po2 implies that po2 > pcrit
p_po2_exceeds_pcrit <- mapply(FUN = calc_p_po2, temp_po2_grid$temp, temp_po2_grid$po2, taxa.name = taxa.2.use, w = W)
# put in a data frame
smoothing_df <- data.frame(temp = temp_po2_grid$temp,
                           po2 = temp_po2_grid$po2,
                           p_po2 = p_po2_exceeds_pcrit)
# 2-D smoothing via a loess smoother
model <- loess(p_po2 ~ temp * po2, data = smoothing_df)
# create a prediction data frame based on observations of temp and po2
predict_df <- data.frame(temp = dat$temp,
                         po2 = dat$po2)
# use predict to get smoothed value for probability that po2 > pcrit
predict_p_pcrit <- predict(model, newdata = predict_df)
predict_p_pcrit <- predict(model, newdata = smoothing_df)
smoothing_df$predict_p_crit <- predict_p_pcrit
breaks.2.plot <- c(0.05, 0.25, 0.5, 0.75, 0.95)

ggplot(smoothing_df, aes(x= po2, y = temp, z = predict_p_crit)) + 
  geom_contour( col = "black", breaks = breaks.2.plot, linewidth = 1 ) +
  metR::geom_text_contour(aes(z = predict_p_crit), breaks = breaks.2.plot, stroke = 0.5, skip = 0) +
  xlim(0, 15) + 
  ylim(2, 15)
# add this to the dataframe 'dat'


dat$p_pcrit <- predict_p_pcrit
# Set all rows with po2 >-15 to 1
dat$p_pcrit[dat$po2>=15] = 1


# Plot results
ggplot(data = dat, aes(x = po2, y=-depth, col = p_pcrit)) +
  facet_wrap("region", ncol=3, labeller = labeller(region=labs))+
  geom_point(alpha = 0.75)+
  labs(col = bquote(P(pO[2]>p["crit"]))) +
  ylab("Depth(m) ") +
  xlim(c(0, 25)) + 
  xlab(bquote(pO[2]~"(kPa)")) +
  theme(panel.spacing = unit(2, "lines")) +
  scale_colour_viridis_c(option = "plasma", direction = 1, limits = c(0,1),
                         oob = scales::squish)  

savefilename <- paste0("figures/", 
                       taxa.2.use, 
                      "prob_pcrit.png")
ggsave(filename = savefilename,
       plot = last_plot(),
       width = 8,
       height = 3.5,
       units = "in")

