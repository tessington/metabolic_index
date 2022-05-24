# calc o2 solubility, relies on o2 in umol/kg
gsw_O2sol_SP_pt <- function(sal,pt) {
  x = sal
  pt68 = pt*1.00024
  y = log((298.15 - pt68)/(273.15 + pt68))
  
  a0 =  5.80871
  a1 =  3.20291
  a2 =  4.17887
  a3 =  5.10006
  a4 = -9.86643e-2
  a5 =  3.80369
  b0 = -7.01577e-3
  b1 = -7.70028e-3
  b2 = -1.13864e-2
  b3 = -9.51519e-3
  c0 = -2.75915e-7
  
  O2sol = exp(a0 + y*(a1 + y*(a2 + y*(a3 + y*(a4 + a5*y)))) + x*(b0 + y*(b1 + y*(b2 + b3*y)) + c0*x))
  return(O2sol)
}

load_all_hauls <- function() {
  haul = nwfscSurvey::PullHaul.fn(SurveyName = "NWFSC.Combo")
  haul <- plyr::rename(haul, replace=c("salinity_at_gear_psu_der" = "sal", 
                                       "temperature_at_gear_c_der" = "temp", 
                                       "o2_at_gear_ml_per_l_der" = "o2",
                                       "depth_hi_prec_m" = "depth"))
                               
  # read in the grid cell data from the survey design
  grid_cells = readxl::read_excel("data/data/Selection Set 2018 with Cell Corners.xlsx")
  grid_cells = dplyr::mutate(grid_cells,
                             depth_min = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[1]),
                             depth_max = as.numeric(unlist(strsplit(grid_cells$Depth.Range,"-"))[2]))
  
  # convert grid_cells to sp object
  grid = SpatialPoints(cbind(grid_cells$Cent.Long,grid_cells$Cent.Lat),
                       proj4string = CRS("+proj=longlat +datum=WGS84"))
  r = raster::rasterize(x=grid, y = raster(nrow=length(unique(grid_cells$Cent.Lat)),
                                           ncol=length(unique(grid_cells$Cent.Long))))
  rasterToPoints(r)
  
  raster = aggregate(r, fact = 2)
  raster = projectRaster(raster, crs = "+proj=tmerc +lat_0=31.96 +lon_0=-121.6 +k=1 +x_0=390000 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  # create matrix of point data with coordinates and depth from raster
  grid = as.data.frame(rasterToPoints(raster))
  
  # Figure out the grid cell corresponding to each tow location
  haul$Cent.Lat = NA
  haul$Cent.Lon = NA
  haul$Cent.ID = NA
  for(i in 1:nrow(haul)) {
    indx = which(grid_cells$NW.LAT > haul$latitude_dd[i] &
                   grid_cells$SW.LAT < haul$latitude_dd[i] &
                   grid_cells$NW.LON < haul$longitude_dd[i] &
                   grid_cells$NE.LON > haul$longitude_dd[i])
    if(length(indx) > 0) {
      haul$Cent.ID[i] = grid_cells$Cent.ID[indx]
      haul$Cent.Lat[i] = grid_cells$Cent.Lat[indx]
      haul$Cent.Lon[i] = grid_cells$Cent.Long[indx]
    }
  }
  
  # project lat/lon to UTM, after removing missing values and unsatisfactory hauls
  haul = haul %>% filter(!is.na(Cent.Lon), performance == "Satisfactory")
  
  haul_trans = haul
  coordinates(haul_trans) <- c("Cent.Lon", "Cent.Lat")
  proj4string(haul_trans) <- CRS("+proj=longlat +datum=WGS84")
  newproj = paste("+proj=utm +zone=10 ellps=WGS84")
  haul_trans <- spTransform(haul_trans, CRS(newproj))
  haul_trans = as.data.frame(haul_trans)
  haul_trans$Cent.Lon = haul_trans$Cent.Lon/10000
  haul_trans$Cent.Lat = haul_trans$Cent.Lat/10000
  haul_trans$year = as.numeric(substr(haul_trans$date_yyyymmdd,1,4))
  
  haul$X = haul_trans$Cent.Lon
  haul$Y = haul_trans$Cent.Lat
  haul$year = haul_trans$year
  #haul$year_centered = haul$year - mean(unique(haul$year))
  
  return(haul)
  
   
}

calc_po2_mi <- function(dat) {
  #O2 from trawl data is in ml/l - may need to be converted to umol/kg
  gas_const = 8.31
  partial_molar_vol = 0.000032
  kelvin = 273.15
  boltz = 0.000086173324
  
  #calculate percent saturation for O2 - assumes  units of mL O2/L
  # Input:       S = Salinity (pss-78)
  #              T = Temp (deg C) ! use potential temp
  #depth is in meters
  #[umole/kg] = [ml/L]*44660/(sigmatheta(P=0,theta,S) + 1000)
  dat$SA = gsw_SA_from_SP(dat$sal,dat$depth,dat$longitude_dd,dat$latitude_dd) #absolute salinity for pot T calc
  dat$pt = gsw_pt_from_t(dat$SA,dat$temp,dat$depth) #potential temp at a particular depth
  dat$CT = gsw_CT_from_t(dat$SA,dat$temp,dat$depth) #conservative temp
  dat$sigma0 = gsw_sigma0(dat$SA,dat$CT)
  dat$o2_umolkg = dat$o2*44660/(dat$sigma0+1000)
  
  
  dat$O2_Sat0 = gsw_O2sol_SP_pt(dat$sal,dat$pt)
  
  #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
  dat$press = exp(dat$depth*10000*partial_molar_vol/gas_const/(dat$temp+kelvin))
  dat$O2_satdepth = dat$O2_Sat0*dat$press
  
  #solubility at p=0
  dat$sol0 = dat$O2_Sat0/0.209
  dat$sol_Dep = dat$sol0*dat$press
  dat$po2 = dat$o2_umolkg/dat$sol_Dep
  
  # species-specific parameters
  Ao = 1.16625e-13
  Eo = 0.8736 # from cod, 0.8736.  Make it one half or double
  B = 1200 # size in grams, roughly average (initial calculations used 10000g)
  n = -0.208 # borrowed from cod 
  
  dat$mi = B^n*Ao*dat$po2/exp(-1*Eo/(boltz*(dat$temp+kelvin))) 
  return(dat)
}

load_data <- function(fit.model= F, spc, constrain_latitude = F) {
  
  dat <- readRDS("survey_data/joined_nwfsc_data.rds")
  
  # analyze sablefish for years and hauls with adequate oxygen and temperature data, within range of occurrence
  dat = dplyr::filter(dat, species == spc, year%in%seq(2010,2015))
  if (constrain_latitude) dat <- dplyr::filter(dat, latitude_dd >=43)

  
  
  # get julian day
  dat$julian_day <- rep(NA, nrow(dat))
  for (i in 1:nrow(dat)) dat$julian_day[i] <- as.POSIXlt(dat$date[i], format = "%Y-%b-%d")$yday

  
  # create temporary data file, matching J-SCOPE extent, for model fitting
  if(fit.model) {
    # constraint to J-SCOPE extent
    # UTM transformation
    dat_ll = dat
    coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
    proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
    # convert to utm with spTransform
    dat_utm = spTransform(dat_ll, 
                          CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
    # convert back from sp object to data frame
    dat$X <- as.data.frame(dat_utm)$longitude_dd
    dat$Y <- as.data.frame(dat_utm)$latitude_dd
    # scale depth and julian dat
    dat$log_depth_scaled <- scale(log(dat$depth))
    dat$log_depth_scaled2 <- dat$log_depth_scaled^2
    dat$jday_scaled <- scale(dat$julian_day)
    dat$jday_scaled2 <- dat$jday_scaled^2
    # create temporary data file for fitting
    fit_dat <- dplyr::filter(dat, !is.na(o2))
    
    
    c_spde <-make_mesh(data = fit_dat, xy_cols = c("X", "Y"), n_knots = 250) # choose # knots
    # fit dissolved oxygen 
    o2_model <-  sdmTMB(formula = o2 ~ -1 + log_depth_scaled + log_depth_scaled2 
                        +  as.factor(year) + jday_scaled + jday_scaled2,
                        data = fit_dat,
                        time = "year", spde = c_spde, anisotropy = TRUE,
                        silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                        control = sdmTMBcontrol(step.min = 0.01, step.max = 1))
    # get predictions
    pred_o2 <- predict(o2_model,
                                 newdata = dat,
                                 return_tmb_object = F)
    #impute missing values
    index <- which(is.na(dat$o2))
    dat$o2[index] <- pred_o2$est[index]
    
    # impute salinity, following same steps
    sal_model <- sdmTMB(formula = o2 ~ -1 + log_depth_scaled + log_depth_scaled2 
                         +  as.factor(year) + jday_scaled + jday_scaled2,
                         data = fit_dat,
                         time = "year", spde = c_spde, anisotropy = TRUE,
                         silent = TRUE, spatial_trend = FALSE, spatial_only = FALSE,
                         control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
    )
  pred_sal <- predict(sal_model,
                          newdat = dat,
                          return_tmb_object = F)
  index <- which(is.na(dat$sal))
  dat$sal[index] <- pred_sal$est[index]
  }
  
  # compute metabolic index (mi) --------------------------------------------
  # converted from Halle Berger matlab script
  
  #O2 from trawl data is in ml/l 
  # just in case, remove any missing or nonsense values from sensors
  dat <- dplyr::filter(dat, !is.na(o2), !is.na(sal), !is.na(temp), is.finite(sal))
  
  gas_const = 8.31
  partial_molar_vol = 0.000032
  kelvin = 273.15
  boltz = 0.000086173324
  
  #calculate percent saturation for O2 - assumes  units of mL O2/L
  # Input:       S = Salinity (pss-78)
  #              T = Temp (deg C) ! use potential temp
  #depth is in meters
  #[umole/kg] = [ml/L]*44660/(sigmatheta(P=0,theta,S) + 1000)
  dat$SA = gsw_SA_from_SP(dat$sal,dat$depth,dat$longitude_dd,dat$latitude_dd) #absolute salinity for pot T calc
  dat$pt = gsw_pt_from_t(dat$SA,dat$temp,dat$depth) #potential temp at a particular depth
  dat$CT = gsw_CT_from_t(dat$SA,dat$temp,dat$depth) #conservative temp
  dat$sigma0 = gsw_sigma0(dat$SA,dat$CT)
  dat$o2_umolkg = dat$o2*44660/(dat$sigma0+1000)
 
  
  dat$O2_Sat0 = gsw_O2sol_SP_pt(dat$sal,dat$pt)
  
  #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
  dat$press = exp(dat$depth*10000*partial_molar_vol/gas_const/(dat$temp+kelvin))
  dat$O2_satdepth = dat$O2_Sat0*dat$press
  
  #solubility at p=0
  dat$sol0 = dat$O2_Sat0/0.209
  dat$sol_Dep = dat$sol0*dat$press
  dat$po2 = dat$o2_umolkg/dat$sol_Dep
  
  # species-specific parameters, Ao doesn't really matter because it is a constant
  Ao = 1.16625e-13
  Eo = 0.8736 # from cod, 0.8736.  Make it one half or double
  # the below also doesn't matter because B^n is a constant and the index is scaled for fitting.  
  B = 1200 # size in grams, roughly average (initial calculations used 10000g)
  n = -0.208 # borrowed from cod 
  
  dat$mi = B^n*Ao*dat$po2/exp(-1*Eo/(boltz*(dat$temp+kelvin)))
  dat <- dplyr::filter(dat, !is.na(temp), !is.na(mi))
  
  # prepare data and models -------------------------------------------------
  
  dat <- dplyr::select(dat, trawl_id, species, year, longitude_dd, latitude_dd, cpue_kg_km2,
                o2, temp, depth, mi, po2, julian_day, pass)
  
  
  # UTM transformation
  dat_ll = dat
  coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = spTransform(dat_ll, 
                        CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, longitude = longitude_dd, 
                      latitude = latitude_dd)
  return(dat)
}

load_data_jscope <- function(spc, years) {
  
  #load trawl data
  hauldat <- readRDS("survey_data/joined_jscope_trawl_2003_2018.rds")
  
  hauldat <- dplyr::filter(hauldat, year %in% years)
  catchdat <- readRDS("survey_data/joined_nwfsc_data.rds")
  catchdat <- dplyr::select(catchdat, c(trawl_id, project, year, pass, vessel, tow, date,
                                        longitude_dd, latitude_dd, area_swept_ha, cpue_kg_km2,
                                        species,performance, survey, depth)
  )
  
  hauldat <- dplyr::select(hauldat, c(trawl_id, o2_model, temp_model, sal_model))
                           
 
  # analyze sablefish for years and hauls with adequate oxygen and temperature data, within range of occurrence
  catchdat = dplyr::filter(catchdat, species == spc, year%in%years)
  # merge with new jscope
  dat <- inner_join(catchdat, hauldat, by = "trawl_id")
  dat <- dplyr::rename(dat, o2 = o2_model, temp = temp_model, sal = sal_model)  
  
  # remove nonsense rows
  dat <- dat %>%
    filter(!is.na(o2), o2>0, !is.na(sal), !is.infinite(sal), !is.na(temp))
  
  # get julian day
  dat$julian_day <- rep(NA, nrow(dat))
  for (i in 1:nrow(dat)) dat$julian_day[i] <- as.POSIXlt(dat$date[i], format = "%Y-%b-%d")$yday
  

  
  # get Metabolic Indiex
  # converted from Halle Berger matlab script
  #O2 from trawl data is in ml/l - may need to be converted to umol/kg
  gas_const = 8.31
  partial_molar_vol = 0.000032
  kelvin = 273.15
  boltz = 0.000086173324
  
  #calculate percent saturation for O2 - assumes  units of mL O2/L
  # Input:       S = Salinity (pss-78)
  #              T = Temp (deg C) ! use potential temp
  #depth is in meters
  #[umole/kg] = [ml/L]*44660/(sigmatheta(P=0,theta,S) + 1000)
  dat$SA = gsw_SA_from_SP(dat$sal,dat$depth,dat$longitude_dd,dat$latitude_dd) #absolute salinity for pot T calc
  dat$pt = gsw_pt_from_t(dat$SA,dat$temp,dat$depth) #potential temp at a particular depth
  dat$CT = gsw_CT_from_t(dat$SA,dat$temp,dat$depth) #conservative temp
  dat$sigma0 = gsw_sigma0(dat$SA,dat$CT)
  dat$o2_umolkg = dat$o2*44660/(dat$sigma0+1000)
 
  
  dat$O2_Sat0 = gsw_O2sol_SP_pt(dat$sal,dat$pt)
  
  #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
  dat$press = exp(dat$depth*10000*partial_molar_vol/gas_const/(dat$temp+kelvin))
  dat$O2_satdepth = dat$O2_Sat0*dat$press
  
  #solubility at p=0
  dat$sol0 = dat$O2_Sat0/0.209
  dat$sol_Dep = dat$sol0*dat$press
  dat$po2 = dat$o2_umolkg/dat$sol_Dep
  
  # species-specific parameters
  Ao = 1.16625e-13
  Eo = 0.8736 # from cod, 0.8736.  Make it one half or double
  B = 1200 # size in grams, roughly average (initial calculations used 10000g)
  n = -0.208 # borrowed from cod 
  
  dat$mi = B^n*Ao*dat$po2/exp(-1*Eo/(boltz*(dat$temp+kelvin)))
  dat$longitude <- dat$longitude / 10
  dat$latitude <- dat$latitude / 10
  # prepare data and models -------------------------------------------------
  
  dat <- dplyr::select(dat, species, trawl_id, year, longitude_dd, latitude_dd, cpue_kg_km2,
                       o2, temp, depth, mi, po2, julian_day, pass)
  
  
  # UTM transformation
  dat_ll = dat
  coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
  proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = spTransform(dat_ll, 
                        CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat = as.data.frame(dat_utm)
  dat = dplyr::rename(dat, longitude = longitude_dd, 
                      latitude = latitude_dd)
  return(dat)
}
load_bc_data<- function(species.name = "sablefish", survey.name = "SYN WCVI") {
  env_data <- readRDS("survey_data/bc-synoptic-env.rds")
  trawl_data <- readRDS("survey_data/bc-synoptic-trawls.rds")
  dat_all <- left_join(trawl_data, env_data, by = "fishing_event_id")
  
  # renaming to match wc data
  dat_all <- rename(dat_all, o2 = do, temp = temperature, cpue_kg_km2 = density_kgpm2, sal = salinity, 
                    depth = depth_m, longitude_dd = longitude, latitude_dd = latitude)
  
  
    # analyze years and hauls with adequate oxygen and temperature data, within range of occurrence
    dat = filter(dat_all, species == species.name, survey == survey.name,
                 !is.na(temp), !is.na(o2), !is.na(sal), !is.na(depth),
                 latitude_dd > min(latitude_dd[which(cpue_kg_km2>0)]),
                 latitude_dd <= max(latitude_dd[which(cpue_kg_km2>0)]),
                 longitude_dd > min(longitude_dd[which(cpue_kg_km2>0)]),
                 longitude_dd < max(longitude_dd[which(cpue_kg_km2>0)]),
                 year != c(2016))
    gas_const = 8.31
    partial_molar_vol = 0.000032
    kelvin = 273.15
    boltz = 0.000086173324
    
    #calculate percent saturation for O2 - assumes  units of mL O2/L
    # Input:       S = Salinity (pss-78)
    #              T = Temp (deg C) ! use potential temp
    #depth is in meters
    #[umole/kg] = [ml/L]*44660/(sigmatheta(P=0,theta,S) + 1000)
    dat$SA = gsw_SA_from_SP(dat$sal,dat$depth,dat$longitude_dd,dat$latitude_dd) #absolute salinity for pot T calc
    dat$pt = gsw_pt_from_t(dat$SA,dat$temp,dat$depth) #potential temp at a particular depth
    dat$CT = gsw_CT_from_t(dat$SA,dat$temp,dat$depth) #conservative temp
    dat$sigma0 = gsw_sigma0(dat$SA,dat$CT)
    dat$o2_umolkg = dat$o2*44660/(dat$sigma0+1000)
    # calc o2 solubility, relies on o2 in umol/kg
  
    
    dat$O2_Sat0 = gsw_O2sol_SP_pt(dat$sal,dat$pt)
    
    #= o2satv2a(sal,pt) #uses practical salinity and potential temp - solubity at p =1 atm
    dat$press = exp(dat$depth*10000*partial_molar_vol/gas_const/(dat$temp+kelvin))
    dat$O2_satdepth = dat$O2_Sat0*dat$press
    
    #solubility at p=0
    dat$sol0 = dat$O2_Sat0/0.209
    dat$sol_Dep = dat$sol0*dat$press
    dat$po2 = dat$o2_umolkg/dat$sol_Dep
    
    # species-specific parameters
    if(species.name == "sablefish"){
      Ao = 1.16625e-13 # sablefish specific
      B = 1200 # size in grams, average in survey
      Eo = 0.8736  # from cod, 0.8736.  Make it one half or double
    } else{
      Ao = 3.11E-14 # cod
      B = 1250 # eyeball estimate from Anderson et al. 2019 report (varies by survey in reality)
      Eo = 0.8736
    }
    n = -0.208 # borrowed from cod 
    
    dat$mi = B^n*Ao*dat$po2/exp(-1*Eo/(boltz*(dat$temp+kelvin)))
    
    # prepare data and models -------------------------------------------------
   
    dat <- dplyr::select(dat, species, year, longitude_dd, latitude_dd, cpue_kg_km2,
                  o2, temp, depth, mi, po2)
    dat_ll = dat
    coordinates(dat_ll) <- c("longitude_dd", "latitude_dd")
    proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
    # convert to utm with spTransform
    dat_utm = spTransform(dat_ll, 
                          CRS("+proj=utm +zone=9 +datum=WGS84 +units=km")) # check on zone
    # convert back from sp object to data frame
    dat = as.data.frame(dat_utm)
    dat = dplyr::rename(dat, longitude = longitude_dd, 
                        latitude = latitude_dd)
 }

back.convert <- function(x, mean_orig, sd_orig) (x* sd_orig+mean_orig)

plot_map <- function(dat, column = "est") {
  ggplot(dat, aes_string("longitude", "latitude", fill = column)) +
    geom_tile() +
    facet_wrap(~year) +
    coord_fixed()
}

get_models_compare <- function() {
  formula <- c("log_depth_scaled + log_depth_scaled2 + as.factor(year)", 
                                      "log_depth_scaled + log_depth_scaled2 + as.factor(year) + temp_trawl", 
                                      "log_depth_scaled + log_depth_scaled2 + as.factor(year) + po2_trawl",
                                      "log_depth_scaled + log_depth_scaled2 + as.factor(year) + mi_trawl",
                                      "log_depth_scaled + log_depth_scaled2 + as.factor(year) + temp_trawl + po2_trawl",
                                      "log_depth_scaled + log_depth_scaled2 + as.factor(year) + temp_trawl + po2_trawl + temp_trawl * po2_trawl",
                                      "log_depth_scaled + log_depth_scaled2 + as.factor(year) + breakpt(po2_trawl)",
                                      "log_depth_scaled + log_depth_scaled2 + as.factor(year) + breakpt(po2_trawl) + temp_trawl",
                                      "log_depth_scaled + log_depth_scaled2 + as.factor(year) + breakpt(mi_trawl)",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year)", 
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + temp_jscope", 
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + po2_jscope",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + mi_jscope",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + temp_jscope + po2_jscope",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + temp_jscope + po2_jscope + temp_jscope * po2_jscope",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + breakpt(po2_jscope)",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + breakpt(po2_jscope) + temp_jscope",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + breakpt(mi_jscope)"
                                      )
}
 

get_models <- function() {
  formula <- c("log_depth_scaled + log_depth_scaled2 + as.factor(year)", 
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + temp", 
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + po2",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + mi",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + temp + po2",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + temp + po2 + temp * po2",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + breakpt(po2)",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + breakpt(po2) + temp",
               "log_depth_scaled + log_depth_scaled2 + as.factor(year) + breakpt(mi)"
  )
}                     
get_bp_parameters <- function(m) {
  fixed_effects <- m$model$par
  bp_ind <- grep("b_threshold", names(fixed_effects))
  bp <- fixed_effects[bp_ind[2]]
  slope <- fixed_effects[bp_ind[1]]
  return(c(bp, slope))
}

get_partial_resisuals<- function(m, dat) {
  # get model residuals
  fitted <- predict(m)
  residual <- fitted$cpue_kg_km2 - exp(fitted$est)
  bp_pars <- get_bp_parameters(m)
  # apply breakpoint function
  for (i in 1:nrow(dat)) dat$bp[i] <- ifelse(dat$po2[i]<bp_pars[1], dat$po2[i]*bp_pars[2],bp_pars[1]*bp_pars[2])
  partial_residuals <- residual * exp(dat$bp)
  return(partual_residuals)
}

