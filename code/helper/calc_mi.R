# mock code for calculating 95 percent CI of metabolic index


# Code to generate MI estimate (confidence interval) given pO2, temp, and other info
calc_mi <- function(pO2, logw, inv_temperature, betas, var_covar, method = "smr", confidence_level = 0.95) {
  # code to calculate 95 % confidence interval of metabolic index given pO2 and temperature
  # logw is log(w / wref)
  # inv_temperature is (1 / kb) ( )1/T - 1 / Tref), where t is in kelvin and kb is boltzmann's constant
  if(method == "smr") x_predict <- as.vector(c(-1, logw, inv_temperature, -1))
  if(!method == "smr") x_predict <- c(-1, logw, inv_temperature, 0)
  
  log_mi_predict <- x_predict %*% betas + log(pO2)
  log_mi_se <-sqrt( t(x_predict) %*% var_covar %*% x_predict )
  log_lower_bound <- log_mi_predict - log_mi_se * qnorm(0.5 + confidence_level / 2) 
  log_upper_bound <- log_mi_predict + log_mi_se * qnorm(0.5 + confidence_level / 2)
  return(list ( mi = exp(log_mi_predict), lower_bound = exp(log_lower_bound), upper_bound = exp(log_upper_bound) ) )
}

betas <- taxa_estimates$parameters[,1]
calc_mi(5, log(500 / 10), inv_temperature = 0.75,  betas = betas, var_covar = taxa_estimates$var_covar, method = "smr", confidence_level = 0.95)
