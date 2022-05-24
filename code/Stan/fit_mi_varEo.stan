//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower =0> nspc; // n is the number of species in the data set
  vector[N] y;
  vector[N] logb;
  vector[N] inv_temp;
  int spc[N];
}

// The parameters accepted by the model.
parameters {
  real logmu_Ao;
  real mu_Eo;
  real mu_n;
  real<lower = 0> logsigma_Ao;
  real<lower = 0> sigma_Eo;
  real<lower = 0> sigma_n;
  vector[nspc] logAo_raw;
  vector[nspc] Eo_raw;
  vector[nspc] n_raw;
  real logsigma;
}

transformed parameters {
 real logAo[nspc];
 real Eo[nspc];
 real n[nspc];
 real sigma;
 vector[N] log_lik;
 vector[N] mu;
 
 sigma = exp(logsigma);
 // note, make Ao normally distributed in log space.  This is needed
 // to generate positive values of Ao, and to generate a posterior that can
 // be used in TMB later
 for (i in 1:nspc) logAo[i] = logmu_Ao  +  logAo_raw[i] * logsigma_Ao;
 for (i in 1:nspc) Eo[i] = mu_Eo  +  Eo_raw[i] * sigma_Eo;
 for (i in 1:nspc) n[i] = mu_n  +  n_raw[i] * sigma_n;
 for (i in 1:N) mu[i] = logAo[spc[i]] + Eo[spc[i]] * inv_temp[i] + n[spc[i]] * logb[i];
 for (i in 1:N) log_lik[i] = normal_lpdf(y[i]|mu[i], sigma);

}



// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  // Priors
  for (i in 1:nspc) logAo_raw[i] ~ normal(0,1);
  for (i in 1:nspc) Eo_raw[i] ~ normal(0,1);
  for (i in 1:nspc) n_raw[i] ~ normal(0,1);
  logsigma_Ao ~ cauchy(0, 10);
  sigma_Eo ~ cauchy(0,1);
  sigma_n ~ cauchy(0,1);
  mu_n ~ cauchy(0, 1); 
  logsigma ~ cauchy(0, 10);
  logmu_Ao ~ cauchy(0, 10);
  y ~ normal(mu, sigma);
}

