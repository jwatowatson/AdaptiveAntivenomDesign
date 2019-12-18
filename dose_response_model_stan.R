library(rstan)
rstan_options(auto_write = TRUE)

mod = "

data {
  int<lower=0> N;
  int<lower=0,upper=1> y_tox[N];
  int<lower=0,upper=1> y_eff[N];
  vector[N] dose;
  real beta_tox;
  real beta_tox_sd;
  real alpha_tox;
  real alpha_tox_sd;
  real<lower=0> mu_antivenom_mean;
  real<lower=0> mu_antivenom_sd;
  real<lower=0> sd_antivenom_mean;
  real<lower=0> sd_antivenom_sd;
}

parameters {
  real beta_t;
  real alpha_t;
  real<lower=0> mu_antivenom;
  real<lower=0> sd_antivenom;
}

model {
  // priors
  beta_t ~ normal(beta_tox,beta_tox_sd);
  alpha_t ~ normal(alpha_tox,alpha_tox_sd);
  mu_antivenom ~ normal(mu_antivenom_mean, mu_antivenom_sd);
  sd_antivenom ~ normal(sd_antivenom_mean, sd_antivenom_sd);
  
  // likelihood
  y_tox ~ bernoulli_logit(alpha_t + beta_t * log2(dose));
  y_eff ~ bernoulli(Phi_approx( (dose - mu_antivenom)/sd_antivenom));
}
"

dose_response_model = stan_model(model_code = mod)
