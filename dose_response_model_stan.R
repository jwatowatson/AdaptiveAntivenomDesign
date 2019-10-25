library(rstan)

mod = "

data {
  int<lower=1> N;
  int<lower=0,upper=1> y_tox[N];
  int<lower=0,upper=1> y_eff[N];
  vector[N] x;
  real beta_tox;
  real beta_tox_sd;
  real alpha_tox;
  real alpha_tox_sd;
  real beta_eff;
  real beta_eff_sd;
  real alpha_eff;
  real alpha_eff_sd;
}

parameters {
  real beta_t;
  real alpha_t;
  real beta_e;
  real alpha_e;
}

model {
  // priors
  beta_t ~ normal(beta_tox,beta_tox_sd);
  alpha_t ~ normal(alpha_tox,alpha_tox_sd);
  beta_e ~ normal(beta_eff,beta_eff_sd);
  alpha_e ~ normal(alpha_eff,alpha_eff_sd);
  
  // likelihood
  y_tox ~ bernoulli_logit(alpha_t + beta_t * x);
  y_eff ~ bernoulli_logit(alpha_e + beta_e * x);
}
"

dose_response_model = stan_model(model_code = mod)
