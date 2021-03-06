dose_response = function(x, alpha_val, beta_val, y_star){
  # x is on the log2 scale
  z = inv.logit(alpha_val + beta_val*x) - y_star
  return(z)
}

Estimate_log2_Vstar = function(model_params, true_model, MTT, TEL){

  if(is.null(model_params) & is.null(true_model)) {
    stop('One of model_params or true_model has to be specified (Estimate_log2_Vstar)')
  }
  if(!is.null(model_params) & !is.null(true_model)) {
    stop('Ambiguous whether true model or model_params')
  }
  if(is.null(model_params) & !is.null(true_model)) { # Use the mis-specified model
    TED = uniroot.all(function(x,TEL) true_model$eff(x)-TEL, lower = 0, upper = 1000,TEL=TEL)
    log2MTD = log2(uniroot.all(function(x,MTT) true_model$tox(x)-MTT, lower = 0, upper = 1000,MTT=MTT))
  }
  if(!is.null(model_params) & is.null(true_model)) {
    TED = qnorm(p = TEL, mean = model_params$mu_antivenom, sd = model_params$sd_antivenom)
    log2MTD = uniroot.all(dose_response, lower = -10, upper = 100,
                          alpha_val = model_params$alpha_tox,
                          y_star = MTT,
                          beta_val = model_params$beta_tox)
  }


  log2Vstar = min(log2(TED), log2MTD)
  out = list(TED=log2(TED), MTD=log2MTD, Vstar=log2Vstar)
  return(out)
}

# This function does the Bayesian posterior estimation
estimate_posterior_params = function(Trial_data,
                                     prior_model_params,
                                     current_model_params,
                                     use_SoC_data,
                                     iter=2*10^3){
  if(!use_SoC_data){
    if(sum(Trial_data$Randomisation_code==1) < 2){
      return(prior_model_params)
    } else {
      Trial_data = dplyr::filter(Trial_data, Randomisation_code==1)
    }
  }
  # sample from the posterior
  my_fit = rstan::optimizing(dose_response_model,
                    data = list(N = nrow(Trial_data),
                                y_tox = as.integer(Trial_data$tox_outcome),
                                y_eff = as.integer(Trial_data$eff_outcome),
                                dose = 2^Trial_data$log2_dose,
                                beta_tox_sd = prior_model_params$beta_tox_sd,
                                beta_tox = prior_model_params$beta_tox,
                                alpha_tox = prior_model_params$alpha_tox,
                                alpha_tox_sd = prior_model_params$alpha_tox_sd,
                                mu_antivenom_mean = prior_model_params$mu_antivenom,
                                mu_antivenom_sd = prior_model_params$mu_antivenom_sd,
                                sd_antivenom_mean = prior_model_params$sd_antivenom,
                                sd_antivenom_sd = prior_model_params$sd_antivenom_sd),
                   # iter = iter, chains = 1, verbose = F, refresh = -1,
                    init = list(chain1=current_model_params))

  # extract mean parameter values
  posterior_model_params = list(beta_tox = my_fit$par['beta_t'],
                                alpha_tox=  my_fit$par['alpha_t'],
                                mu_antivenom = my_fit$par['mu_antivenom'],
                                sd_antivenom = my_fit$par['sd_antivenom'])
  return(posterior_model_params)
}

# The escalation and de-escalation rules for the rule-based design
escalation_rule = function(Trial_data, current_dose, max_increment,
                           MTT, TEL, epsilon=0.01){

  ind = Trial_data$log2_dose == log2(current_dose)
  eff_current_dose = mean(Trial_data$eff_outcome[ind])
  tox_current_dose = mean(Trial_data$tox_outcome[ind])

  N_dose = sum(ind)
  if(N_dose == 0) return(current_dose)
  tox_decision = NA
  eff_decision = NA

  if(N_dose < 20){
    if(tox_current_dose >= 1/3){
      decision = -1
    }
    if(tox_current_dose > 0 & tox_current_dose < 1/3){
      decision = 0
    }
    if(tox_current_dose == 0){
      if(eff_current_dose == 1) {
        decision = 0
      } else {
        decision = 1
      }
    }

  } else {
    if(tox_current_dose > MTT ){
      decision = -1
    } else {
      if(eff_current_dose >= TEL+epsilon ) {
        decision = -1
      }
      if(eff_current_dose <= TEL-epsilon){
        decision = 1
      }
      if(eff_current_dose < TEL + epsilon & eff_current_dose > TEL - epsilon){
        decision = 0
      }
    }
  }

  next_dose = current_dose + decision * max_increment
  return(next_dose)
}
