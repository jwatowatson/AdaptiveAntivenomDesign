dose_response = function(x, alpha_val, beta_val, y_star){
  # x is on the log2 scale
  z = inv.logit(alpha_val + beta_val*x) - y_star
  return(z)
}

Estimate_log2_Vstar = function(model_params, MTT, TEL){
  TED = uniroot.all(dose_response, lower = -10, upper = 100, 
                    alpha_val = model_params$alpha_eff,
                    y_star = TEL, 
                    beta_val = model_params$beta_eff)
  
  MTD = uniroot.all(dose_response, lower = -10, upper = 100, 
                    alpha_val = model_params$alpha_tox,
                    y_star = MTT, 
                    beta_val = model_params$beta_tox)
  
  Vstar = min(TED, MTD)
  out = list(TED=TED, MTD=MTD, Vstar=Vstar)
  return(out)
}

# This function does the Bayesian posterior estimation
estimate_posterior_params = function(Trial_data=Trial_data,
                                     prior_model_params=prior_model_params,iter=10^3){
  
  my_fit = sampling(dose_response_model,
                    data = list(N = nrow(Trial_data),
                                y_tox = as.integer(Trial_data$tox_outcome),
                                y_eff = as.integer(Trial_data$eff_outcome),
                                x = Trial_data$log2_dose,
                                beta_tox_sd = prior_model_params$beta_tox_sd,
                                beta_tox = prior_model_params$beta_tox,
                                alpha_tox = prior_model_params$alpha_tox,
                                alpha_tox_sd = prior_model_params$alpha_tox_sd,
                                beta_eff_sd = prior_model_params$beta_eff_sd,
                                beta_eff = prior_model_params$beta_eff,
                                alpha_eff = prior_model_params$alpha_eff,
                                alpha_eff_sd = prior_model_params$alpha_eff_sd),
                    iter = iter, 
                    chains = 1, verbose = F, refresh = -1)
  
  model_params = summary(my_fit, probs = NA)$summary
  
  posterior_model_params = list(beta_tox = model_params['beta_t','mean'],
                                alpha_tox=  model_params['alpha_t','mean'],
                                beta_eff = model_params['beta_e','mean'], 
                                alpha_eff = model_params['alpha_e','mean'])
  return(posterior_model_params)
}

# The escalation and de-escalation rules for the rule-based design
escalation_rule = function(Trial_data, current_dose, max_increment, MTT, TEL, epsilon=0.01){
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