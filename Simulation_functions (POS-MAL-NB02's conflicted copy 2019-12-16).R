# This generates outcomes from the simulation truth
generate_sim_truth_outcomes = function(log2_doses, 
                                       model_params, 
                                       true_model){
  
  if(is.null(model_params_true) & is.null(true_model)) {
    stop('One of model_params_true or true_model has to be specified')
  }
  
  N_patients = length(log2_doses)
  
  if(!is.null(model_params)){
    # well-specified logistic model
    ptox = inv.logit(model_params$alpha_tox + model_params$beta_tox * log2_doses)
    peff = pnorm(log2_doses/model_params$beta_venom, 
                 mean = model_params$mu_venom, 
                 sd = model_params$sd_venom)
  } else {
    # mis-specified linear model
    ptox = true_model$tox(2^log2_doses)
    peff = true_model$eff(2^log2_doses)
  }
  
  # sample outcomes based on probabilities
  tox_outcome = eff_outcome = array(dim=N_patients)
  for(i in 1:N_patients){
    tox_outcome[i] = sample(x = 0:1,size = 1, 
                            prob = c(1-ptox[i], ptox[i]), 
                            replace = T)
    eff_outcome[i] = sample(x = 0:1,size = 1, 
                            prob = c(1-peff[i], peff[i]), 
                            replace = T)
  }
  out = data.frame(tox_outcome=tox_outcome,
                   eff_outcome=eff_outcome)
  return(out)
}

# generate outcomes under the simulation model truth
Generate_outcomes = function(model_params, 
                             true_model,
                             log2_dose, 
                             N_patients,  
                             p = p, 
                             SoC = SoC){
  
  # generate randomisation code
  Randomisation_code = sample(1:2, 
                              size = N_patients, 
                              replace = T, 
                              prob = c(1-p,p))
  
  # determine administered doses based on randomisation codes
  log2_doses_administered = c(log2_dose, log2(SoC))[Randomisation_code]
  
  # generate outcomes 
  out = generate_sim_truth_outcomes(log2_doses = log2_doses_administered, 
                                    model_params = model_params,
                                    true_model = true_model)
  
  return(data.frame(tox_outcome = out$tox_outcome, 
                    eff_outcome = out$eff_outcome, 
                    log2_dose = log2_doses_administered,
                    Randomisation_code = Randomisation_code))
}


# Run a single trial under the model-based design
run_trial = function(model_params_true, 
                     true_model,
                     prior_model_params, 
                     N_max, 
                     max_increment, MTT, TEL, 
                     N_batch, p, starting_dose, SoC,
                     use_SoC_data){
  
  # Set up parameters and choose the optimal dose based on the prior
  N_current = 0
  out = Estimate_log2_Vstar(model_params = prior_model_params,
                            true_model = NULL,
                            MTT = MTT, TEL = TEL)
  
  current_dose = starting_dose # current_dose is the current model_based dose
  # simulate outcomes for the 1st set of patients
  Trial_data = Generate_outcomes(model_params = model_params_true, 
                                 true_model = true_model,
                                 log2_dose = log2(current_dose),
                                 N_patients = N_batch, 
                                 p = p, 
                                 SoC = SoC)
  Trial_data$optimal_dose = 2^out$Vstar
  Trial_data$MTD = 2^out$MTD
  Trial_data$TED = 2^out$TED
  N_current = N_current + N_batch
  
  # iteratively simulate batches of patients until reach max sample size
  while(N_current < N_max){
    # Estimate Bayesian MAP
    posterior_model_params = estimate_posterior_params(Trial_data=Trial_data,
                                                       prior_model_params=prior_model_params,
                                                       use_SoC_data=use_SoC_data)
    out = Estimate_log2_Vstar(model_params = posterior_model_params,
                              true_model = NULL,
                              MTT = MTT, TEL = TEL)
    # update dose according to model
    current_dose = min(max(2^Trial_data$log2_dose) + max_increment, 
                       round(2^out$Vstar))
    if(current_dose<1) {
      print('Estimated optimal dose less than 1 vial')
      current_dose=1
    }
    # simulate outcomes 
    Next_batch = Generate_outcomes(model_params = model_params_true, 
                                   true_model = true_model,
                                   log2_dose = log2(current_dose),
                                   N_patients = N_batch, 
                                   p = p, 
                                   SoC = SoC)
    # Record the estimated optimal dose, estimated MTD, estimated TED
    Next_batch$optimal_dose = 2^out$Vstar
    Next_batch$MTD = 2^out$MTD
    Next_batch$TED = 2^out$TED
    
    Trial_data = rbind(Trial_data,Next_batch)
    N_current = N_current + N_batch
  }
  return(Trial_data)
}

## Run a single trial under the rule-based design (rule_based type trial)
run_3plus3_trial = function(model_params_true, 
                            true_model, 
                            N_max, 
                            starting_dose, 
                            N_batch, 
                            max_increment,
                            MTT, TEL, epsilon, p, SoC){
  
  # Set up parameters and choose the optimal dose based on the prior
  N_current = 0
  current_dose = starting_dose
  
  # iteratively simulate batches of patients until reach max sample size
  while(N_current < N_max){
    # simulate outcomes for the 1st set of patients
    Next_batch = Generate_outcomes(model_params = model_params_true, 
                                   true_model = true_model,
                                   log2_dose = log2(current_dose),
                                   N_patients = N_batch, 
                                   p = p, 
                                   SoC = SoC)
    if(N_current == 0) {
      Trial_data = Next_batch
    } else {
      Trial_data = rbind(Trial_data,Next_batch)
    }
    N_current = N_current + N_batch
    
    current_dose = escalation_rule(Trial_data = Trial_data, 
                                   current_dose = current_dose, 
                                   max_increment = max_increment,
                                   MTT = MTT, TEL = TEL, epsilon = epsilon)
    if(current_dose<1) {
      print('Estimated optimal dose less than 1 vial')
      current_dose=1
    }
  }
  return(Trial_data)
}


# Wrapper function that runs an adaptive design trial (model or rule based)
Full_Simulation = function(model_params_true, # parameters for sim truth when well-specified
                           true_model, # sim truth when mis-specified
                           prior_model_params, # parameters for the prior distribution
                           N_trials, # number of simulated trials
                           TEL,  # target efficacy level
                           MTT,  # maximum tolerated toxicity
                           N_max, # max number of patients per trial
                           N_batch, # number of patients until re-assess optimal dose
                           max_increment, # maximum dose increment
                           Randomisation_p_SOC, # randomisation ratio to SoC
                           sim_title='',
                           FORCE_RERUN, 
                           N_cores=NA, 
                           design_type,
                           starting_dose, 
                           SoC, 
                           epsilon = 0.01,
                           use_SoC_data = T # only for the model-based: specified if the model update uses the SoC data
){
  if(is.null(model_params_true) & is.null(true_model)) {
    stop('One of model_params_true or true_model has to be specified')
  }
  if(!is.null(model_params_true) & is.null(true_model)) { 
    specification = 'well_specified'
  } 
  if(is.null(model_params_true) & !is.null(true_model)) {
    specification = 'mis_specified'
  }
  if(! design_type %in% c('model_based','rule_based')) {
    stop('Unknown design_type. This has to be = model_based / rule_based')
  }
  if(use_SoC_data) {
    SoC_use = 'All_data'
  } else {
    SoC_use = 'Adaptive_data'
  }
  
  ##*** Run the simulated trial ***
  # We don't rerun unless FORCE_RERUN is TRUE
  f_name = paste(sim_title,'_',design_type,'_',specification, '_', SoC_use,'.RData', sep = '')
  print(f_name)
  if(FORCE_RERUN | (!f_name %in% list.files(path = 'SimulationOutputs/')) ){
    
    if(is.na(N_cores)) N_cores = detectCores()-1
    registerDoParallel(N_cores)
    Summary_trials = foreach(s=1:N_trials, .combine = cbind) %dopar% {
      
      if(design_type == 'model_based'){
        Trial_data = run_trial(model_params_true = model_params_true,
                               true_model = true_model,
                               N_max = N_max, 
                               prior_model_params = prior_model_params, 
                               N_batch = N_batch,
                               max_increment = max_increment,
                               MTT = MTT, 
                               TEL = TEL,
                               starting_dose = starting_dose, 
                               p = Randomisation_p_SOC,
                               SoC = SoC,
                               use_SoC_data = use_SoC_data)
      } 
      if(design_type == 'rule_based'){
        Trial_data = run_3plus3_trial(model_params_true = model_params_true, 
                                      true_model = true_model,
                                      N_max = N_max, 
                                      starting_dose = starting_dose, 
                                      N_batch = N_batch, 
                                      max_increment = max_increment,
                                      MTT = MTT, 
                                      TEL = TEL, 
                                      epsilon=epsilon,
                                      p = Randomisation_p_SOC,
                                      SoC = SoC)
      }
      ind_NA = Trial_data$Randomisation_code==2
      Trial_data$log2_dose[ind_NA] = NA
      out = cbind(opt_dose = Trial_data$optimal_dose, 
                  MTD = Trial_data$MTD,
                  TED = Trial_data$TED,
                  assigned_dose = 2^Trial_data$log2_dose)
      out
    }
    save(Summary_trials, file = paste('SimulationOutputs/',f_name,sep=''))
  } else {
    load(paste('SimulationOutputs/',f_name,sep=''))
  }
  
}
