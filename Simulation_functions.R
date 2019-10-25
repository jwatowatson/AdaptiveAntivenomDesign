generate_sim_truth_outcomes = function(doses_administered, model_params, true_model){
  
  N_patients = length(doses_administered)
  # options for true_model:
  # logistic: then model_params need to be 
  ptox = inv.logit(model_params$alpha_tox + model_params$beta_tox * doses_administered)
  peff = inv.logit(model_params$alpha_eff + model_params$beta_eff * doses_administered)
  
  tox_outcome = eff_outcome = array(dim=N_patients)
  for(i in 1:N_patients){
    tox_outcome[i] = sample(x = 0:1,size = 1, 
                            prob = c(1-ptox[i], ptox[i]), 
                            replace = T)
    eff_outcome[i] = sample(x = 0:1,size = 1, 
                            prob = c(1-peff[i], peff[i]), 
                            replace = T)
  }
  out = data.frame(tox_outcome=tox_outcome,eff_outcome=eff_outcome)
  return(out)
}


# generate outcomes under the simulation model truth
Generate_outcomes = function(model_params, true_model='logistic',
                             log2_dose, N_patients=1,  
                             p = p, SoC = SoC){
  
  # generate randomisation code
  Randomisation_code = sample(1:2, 
                              size = N_patients, 
                              replace = T, 
                              prob = c(1-p,p))
  doses_administered = c(log2_dose, log2(SoC))[Randomisation_code]
  
  out = generate_sim_truth_outcomes(doses_administered, model_params, true_model)
  
  return(data.frame(tox_outcome = out$tox_outcome, 
                    eff_outcome = out$eff_outcome, 
                    log2_dose = doses_administered,
                    Randomisation_code = Randomisation_code))
}


# Run a single trial under the model-based design
run_trial = function(model_params_true, prior_model_params, N_max=200, 
                     max_increment, MTT, TEL, N_batch, p, starting_dose, SoC){
  
  require(doParallel)
  # Set up parameters and choose the optimal dose based on the prior
  N_current = 0
  out = Estimate_log2_Vstar(model_params = prior_model_params,
                                       MTT = MTT, TEL = TEL)
  
  current_dose = starting_dose # current_dose is the current adaptive dose
  # simulate outcomes for the 1st set of patients
  Trial_data = Generate_outcomes(model_params = model_params_true, 
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
    posterior_model_params  = estimate_posterior_params(Trial_data=Trial_data,
                                                        prior_model_params=prior_model_params)
    out = Estimate_log2_Vstar(model_params = posterior_model_params,
                                         MTT = MTT, TEL = TEL)
    # update dose according to model
    current_dose = min(max(2^Trial_data$log2_dose) + max_increment, 
                       round(2^out$Vstar))
    # simulate outcomes 
    Next_batch = Generate_outcomes(model_params = model_params_true, 
                                   log2_dose = log2(current_dose),
                                   N_patients = N_batch, 
                                   p = p, 
                                   SoC = SoC)
  
    Next_batch$optimal_dose = 2^out$Vstar
    Next_batch$MTD = 2^out$MTD
    Next_batch$TED = 2^out$TED
    
    Trial_data = rbind(Trial_data,Next_batch)
    N_current = N_current + N_batch
  }
  return(Trial_data)
}

## Run a single trial under the rule-based design (3+3 type trial)
run_3plus3_trial = function(model_params_true, N_max = 200, 
                            starting_dose=8, N_batch=3, max_increment,
                            MTT, TEL, epsilon, p){
  
  # Set up parameters and choose the optimal dose based on the prior
  N_current = 0
  current_dose = starting_dose
  
  # iteratively simulate batches of patients until reach max sample size
  while(N_current < N_max){
    # simulate outcomes for the 1st set of patients
    Next_batch = Generate_outcomes(model_params = model_params_true, 
                                   log2_dose = log2(current_dose),
                                   N_patients = N_batch, 
                                   p = p, 
                                   SoC = NA)
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
    
    
  }
  return(Trial_data)
}


# Wrapper function that plots the prior and runs an adaptive design trial
Full_Simulation = function(model_params_true, 
                           prior_model_params,
                           N_trials = 20, # number of simulated trials
                           TEL = 0.95,  # target efficacy level
                           MTT = 0.05,  # maximum tolerated toxicity
                           N_max = 200, # max number of patients per trial
                           N_batch = 1, # number of patients until re-assess optimal dose
                           max_increment = 1,
                           Randomisation_p_SOC,
                           sim_title='',
                           FORCE_RERUN, 
                           N_cores=NA, 
                           individ_plots=F,
                           design_type='Adaptive',
                           starting_dose, 
                           SoC, 
                           epsilon = 0.01
){
  
  f_name = paste(sim_title,'_',design_type, '.RData', sep = '')
  if(! design_type %in% c('Adaptive','3+3')) stop('Unknown design_type')
  
  plot_prior_versus_truth(model_params_true = model_params_true,
                          prior_model_params = prior_model_params,
                          individ_plots = individ_plots,
                          sim_title = sim_title)
  
  ##*** Run the simulated trial ***
  if(FORCE_RERUN | (!f_name %in% list.files(path = 'SimulationOutputs/')) ){
    
    if(is.na(N_cores)) N_cores = detectCores()-1
    registerDoParallel(N_cores)
    Summary_trials = foreach(s=1:N_trials, .combine = cbind) %dopar% {
      
      if(design_type == 'Adaptive'){
        Trial_data = run_trial(model_params_true = model_params_true,
                               N_max = N_max, 
                               prior_model_params = prior_model_params, 
                               N_batch = N_batch,
                               max_increment = max_increment,
                               MTT = MTT, 
                               TEL = TEL,
                               starting_dose = starting_dose, 
                               p = Randomisation_p_SOC,
                               SoC = SoC)
        ind_NA = Trial_data$Randomisation_code==2
        Trial_data$log2_dose[ind_NA] = NA
        out = cbind(opt_dose = Trial_data$optimal_dose, 
                    MTD = Trial_data$MTD,
                    TED = Trial_data$TED,
                    assigned_dose = 2^Trial_data$log2_dose)
      } 
      if(design_type == '3+3'){
        Trial_data = run_3plus3_trial(model_params_true = model_params_true, N_max = N_max, 
                                      starting_dose = starting_dose, 
                                      N_batch = N_batch, 
                                      max_increment = max_increment,
                                      MTT = MTT, 
                                      TEL = TEL, epsilon=epsilon)
        out = cbind(opt_dose = NA, assigned_dose = 2^Trial_data$log2_dose)
      }
      out
    }
    save(Summary_trials, file = paste('SimulationOutputs/',f_name,sep=''))
  } else {
    load(paste('SimulationOutputs/',f_name,sep=''))
  }
  
  out=Estimate_log2_Vstar(model_params = model_params_true,MTT = .05,TEL = .95)
  print('done the trial simulation, now plotting results')
  ## *** Plot the estimated optimal dose ***
  # This is only for the Adaptive type design
  if(design_type == 'Adaptive'){
    cols_opt_dose = grep('opt_dose', colnames(Summary_trials))
    Summary_optimal_dose = Summary_trials[, cols_opt_dose]
    q_upper = apply(Summary_optimal_dose, 1, quantile, prob=0.05, na.rm=T)
    q_lower = apply(Summary_optimal_dose, 1, quantile, prob=0.95, na.rm=T)
    plot(1:nrow(Summary_optimal_dose), Summary_optimal_dose[,1],
         type='n', main = 'The estimated optimal doses', ylim=range(c(q_upper,q_lower)),
         xlab = 'Patient recruitment index', ylab = 'Estimated optimal dose')
    
    
    abline(h = 2^out$Vstar,col='red',lwd=2)
    lines(1:nrow(Summary_optimal_dose), apply(Summary_optimal_dose, 1, mean, na.rm=T),lwd=3)
    lines(1:nrow(Summary_optimal_dose), q_upper, lwd=1, lty=2)
    lines(1:nrow(Summary_optimal_dose), q_lower, lwd=1, lty=2)
    
    legend('topright', legend = c('mean value','5&95 quantiles'), 
           lwd = c(2,1), lty=c(1,2), inset=0.01, bty='y', bg = 'white')
  }
  
  ## *** Plot the assigned dose ***
  cols_assigned_dose = grep('assigned_dose', colnames(Summary_trials))
  Summary_assigned_dose = Summary_trials[, cols_assigned_dose]
  ind = complete.cases(Summary_assigned_dose[,1])
  q_upper = apply(Summary_assigned_dose, 1, quantile, prob=0.05, na.rm=T)
  q_lower = apply(Summary_assigned_dose, 1, quantile, prob=0.95, na.rm=T)
  
  if(design_type=='Adaptive') design_main = 'The assigned doses for adaptive arm'
  if(design_type=='3+3') design_main = 'The assigned dose'
  
  plot(which(ind), Summary_assigned_dose[ind,1], col = adjustcolor('grey',alpha.f = .5),
       type='n', main = design_main, 
       ylim=range(c(q_upper,q_lower)),xlim=c(1,nrow(Summary_assigned_dose)),
       xlab = 'Patient recruitment index', ylab = 'Assigned dose')
  
  abline(h = 2^out$Vstar,col='red',lwd=3)
  lines(1:nrow(Summary_assigned_dose), apply(Summary_assigned_dose, 1, median, na.rm=T),lwd=2)
  lines(1:nrow(Summary_assigned_dose), q_upper, lwd=1, lty=2)
  lines(1:nrow(Summary_assigned_dose), q_lower, lwd=1, lty=2)
  
  legend('topright', legend = c('median value','5&95 quantiles'), 
         lwd = c(2,1), lty=c(1,2), inset=0.01, bty='y', bg = 'white')
  
  # Histogram of final assigned dose
  last_dose = round(Summary_assigned_dose[nrow(Summary_assigned_dose),])
  barplot(table(last_dose), xlab = 'Dose', ylab = 'Number of trials', 
          main = 'Final estimated optimal dose')
  
}
