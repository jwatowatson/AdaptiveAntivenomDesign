run_all_scenarios = function(sim_title,
                             model_params_true,
                             prior_model_params,
                             true_model,
                             N_trials,
                             MTT, TEL, N_max, N_batch,
                             max_increment, 
                             Randomisation_p_SOC, 
                             FORCE_RERUN, 
                             N_cores, 
                             starting_dose, 
                             SoC,
                             use_SoC_data){
  
  #************** model_based simulations ************
  #************** Well specified ******************
  tic()
  writeLines(paste(sim_title, ', model based design, well-specified, ...'))
  Full_Simulation(model_params_true = model_params_true,
                  true_model = NULL,
                  prior_model_params = prior_model_params,
                  N_trials = N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max,
                  N_batch = N_batch,
                  max_increment = max_increment,
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = sim_title,
                  FORCE_RERUN=FORCE_RERUN,
                  N_cores = N_cores,
                  design_type = 'model_based',
                  starting_dose = starting_dose,
                  SoC = SoC,
                  use_SoC_data = use_SoC_data)
  toc()
  
  #************** Mis-specified ******************
  tic()
  writeLines(paste(sim_title, ', model based design, mis-specified, ...'))
  Full_Simulation(model_params_true = NULL,
                  true_model = true_model,
                  prior_model_params = prior_model_params,
                  N_trials = N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = sim_title,
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  design_type = 'model_based',
                  starting_dose = starting_dose,
                  SoC = SoC,
                  use_SoC_data = use_SoC_data)
  toc()
  
  #************** rule_based simulations ************
  #************** Well specified ******************
  tic()
  writeLines(paste(sim_title, ', rule based design, well-specified, ...'))
  Full_Simulation(model_params_true = model_params_true,
                  true_model = NULL,
                  prior_model_params = NA,
                  N_trials = 10*N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = sim_title,
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  design_type = 'rule_based',
                  starting_dose = starting_dose,
                  SoC = SoC)
  toc()
  
  #************** Mis-specified ******************
  tic()
  writeLines(paste(sim_title, ', rule based design, mis-specified, ...'))
  Full_Simulation(model_params_true = NULL,
                  true_model = true_model,
                  prior_model_params = NA,
                  N_trials = 10*N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = sim_title,
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  design_type = 'rule_based',
                  starting_dose = starting_dose,
                  SoC = SoC)
  toc()
  
}