---
title: "Simulation study comparing rule-based and model-based adaptive designs"
author: "James Watson"
date: "October 10, 2019"
output:
  pdf_document: 
    fig_caption: yes
    keep_md: yes
---



## Setup model parameters


```r
source('Core_functions.R')
source('Simulation_functions.R')
source('plotting_functions.R')
source('wrapper_function_all_scenarios.R')
TEL = 0.95  # Target Efficacy Level
MTT = 0.05  # Maximum Tolerated Toxicity
N_max = 250 # Maximum number of patients recruited
max_increment = 1 # Maximum dose increment to doses previously unseen
N_trials = 2000
N_batch = 3
Randomisation_p_SOC = 0.2 # proportion randomised to the standard of care dose
starting_dose = 12 # starting dose in adaptive arm
SoC = 8 # Standard of Care

# function that solves for the beta value based on interpretable parameters
solve_beta = function(alpha_val, v_star, y_star){
  beta_val = ( logit(y_star) - alpha_val ) / log2(v_star)
  return(beta_val)
}
```

## Priors


```r
#******* Prior point estimates *********
Prior_TED = 12;                   # prior estimate of the Target Efficacious Dose
Prior_alpha_eff = logit(1/10)       # prior estimate of the efficacy with one vial
Prior_beta_eff = solve_beta(alpha_val = Prior_alpha_eff, v_star = Prior_TED, y_star = TEL)

Prior_MTD = 32;                    # prior estimate of the Maximum Tolerated Dose
Prior_alpha_tox = logit(1/1000)    # prior estimate of the toxicity after 1 vial
Prior_beta_tox = solve_beta(alpha_val = Prior_alpha_tox, v_star = Prior_MTD, y_star = MTT)

#******* Prior uncertainty estimates *******
prior_model_params = list(beta_tox = Prior_beta_tox,
                          beta_tox_sd = .05,
                          alpha_tox=Prior_alpha_tox,
                          alpha_tox_sd = 2,
                          beta_eff=Prior_beta_eff, 
                          beta_eff_sd = .05, 
                          alpha_eff=Prior_alpha_eff,
                          alpha_eff_sd = 2)
```

## Mis-specified dose response

This function computes linear dose-response curves based on the parameters (MTT, TEL etc).


```r
make_linear_dose_response = function(MTT, true_MTD, tox_zero=1,
                                     TEL, true_TED, eff_zero=1){
  xs = seq(0, 1000, by=.1)
  tox_MTD = true_MTD
  slope_tox = MTT/(tox_MTD - tox_zero)
  intercept_tox = -slope_tox * tox_zero
  ys = slope_tox*xs + intercept_tox
  ys[ys>1]=1; ys[ys<0]=0
  f_true_tox = approxfun(x = xs, y = ys, rule = 2)
  
  eff_TED = true_TED
  slope_eff = TEL/(eff_TED - eff_zero)
  intercept_eff = -slope_eff * eff_zero
  ys = slope_eff*xs + intercept_eff
  ys[ys>1]=1; ys[ys<0]=0
  f_true_eff = approxfun(x = xs, ys, rule = 2)
  true_model = list(tox = f_true_tox, eff = f_true_eff)
  
  return(true_model)
}
```

## Simulation 1


```r
true_alpha_eff = logit(1/50)
true_TED = 20  # 
true_beta_eff = solve_beta(alpha_val = true_alpha_eff, v_star = true_TED, y_star = TEL)


true_alpha_tox = logit(1/500)  # toxicity at 1 vial
true_MTD = 8  # simulation truth for the MTD 
true_beta_tox = solve_beta(alpha_val = true_alpha_tox, v_star = true_MTD, y_star = MTT)

# well-specified simulation truth
model_params_true = list(beta_tox = true_beta_tox,
                         alpha_tox=true_alpha_tox,
                         beta_eff=true_beta_eff, 
                         alpha_eff=true_alpha_eff)

# mis-specified simulation truth
true_model = make_linear_dose_response(MTT, true_MTD, tox_zero=1,
                                       TEL, true_TED, eff_zero=1)

run_all_scenarios(model_params_true = model_params_true,
                  true_model = true_model,
                  prior_model_params = prior_model_params,
                  N_trials = N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = 'Simulation scenario 1',
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  starting_dose = starting_dose,
                  SoC = SoC, use_SoC_data = T)
```

```
## Simulation scenario 1 , model based design, well-specified, ...
## [1] "Simulation scenario 1_model_based_well_specified_All_data.RData"
## 0.108 sec elapsed
## Simulation scenario 1 , model based design, mis-specified, ...
## [1] "Simulation scenario 1_model_based_mis_specified_All_data.RData"
## 0.055 sec elapsed
## Simulation scenario 1 , rule based design, well-specified, ...
## [1] "Simulation scenario 1_rule_based_well_specified_All_data.RData"
## 0.229 sec elapsed
## Simulation scenario 1 , rule based design, mis-specified, ...
## [1] "Simulation scenario 1_rule_based_mis_specified_All_data.RData"
## 0.188 sec elapsed
```

```r
# Comparison in well-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 1',
                      model_params_true = model_params_true,
                      true_model = NULL,
                      prior_model_params = prior_model_params,use_SoC_data = T)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_1-1.png)<!-- --> 

```
## For the rule-based design, 13% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 27% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
# Comparison in mis-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 1',
                      model_params_true = NULL,
                      true_model = true_model,
                      prior_model_params = prior_model_params,use_SoC_data = T)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_1-2.png)<!-- --> 

```
## For the rule-based design, 5% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 22% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
run_all_scenarios(model_params_true = model_params_true,
                  true_model = true_model,
                  prior_model_params = prior_model_params,
                  N_trials = N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = 'Simulation scenario 1',
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  starting_dose = starting_dose,
                  SoC = SoC, use_SoC_data = F)
```

```
## Simulation scenario 1 , model based design, well-specified, ...
## [1] "Simulation scenario 1_model_based_well_specified_Adaptive_data.RData"
## 0.059 sec elapsed
## Simulation scenario 1 , model based design, mis-specified, ...
## [1] "Simulation scenario 1_model_based_mis_specified_Adaptive_data.RData"
## 0.061 sec elapsed
## Simulation scenario 1 , rule based design, well-specified, ...
## [1] "Simulation scenario 1_rule_based_well_specified_All_data.RData"
## 0.101 sec elapsed
## Simulation scenario 1 , rule based design, mis-specified, ...
## [1] "Simulation scenario 1_rule_based_mis_specified_All_data.RData"
## 0.102 sec elapsed
```

```r
# Comparison in well-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 1',
                      model_params_true = model_params_true,
                      true_model = NULL,
                      prior_model_params = prior_model_params,use_SoC_data = F)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_1-3.png)<!-- --> 

```
## For the rule-based design, 13% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 24% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
# Comparison in mis-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 1',
                      model_params_true = NULL,
                      true_model = true_model,
                      prior_model_params = prior_model_params,use_SoC_data = F)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_1-4.png)<!-- --> 

```
## For the rule-based design, 5% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 20% of trials give patient 252 a dose within +/-10% of the true optimal dose
```


## Simulation 2


```r
true_alpha_eff = logit(1/20)
true_TED = 8  # simulation truth for the MED 
true_beta_eff = solve_beta(alpha_val = true_alpha_eff, v_star = true_TED, y_star = TEL)

true_alpha_tox = logit(1/500)  # toxicity at 1 vial
true_MTD = 20  # simulation truth for the MTD 
true_beta_tox = solve_beta(alpha_val = true_alpha_tox, v_star = true_MTD, y_star = MTT)


# well-specified simulation truth
model_params_true = list(beta_tox = true_beta_tox,
                         alpha_tox = true_alpha_tox,
                         beta_eff = true_beta_eff, 
                         alpha_eff = true_alpha_eff)

# mis-specified simulation truth
true_model = make_linear_dose_response(MTT, true_MTD, tox_zero=1,
                                       TEL, true_TED, eff_zero=1)

run_all_scenarios(model_params_true = model_params_true,
                  true_model = true_model,
                  prior_model_params = prior_model_params,
                  N_trials = N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = 'Simulation scenario 2',
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  starting_dose = starting_dose,
                  SoC = SoC,
                  use_SoC_data = T)
```

```
## Simulation scenario 2 , model based design, well-specified, ...
## [1] "Simulation scenario 2_model_based_well_specified_All_data.RData"
## 0.053 sec elapsed
## Simulation scenario 2 , model based design, mis-specified, ...
## [1] "Simulation scenario 2_model_based_mis_specified_All_data.RData"
## 0.053 sec elapsed
## Simulation scenario 2 , rule based design, well-specified, ...
## [1] "Simulation scenario 2_rule_based_well_specified_All_data.RData"
## 0.132 sec elapsed
## Simulation scenario 2 , rule based design, mis-specified, ...
## [1] "Simulation scenario 2_rule_based_mis_specified_All_data.RData"
## 0.118 sec elapsed
```

```r
# Comparison in well-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 2',
                      model_params_true = model_params_true,
                      true_model = NULL,
                      prior_model_params = prior_model_params,use_SoC_data = T)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_2-1.png)<!-- --> 

```
## For the rule-based design, 32% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 42% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
# Comparison in mis-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 2',
                      model_params_true = NULL,
                      true_model = true_model,
                      prior_model_params = prior_model_params,use_SoC_data = T)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_2-2.png)<!-- --> 

```
## For the rule-based design, 62% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 75% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
run_all_scenarios(model_params_true = model_params_true,
                  true_model = true_model,
                  prior_model_params = prior_model_params,
                  N_trials = N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = 'Simulation scenario 2',
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  starting_dose = starting_dose,
                  SoC = SoC,
                  use_SoC_data = F)
```

```
## Simulation scenario 2 , model based design, well-specified, ...
## [1] "Simulation scenario 2_model_based_well_specified_Adaptive_data.RData"
## 0.053 sec elapsed
## Simulation scenario 2 , model based design, mis-specified, ...
## [1] "Simulation scenario 2_model_based_mis_specified_Adaptive_data.RData"
## 0.053 sec elapsed
## Simulation scenario 2 , rule based design, well-specified, ...
## [1] "Simulation scenario 2_rule_based_well_specified_All_data.RData"
## 0.104 sec elapsed
## Simulation scenario 2 , rule based design, mis-specified, ...
## [1] "Simulation scenario 2_rule_based_mis_specified_All_data.RData"
## 0.101 sec elapsed
```

```r
# Comparison in well-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 2',
                      model_params_true = model_params_true,
                      true_model = NULL,
                      prior_model_params = prior_model_params,use_SoC_data = F)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_2-3.png)<!-- --> 

```
## For the rule-based design, 32% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 41% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
# Comparison in mis-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 2',
                      model_params_true = NULL,
                      true_model = true_model,
                      prior_model_params = prior_model_params,use_SoC_data = F)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_2-4.png)<!-- --> 

```
## For the rule-based design, 62% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 75% of trials give patient 252 a dose within +/-10% of the true optimal dose
```


## Simulation 3


```r
true_alpha_eff = logit(1/100)
true_TED = 60  # simulation truth for the MED 
true_beta_eff = solve_beta(alpha_val = true_alpha_eff, v_star = true_TED, y_star = TEL)


true_alpha_tox = logit(1/1000)  # toxicity at 1 vial
true_MTD = 30  # simulation truth for the MTD 
true_beta_tox = solve_beta(alpha_val = true_alpha_tox, v_star = true_MTD, y_star = MTT)


# well-specified simulation truth
model_params_true = list(beta_tox = true_beta_tox,
                         alpha_tox = true_alpha_tox,
                         beta_eff = true_beta_eff, 
                         alpha_eff = true_alpha_eff)

# mis-specified simulation truth
true_model = make_linear_dose_response(MTT, true_MTD, tox_zero=1,
                                       TEL, true_TED, eff_zero=1)

run_all_scenarios(model_params_true = model_params_true,
                  true_model = true_model,
                  prior_model_params = prior_model_params,
                  N_trials = N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = 'Simulation scenario 3',
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  starting_dose = starting_dose,
                  SoC = SoC, use_SoC_data = T)
```

```
## Simulation scenario 3 , model based design, well-specified, ...
## [1] "Simulation scenario 3_model_based_well_specified_All_data.RData"
## 0.055 sec elapsed
## Simulation scenario 3 , model based design, mis-specified, ...
## [1] "Simulation scenario 3_model_based_mis_specified_All_data.RData"
## 0.059 sec elapsed
## Simulation scenario 3 , rule based design, well-specified, ...
## [1] "Simulation scenario 3_rule_based_well_specified_All_data.RData"
## 0.118 sec elapsed
## Simulation scenario 3 , rule based design, mis-specified, ...
## [1] "Simulation scenario 3_rule_based_mis_specified_All_data.RData"
## 0.119 sec elapsed
```

```r
# Comparison in well-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 3',
                      model_params_true = model_params_true,
                      true_model = NULL,
                      prior_model_params = prior_model_params,use_SoC_data = T)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_3-1.png)<!-- --> 

```
## For the rule-based design, 38% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 28% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
# Comparison in mis-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 3',
                      model_params_true = NULL,
                      true_model = true_model,
                      prior_model_params = prior_model_params, use_SoC_data = T)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_3-2.png)<!-- --> 

```
## For the rule-based design, 28% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 27% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
run_all_scenarios(model_params_true = model_params_true,
                  true_model = true_model,
                  prior_model_params = prior_model_params,
                  N_trials = N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = 'Simulation scenario 3',
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  starting_dose = starting_dose,
                  SoC = SoC, use_SoC_data = F)
```

```
## Simulation scenario 3 , model based design, well-specified, ...
## [1] "Simulation scenario 3_model_based_well_specified_Adaptive_data.RData"
## 0.055 sec elapsed
## Simulation scenario 3 , model based design, mis-specified, ...
## [1] "Simulation scenario 3_model_based_mis_specified_Adaptive_data.RData"
## 0.055 sec elapsed
## Simulation scenario 3 , rule based design, well-specified, ...
## [1] "Simulation scenario 3_rule_based_well_specified_All_data.RData"
## 0.133 sec elapsed
## Simulation scenario 3 , rule based design, mis-specified, ...
## [1] "Simulation scenario 3_rule_based_mis_specified_All_data.RData"
## 0.121 sec elapsed
```

```r
# Comparison in well-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 3',
                      model_params_true = model_params_true,
                      true_model = NULL,
                      prior_model_params = prior_model_params,use_SoC_data = F)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_3-3.png)<!-- --> 

```
## For the rule-based design, 38% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 27% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
# Comparison in mis-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 3',
                      model_params_true = NULL,
                      true_model = true_model,
                      prior_model_params = prior_model_params, use_SoC_data = F)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_3-4.png)<!-- --> 

```
## For the rule-based design, 28% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 25% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

## Simulation 4


```r
true_alpha_eff = logit(1/100)
true_TED = 30  # simulation truth for the MED 
true_beta_eff = solve_beta(alpha_val = true_alpha_eff, v_star = true_TED, y_star = TEL)


true_alpha_tox = logit(1/1000)  # toxicity at 1 vial
true_MTD = 60  # simulation truth for the MTD 
true_beta_tox = solve_beta(alpha_val = true_alpha_tox, v_star = true_MTD, y_star = MTT)

# well-specified simulation truth
model_params_true = list(beta_tox = true_beta_tox,
                         alpha_tox = true_alpha_tox,
                         beta_eff = true_beta_eff, 
                         alpha_eff = true_alpha_eff)

# mis-specified simulation truth
# mis-specified simulation truth
true_model = make_linear_dose_response(MTT, true_MTD, tox_zero=1,
                                       TEL, true_TED, eff_zero=1)

run_all_scenarios(model_params_true = model_params_true,
                  true_model = true_model,
                  prior_model_params = prior_model_params,
                  N_trials = N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = 'Simulation scenario 4',
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  starting_dose = starting_dose,
                  SoC = SoC,
                  use_SoC_data = T)
```

```
## Simulation scenario 4 , model based design, well-specified, ...
## [1] "Simulation scenario 4_model_based_well_specified_All_data.RData"
## 0.064 sec elapsed
## Simulation scenario 4 , model based design, mis-specified, ...
## [1] "Simulation scenario 4_model_based_mis_specified_All_data.RData"
## 0.066 sec elapsed
## Simulation scenario 4 , rule based design, well-specified, ...
## [1] "Simulation scenario 4_rule_based_well_specified_All_data.RData"
## 0.131 sec elapsed
## Simulation scenario 4 , rule based design, mis-specified, ...
## [1] "Simulation scenario 4_rule_based_mis_specified_All_data.RData"
## 0.119 sec elapsed
```

```r
# Comparison in well-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 4',
                      model_params_true = model_params_true,
                      true_model = NULL,
                      prior_model_params = prior_model_params, use_SoC_data = T)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_4-1.png)<!-- --> 

```
## For the rule-based design, 13% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 65% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
# Comparison in mis-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 4',
                      model_params_true = NULL,
                      true_model = true_model,
                      prior_model_params = prior_model_params, use_SoC_data = T)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_4-2.png)<!-- --> 

```
## For the rule-based design, 88% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 3% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
run_all_scenarios(model_params_true = model_params_true,
                  true_model = true_model,
                  prior_model_params = prior_model_params,
                  N_trials = N_trials,
                  MTT = MTT, TEL = TEL,
                  N_max = N_max, 
                  N_batch = N_batch,
                  max_increment = max_increment, 
                  Randomisation_p_SOC = Randomisation_p_SOC,
                  sim_title = 'Simulation scenario 4',
                  FORCE_RERUN=FORCE_RERUN, 
                  N_cores = N_cores, 
                  starting_dose = starting_dose,
                  SoC = SoC,
                  use_SoC_data = F)
```

```
## Simulation scenario 4 , model based design, well-specified, ...
## [1] "Simulation scenario 4_model_based_well_specified_Adaptive_data.RData"
## 0.053 sec elapsed
## Simulation scenario 4 , model based design, mis-specified, ...
## [1] "Simulation scenario 4_model_based_mis_specified_Adaptive_data.RData"
## 0.055 sec elapsed
## Simulation scenario 4 , rule based design, well-specified, ...
## [1] "Simulation scenario 4_rule_based_well_specified_All_data.RData"
## 0.103 sec elapsed
## Simulation scenario 4 , rule based design, mis-specified, ...
## [1] "Simulation scenario 4_rule_based_mis_specified_All_data.RData"
## 0.104 sec elapsed
```

```r
# Comparison in well-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 4',
                      model_params_true = model_params_true,
                      true_model = NULL,
                      prior_model_params = prior_model_params, use_SoC_data = F)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_4-3.png)<!-- --> 

```
## For the rule-based design, 13% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 56% of trials give patient 252 a dose within +/-10% of the true optimal dose
```

```r
# Comparison in mis-specified case
compare_rule_vs_model(sim_title = 'Simulation scenario 4',
                      model_params_true = NULL,
                      true_model = true_model,
                      prior_model_params = prior_model_params, use_SoC_data = F)
```

![](Sample_Size_Simulation_Study_files/figure-latex/simulation_scenario_4-4.png)<!-- --> 

```
## For the rule-based design, 88% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 19% of trials give patient 252 a dose within +/-10% of the true optimal dose
```


# Idiosyncratic toxicity

This simulation assumes a dose-independent probability of having a toxic event. 
We simulate a trial whereby the overall toxicity is 15%


```r
model_params_true = NULL
tic()
Full_Simulation(model_params_true = NULL,
                true_model = true_model,
                prior_model_params = prior_model_params,
                N_trials = N_trials,
                MTT = MTT, TEL = TEL,
                N_max = N_max, 
                N_batch = N_batch,
                max_increment = max_increment, 
                Randomisation_p_SOC = Randomisation_p_SOC,
                sim_title = 'Simulation scenario 5',
                FORCE_RERUN=FORCE_RERUN, 
                N_cores = N_cores, 
                design_type = 'model_based',
                starting_dose = starting_dose,
                SoC = SoC,use_SoC_data = T)
```

```
## [1] "Simulation scenario 5_model_based_mis_specified_All_data.RData"
```

```r
Full_Simulation(model_params_true = NULL,
                true_model = true_model,
                prior_model_params = prior_model_params,
                N_trials = 10*N_trials,
                MTT = MTT, TEL = TEL,
                N_max = N_max, 
                N_batch = N_batch,
                max_increment = max_increment, 
                Randomisation_p_SOC = Randomisation_p_SOC,
                sim_title = 'Simulation scenario 5',
                FORCE_RERUN=FORCE_RERUN, 
                N_cores = N_cores, 
                design_type = 'rule_based',
                starting_dose = starting_dose,
                SoC = SoC, use_SoC_data = T)
```

```
## [1] "Simulation scenario 5_rule_based_mis_specified_All_data.RData"
```

```r
toc()
```

```
## 205.767 sec elapsed
```

```r
# # Comparison in mis-specified idiosyncratic case
compare_rule_vs_model(sim_title = 'Simulation scenario 5',
                      model_params_true = NULL,
                      true_model = true_model,
                      prior_model_params = prior_model_params,use_SoC_data = T,
                      idiosyncratic = T,plot_vstar = F)
```

![](Sample_Size_Simulation_Study_files/figure-latex/idiosyncratic-1.png)<!-- --> 

```
## For the rule-based design, 88% of trials give patient 252 a dose within +/-10% of the true optimal dose
## For the model-based design, 0% of trials give patient 252 a dose within +/-10% of the true optimal dose
```


