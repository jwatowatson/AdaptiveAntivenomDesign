plot_prior_versus_truth = function(model_params_true,
                                   true_model,
                                   prior_model_params,
                                   sim_title, plot_vstar = T, MTT, TEL){
  
  if(is.null(model_params_true) & is.null(true_model)) {
    stop('One of model_params_true or true_model has to be specified (plot_prior_versus_truth)')
  }
  if(is.null(model_params_true)) { 
    mis_specified=TRUE
  } else {
    mis_specified=FALSE
  }
  
  # ** Set up parameters **
  K = 1000
  xs = seq(1,5, length.out = K)
  N_samples = 5000
  ys_tox = array(dim = c(N_samples, K))
  ys_eff = array(dim = c(N_samples, K))
  mypallete = brewer.pal(8,'Dark2')
  
  ## *** Plot the prior versus simulation truth ***
  if(mis_specified){
    ys_truth_eff=true_model$eff(2^xs)
    ys_truth_tox=true_model$tox(2^xs)
  } else {
    ys_truth_eff=pnorm(q = 2^xs, mean = model_params_true$mu_antivenom, 
                       sd = model_params_true$sd_antivenom)
    ys_truth_tox=inv.logit(model_params_true$alpha_tox + model_params_true$beta_tox*xs)
  }
  plot(10 * 2^xs, 100*ys_truth_eff, type='l',xlab = 'Dose (mL)', ylab = 'Probability of outcome (%)', 
        ylim=c(0,100),lwd=2, col=mypallete[3], yaxt='n')
  axis(2, at = c(5,20,50,80,95))
  
  lines(10 * 2^xs, 100*ys_truth_tox,lwd=2,col=mypallete[1])
  abline(h = 100*c(TEL, MTT))
  # plot prior guess
  lines(10 * 2^xs, 100*inv.logit(Prior_alpha_tox + Prior_beta_tox*xs),
        col=mypallete[5],lty=3,lwd=2)
  for(i in 1:N_samples){
    alpha_tox = rnorm(1, prior_model_params$alpha_tox, prior_model_params$alpha_tox_sd)
    beta_tox = rnorm(1, prior_model_params$beta_tox, prior_model_params$beta_tox_sd)
    ys_tox[i,] = inv.logit(alpha_tox + beta_tox*xs)
  }
  xs=10*2^xs
  polygon(x = c(xs,rev(xs)), y = 100*c(apply(ys_tox,2,quantile,prob=0.025),
                                       rev(apply(ys_tox,2,quantile,prob=0.975))), 
          col = adjustcolor(mypallete[5],alpha.f = .2), border = NA)
  
  lines(xs, 100*pnorm(q = xs/10, mean = prior_model_params$mu_antivenom, 
                  sd = prior_model_params$sd_antivenom),lwd=2,
        col=mypallete[8],lty=3)
  for(i in 1:N_samples){
    mymu = rnorm(1, prior_model_params$mu_antivenom, prior_model_params$mu_antivenom_sd)
    mysd = rtruncnorm(1, prior_model_params$sd_antivenom, prior_model_params$sd_antivenom_sd,a = 0)
    ys_eff[i,] = pnorm(q = xs/10, mean = mymu, sd = mysd)
  }
  polygon(x = c(xs,rev(xs)), y = 100*c(apply(ys_eff,2,quantile,prob=0.025,na.rm=T),
                                       rev(apply(ys_eff,2,quantile,prob=0.975,na.rm=T))), 
          col = adjustcolor(mypallete[8],alpha.f = .2), border = NA)
  out = Estimate_log2_Vstar(model_params = model_params_true, 
                            true_model = true_model, 
                            MTT = MTT, TEL = TEL)
  out_prior = Estimate_log2_Vstar(model_params = prior_model_params, 
                                  true_model = NULL, 
                                  MTT = MTT, TEL = TEL)
  if(plot_vstar){
    abline(v = 10 * 2^out$Vstar, col='red',lwd=3)
    abline(v = 10 * 2^out_prior$Vstar, col='red',lwd=2,lty=3)
  }
  
  
  legend('topleft',col= c(mypallete[c(3,8,1,5)],'red','red'),bg = 'white',bty='y',cex = .8,
         legend = c('True efficacy','Prior efficacy estimate',
                    'True toxicity', 'Prior toxicity estimate', 
                    'True optimal dose', 'Prior optimal dose'),
         lwd=2,lty=c(1,3,1,3,1,3), inset = 0.02)
  
}


compare_rule_vs_model = function(sim_title,
                                 true_model,
                                 model_params_true,
                                 prior_model_params, 
                                 use_SoC_data,
                                 idiosyncratic = F, 
                                 plot_vstar = T, MTT, TEL){
  
  if(is.null(model_params_true) & is.null(true_model)) {
    stop('One of model_params_true or true_model has to be specified (compare_rule_vs_model)')
  }
  if(is.null(model_params_true) & !is.null(true_model)) { 
    specification_type='mis_specified'
  } 
  if(!is.null(model_params_true) & is.null(true_model)) {
    specification_type='well_specified'
  }
  if(!is.null(model_params_true) & !is.null(true_model)) {
    stop('ambiguous specification')
  }
  if(use_SoC_data) {
    SoC_use = 'All_data'
  } else {
    SoC_use = 'Adaptive_data'
  }
  
  # load saved data
  f_name_rule = paste('SimulationOutputs/', sim_title, '_', 'rule_based', '_', specification_type, '_', 'All_data', '.RData', sep = '')
  f_name_model = paste('SimulationOutputs/',sim_title,'_','model_based', '_', specification_type, '_', SoC_use, '.RData', sep = '')
  load(f_name_model)
  Summary_model = Summary_trials
  load(f_name_rule)
  Summary_rule = Summary_trials
  rm(Summary_trials)
  
  # uncertainty interval quantiles (upper and lower)
  qu=.95
  ql=.05
  
  par(mfrow=c(2,2), las=1, bty='n', family = 'serif', cex.lab=1.5, cex.axis=1.5, mar=c(5,5,4,2))
  # ******* Plot the prior versus the simulation truth *********
  plot_prior_versus_truth(model_params_true = model_params_true,
                          true_model = true_model,
                          prior_model_params = prior_model_params,
                          sim_title = sim_title, plot_vstar = plot_vstar,
                          MTT = MTT, TEL = TEL)
  mtext(text = 'a',side = 3,line = 2,at = 0,cex = 2)
  
  # ******* Plot the estimated optimal dose from the model-based design *******
  mypallete = brewer.pal(8,'Dark2')
  cols_MTD = grep('MTD', colnames(Summary_model))
  cols_TED = grep('TED', colnames(Summary_model))
  Summary_MTD = Summary_model[, cols_MTD]
  Summary_TED = Summary_model[, cols_TED]
  q_upperMTD = log2(apply(Summary_MTD, 1, quantile, prob=ql, na.rm=T))
  q_lowerMTD = log2(apply(Summary_MTD, 1, quantile, prob=qu, na.rm=T))
  q_upperTED = log2(apply(Summary_TED, 1, quantile, prob=ql, na.rm=T))
  q_lowerTED = log2(apply(Summary_TED, 1, quantile, prob=qu, na.rm=T))
  plot(1:nrow(Summary_MTD), log2(Summary_MTD[,1]),
       type='n', main = '', yaxt='n',
       ylim= range(c(q_upperMTD,q_lowerMTD,q_upperTED,q_lowerTED)),
       xlab = 'Number of patients enrolled', ylab = 'Dose (mL)')
  my_ys = log2(c(10,100, 200, 400, 800, 1200)/10)
  axis(2, at = my_ys, labels = 10 * 2^(my_ys))
  out = Estimate_log2_Vstar(model_params = model_params_true,
                            true_model = true_model,
                            MTT = MTT,TEL = TEL)
  vstar = 10 * 2^out$Vstar
  
  if(plot_vstar) abline(h = out$MTD,col=mypallete[1], lwd=2)
  abline(h = out$TED, col = mypallete[3], lwd=2)
  
  polygon(x = c(1:nrow(Summary_MTD),rev(1:nrow(Summary_MTD))), 
          y = c(q_upperMTD, rev(q_lowerMTD)), 
          col = adjustcolor(mypallete[5],alpha.f = .2), border = NA)
  lines(1:nrow(Summary_MTD), log2(apply(Summary_MTD, 1, mean, na.rm=T)),
        lwd=3,col=mypallete[5],lty=3)
  
  polygon(x = c(1:nrow(Summary_TED),rev(1:nrow(Summary_TED))), 
          y = c(q_upperTED, rev(q_lowerTED)), 
          col = adjustcolor(mypallete[8],alpha.f = .2), border = NA)
  lines(1:nrow(Summary_TED), log2(apply(Summary_TED, 1, mean, na.rm=T)),
        lwd=3,col=mypallete[8],lty=3)
  
  legend('topright', legend = c('true MTD','true TED', 'estimated MTD','estimated TED'), 
         col=mypallete[c(1,3,5,8)], lty=c(1,1,3,3),
         lwd = c(2,2,3,3), inset=0.01, bty='y', bg = 'white')
  mtext(text = 'b',side = 3,line = 2,at = 0,cex = 2)
  
  # ******* Plot the assigned doses *******
  mypallete = brewer.pal(9,'Set1')[c(1,2,7)]
  
  cols_assigned_dose_model = grep('assigned_dose', colnames(Summary_model))
  cols_assigned_dose_rule = grep('assigned_dose', colnames(Summary_rule))
  Summary_assigned_model = 10*Summary_model[, cols_assigned_dose_model]
  Summary_assigned_rule = 10*Summary_rule[, cols_assigned_dose_rule]
  
  q_upper_model = apply(Summary_assigned_model, 1, quantile, prob=ql, na.rm=T)
  q_lower_model = apply(Summary_assigned_model, 1, quantile, prob=qu, na.rm=T)
  q_upper_rule = apply(Summary_assigned_rule, 1, quantile, prob=ql, na.rm=T)
  q_lower_rule = apply(Summary_assigned_rule, 1, quantile, prob=qu, na.rm=T)
  
  plot(NA, NA, col = adjustcolor('grey',alpha.f = .5),
       type='n', main = '', 
       ylim=range(c(q_upper_rule,q_lower_rule,q_upper_model,q_lower_model)),
       xlim=c(1,nrow(Summary_assigned_model)),
       xlab = 'Number of patients enrolled', ylab = 'Dose (mL)')
  mtext(text = 'c',side = 3,line = 2,at = 0,cex = 2)
  
  abline(h = vstar,col=mypallete[1],lwd=3)
  
  polygon(x = c(1:nrow(Summary_assigned_model),rev(1:nrow(Summary_assigned_model))), 
          y = c(q_upper_model, rev(q_lower_model)), 
          col = adjustcolor(mypallete[2],alpha.f = .2), border = NA)
  polygon(x = c(1:nrow(Summary_assigned_model),rev(1:nrow(Summary_assigned_model))), 
          y = c(q_upper_rule, rev(q_lower_rule)), 
          col = adjustcolor(mypallete[3],alpha.f = .2), border = NA)
  
  lines(1:nrow(Summary_assigned_model), apply(Summary_assigned_model, 1, mean, na.rm=T),
        col=mypallete[2], lwd=2)
  
  lines(1:nrow(Summary_assigned_rule), apply(Summary_assigned_rule, 1, mean, na.rm=T),
        col=mypallete[3], lwd=2)
  
  if(use_SoC_data) {
    SoC_legend = ' (all data)'
  } else {
    SoC_legend = ' (adaptive data)'
  }
  legend('topright', legend = c(paste('Model based',SoC_legend),'Rule based'), 
         title = 'Assigned dose',
         lwd = 2, lty=1, inset=0.01, bty='n', bg = 'white', col=mypallete[c(2,3)])
  
  # ******* Compare final doses: histogram *******
  # Histogram of final assigned dose
  Summary_assigned_model[nrow(Summary_assigned_model), ] = apply(Summary_assigned_model,2,function(x){
    y= rev(x[complete.cases(x)])
    return(y[1])
  } )
  Summary_assigned_rule[nrow(Summary_assigned_rule), ] = apply(Summary_assigned_rule,2,function(x){
    y= rev(x[complete.cases(x)])
    return(y[1])
  } )
  last_dose_model = Summary_assigned_model[nrow(Summary_assigned_model),]
  last_dose_rule = Summary_assigned_rule[nrow(Summary_assigned_rule),]
  min_v = min(c(last_dose_model,last_dose_rule))
  max_v = max(c(last_dose_model,last_dose_rule))
  my_breaks = seq(min_v-0.5,max_v+0.5, length.out = (max_v-min_v+2)/10)
  h1 = hist(x = last_dose_model, 
            breaks = my_breaks, plot = F)
  h2 = hist(x = last_dose_rule, 
            breaks = my_breaks,plot = F)
  
  hist(x = last_dose_model, 
       xlab = 'Dose assigned to last patient (mL)', 
       ylab = '',  main = '', freq = F,
       col=adjustcolor(col = mypallete[2],alpha.f = .3),yaxt='n',
       breaks = my_breaks, ylim=c(0,max(c(h1$density,h2$density))))
  hist(x = last_dose_rule, add=T, freq = F,
       col=adjustcolor(col = mypallete[3],alpha.f = .3),
       breaks = my_breaks)
  abline(v=vstar, col=mypallete[1],lwd=3)
  mtext(text = 'd',side = 3,line = 2,at = min_v,cex = 2)
  
  legend('topright', legend = c('Model based','Rule based'), title = 'Final dose',
        inset=0.01, bty='n', bg = 'white', 
        fill= adjustcolor(mypallete[c(2,3)],alpha.f = .3))
  
  within10percent_rule = round(100*mean(last_dose_rule >= vstar-.1*vstar & last_dose_rule <= vstar+.1*vstar))
  within10percent_model = round(100*mean(last_dose_model >= vstar-.1*vstar & last_dose_model <= vstar+.1*vstar))
  writeLines(sprintf('For the rule-based design, %s%% of trials give patient %s a dose within +/-10%% of the true optimal dose',
                     within10percent_rule, nrow(Summary_assigned_model)))
  writeLines(sprintf('For the model-based design, %s%% of trials give patient %s a dose within +/-10%% of the true optimal dose',
                     within10percent_model, nrow(Summary_assigned_model)))
}