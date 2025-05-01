
# Funcion para simular ensayo clinico frecuentista. Aqui es donde se van a diferenciar
# las diferentes distribuciones (weibull, piecewise exponential)
# De momento solo Weibull (ahora no tiene mucho sentido pero sera util para diferenciar)

sim_trials <- function(n_sim = 1, 
                       analysis,
                           sample_size, 
                           ratio,
                           rand_type,
                           block_size = FALSE,
                           Tmax, 
                           #dist_type = "weibull",
                           scenarios_eff,
                           shape_parameter = 1,
                           censor,
                           alpha = 0.05,
                           test.type = 2,
                           IA = NULL,
                           method_IA = NULL,
                           n_exp_events = NULL,
                           seed = 1,
                           HR_1 = FALSE,
                           Plot_Power = FALSE,
                           Plot_Power_scenarios = FALSE,
                           plot_pvalues = FALSE,
                           Plot_Control_Scenarios = FALSE,
                           desired_HRs = NULL,
                           prior = NULL,
                           modelo = NULL,
                           modelo_bayes_test = NULL,
                           prior_type = NULL,
                           P_HR_data_Boundary = NULL,
                           prior_gamma = NULL,
                           w = NULL,
                           external_data = NULL, 
                           a0 = NULL, 
                           Dynamic_Borrowing_PP = FALSE){
  start<- Sys.time()
  
  set.seed(seed)
  

  if(analysis == "freq"){
    
  res <- sim_freq(n_sim = n_sim, 
                     sample_size = sample_size, 
                     ratio = ratio,
                     rand_type = rand_type,
                     Tmax = Tmax, 
                     scenarios_eff = scenarios_eff,
                     shape_parameter = shape_parameter,
                     censor = censor,
                     test.type = test.type,
                     alpha = alpha,
                     method_IA = method_IA,
                     IA = IA,
                     n_exp_events = n_exp_events,
                     HR_1 = HR_1,
                     Plot_Power = Plot_Power,
                     Plot_Power_scenarios = Plot_Power_scenarios,
                     plot_pvalues = plot_pvalues,
                     Plot_Control_Scenarios = Plot_Control_Scenarios,
                     desired_HRs = desired_HRs,
                     prior_gamma)
  
  # Esto es para calcular el tiempo de calculo
  
  end <- Sys.time()
  time <- end - start
  seconds <- as.numeric(time, units = "secs")
  time <- sprintf("%02d:%02d:%02d", as.integer(seconds %/% 3600),
                  as.integer((seconds %% 3600) %/% 60),
                  as.integer(seconds %% 60))
  cat("Tiempo de duracion:", time)
  
  return(res)
  
  }
  
  if(analysis == "bayes"){
    
    res <- sim_bayes(n_sim = n_sim, 
                       sample_size = sample_size, 
                       ratio = ratio,
                       rand_type = rand_type,
                       Tmax = Tmax, 
                       scenarios_eff = scenarios_eff,
                       shape_parameter = shape_parameter,
                       censor = censor,
                       test.type = test.type,
                       alpha = alpha,
                       method_IA = method_IA,
                       IA = IA,
                       n_exp_events = n_exp_events,
                       HR_1 = HR_1,
                       Plot_Power = Plot_Power,
                       prior = prior,
                       modelo = modelo,
                       modelo_bayes_test = modelo_bayes_test,
                       prior_type = prior_type,
                       P_HR_data_Boundary = P_HR_data_Boundary,
                       prior_gamma = prior_gamma,
                       Plot_Power_scenarios = Plot_Power_scenarios,
                       desired_HRs = desired_HRs,
                       plot_pvalues = plot_pvalues,
                       Plot_Control_Scenarios = Plot_Control_Scenarios,
                       w = w,
                       external_data = external_data,
                       a0 = a0,
                       Dynamic_Borrowing_PP = Dynamic_Borrowing_PP)
    
    # Esto es para calcular el tiempo de calculo
    
    end <- Sys.time()
    time <- end - start
    seconds <- as.numeric(time, units = "secs")
    time <- sprintf("%02d:%02d:%02d", as.integer(seconds %/% 3600),
                    as.integer((seconds %% 3600) %/% 60),
                    as.integer(seconds %% 60))
    cat("Tiempo de duracion:", time)
    
    return(res)
    
  }
}