
# Funcion para generar datos de supervivencia

surv_analysis_bayes <- function(sample_size = sample_size[i], 
                          ratio = ratio,
                          rand_type = rand_type,
                          block_size = FALSE,
                          Tmax = Tmax, 
                          parameters = parameters,
                          shape_parameter = shape_parameter,
                          scenario_control = scenario_control[j],
                          scenario_treatment = scenario_treatment[j],
                          censor = censor,
                          alpha = alpha,
                          test.type = test.type,
                          IA = NULL, 
                          method_IA = NULL,
                          n_exp_events = n_exp_events,
                          prior,
                          modelo = modelo,
                          modelo_bayes_test = modelo_bayes_test,
                          prior_type,
                          P_HR_data_Boundary,
                          prior_gamma,
                          w = NULL,
                          external_data = NULL, 
                          a0 = NULL, 
                          Dynamic_Borrowing_PP = FALSE){
  
  
  ############ 1) Aleatorizacion  ############ 
  #set.seed(seed)
  pat <- rand(sample_size, ratio, rand_type = rand_type)
  
  
  ######## 2) Generar tiempos supervivencia brazos #############
  
  #set.seed(seed)
  
  simul_data <- gen_surv_data(pat, parameters,
                              shape_parameter,
                              censor, Tmax)

  ######## 3) Fit regresion de Cox y ver si se rechaza la H0 (INCLUYE IAs) #############

  
  analysis <- analysis_surv_data_bayes( simul_data = simul_data, 
                                        alpha = alpha,
                                        test.type = test.type, 
                                        method_IA = method_IA, 
                                        n_exp_events = n_exp_events,
                                        prior = prior, 
                                        prior_type = prior_type, 
                                        modelo = modelo, 
                                        modelo_bayes_test = modelo_bayes_test,
                                        IA = IA, 
                                        P_HR_data_Boundary = P_HR_data_Boundary,
                                        shape_parameter = shape_parameter,
                                        prior_gamma = prior_gamma,
                                        parameters = parameters,
                                        w = w,
                                        external_data = external_data,
                                        a0 = a0,
                                        Dynamic_Borrowing_PP = Dynamic_Borrowing_PP)

  return(analysis)
  
}