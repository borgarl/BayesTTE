
# Funcion para generar datos de supervivencia

surv_analysis_freq <- function(sample_size = sample_size[i], 
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
                          n_exp_events = n_exp_events){
  
  
  ############ 1) Aleatorizacion  ############ 
  
  # Compartimos freq con bayes ya que está igual
  #set.seed(seed)
  pat <- rand(sample_size, ratio, rand_type = rand_type)
  
  
  ######## 2) Generar tiempos supervivencia brazos #############
  
  # Compartimos freq con bayes ya que está igual
  #set.seed(seed)
  
  simul_data <- gen_surv_data(pat, parameters,
                              shape_parameter,
                              censor, Tmax)

 
  ######## 3) Fit regresion de Cox y ver si se rechaza la H0 (INCLUYE IAs) #############


  analysis <- as.data.frame(analysis_surv_data_freq(simul_data, alpha = alpha, 
                                             test.type = test.type,
                                             IA = IA, method_IA = method_IA,
                                             sample_size = sample_size,
                                             n_exp_events = n_exp_events,
                                             parameters = parameters,
                                             shape_parameter = shape_parameter))
 
  return(analysis)
  
}