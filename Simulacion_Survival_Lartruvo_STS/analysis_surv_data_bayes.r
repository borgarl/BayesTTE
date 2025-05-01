
analysis_surv_data_bayes <- function(simul_data, alpha, test.type, prior, prior_type, 
                                                modelo = modelo, IA = NULL, method_IA = NULL,
                                                P_HR_data_Boundary = 1, n_exp_events = NULL, modelo_bayes_test,
                                                shape_parameter, prior_gamma, parameters, w = NULL,
                                                external_data = NULL, a0 = NULL, Dynamic_Borrowing_PP = FALSE){
 
  
  
  
  # Mirar la carpeta mixture prior para mirar referencias.
  
  True_HR <- (parameters[1] / parameters[2])^shape_parameter
  res_df <- as.data.frame(matrix(0, nrow = length(prior_type), ncol = 15)) 
  
  colnames(res_df) <- c("Model", "Prior", "Mean_Posterior", "Lower_Cred_Int", "Upper_Cred_Int",
                        "Median_Posterior", "Reject", "1-post_prob", "reject_Cred_Int","Test", "P_HR_Data", 
                        "Coverage_Probability", "MSE", "WAIC_Est", "WAIC_SE")
  
  lista_prior <- vector(mode = "list")
  
  
  if(modelo_bayes_test == "weibull") {
    
    
    sample_size <- dim(simul_data)[1]
    
    if(!is.null(IA)){ boundaries <- get_boundaries_IA(IA, alpha, method_IA = method_IA,
                                                      test.type = test.type,
                                                      sample_size, n_exp_events) }else{boundaries = alpha}
    
    
    n_analyses <- ifelse(is.null(IA), 1, length(IA) + 1)
    sum_events <- sum(simul_data$status)
    
    if(!is.null(w)){ lista_prior$w <- w } 
    
    if(!is.null(a0)){
      
      lista_prior <- list(N_external = dim(external_data)[1],
                          t_external = as.integer(external_data$arm == 2), # solo para brazo control
                          y_external = external_data$time,
                          v_external = external_data$status,
                          a0 = a0)

    } 

    res <- NULL
    reject <- 0
    reject_Cred_Int <- 0
    Test <- "Not significant"
    
    # Pasamos inputs para el modelo Bayesiano
    
    lista_prior$alpha_shape0 <- 0.1
    lista_prior$beta_shape0 <- 4
    lista_prior$alpha_scale0_control = prior_gamma[1,1]
    lista_prior$beta_scale0_control = prior_gamma[1,2]
    lista_prior$alpha_scale0_treatment = prior_gamma[2,1]
    lista_prior$beta_scale0_treatment = prior_gamma[2,2]
    
    for(j in 1:n_analyses) {
      
      if (is.null(IA) || j == n_analyses) {
        IA_data <- simul_data
      } else {
        events_IA <- round(IA[j] * sum_events)
        IA_data <- simul_data[cumsum(simul_data$status) <= events_IA,]
      }
      
      lista_prior$N <- dim(IA_data)[1]
      lista_prior$t <- as.integer(IA_data$arm == 2)
      lista_prior$y <- IA_data$time
      lista_prior$v <- IA_data$status
      
      if (Dynamic_Borrowing_PP == TRUE){
        
        IA_data_control <- subset(IA_data, arm == 1)
        IA_data_control$dataset <- 'current'

        combined_data <- rbind(IA_data_control, external_data)
        combined_data$status <- as.numeric(combined_data$status == "TRUE")
        combined_data$dataset_bin <- as.integer(combined_data$dataset == "current")
        
        # Hacemos un cox para evaluar la similitud entre los 2 datasets
        cox_model <- coxph(Surv(time, status) ~ as.factor(dataset), data = combined_data)
        
        # Estas son las condiciones que vamos a considerar para el préstamo de evidencia en función
        # de cuánto se parezcan los datos del ensayo y los del histórico.
        
        # De momento lo dejo incrustado en la función principal aunque más adelante se podría poner como input.
        
        get_borrowing_strength <- function(HR) {
          if (HR >= 0.7 && HR <= 0.75) {
            return(0.05)
          } else if (HR > 0.75 && HR < 0.8) {
            return(0.1)
          } else if (HR >= 0.8 && HR <= 0.85) {
            return(0.2)
          } else if (HR > 0.85 && HR < 0.9) {
            return(0.3)
          }  else if (HR >= 0.9 && HR <= 0.95) {
            return(0.4)
          } else if (HR > 0.95 && HR < 1.05) {
            return(0.5)
          } else if (HR >= 1.05 && HR <= 1.1) {
            return(0.4)
          } else if (HR > 1.1 && HR < 1.15) {
            return(0.3)
          }else if (HR > 1.15 && HR < 1.2) {
            return(0.2)
          }  else if (HR >= 1.2 && HR <= 1.25) {
            return(0.1)
          } else if (HR > 1.25 && HR < 1.3) {
            return(0.05)
          } else {
            return(0)
          }
        }
        
        HR <- exp(coef(cox_model)['as.factor(dataset)external'])
        
        # Metemos el valor escogido en función de la condición para la power prior
        
        if(!is.null(a0)){
          
         lista_prior$a0 <- get_borrowing_strength(HR)
            
        }else{
        
          lista_prior$N_external <- dim(external_data)[1]
          lista_prior$t_external = as.integer(external_data$arm == 2) # solo para brazo control
          lista_prior$y_external = external_data$time
          lista_prior$v_external = external_data$status
          lista_prior$a0 = get_borrowing_strength(HR)
          
        }
        #lista_prior$a0 <- get_borrowing_strength(HR)
        
        # Esto es para sacar el HR entre los dos datasets pero con el modelo bayesiano que he hecho. 
        # Como he dicho en el informe, no lo hago aquí por temas de tiempo de computación (sería el doble)
        # Para tener el mismo que resultado que usando el frecuentista (comprobado).
        
        # stan_data <- list(
        #   N = nrow(combined_data),
        #   status = as.integer(combined_data$status),
        #   time = combined_data$time,
        #   dataset =  combined_data$dataset_bin
        # )
        # 
        # fit_CPP <- sampling(stan_model, data = stan_data, iter = 2000, chains = 3)
        # print(fit_CPP, digits = 6)
        
      }
      
      
      f <- rstan::sampling(
        object  = modelo,
        data    = lista_prior,
        chains  = 3,       
        thin    = 3,
        iter    = 3400,   
        warmup  = 1900,
        refresh = 0,
        seed = 12
      )
      
      print(f, pars = c("shape", "scale_control", "scale_treatment", "HR", "lp__"), digits = 3, probs = c(0.025, 0.5, 0.975))
      
      #print(f, pars = c("used_a0","shape", "scale_control", "scale_treatment", "HR", "lp__"), digits = 3, probs = c(0.025, 0.5, 0.975))
      
      
      # Esta funcion sera importante para presentar resultados pero para simulacion no es necesaria
      # 
      # plots_bayes_sep <- plot_bayes(f, lista_prior, shape_parameter, parameters)
      # 
      # plots_bayes <- plot_bayes(f, lista_prior, shape_parameter, parameters)
      # traceplot <- plots_bayes[[1]]
      # plot_posterior_shape <- plots_bayes[[2]]
      # plot_posterior_scale <- plots_bayes[[3]]
      # plot_posterior_slope <- plots_bayes[[4]]
      # plot_dist_prior_like_posterior <- plots_bayes[[5]] # no tiene sentido si solo dibujo la posterior de la slope en nada se parece a la conjunta. Mirar parece que esta mal.
      
      cuantiles <- c((alpha/test.type), 1-(alpha/test.type))
      post_shape = extract(f, "shape")$shape
      post_scale_control = extract(f, "scale_control")$scale_control
      post_scale_experimental = extract(f, "scale_treatment")$scale_treatment

      post_hr = (extract(f, "HR")$HR)
      cred_int <- quantile(post_hr, c(cuantiles[1], cuantiles[2])) # Ponemos todos los resultados ya en esos terminos
      
      # Calcular P(HR < P_HR_data_Boundary | data).
      # Esto no define el criterio de parada, solo si esperamos un efecto, cual es la P de la dist a posteriori.
      # Se hace esto para el IA pero siempre HR<1 (1 = no eficacia)
      P_HR_data <- mean(post_hr < P_HR_data_Boundary)
      # Parar el estudio por eficacia si P(HR < 1 | data) > 0.995 por ejemplo (1-BOUNDARY)
      # Desde un punto de vista Bayesiano no tiene sentido ajustar por el nº de veces que ves los datos
      # Aun asi dejo la opcion de ajustar la P posterior condicionado a los datos de que se pare por eficacia
      
      # Tambien dejo una opcion donde no se ajuste (como tiene que ser) con la opcion method_IA == "Bayes"
      
      # Aunque esta condicion la he hecho para evaluar cada IA (es mas facil tambien ajustar por alfa),
      # lo dejo comentado y evaluo cada IA en funcion de si el intervalo de credibilidad al 95% contiene el 1 o no
      
      Stop_efficacy <- mean(post_hr < 1) 
      one_minus_post_prob <- 1-Stop_efficacy
      
      # Esto discrimina mejor por lo que he visto en la comparacion de los p-values. Ahora me fio mas de esta medida ya que
      # el bayesiano no informativo es igual que el frecuentista.
      
      # La unica razon por la que no pongo el break aqui (solo afecta si hay IAs) es para poder tener el "pvalor" bayesiano
      # Y para ajustar de alguna manera el control del "ET1" en caso de IAs.

      if ((cred_int[1] < 1 && cred_int[2] > 1) || (cred_int[1] > 1 && cred_int[2] > 1)) {
        reject_Cred_Int <- 0
      } else {
        reject_Cred_Int <- 1
      }
      
      if(Stop_efficacy > 1-boundaries[j]/test.type) {
        reject <- 1
        Test <- ifelse(j < n_analyses, paste("IA", j), "FA")
        break
      } else {
        reject <- 0
      }
      
    }
    
    Estimate <- c(mean(post_hr), cred_int[1], cred_int[2], median(post_hr))
    res <- c(Estimate, reject, reject_Cred_Int, one_minus_post_prob, Test, P_HR_data) 
    # Esto no es el coverage probability si no otra manera de ver si es significativo
    # tengo que cambiarlo porque son las veces que el verdadera HR esta contenido en el intervalo
    # y no si contiene el 1 o no.
    Coverage_Probability <- ifelse((cred_int[1] > 1) || (cred_int[2] < 1), 1, 0)
    MSE <- (as.numeric(res[1]) - True_HR)^2
    log_lik <- extract(f, "log_lik")$log_lik
    waic <- waic(log_lik)
    
    res_df[1] <- "Weibull"
    res_df[2] <- prior_type
    res_df[3] <- Estimate[1]
    res_df[4] <- Estimate[2]
    res_df[5] <- Estimate[3]
    res_df[6] <- Estimate[4]
    res_df[7] <- reject
    res_df[8] <- one_minus_post_prob
    res_df[9] <- reject_Cred_Int
    res_df[10] <- Test
    res_df[11] <- P_HR_data
    res_df[12] <- Coverage_Probability
    res_df[13] <- MSE
    res_df[14] <- waic$estimates["waic", "Estimate"]
    res_df[15] <- waic$estimates["waic", "SE"]
    
  }
  
  return(res_df) 

}


  